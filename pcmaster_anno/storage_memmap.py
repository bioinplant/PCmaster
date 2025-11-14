# pcmaster_anno/storage_memmap.py
import os
import json
from typing import Tuple, Dict, Optional, List
import numpy as np
import pandas as pd


def _bool_series_from_var(var, key: str):
    """Return boolean mask for var[key] indicating inclusion.
    Accept 'yes','true','1','y' (case insensitive) or boolean True.
    """
    if key not in var.columns:
        return None
    s = var[key].astype(str).fillna("").str.lower()
    mask = s.isin(["yes", "true", "1", "y", "t"])
    # also treat literal 'True'/'False' converted to str above
    return mask.values


def _pad_or_truncate_matrix(mat: np.ndarray, target_cols: int, dtype=np.float32):
    n, c = mat.shape
    if c == target_cols:
        return mat.astype(dtype)
    if c < target_cols:
        pad = np.zeros((n, target_cols - c), dtype=dtype)
        return np.concatenate([mat.astype(dtype), pad], axis=1)
    # truncate (take first columns)
    return mat[:, :target_cols].astype(dtype)


def save_adata_to_memmap(
    adata,
    out_dir: str,
    rc_var_key: str = "anno_genes_rc",
    fc_var_key: str = "anno_genes_fc_hp_sa",
    target_rc: Optional[int] = None,
    target_fc: int = 2000,
    dtype=np.float32,
) -> Tuple[Dict, Dict]:
    """
    Create two memmap datasets:
      - features_rc.dat: used by PCMA_MYGO_RC (image-like, default target_rc = 224*224)
      - features_fc.dat: used by PCMA_MYGO_FC / SA / HP (vector-like, target_fc default 2000)
    Also save:
      - metadata.csv (global)
      - metadata.train.csv / metadata.valid.csv / metadata.test.csv (if adata.obs['tvt'] exists)
      - mappings.json (contains int2obs, obs2int, gene lists used)

    Returns: info dict, mappings dict
    """
    os.makedirs(out_dir, exist_ok=True)

    obs = getattr(adata, "obs", None)
    if obs is None:
        raise ValueError("adata must have .obs")

    obs = obs.copy()
    if "cellname" not in obs.columns:
        if hasattr(adata, "obs_names"):
            obs["cellname"] = list(adata.obs_names)
        else:
            obs["cellname"] = [str(i) for i in range(len(obs))]

    if "celltype" not in obs.columns:
        obs["celltype"] = obs.get("cell_type", "unknown")

    # global metadata
    meta_df = pd.DataFrame(
        {
            "cellname": list(obs["cellname"].astype(str)),
            "celltype_orig": list(obs["celltype"].astype(str)),
        }
    )

    # preserve tvt if present
    if "tvt" in obs.columns:
        meta_df["tvt"] = list(obs["tvt"].astype(str))

    unique = pd.Series(meta_df["celltype_orig"]).fillna("unknown").astype(str).unique().tolist()
    map_obs2int = {lab: idx for idx, lab in enumerate(unique)}
    map_int2obs = {v: k for k, v in map_obs2int.items()}
    meta_df["celltype_int"] = meta_df["celltype_orig"].map(map_obs2int)

    # X -> dense numpy
    X = adata.X
    try:
        Xmat = np.asarray(X.todense())
    except Exception:
        Xmat = np.asarray(X)

    n, total_feat = Xmat.shape

    # select rc and fc gene indices using var annotations (fallback behavior)
    var = getattr(adata, "var", None)
    rc_mask = None
    fc_mask = None
    rc_genes: List[str] = []
    fc_genes: List[str] = []
    if var is not None:
        rc_mask = _bool_series_from_var(var, rc_var_key)
        fc_mask = _bool_series_from_var(var, fc_var_key)
        # if mask not None but all False, treat as None
        if rc_mask is not None and rc_mask.sum() == 0:
            rc_mask = None
        if fc_mask is not None and fc_mask.sum() == 0:
            fc_mask = None

    var_names = list(getattr(adata, "var_names", list(range(total_feat))))
    if rc_mask is not None:
        rc_idx = np.where(rc_mask)[0].tolist()
        rc_genes = [var_names[i] for i in rc_idx]
        X_rc = Xmat[:, rc_idx]
    else:
        # fallback: use all genes for RC if user didn't annotate
        X_rc = Xmat.copy()
        rc_genes = var_names.copy()

    if fc_mask is not None:
        fc_idx = np.where(fc_mask)[0].tolist()
        fc_genes = [var_names[i] for i in fc_idx]
        X_fc = Xmat[:, fc_idx]
    else:
        # fallback: use all genes for FC then we'll truncate/pad to target_fc
        X_fc = Xmat.copy()
        fc_genes = var_names.copy()

    # set default rc target if None
    if target_rc is None:
        target_rc = 224 * 224

    # pad/truncate
    X_rc_prepared = _pad_or_truncate_matrix(X_rc, target_rc, dtype=dtype)
    X_fc_prepared = _pad_or_truncate_matrix(X_fc, target_fc, dtype=dtype)

    # write memmaps
    rc_path = os.path.join(out_dir, "features_rc.dat")
    fp_rc = np.memmap(rc_path, dtype=dtype, mode="w+", shape=(n, target_rc))
    fp_rc[:] = X_rc_prepared[:]
    del fp_rc

    fc_path = os.path.join(out_dir, "features_fc.dat")
    fp_fc = np.memmap(fc_path, dtype=dtype, mode="w+", shape=(n, target_fc))
    fp_fc[:] = X_fc_prepared[:]
    del fp_fc

    # save global metadata and splits if tvt present
    meta_csv = os.path.join(out_dir, "metadata.csv")
    meta_df.to_csv(meta_csv, index=False)

    train_csv = valid_csv = test_csv = None
    if "tvt" in meta_df.columns:
        train_df = meta_df[meta_df["tvt"] == "train"].reset_index(drop=True)
        valid_df = meta_df[meta_df["tvt"] == "valid"].reset_index(drop=True)
        test_df = meta_df[meta_df["tvt"] == "test"].reset_index(drop=True)

        train_csv = os.path.join(out_dir, "metadata.train.csv")
        valid_csv = os.path.join(out_dir, "metadata.valid.csv")
        test_csv = os.path.join(out_dir, "metadata.test.csv")

        train_df.to_csv(train_csv, index=False)
        valid_df.to_csv(valid_csv, index=False)
        test_df.to_csv(test_csv, index=False)

    mappings = {
        "obs2int": map_obs2int,
        "int2obs": map_int2obs,
        "rc_genes": rc_genes,
        "fc_genes": fc_genes,
        "rc_memmap": os.path.basename(rc_path),
        "fc_memmap": os.path.basename(fc_path),
    }
    mappings_path = os.path.join(out_dir, "mappings.json")
    with open(mappings_path, "w", encoding="utf-8") as f:
        json.dump(mappings, f, ensure_ascii=False, indent=2)

    info = {
        "memmap_rc_path": rc_path,
        "memmap_fc_path": fc_path,
        "metadata_csv": meta_csv,
        "metadata_train_csv": train_csv,
        "metadata_valid_csv": valid_csv,
        "metadata_test_csv": test_csv,
        "mappings_json": mappings_path,
        "n_samples": n,
        "n_features_total": total_feat,
        "n_features_rc": X_rc_prepared.shape[1],
        "n_features_fc": X_fc_prepared.shape[1],
    }
    return info, mappings


def load_memmap_info(out_dir: str) -> Dict:
    """Return dict with memmap paths, metadata and mapping paths + inferred n_features"""
    rc = os.path.join(out_dir, "features_rc.dat")
    fc = os.path.join(out_dir, "features_fc.dat")
    meta_csv = os.path.join(out_dir, "metadata.csv")
    mappings_json = os.path.join(out_dir, "mappings.json")
    if not os.path.exists(rc):
        raise FileNotFoundError(rc)
    if not os.path.exists(fc):
        raise FileNotFoundError(fc)
    if not os.path.exists(meta_csv):
        raise FileNotFoundError(meta_csv)
    meta = pd.read_csv(meta_csv)
    n = meta.shape[0]
    # import os as _os
    # filesize_rc = _os.path.getsize(rc)
    # filesize_fc = _os.path.getsize(fc)
    filesize_rc = os.path.getsize(rc)
    filesize_fc = os.path.getsize(fc)
    dtype = np.float32
    feat_rc = filesize_rc // (np.dtype(dtype).itemsize * n)
    feat_fc = filesize_fc // (np.dtype(dtype).itemsize * n)
    return {
        "memmap_rc_path": rc,
        "memmap_fc_path": fc,
        "metadata_csv": meta_csv,
        "metadata_train_csv": os.path.join(out_dir, "metadata.train.csv") if os.path.exists(os.path.join(out_dir, "metadata.train.csv")) else None,
        "metadata_valid_csv": os.path.join(out_dir, "metadata.valid.csv") if os.path.exists(os.path.join(out_dir, "metadata.valid.csv")) else None,
        "metadata_test_csv": os.path.join(out_dir, "metadata.test.csv") if os.path.exists(os.path.join(out_dir, "metadata.test.csv")) else None,
        "mappings_json": mappings_json,
        "n_samples": n,
        "n_features_rc": int(feat_rc),
        "n_features_fc": int(feat_fc),
    }
