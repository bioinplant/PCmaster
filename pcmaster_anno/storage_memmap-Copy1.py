# pcmaster_anno/storage_memmap.py
import os
import json
from typing import Tuple, Dict, Optional
import numpy as np
import pandas as pd


def save_adata_to_memmap(adata, out_dir: str, ensure_length: Optional[int]=None, dtype=np.float32):
    """
    Convert an AnnData-like object to on-disk memmap storage.

    Produces files under out_dir:
      - features.dat (numpy memmap binary)
      - metadata.csv (columns: cellname, celltype_orig, celltype_int)
      - mappings.json ({"obs2int":..., "int2obs":...})

    Returns: info dict and mapping dict
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

    # X -> dense numpy
    X = adata.X
    try:
        Xmat = np.asarray(X.todense())
    except Exception:
        Xmat = np.asarray(X)

    n, feat = Xmat.shape
    if ensure_length is not None and feat < ensure_length:
        pad = np.zeros((n, ensure_length - feat), dtype=dtype)
        Xmat = np.concatenate([Xmat.astype(dtype), pad], axis=1)
        feat = ensure_length
    else:
        Xmat = Xmat.astype(dtype)

    memmap_path = os.path.join(out_dir, "features.dat")
    fp = np.memmap(memmap_path, dtype=dtype, mode='w+', shape=(n, feat))
    fp[:] = Xmat[:]
    del fp

    meta_df = pd.DataFrame({
        "cellname": list(obs["cellname"].astype(str)),
        "celltype_orig": list(obs["celltype"].astype(str))
    })
    unique = pd.Series(meta_df["celltype_orig"]).fillna("unknown").astype(str).unique().tolist()
    map_obs2int = {lab: idx for idx, lab in enumerate(unique)}
    map_int2obs = {v: k for k, v in map_obs2int.items()}
    meta_df["celltype_int"] = meta_df["celltype_orig"].map(map_obs2int)

    meta_csv = os.path.join(out_dir, "metadata.csv")
    meta_df.to_csv(meta_csv, index=False)

    mappings = {"obs2int": map_obs2int, "int2obs": map_int2obs}
    with open(os.path.join(out_dir, "mappings.json"), "w", encoding="utf-8") as f:
        json.dump(mappings, f, ensure_ascii=False, indent=2)

    info = {
        "memmap_path": memmap_path,
        "metadata_csv": meta_csv,
        "mappings_json": os.path.join(out_dir, "mappings.json"),
        "n_samples": n,
        "n_features": feat
    }
    return info, mappings


def load_memmap_info(out_dir: str) -> Dict:
    """Return dict with memmap_path, metadata and mapping paths + inferred n_features"""
    memmap_path = os.path.join(out_dir, "features.dat")
    meta_csv = os.path.join(out_dir, "metadata.csv")
    mappings_json = os.path.join(out_dir, "mappings.json")
    if not os.path.exists(memmap_path):
        raise FileNotFoundError(memmap_path)
    if not os.path.exists(meta_csv):
        raise FileNotFoundError(meta_csv)
    meta = pd.read_csv(meta_csv)
    n = meta.shape[0]
    import os as _os
    filesize = _os.path.getsize(memmap_path)
    dtype = np.float32
    feat = filesize // (np.dtype(dtype).itemsize * n)
    return {"memmap_path": memmap_path, "metadata_csv": meta_csv, "mappings_json": mappings_json, "n_samples": n, "n_features": int(feat)}

