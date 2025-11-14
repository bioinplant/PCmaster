# pcmaster_anno/predict_v2.py
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
import joblib

from .memmap_dataset import MemmapDataset
from torch.utils.data import DataLoader


def predict_torch_model(model_or_path, memmap_path, metadata_csv, mode="rc", image_size=None, batch_size=128, device="cpu"):
    from pcmaster_anno.models_small import PCMA_MYGO_RC, PCMA_MYGO_FC, PCMA_MYGO_SA

    # If a string path is passed, auto-detect which model to build
    if isinstance(model_or_path, str):
        path = model_or_path
        # crude detection: RC uses square image size, FC/SA uses vector
        if mode == "rc":
            # need num_classes from metadata
            meta = pd.read_csv(metadata_csv)
            n_classes = meta["celltype_int"].nunique()
            model = PCMA_MYGO_RC(num_classes=n_classes, variant='resnet3')
        elif mode == "fc":
            meta = pd.read_csv(metadata_csv)
            n_classes = meta["celltype_int"].nunique()
            model = PCMA_MYGO_FC(feature_len=None, num_classes=n_classes)
        elif mode == "sa":
            meta = pd.read_csv(metadata_csv)
            n_classes = meta["celltype_int"].nunique()
            model = PCMA_MYGO_SA(feature_len=None, num_classes=n_classes)
        else:
            raise ValueError("Unsupported mode: " + mode)
        import torch
        model.load_state_dict(torch.load(path, map_location="cpu"))
    else:
        model = model_or_path

    ds = MemmapDataset(memmap_path, metadata_csv, mode=mode, image_size=image_size)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False, num_workers=0)
    model = model.to(device)
    model.eval()
    all_probs = []
    all_cellnames = []
    all_true = []
    with torch.no_grad():
        for xb, yb, cellnames in loader:
            xb = xb.to(device)
            out = model(xb)
            probs = F.softmax(out, dim=1).cpu().numpy()
            all_probs.append(probs)
            all_cellnames.extend([str(x) for x in cellnames])
            all_true.extend([int(x) for x in yb.numpy()])
    probs_arr = np.vstack(all_probs)
    n_classes = probs_arr.shape[1]
    prob_cols = [f"prob_{i}" for i in range(n_classes)]
    out_df = pd.DataFrame(probs_arr, columns=prob_cols)
    out_df.insert(0, "orig_celltype", all_true)
    out_df.insert(0, "cellname", all_cellnames)
    out_df["pred_label"] = out_df[prob_cols].values.argmax(axis=1).astype(int)
    return out_df


def predict_sgd_svm(model_path, memmap_path, metadata_csv, scaler_path=None, batch_size=2048):
    import os
    clf = joblib.load(model_path)
    sub_meta = pd.read_csv(metadata_csv)
    parent = os.path.dirname(metadata_csv)
    global_meta = pd.read_csv(os.path.join(parent, "metadata.csv"))
    total_n = len(global_meta)
    filesize = os.path.getsize(memmap_path)
    feat = filesize // (np.dtype(np.float32).itemsize * total_n)
    mm = np.memmap(memmap_path, dtype=np.float32, mode='r', shape=(total_n, feat))

    idx_map = {str(r): i for i, r in enumerate(global_meta['cellname'].tolist())}
    indices = [idx_map.get(str(cn)) for cn in sub_meta['cellname'].tolist()]
    indices = [i for i in indices if i is not None]

    scaler = None
    if scaler_path and os.path.exists(scaler_path):
        scaler = joblib.load(scaler_path)

    results = []
    # iterate over indices in batches
    for start in range(0, len(indices), batch_size):
        batch_idx = indices[start:start+batch_size]
        Xb = mm[batch_idx].astype(np.float32)
        if scaler is not None:
            Xb = scaler.transform(Xb)
        # compute probabilities/decision_function
        if hasattr(clf, 'predict_proba'):
            probs = clf.predict_proba(Xb)
        else:
            dec = clf.decision_function(Xb)
            if dec.ndim == 1:
                p1 = 1.0 / (1.0 + np.exp(-dec))
                probs = np.vstack([1-p1, p1]).T
            else:
                ex = np.exp(dec - dec.max(axis=1, keepdims=True))
                probs = ex / ex.sum(axis=1, keepdims=True)
        # map back to cellnames and true labels (from subset)
        for i_local, gi in enumerate(batch_idx):
            cname = sub_meta['cellname'].iloc[start + i_local]
            true = int(sub_meta['celltype_int'].iloc[start + i_local])
            results.append((cname, true, probs[i_local]))
    n_classes = results[0][2].shape[0]
    cols = [f"prob_{i}" for i in range(n_classes)]
    rows = []
    for cname, t, p in results:
        rows.append([cname, t] + p.tolist())
    df = pd.DataFrame(rows, columns=["cellname", "orig_celltype"] + cols)
    df['pred_label'] = df[[f"prob_{i}" for i in range(n_classes)]].values.argmax(axis=1).astype(int)
    return df



