# pcmaster_anno/predict_v2.py
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
import joblib

from .memmap_dataset import MemmapDataset
from torch.utils.data import DataLoader


def predict_torch(model, memmap_path, metadata_csv, mode="rc", image_size=None, batch_size=128, device="cpu"):
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
    clf = joblib.load(model_path)
    meta = pd.read_csv(metadata_csv)
    n = len(meta)
    mm = np.memmap(memmap_path, dtype=np.float32, mode='r', shape=(n, -1))
    scaler = None
    if scaler_path and os.path.exists(scaler_path):
        scaler = joblib.load(scaler_path)
    results = []
    for start in range(0, n, batch_size):
        end = min(start+batch_size, n)
        Xb = mm[start:end].astype(np.float32)
        if scaler is not None:
            Xb = scaler.transform(Xb)
        cellnames = meta['cellname'].iloc[start:end].tolist()
        true = meta['celltype_int'].iloc[start:end].tolist()
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
        for i, cname in enumerate(cellnames):
            results.append((cname, true[i], probs[i]))
    n_classes = results[0][2].shape[0]
    cols = [f"prob_{i}" for i in range(n_classes)]
    rows = []
    for cname, t, p in results:
        rows.append([cname, t] + p.tolist())
    df = pd.DataFrame(rows, columns=["cellname", "orig_celltype"] + cols)
    df['pred_label'] = df[[f"prob_{i}" for i in range(n_classes)]].values.argmax(axis=1).astype(int)
    return df


