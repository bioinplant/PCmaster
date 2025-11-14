# pcmaster_anno/train_v2.py
import os
import torch
from torch.utils.data import DataLoader
import torch.nn.functional as F
from sklearn.linear_model import SGDClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import StandardScaler
import numpy as np
import joblib
import pandas as pd

from .memmap_dataset import MemmapDataset
from .models_small import PCMA_MYGO_FC, PCMA_MYGO_RC, PCMA_MYGO_SA
from .utils import save_json, load_json, save_torch_state


def train_torch_model(model: torch.nn.Module,
                      memmap_path: str,
                      metadata_csv: str,
                      mode: str,
                      image_size: int = None,
                      device: str = "cpu",
                      epochs: int = 10,
                      batch_size: int = 64,
                      lr: float = 1e-3,
                      out_model_path: str = "model.pth"):
    ds = MemmapDataset(memmap_path, metadata_csv, mode=mode, image_size=image_size)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=True, num_workers=0)
    model = model.to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    for epoch in range(epochs):
        model.train()
        total_loss = 0.0
        for xb, yb, _ in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            out = model(xb)
            loss = F.cross_entropy(out, yb)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total_loss += loss.item() * xb.size(0)
        avg = total_loss / len(ds)
        print(f"Epoch {epoch+1}/{epochs} loss={avg:.6f}")
    save_torch_state(model, out_model_path)
    print("Saved to", out_model_path)
    return out_model_path


def train_sgd_svm_with_scaler_and_calibration(memmap_path: str,
                                              metadata_csv: str,
                                              out_model_path: str,
                                              scaler_path: str = None,
                                              calibrated_out_path: str = None,
                                              batch_size: int = 2048,
                                              epochs: int = 3,
                                              learning_rate: float = 0.01,
                                              l2_reg: float = 1e-4,
                                              use_hinge: bool = True,
                                              calibrate_on_valid: bool = True,
                                              valid_metadata_csv: str = None):
    """
    Train SGDClassifier (linear SVM-like) with an online StandardScaler and optional calibration on validation set.

    - Saves scaler (joblib) and final model (joblib). If calibrate_on_valid=True and valid_metadata_csv provided,
      performs CalibratedClassifierCV(cv='prefit', method='sigmoid') on validation set.
    """
    meta = pd.read_csv(metadata_csv)
    classes = sorted(meta["celltype_int"].unique().tolist())
    n = len(meta)
    mm = np.memmap(memmap_path, dtype=np.float32, mode='r', shape=(n, -1))

    scaler = StandardScaler()
    # partial fit scaler on training memmap (in chunks)
    for start in range(0, n, batch_size):
        end = min(start+batch_size, n)
        Xb = mm[start:end].astype(np.float32)
        scaler.partial_fit(Xb)
    if scaler_path:
        joblib.dump(scaler, scaler_path)
        print("Saved scaler to", scaler_path)

    clf = SGDClassifier(loss='hinge' if use_hinge else 'log', penalty='l2', alpha=l2_reg, learning_rate='optimal', eta0=learning_rate)
    # partial_fit
    for epoch in range(epochs):
        idxs = np.random.permutation(n)
        for start in range(0, n, batch_size):
            batch_idx = idxs[start:start+batch_size]
            Xb = mm[batch_idx].astype(np.float32)
            Xb = scaler.transform(Xb)
            yb = meta.iloc[batch_idx]["celltype_int"].values.astype(np.int64)
            if epoch == 0 and start == 0:
                clf.partial_fit(Xb, yb, classes=classes)
            else:
                clf.partial_fit(Xb, yb)
        print(f"Epoch {epoch+1}/{epochs} done for SGD SVM")
    joblib.dump(clf, out_model_path)
    print("Saved base SGD model to", out_model_path)

    # calibration on validation set if requested
    if calibrate_on_valid and valid_metadata_csv is not None and calibrated_out_path is not None:
        vmeta = pd.read_csv(valid_metadata_csv)
        vn = len(vmeta)
        vmm = np.memmap(memmap_path, dtype=np.float32, mode='r', shape=(n, -1))
        # collect validation data from metadata: find indices of validation rows in meta
        # assume validation rows are contiguous or we find them by joining on cellname
        train_meta = meta.copy()
        valid_meta = vmeta.copy()
        # build index mapping from cellname -> index in main meta
        idx_map = {str(r): i for i, r in enumerate(train_meta['cellname'].tolist())}
        valid_idxs = [idx_map.get(str(cn)) for cn in valid_meta['cellname'].tolist()]
        valid_idxs = [i for i in valid_idxs if i is not None]
        if len(valid_idxs) == 0:
            print("No valid indices found for calibration. Skipping calibration.")
            return out_model_path
        Xv = vmm[valid_idxs].astype(np.float32)
        Xv = scaler.transform(Xv)
        yv = train_meta.iloc[valid_idxs]['celltype_int'].values.astype(int)
        # calibrated classifier
        calibrator = CalibratedClassifierCV(base_estimator=clf, method='sigmoid', cv='prefit')
        calibrator.fit(Xv, yv)
        joblib.dump(calibrator, calibrated_out_path)
        print("Saved calibrated classifier to", calibrated_out_path)
        return calibrated_out_path
    return out_model_path

