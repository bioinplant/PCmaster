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
from .config import PCMAConfig
from sklearn.base import clone
from .config import PCMAConfig

def _evaluate_model(model, loader, device):
    model.eval()
    total_loss = 0.0
    correct = 0
    total = 0
    with torch.no_grad():
        for xb, yb, _ in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            out = model(xb)
            loss = F.cross_entropy(out, yb, reduction='sum')
            total_loss += loss.item()
            preds = out.argmax(dim=1)
            correct += (preds == yb).sum().item()
            total += xb.size(0)
    avg_loss = total_loss / max(1, total)
    acc = correct / max(1, total)
    return avg_loss, acc

def train_torch_model(model: torch.nn.Module,
                      memmap_path: str,
                      metadata_csv: str,
                      mode: str,
                      image_size: int = None,
                      device: str = "cpu",
                      epochs: int = 10,
                      batch_size: int = 64,
                      lr: float = 1e-3,
                      out_model_path: str = "model.pth",
                      cfg: PCMAConfig = None,
                      valid_metadata_csv: str = None):
    """
    Train a torch model with optional validation, LR scheduler and early stopping.
    - cfg: PCMAConfig object (if None default values used)
    - valid_metadata_csv: path to validation subset CSV; if provided we'll evaluate the model on it after each epoch
    """
    if cfg is None:
        cfg = PCMAConfig()

    device = device if device is not None else ("cuda" if torch.cuda.is_available() else "cpu")
    ds = MemmapDataset(memmap_path, metadata_csv, mode=mode, image_size=image_size)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=True, num_workers=cfg.num_workers, pin_memory=(device != "cpu" and torch.cuda.is_available()))
    model = model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)

    # Use OneCycleLR for smoother convergence
    scheduler = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=lr, steps_per_epoch=len(loader), epochs=epochs
    )

    best_val_loss = float('inf')
    best_val_acc = 0.0
    no_improve = 0

    # optional validation loader
    val_loader = None
    if valid_metadata_csv is not None:
        val_ds = MemmapDataset(memmap_path, valid_metadata_csv, mode=mode, image_size=image_size)
        val_loader = DataLoader(
            val_ds,
            batch_size=batch_size,
            shuffle=False,
            num_workers=cfg.num_workers,
            pin_memory=(device != "cpu" and torch.cuda.is_available())
        )

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
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            opt.step()
            scheduler.step()
            total_loss += loss.item() * xb.size(0)
        train_loss = total_loss / max(1, len(ds))

        if val_loader is not None:
            val_loss, val_acc = _evaluate_model(model, val_loader, device)
            print(f"Epoch {epoch+1}/{epochs} train_loss={train_loss:.6f} val_loss={val_loss:.6f} val_acc={val_acc:.4f}")
            if val_acc > best_val_acc + 1e-2:
                best_val_acc = val_acc
                no_improve = 0
                save_torch_state(model, out_model_path)
                print(f"  -> New best model saved to {out_model_path}")
            else:
                no_improve += 1
                if no_improve >= cfg.earlystop_patience:
                    print(f"Early stopping: no improvement for {no_improve} epochs")
                    break
        else:
            print(f"Epoch {epoch+1}/{epochs} train_loss={train_loss:.6f}")
            save_torch_state(model, out_model_path)

    if val_loader is None:
        save_torch_state(model, out_model_path)
        print("Saved final model to", out_model_path)
    return out_model_path


def train_sgd_svm_with_scaler_and_calibration(
        memmap_path: str,
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
        valid_metadata_csv: str = None,
        cfg: PCMAConfig = None):
    """
    Train an SGD SVM on memmap data, with optional StandardScaler and probability calibration.
    """
    if cfg is None:
        cfg = PCMAConfig()

    meta_train = pd.read_csv(metadata_csv)
    parent = os.path.dirname(metadata_csv)
    global_meta = pd.read_csv(os.path.join(parent, "metadata.csv"))
    total_n = len(global_meta)

    filesize = os.path.getsize(memmap_path)
    feat = filesize // (np.dtype(np.float32).itemsize * total_n)
    mm = np.memmap(memmap_path, dtype=np.float32, mode='r', shape=(total_n, feat))

    idx_map = {str(r): i for i, r in enumerate(global_meta['cellname'].tolist())}
    train_indices = [idx_map.get(str(cn)) for cn in meta_train['cellname'].tolist()]
    train_indices = [i for i in train_indices if i is not None]

    classes = np.sort(meta_train["celltype_int"].unique())

    # fit scaler incrementally
    scaler = StandardScaler()
    for start in range(0, len(train_indices), batch_size):
        batch_idx = train_indices[start:start+batch_size]
        Xb = mm[batch_idx].astype(np.float32)
        scaler.partial_fit(Xb)
    if scaler_path:
        joblib.dump(scaler, scaler_path)
        print(f"Saved scaler to {scaler_path}")

    clf = SGDClassifier(
        loss='hinge' if use_hinge else 'log',
        penalty='l2',
        alpha=l2_reg,
        learning_rate='optimal',
        eta0=learning_rate,
        max_iter=1,
        warm_start=True
    )

    for epoch in range(epochs):
        perm = np.random.permutation(len(train_indices))
        for start in range(0, len(train_indices), batch_size):
            batch_idx_local = perm[start:start+batch_size]
            batch_idx = [train_indices[i] for i in batch_idx_local]
            Xb = mm[batch_idx].astype(np.float32)
            Xb = scaler.transform(Xb)
            yb = meta_train.iloc[batch_idx_local]["celltype_int"].values.astype(np.int64)
            if epoch == 0 and start == 0:
                clf.partial_fit(Xb, yb, classes=classes)
            else:
                clf.partial_fit(Xb, yb)
        print(f"Epoch {epoch+1}/{epochs} done for SGD SVM")

    joblib.dump(clf, out_model_path)
    print(f"Saved base SGD model to {out_model_path}")

    return out_model_path
