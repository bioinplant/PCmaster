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
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    # scheduler reduces lr on plateau of validation loss
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, mode='min', factor=0.5, patience=cfg.lr_patience, min_lr=cfg.min_lr)

    best_val_loss = float('inf')
    no_improve = 0

    # optional validation loader (constructed once)
    val_loader = None
    if valid_metadata_csv is not None:
        val_ds = MemmapDataset(memmap_path, valid_metadata_csv, mode=mode, image_size=image_size)
        val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, num_workers=cfg.num_workers, pin_memory=(device != "cpu" and torch.cuda.is_available()))

    for epoch in range(epochs):
        model.train()
        total_loss = 0.0
        for xb, yb, _ in loader:
            # if you want feature scaling for FC/SA, you can add scaler logic here (apply before .to(device))
            xb = xb.to(device)
            yb = yb.to(device)
            out = model(xb)
            loss = F.cross_entropy(out, yb)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total_loss += loss.item() * xb.size(0)
        train_loss = total_loss / max(1, len(ds))
        # validation
        if val_loader is not None:
            val_loss, val_acc = _evaluate_model(model, val_loader, device)
            print(f"Epoch {epoch+1}/{epochs} train_loss={train_loss:.6f} val_loss={val_loss:.6f} val_acc={val_acc:.4f}")
            scheduler.step(val_loss)
            # early stopping / save best
            if val_loss < best_val_loss - 1e-8:
                best_val_loss = val_loss
                no_improve = 0
                # save current best
                save_torch_state(model, out_model_path)
                print(f"  -> New best model saved to {out_model_path}")
            else:
                no_improve += 1
                if no_improve >= cfg.earlystop_patience:
                    print(f"Early stopping: no improvement for {no_improve} epochs (patience={cfg.earlystop_patience})")
                    break
        else:
            # no validation provided, just report train loss and save final model
            print(f"Epoch {epoch+1}/{epochs} train_loss={train_loss:.6f}")
            save_torch_state(model, out_model_path)
    # If validation was used we already saved best model; if not, ensure final saved
    if val_loader is None:
        save_torch_state(model, out_model_path)
        print("Saved to", out_model_path)
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

    # calibration
    if calibrate_on_valid and valid_metadata_csv and calibrated_out_path:
        valid_meta = pd.read_csv(valid_metadata_csv)
        valid_idxs = [idx_map.get(str(cn)) for cn in valid_meta['cellname'].tolist()]
        valid_idxs = [i for i in valid_idxs if i is not None]
        if len(valid_idxs) == 0:
            print("⚠️ No valid indices found for calibration. Skipping calibration.")
            return out_model_path
        Xv = mm[valid_idxs].astype(np.float32)
        Xv = scaler.transform(Xv)
        yv = valid_meta['celltype_int'].values.astype(int)

        # clone to avoid corrupting clf
        try:
            calibrator = CalibratedClassifierCV(base_estimator=clone(clf), method='sigmoid', cv='prefit')
        except TypeError:
            calibrator = CalibratedClassifierCV(estimator=clone(clf), method='sigmoid', cv='prefit')

        calibrator.fit(Xv, yv)
        joblib.dump(calibrator, calibrated_out_path)
        print(f"Saved calibrated classifier to {calibrated_out_path}")
        return calibrated_out_path

    return out_model_path
# train_sgd_svm_with_scaler_and_calibration error
# HP model (SVM with calibration)
# svm_model_path = os.path.join(out_dir, "model_svm.pkl")
# scaler_path = os.path.join(out_dir, "scaler.pkl")
# svm_calibrated_path = os.path.join(out_dir, "model_svm_calibrated.pkl")

# train_v2.train_sgd_svm_with_scaler_and_calibration(
#     memmap_fc,
#     train_csv,
#     out_model_path=svm_model_path,
#     scaler_path=scaler_path,
#     calibrated_out_path=svm_calibrated_path,
#     valid_metadata_csv=valid_csv,
#     epochs=3
# )

# df_svm = predict_v2.predict_sgd_svm(svm_calibrated_path, memmap_fc, valid_csv, scaler_path=scaler_path)
# print("SVM results head:\n", df_svm.head())
# Saved scaler to ./my_dataset_root_bak1/scaler.pkl
# Epoch 1/3 done for SGD SVM
# Epoch 2/3 done for SGD SVM
# Epoch 3/3 done for SGD SVM
# Saved base SGD model to ./my_dataset_root_bak1/model_svm.pkl
# D:\aaa-new\setups\conda\conda-files\envs\PCmaster_anno_win\lib\site-packages\sklearn\calibration.py:321: FutureWarning: `base_estimator` was renamed to `estimator` in version 1.2 and will be removed in 1.4.
#   warnings.warn(
# ---------------------------------------------------------------------------
# TypeError                                 Traceback (most recent call last)
# Cell In[12], line 6
#       3 scaler_path = os.path.join(out_dir, "scaler.pkl")
#       4 svm_calibrated_path = os.path.join(out_dir, "model_svm_calibrated.pkl")
# ----> 6 train_v2.train_sgd_svm_with_scaler_and_calibration(
#       7     memmap_fc,
#       8     train_csv,
#       9     out_model_path=svm_model_path,
#      10     scaler_path=scaler_path,
#      11     calibrated_out_path=svm_calibrated_path,
#      12     valid_metadata_csv=valid_csv,
#      13     epochs=3
#      14 )
#      16 df_svm = predict_v2.predict_sgd_svm(svm_calibrated_path, memmap_fc, valid_csv, scaler_path=scaler_path)
#      17 print("SVM results head:\n", df_svm.head())

# File D:\pcmaster_anno\train_v2.py:207, in train_sgd_svm_with_scaler_and_calibration(memmap_path, metadata_csv, out_model_path, scaler_path, calibrated_out_path, batch_size, epochs, learning_rate, l2_reg, use_hinge, calibrate_on_valid, valid_metadata_csv, cfg)
#     204 except TypeError:
#     205     calibrator = CalibratedClassifierCV(estimator=clone(clf), method='sigmoid', cv='prefit')
# --> 207 calibrator.fit(Xv, yv)
#     208 joblib.dump(calibrator, calibrated_out_path)
#     209 print(f"Saved calibrated classifier to {calibrated_out_path}")

# File D:\aaa-new\setups\conda\conda-files\envs\PCmaster_anno_win\lib\site-packages\sklearn\calibration.py:338, in CalibratedClassifierCV.fit(self, X, y, sample_weight, **fit_params)
#     335 self.calibrated_classifiers_ = []
#     336 if self.cv == "prefit":
#     337     # `classes_` should be consistent with that of estimator
# --> 338     check_is_fitted(self.estimator, attributes=["classes_"])
#     339     self.classes_ = self.estimator.classes_
#     341     pred_method, method_name = _get_prediction_method(estimator)

# File D:\aaa-new\setups\conda\conda-files\envs\PCmaster_anno_win\lib\site-packages\sklearn\utils\validation.py:1376, in check_is_fitted(estimator, attributes, msg, all_or_any)
#    1370     msg = (
#    1371         "This %(name)s instance is not fitted yet. Call 'fit' with "
#    1372         "appropriate arguments before using this estimator."
#    1373     )
#    1375 if not hasattr(estimator, "fit"):
# -> 1376     raise TypeError("%s is not an estimator instance." % (estimator))
#    1378 if attributes is not None:
#    1379     if not isinstance(attributes, (list, tuple)):

# TypeError: None is not an estimator instance.

# part of ref for origin simple PCMA_MYGO_HP svm
# def train_PCMA_MYGO_HP_with_unknown(
#     train_data=None,train_label=None,valid_data=None,valid_label=None,test_data=None,test_label=None,
#     the_threshold=None,
#     mapping_1=None,mapping_2=None,mapping_1_test=None,mapping_2_test=None,
# ):
#     if the_threshold is None:
#         the_threshold = 0.7
#     Classifier = LinearSVC(C=1.2)
#     clf = CalibratedClassifierCV(Classifier)
#     start_time = time.time()
#     clf.fit(train_data, train_label)
    
#     del train_data,train_label,valid_data,valid_label
#     gc.collect()
#     print_memory_used()
#     # print(f'num_classes = {num_classes}')
    
#     model = clf
#     end_time = time.time()
#     execution_time = end_time - start_time
#     print(f"time cost：{execution_time} s")
#     start_time = time.time()
    
#     predicted = clf.predict(test_data)
#     prob = np.max(clf.predict_proba(test_data), axis = 1)