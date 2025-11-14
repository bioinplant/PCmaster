# pcmaster_anno/memmap_dataset.py
import math
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import torch
import os


class MemmapDataset(Dataset):
    """Dataset reading features from a numpy memmap and a subset metadata CSV.

    memmap_path: path to memmap file (must store n_rows x feat_cols where n_rows == len(metadata.csv))
    metadata_csv: path to the subset metadata file you want to iterate (e.g., metadata.train.csv).
      The parent folder must contain the global metadata.csv (all rows).
    mode: 'fc'|'rc'|'sa'|'hp'
    image_size: side length for rc mode (if None inferred from sqrt)
    """
    def __init__(self, memmap_path: str, metadata_csv: str, mode: str = "rc", image_size: int = None, dtype=np.float32):
        self.memmap_path = memmap_path
        self.sub_meta = pd.read_csv(metadata_csv)
        self.dtype = dtype
        self.mode = mode
        self.image_size = image_size

        # load global metadata (must be in same folder)
        parent = os.path.dirname(metadata_csv)
        global_meta_path = os.path.join(parent, "metadata.csv")
        if not os.path.exists(global_meta_path):
            raise FileNotFoundError(f"Global metadata.csv not found in {parent}; required to map subset indices.")
        self.global_meta = pd.read_csv(global_meta_path)
        self.total_n = len(self.global_meta)

        # compute feature dim from file size
        filesize = os.path.getsize(memmap_path)
        feat = filesize // (np.dtype(self.dtype).itemsize * self.total_n)
        self.feat = int(feat)

        # build mapping from cellname -> global index
        name_to_idx = {str(c): i for i, c in enumerate(self.global_meta["cellname"].tolist())}
        self.indices = []
        for cn in self.sub_meta["cellname"].tolist():
            idx = name_to_idx.get(str(cn))
            if idx is None:
                raise ValueError(f"Cellname {cn} from subset metadata not found in global metadata.")
            self.indices.append(int(idx))

        self._memmap = None

    @property
    def memmap(self):
        # open memmap with global shape (total_n, feat); we will index into it using global indices
        if self._memmap is None:
            self._memmap = np.memmap(self.memmap_path, dtype=self.dtype, mode="r", shape=(self.total_n, self.feat))
        return self._memmap

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        global_idx = self.indices[idx]
        row = self.memmap[global_idx]
        # label read from subset meta (already mapped to ints by save function)
        label = int(self.sub_meta.iloc[idx]["celltype_int"])
        cellname = str(self.sub_meta.iloc[idx]["cellname"])
        if self.mode in ("fc", "sa", "hp"):
            x = torch.from_numpy(row.astype("float32"))
            return x, label, cellname
        elif self.mode == "rc":
            if self.image_size is None:
                side = int(math.sqrt(self.feat))
            else:
                side = self.image_size
            expected = side * side
            if expected != self.feat:
                if expected > self.feat:
                    pad = np.zeros(expected - self.feat, dtype=self.dtype)
                    arr = np.concatenate([row, pad], axis=0)
                else:
                    arr = row[:expected]
            else:
                arr = row
            img = arr.reshape(side, side)
            img = np.expand_dims(img, 0).astype("float32")
            img = np.repeat(img, 3, axis=0)
            return torch.from_numpy(img), label, cellname
        else:
            raise ValueError("Unknown mode: " + self.mode)
