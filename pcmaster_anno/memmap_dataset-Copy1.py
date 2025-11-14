# pcmaster_anno/memmap_dataset.py
import math
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import torch


class MemmapDataset(Dataset):
    """Dataset reading features from a numpy memmap and metadata CSV.

    mode: 'fc'|'rc'|'sa'|'hp'
    image_size: side length for rc mode (if None inferred from sqrt)
    """
    def __init__(self, memmap_path: str, metadata_csv: str, mode: str = "rc", image_size: int = None, dtype=np.float32):
        self.memmap_path = memmap_path
        self.meta = pd.read_csv(metadata_csv)
        self.n = len(self.meta)
        self.dtype = dtype
        import os
        filesize = os.path.getsize(memmap_path)
        feat = filesize // (np.dtype(self.dtype).itemsize * self.n)
        self.feat = int(feat)
        self._memmap = None
        self.mode = mode
        self.image_size = image_size

    @property
    def memmap(self):
        if self._memmap is None:
            self._memmap = np.memmap(self.memmap_path, dtype=self.dtype, mode='r', shape=(self.n, self.feat))
        return self._memmap

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        row = self.memmap[idx]
        label = int(self.meta.iloc[idx]["celltype_int"])
        cellname = str(self.meta.iloc[idx]["cellname"])
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

            