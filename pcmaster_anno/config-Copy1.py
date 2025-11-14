# pcmaster_anno/config.py
from dataclasses import dataclass
from typing import Tuple

@dataclass
class PCMAConfig:
    # model defaults
    rc_variant: str = "resnet5"          # "resnet3" or "resnet5"
    rc_dropout: float = 0.3

    fc_hidden_dims: Tuple[int,int] = (1024, 256)
    fc_dropout: float = 0.2

    sa_d_model: int = 256
    sa_nhead: int = 4

    # training defaults
    lr: float = 1e-3
    lr_patience: int = 2
    earlystop_patience: int = 3
    min_lr: float = 1e-6

    batch_size: int = 64
    num_workers: int = 0
