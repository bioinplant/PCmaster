# pcmaster_anno/utils.py
import random
import os
import json
import numpy as np
import torch
import torch.backends.cudnn as cudnn

def same_seeds(seed: int = 42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    try:
        torch.cuda.manual_seed_all(seed)
    except Exception:
        pass
    os.environ['PYTHONHASHSEED'] = str(seed)
    try:
        cudnn.deterministic = True
        cudnn.benchmark = False
    except Exception:
        pass


def save_json(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_torch_state(model, path: str):
    torch.save(model.state_dict(), path)


def load_torch_state(model, path: str, map_location='cpu'):
    state = torch.load(path, map_location=map_location)
    if isinstance(state, dict) and any(k.startswith('module.') for k in state.keys()):
        new_state = {k.replace('module.', ''): v for k, v in state.items()}
        model.load_state_dict(new_state)
    else:
        try:
            model.load_state_dict(state)
        except RuntimeError:
            if 'state_dict' in state:
                sd = state['state_dict']
                if any(k.startswith('module.') for k in sd.keys()):
                    new_state = {k.replace('module.', ''): v for k, v in sd.items()}
                    model.load_state_dict(new_state)
                else:
                    model.load_state_dict(sd)
            else:
                raise
    return model

