# pcmaster_anno/models_small.py
import torch.nn as nn
import torch

class PCMA_MYGO_FC(nn.Module):
    def __init__(self, feature_len, hidden_dims=(1024,256), dropout=0.2, num_classes=2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(dropout),
            nn.Linear(feature_len, hidden_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dims[1], num_classes)
        )
    def forward(self, x):
        if x.dim() == 3:
            x = x.view(x.size(0), -1)
        return self.net(x)


class BasicBlock(nn.Module):
    def __init__(self, in_ch, out_ch, stride=1):
        super().__init__()
        self.conv1 = nn.Conv2d(in_ch, out_ch, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_ch)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(out_ch, out_ch, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_ch)
        self.downsample = None
        if stride != 1 or in_ch != out_ch:
            self.downsample = nn.Sequential(
                nn.Conv2d(in_ch, out_ch, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_ch)
            )
    def forward(self, x):
        identity = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.conv2(out)
        out = self.bn2(out)
        if self.downsample is not None:
            identity = self.downsample(identity)
        out += identity
        out = self.relu(out)
        return out


class ResNetSmall(nn.Module):
    def __init__(self, layers, in_ch=3, num_classes=2, base_channels=16, dropout: float = 0.0):
        """
        layers: list like [1] or [1,1] etc.
        dropout: applied before final fc (dropout on features)
        """
        super().__init__()
        self.in_ch = base_channels
        self.conv1 = nn.Conv2d(in_ch, base_channels, kernel_size=7, stride=2, padding=3, bias=False)
        self.bn1 = nn.BatchNorm2d(base_channels)
        self.relu = nn.ReLU(inplace=True)
        self.pool = nn.AdaptiveAvgPool2d((1,1))
        self.layer1 = self._make_layer(base_channels, base_channels, layers[0], stride=1)
        idx = 1
        ch = base_channels
        for ll in layers[1:]:
            setattr(self, f"layer{idx+1}", self._make_layer(ch, ch*2, ll, stride=2))
            ch *= 2
            idx += 1
        final_ch = ch
        self.dropout = nn.Dropout(dropout) if dropout and dropout > 0.0 else nn.Identity()
        self.fc = nn.Linear(final_ch, num_classes)

    def _make_layer(self, in_ch, out_ch, blocks, stride=1):
        layers = []
        layers.append(BasicBlock(in_ch, out_ch, stride))
        for _ in range(1, blocks):
            layers.append(BasicBlock(out_ch, out_ch, stride=1))
        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.layer1(x)
        idx = 2
        while hasattr(self, f"layer{idx}"):
            x = getattr(self, f"layer{idx}")(x)
            idx += 1
        x = self.pool(x)
        x = torch.flatten(x, 1)
        x = self.dropout(x)
        x = self.fc(x)
        return x


def ResNet3(num_classes=2):
    return ResNetSmall(layers=[1], num_classes=num_classes, base_channels=16)

def ResNet5(num_classes=2):
    return ResNetSmall(layers=[1,1], num_classes=num_classes, base_channels=16)



class ResidualBlock1D(nn.Module):
    def __init__(self, in_ch, out_ch, kernel_size=7, stride=1):
        super().__init__()
        padding = kernel_size // 2
        self.conv1 = nn.Conv1d(in_ch, out_ch, kernel_size, stride, padding, bias=False)
        self.bn1 = nn.BatchNorm1d(out_ch)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv1d(out_ch, out_ch, kernel_size, stride=1, padding=padding, bias=False)
        self.bn2 = nn.BatchNorm1d(out_ch)
        self.downsample = None
        if in_ch != out_ch:
            self.downsample = nn.Sequential(
                nn.Conv1d(in_ch, out_ch, kernel_size=1, bias=False),
                nn.BatchNorm1d(out_ch)
            )
    def forward(self, x):
        identity = x
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        if self.downsample is not None:
            identity = self.downsample(identity)
        out += identity
        return self.relu(out)

class PCMA_MYGO_RC_1D(nn.Module):
    def __init__(self, num_classes=2, input_len=10000, base_channels=32, blocks=3, dropout=0.3):
        super().__init__()
        self.proj = nn.Conv1d(1, base_channels, kernel_size=7, stride=2, padding=3, bias=False)
        self.bn = nn.BatchNorm1d(base_channels)
        self.relu = nn.ReLU(inplace=True)

        layers = []
        ch = base_channels
        for _ in range(blocks):
            layers.append(ResidualBlock1D(ch, ch))
        self.res_blocks = nn.Sequential(*layers)

        self.pool = nn.AdaptiveAvgPool1d(1)
        self.drop = nn.Dropout(dropout)
        self.fc = nn.Linear(ch, num_classes)

    def forward(self, x):
        if x.dim() == 2:  # (B, L)
            x = x.unsqueeze(1)
        x = self.relu(self.bn(self.proj(x)))
        x = self.res_blocks(x)
        x = self.pool(x).squeeze(-1)
        x = self.drop(x)
        return self.fc(x)



class PCMA_MYGO_RC(nn.Module):
    def __init__(self, num_classes=2, variant="resnet5", dropout: float = 0.0):
        super().__init__()
        if variant == "resnet3":
            self.net = ResNet3(num_classes=num_classes, )
            # replace with same dropout if desired:
            if dropout and dropout > 0:
                self.net.dropout = nn.Dropout(dropout)
        else:
            # default resnet5
            self.net = ResNet5(num_classes=num_classes)
            if dropout and dropout > 0:
                self.net.dropout = nn.Dropout(dropout)
    def forward(self, x):
        return self.net(x)


class PCMA_MYGO_SA(nn.Module):
    def __init__(self, seq_len, d_model=256, nhead=4, num_classes=2):
        super().__init__()
        self.proj = nn.Linear(seq_len, d_model)
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, batch_first=True)
        self.trans = nn.TransformerEncoder(encoder_layer, num_layers=1)
        self.head = nn.Sequential(nn.LayerNorm(d_model), nn.Linear(d_model, num_classes))
    def forward(self, x):
        if x.dim() == 3:
            x = x.view(x.size(0), -1)
        x = self.proj(x)
        x = x.unsqueeze(1)
        x = self.trans(x)
        x = x.squeeze(1)
        return self.head(x)
    
    