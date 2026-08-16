"""CNN matching the MATLAB Experiment 1 network (Camara et al., IEEE Access, 2026).

The published text cites a 3x3 third convolution (354,278 parameters). The
original MATLAB scripts (main.m / exp1_main_script.m) use 5x5 for that layer,
with no ReLU on the fully connected layers. This implementation follows the
scripts that generated the GISAID results.
"""

from __future__ import annotations

import torch
from torch import nn

from .labels import VARIANT_LABELS


class VariantCNN(nn.Module):
    """Conv 7x7 / 5x5 / 5x5, pool 2x2, FC 64 and 32, dropout 40%."""

    def __init__(self, num_classes: int = len(VARIANT_LABELS)) -> None:
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=7, stride=1, padding=3),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(kernel_size=2, stride=2),
            nn.Conv2d(32, 64, kernel_size=5, stride=1, padding=2),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(kernel_size=2, stride=2),
            nn.Conv2d(64, 64, kernel_size=5, stride=1, padding=2),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(kernel_size=2, stride=2),
        )
        self.classifier = nn.Sequential(
            nn.Flatten(),
            nn.Linear(64 * 8 * 8, 64),
            nn.Dropout(0.4),
            nn.Linear(64, 32),
            nn.Dropout(0.4),
            nn.Linear(32, num_classes),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.classifier(self.features(x))


def count_parameters(model: nn.Module | None = None) -> int:
    if model is None:
        model = VariantCNN()
    return sum(p.numel() for p in model.parameters() if p.requires_grad)
