"""Inference: FASTA -> FCGR -> variant."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.nn.functional import softmax

from .cgr import IMAGE_SIZE, ambiguous_fraction, sequence_to_fcgr
from .fasta import records_from_input
from .labels import (
    MAX_AMBIGUOUS_FRACTION,
    MIN_GENOME_LENGTH,
    VARIANT_LABELS,
    VARIANT_WHO_STATUS,
)
from .model import VariantCNN

DEFAULT_WEIGHTS = Path(__file__).resolve().parent.parent / "models" / "variant_cnn.pt"


@dataclass
class Prediction:
    header: str
    variant: str
    confidence: float
    probabilities: dict[str, float]
    length: int
    ambiguous_fraction: float
    fcgr: np.ndarray
    warning: str | None = None

    @property
    def status(self) -> str:
        return VARIANT_WHO_STATUS[self.variant]


class VariantClassifier:
    def __init__(
        self,
        weights_path: str | Path | None = None,
        device: str | None = None,
        min_confidence: float = 0.70,
    ) -> None:
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.min_confidence = min_confidence
        path = Path(weights_path) if weights_path else DEFAULT_WEIGHTS
        if not path.is_file():
            raise FileNotFoundError(
                f"Weights not found at {path}. Train with `python -m cgr_sar2 train --from-original` "
                "or pass --weights."
            )
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        state = checkpoint["model"] if isinstance(checkpoint, dict) and "model" in checkpoint else checkpoint
        self.labels = list(
            checkpoint.get("labels", VARIANT_LABELS) if isinstance(checkpoint, dict) else VARIANT_LABELS
        )
        self.model = VariantCNN(num_classes=len(self.labels))
        self.model.load_state_dict(state)
        self.model.to(self.device)
        self.model.eval()
        self.pixel_mean = None
        self.pixel_std = None
        if isinstance(checkpoint, dict) and "pixel_mean" in checkpoint:
            self.pixel_mean = np.asarray(checkpoint["pixel_mean"], dtype=np.float32)
            self.pixel_std = np.asarray(checkpoint["pixel_std"], dtype=np.float32)

    def predict_sequence(self, sequence: str, header: str = "") -> Prediction:
        seq = "".join(sequence.split())
        frac_n = ambiguous_fraction(seq)
        warning = None
        if len(seq) < MIN_GENOME_LENGTH:
            warning = (
                f"Short sequence ({len(seq)} nt); the paper uses complete genomes "
                f"(> {MIN_GENOME_LENGTH} bp)."
            )
        elif frac_n > MAX_AMBIGUOUS_FRACTION:
            warning = (
                f"Ambiguous bases ({frac_n:.1%}) exceed the paper threshold "
                f"({MAX_AMBIGUOUS_FRACTION:.0%})."
            )
        fcgr = sequence_to_fcgr(seq)
        image = fcgr
        if self.pixel_mean is not None and self.pixel_std is not None:
            image = (fcgr - self.pixel_mean.squeeze()) / self.pixel_std.squeeze()
        tensor = torch.from_numpy(image.astype(np.float32)).unsqueeze(0).unsqueeze(0).to(self.device)
        with torch.no_grad():
            logits = self.model(tensor)
            probs = softmax(logits, dim=1).squeeze(0).cpu().numpy()
        index = int(np.argmax(probs))
        confidence = float(probs[index])
        variant = self.labels[index]
        if confidence < self.min_confidence:
            extra = (
                f"Low confidence ({confidence:.1%}). The sequence may not belong "
                "to the six modelled variants (Alpha, Beta, Gamma, Delta, Iota, Epsilon)."
            )
            warning = f"{warning} {extra}".strip() if warning else extra
        return Prediction(
            header=header,
            variant=variant,
            confidence=confidence,
            probabilities={name: float(p) for name, p in zip(self.labels, probs)},
            length=len(seq),
            ambiguous_fraction=frac_n,
            fcgr=fcgr,
            warning=warning,
        )

    def predict_fasta(
        self,
        fasta_text: str | None = None,
        fasta_path: str | Path | None = None,
    ) -> list[Prediction]:
        records = records_from_input(fasta_text, fasta_path)
        if not records:
            raise ValueError("No FASTA sequence was provided.")
        return [self.predict_sequence(seq, header=header) for header, seq in records]


def fcgr_preview_rgb(fcgr: np.ndarray, size: int = IMAGE_SIZE) -> np.ndarray:
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("inferno")
    colored = cmap(np.clip(fcgr, 0.0, 1.0))[:, :, :3]
    return (colored * 255).astype(np.uint8)
