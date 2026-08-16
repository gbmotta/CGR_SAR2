"""CGR/FCGR matching the original MATLAB CGR.m + histcounts2 pipeline.

Vertices (CGR.m): A(0,0), C(0,1), G(1,0), T(1,1). Start at (0,0).
Ambiguous bases reset the point to (0,0). All trajectory points (length+1)
are binned with a 64 x 64 histogram, rotated 90 deg CCW, then quantized
to 8 bits: floor((255 / max) * counts).
"""

from __future__ import annotations

import numpy as np

IMAGE_SIZE = 64
K_DEFAULT = 6  # 64 bins; kept for API compatibility

_VX = np.array([0.0, 0.0, 1.0, 1.0], dtype=np.float64)  # A C G T
_VY = np.array([0.0, 1.0, 0.0, 1.0], dtype=np.float64)

_LOOKUP = np.full(256, 255, dtype=np.uint8)
for _ch, _code in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
    _LOOKUP[ord(_ch)] = _code
    _LOOKUP[ord(_ch.lower())] = _code


def encode_bases(sequence: str) -> np.ndarray:
    if not sequence:
        return np.empty(0, dtype=np.uint8)
    raw = np.frombuffer(sequence.encode("ascii", errors="ignore"), dtype=np.uint8)
    return _LOOKUP[raw]


def sequence_to_cgr_points(sequence: str) -> tuple[np.ndarray, np.ndarray]:
    """Return x, y arrays of length n+1, including the (0,0) origin."""
    bases = encode_bases(sequence)
    n = int(bases.size)
    x = np.zeros(n + 1, dtype=np.float64)
    y = np.zeros(n + 1, dtype=np.float64)
    px = 0.0
    py = 0.0
    for i, base in enumerate(bases.tolist(), start=1):
        if base > 3:
            px = 0.0
            py = 0.0
        else:
            px = 0.5 * (px + float(_VX[base]))
            py = 0.5 * (py + float(_VY[base]))
        x[i] = px
        y[i] = py
    return x, y


def sequence_to_fcgr(
    sequence: str,
    k: int = IMAGE_SIZE,
    normalize: bool = True,
) -> np.ndarray:
    """64x64 float32 image. If normalize, values are in [0, 1] (divide by 255)."""
    bins = int(k) if int(k) >= 8 else IMAGE_SIZE
    x, y = sequence_to_cgr_points(sequence)
    # numpy H[i,j] counts x in bin i, y in bin j
    hist, _, _ = np.histogram2d(x, y, bins=bins)
    # Stored MATLAB images match numpy histogram2d(x, y) followed by rot90,
    # not histcounts2's (y, x) layout plus rot90.
    rotated = np.rot90(hist)
    peak = float(rotated.max())
    if peak > 0:
        cgri = np.floor((255.0 / peak) * rotated)
    else:
        cgri = rotated
    if normalize:
        return (cgri / 255.0).astype(np.float32)
    return cgri.astype(np.float32)


def fcgr_to_uint8(fcgr: np.ndarray) -> np.ndarray:
    if fcgr.max() <= 1.0:
        return np.clip(np.round(fcgr * 255.0), 0, 255).astype(np.uint8)
    return np.clip(fcgr, 0, 255).astype(np.uint8)


def ambiguous_fraction(sequence: str) -> float:
    if not sequence:
        return 1.0
    bases = encode_bases(sequence)
    if bases.size == 0:
        return 1.0
    return float(np.mean(bases > 3))
