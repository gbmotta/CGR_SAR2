"""WHO labels for the six variants in Experiment 1."""

from __future__ import annotations

# User-facing order
VARIANT_LABELS: tuple[str, ...] = (
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "Iota",
    "Epsilon",
)

VARIANT_TO_INDEX: dict[str, int] = {name: i for i, name in enumerate(VARIANT_LABELS)}

# LabelIndex in the original MATLAB/Dropbox CSVs (unique() of full GISAID names
# in first-appearance order of table_final).
MATLAB_INDEX_TO_VARIANT: tuple[str, ...] = (
    "Delta",    # 0
    "Alpha",    # 1
    "Beta",     # 2
    "Epsilon",  # 3
    "Iota",     # 4
    "Gamma",    # 5
)

VARIANT_WHO_STATUS: dict[str, str] = {
    "Alpha": "VOC (B.1.1.7)",
    "Beta": "VOC (B.1.351)",
    "Gamma": "VOC (P.1)",
    "Delta": "VOC (B.1.617.2)",
    "Iota": "VOI (B.1.526)",
    "Epsilon": "VOI (B.1.427/B.1.429)",
}

NCBI_LINEAGES: dict[str, tuple[str, ...]] = {
    "Alpha": ("B.1.1.7",),
    "Beta": ("B.1.351",),
    "Gamma": ("P.1",),
    "Delta": ("B.1.617.2",),
    "Iota": ("B.1.526",),
    "Epsilon": ("B.1.427", "B.1.429"),
}

MIN_GENOME_LENGTH = 29_000
MAX_AMBIGUOUS_FRACTION = 0.01


def who_name_from_gisaid(variant: str) -> str:
    """Map a GISAID Variant string to one of the six WHO labels."""
    text = variant.upper()
    for name in VARIANT_LABELS:
        if name.upper() in text:
            return name
    raise ValueError(f"Variante GISAID nao reconhecida: {variant!r}")
