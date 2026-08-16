"""SARS-CoV-2 variant classification with FCGR + CNN."""

from .cgr import K_DEFAULT, IMAGE_SIZE, sequence_to_fcgr, fcgr_to_uint8
from .labels import VARIANT_LABELS, VARIANT_TO_INDEX
from .model import VariantCNN, count_parameters
from .predict import VariantClassifier

__version__ = "1.0.0"
__all__ = [
    "IMAGE_SIZE",
    "K_DEFAULT",
    "VARIANT_LABELS",
    "VARIANT_TO_INDEX",
    "VariantCNN",
    "VariantClassifier",
    "count_parameters",
    "fcgr_to_uint8",
    "sequence_to_fcgr",
]
