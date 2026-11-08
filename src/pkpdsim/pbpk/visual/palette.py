"""The colours of the figures: organs, layers of the body, routes and connections."""

from types import MappingProxyType

__all__ = [
    "BILE",
    "BLOOD",
    "BONE",
    "FAT",
    "FECES",
    "INK",
    "IV",
    "MUSCLE",
    "ORAL",
    "ORGAN_COLORS",
    "PORTAL",
    "RING",
    "SKIN",
    "THORAX",
    "URINE",
    "VAT",
]

INK = "#1b1b1b"
"""The outline of the body and its details."""
SKIN = "#f1d6c2"
FAT = "#f4dc8e"
"""The subcutaneous adipose tissue."""
VAT = "#efcf6a"
"""The visceral adipose tissue, the space of the abdomen between the organs."""
MUSCLE = "#d9907e"
BONE = "#efe6d3"
THORAX = "#e2c3b8"
"""The mediastinum, the space of the thorax between the organs."""
PORTAL = "#6b4fa0"
BILE = "#3f8f55"
URINE = "#c7a200"
FECES = "#7a5a3a"
ORAL = "#d0731f"
IV = "#555555"
BLOOD = "#c0392b"
RING = "#8a8f93"
"""The circle of an organ icon."""

ORGAN_COLORS = MappingProxyType(
    {
        "lu": "#9cc3e6",
        "he": "#c8504b",
        "br": "#c9aee0",
        "li": "#9a5038",
        "ki": "#9d3f5c",
        "st": "#e7a08e",
        "si": "#f0c592",
        "la": "#cf9f6c",
        "sp": "#5e4a9c",
        "pa": "#ecd36f",
        "go": "#e6a3a8",
        "mu": MUSCLE,
        "bo": BONE,
        "ad": FAT,
        "sk": SKIN,
        "gallbladder": "#69b277",
        "bladder": "#efe08a",
        "uterus": "#dc8f8f",
    }
)
"""The colour of every organ of the library and of the further anatomy."""
