"""Fast inter-chain steric clash screening for PDB structures."""

from __future__ import annotations

from pathlib import Path

from pdbcheck.clash import Contact, classify, find_clashes
from pdbcheck.parser import Atom, Structure, parse_pdb

__version__ = "0.1.0"


def check_clashes(path: Path) -> list[Contact]:
    """Read one PDB file and return its inter-chain clashes, worst first."""
    return find_clashes(parse_pdb(path))


__all__ = [
    "Atom",
    "Contact",
    "Structure",
    "__version__",
    "check_clashes",
    "classify",
    "find_clashes",
    "parse_pdb",
]
