"""Write findings as one tidy CSV row per clashing atom pair."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

from pdbcheck.clash import Contact

COLUMNS = [
    "file", "severity", "distance",
    "chain1", "resn1", "resi1", "atom1", "elem1",
    "chain2", "resn2", "resi2", "atom2", "elem2",
    "same_chain",
]


def contact_row(name: str, contact: Contact) -> dict:
    a, b = contact.a, contact.b
    return {
        "file": name,
        "severity": contact.severity,
        "distance": f"{contact.distance:.2f}",
        "chain1": a.chain, "resn1": a.resn, "resi1": a.resi, "atom1": a.name, "elem1": a.element,
        "chain2": b.chain, "resn2": b.resn, "resi2": b.resi, "atom2": b.name, "elem2": b.element,
        "same_chain": "False",
    }


def _bare_row(name: str, severity: str) -> dict:
    row = dict.fromkeys(COLUMNS, "")
    row["file"] = name
    row["severity"] = severity
    return row


def clean_row(name: str) -> dict:
    """A file that was read and held no clashes still earns a row.

    Without it a batch report cannot tell a clean structure apart from one that
    was never looked at.
    """
    return _bare_row(name, "clean")


def error_row(name: str) -> dict:
    """A file that could not be read, so that it is not silently missing."""
    return _bare_row(name, "error")


def write_csv(path: Path, rows: Iterable[dict]) -> None:
    """Replace any report already sitting at path."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
