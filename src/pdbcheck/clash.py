"""Decide which close atom pairs are steric clashes.

Nitrogen, oxygen and sulfur can hydrogen bond; carbon cannot. The same
separation therefore means different things depending on whether a carbon takes
part, which is why the thresholds below are split rather than single-valued.

The two outer values reproduce the usual 0.4 A van der Waals overlap criterion:
two carbons sum to 3.40 A, nitrogen and oxygen to 3.07 A.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree

from pdbcheck.parser import Atom, Structure

HARD = 2.2  # closer than this, no non-bonded pair has a physical explanation
POLAR = 2.6  # N/O/S may hydrogen bond this close, carbon may not
CARBON = 3.0  # a carbon this close to anything is already too tight

# Sugars are joined to the protein, and to each other, by bonds that LINK
# records routinely leave undeclared.
SUGAR = frozenset({
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL",
    "GLC", "A2G", "SIA", "XYS", "BGC",
})
GLYCO_ACCEPTOR = frozenset({("ASN", "ND2"), ("SER", "OG"), ("THR", "OG1")})
GLYCO_BOND = 1.8

SEVERITY_ORDER = ("steric", "tight", "polar", "carbon")


class Contact(NamedTuple):
    severity: str
    distance: float
    a: Atom
    b: Atom


def classify(distance: float, element_a: str, element_b: str) -> str | None:
    """Severity of one non-bonded contact, or None when the distance is healthy."""
    has_carbon = element_a == "C" or element_b == "C"
    if distance < HARD:
        return "steric" if has_carbon else "tight"
    if distance < POLAR:
        return "polar"
    if distance < CARBON and has_carbon:
        return "carbon"
    return None


def is_bonded(a: Atom, b: Atom, links: frozenset, ssbonds: frozenset, distance: float) -> bool:
    """True when a declared or implied covalent bond already explains the distance."""
    key_a = (a.chain, a.resi, a.name)
    key_b = (b.chain, b.resi, b.name)
    if (key_a, key_b) in links or (key_b, key_a) in links:
        return True

    if a.name == "SG" and b.name == "SG":
        return True  # disulfides sit near 2.05 A whether or not SSBOND declares them

    residue_a = (a.chain, a.resi)
    residue_b = (b.chain, b.resi)
    if (residue_a, residue_b) in ssbonds or (residue_b, residue_a) in ssbonds:
        return True

    if distance < GLYCO_BOND:
        for sugar, partner in ((a, b), (b, a)):
            if sugar.resn not in SUGAR:
                continue
            if sugar.name == "C1" and (partner.resn, partner.name) in GLYCO_ACCEPTOR:
                return True
            if partner.resn in SUGAR:
                return True

    return False


def _residue_sort_key(resi: str) -> tuple:
    """Order 9 before 10 before 100A, which plain string order would not."""
    digits = ""
    index = 0
    if resi[:1] == "-":
        digits, index = "-", 1
    while index < len(resi) and resi[index].isdigit():
        digits += resi[index]
        index += 1
    number = int(digits) if digits not in ("", "-") else 0
    return (number, resi[index:])


def _atom_sort_key(atom: Atom) -> tuple:
    return (atom.chain, _residue_sort_key(atom.resi), atom.name)


def find_clashes(structure: Structure) -> list[Contact]:
    """Every inter-chain contact in the structure that no bond accounts for.

    Contacts inside a single chain are left alone in this version. The atoms of
    one chain are held together by bonds this package does not model, so testing
    them would report far more ordinary geometry than genuine error.
    """
    atoms = structure.atoms
    if len(atoms) < 2:
        return []

    coordinates = np.array([(a.x, a.y, a.z) for a in atoms], dtype=float)
    tree = cKDTree(coordinates)

    contacts: list[Contact] = []
    for i, j in tree.query_pairs(CARBON):
        a, b = atoms[i], atoms[j]
        if a.chain == b.chain:
            continue

        distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
        if is_bonded(a, b, structure.links, structure.ssbonds, distance):
            continue

        severity = classify(distance, a.element, b.element)
        if severity is None:
            continue

        if _atom_sort_key(a) > _atom_sort_key(b):
            a, b = b, a  # one fixed order, so a clash always reads the same way
        contacts.append(Contact(severity, round(distance, 2), a, b))

    contacts.sort(key=lambda c: (
        SEVERITY_ORDER.index(c.severity),
        c.distance,
        _atom_sort_key(c.a),
        _atom_sort_key(c.b),
    ))
    return contacts
