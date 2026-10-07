"""Read a PDB file into the atoms and bond records that clash detection needs.

Only the part of the format that can change a verdict is parsed: atom
coordinates, plus the LINK and SSBOND records that declare covalent bonds a
plain distance test would otherwise report as overlaps.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

# Crystallisation additives, buffer components and free ions. None of them are
# part of the structure, and all of them sit close to it by construction.
SOLVENT = frozenset({
    "HOH", "WAT", "SO4", "EDO", "GOL", "PEG",
    "CL", "NA", "MG", "NI", "ZN", "CA",
})

_MIN_COORD_LINE = 54  # an ATOM record is unusable before the z column closes


class Atom(NamedTuple):
    chain: str
    resi: str  # residue number with its insertion code, e.g. 100A
    resn: str
    name: str
    element: str
    x: float
    y: float
    z: float


class Structure(NamedTuple):
    atoms: list[Atom]
    links: frozenset  # declared covalent links, keyed by (chain, resi, atom)
    ssbonds: frozenset  # declared disulfides, keyed by (chain, resi)


def parse_pdb(path: Path) -> Structure:
    """Atoms of the first model, with the connectivity records that explain them."""
    atoms: list[Atom] = []
    links = set()
    ssbonds = set()

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6]

            if record == "ENDMDL":
                break  # every later model is another copy of the same structure

            if record == "SSBOND":
                ssbonds.add((
                    (line[15], line[17:22].strip()),
                    (line[29], line[31:36].strip()),
                ))
            elif record == "LINK  ":
                links.add((
                    (line[21], line[22:27].strip(), line[12:16].strip()),
                    (line[51], line[52:57].strip(), line[42:46].strip()),
                ))
            elif record in ("ATOM  ", "HETATM"):
                atom = _read_atom(line)
                if atom is not None:
                    atoms.append(atom)

    return Structure(atoms, frozenset(links), frozenset(ssbonds))


def _read_atom(line: str) -> Atom | None:
    """One coordinate record, or None when it must not take part in the search."""
    if len(line) < _MIN_COORD_LINE:
        return None

    if line[16] not in " A":
        return None  # keep a single alternate conformation, never two at once

    resn = line[17:20].strip()
    if resn in SOLVENT:
        return None

    element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
    if element == "H":
        return None  # every threshold in this package is a heavy-atom threshold

    try:
        if float(line[54:60]) == 0.0:
            return None  # modelled, but flagged as not actually present
    except ValueError:
        pass  # a missing occupancy column means the atom counts

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    return Atom(
        chain=line[21],
        resi=line[22:27].strip(),
        resn=resn,
        name=line[12:16].strip(),
        element=element,
        x=x,
        y=y,
        z=z,
    )
