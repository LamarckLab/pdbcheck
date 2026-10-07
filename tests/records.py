"""Build minimal PDB records, so every test states the exact geometry it means.

The column offsets here are the ones the format fixes, which is what makes these
helpers worth having: a test that places two atoms 1.5 A apart should fail only
when the detection logic is wrong, never because a field drifted a column.
"""

from __future__ import annotations

from pathlib import Path


def _split_resi(resi) -> tuple[str, str]:
    """Separate a residue number from its insertion code, as the format does."""
    text = str(resi)
    if text and not text[-1].isdigit():
        return text[:-1], text[-1]
    return text, " "


def atom(serial, name, resn, chain, resi, x, y, z, element,
         altloc=" ", occupancy=1.0, record="ATOM  ") -> str:
    number, icode = _split_resi(resi)
    return (
        f"{record:<6.6s}"
        f"{serial:>5d} "
        f"{name:<4.4s}"
        f"{altloc}"
        f"{resn:>3.3s} "
        f"{chain}"
        f"{int(number):>4d}"
        f"{icode}"
        "   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}"
        f"{0.0:6.2f}"
        f"{'':10s}"
        f"{element:>2.2s}"
    )


def link(name1, resn1, chain1, resi1, name2, resn2, chain2, resi2) -> str:
    number1, icode1 = _split_resi(resi1)
    number2, icode2 = _split_resi(resi2)
    return (
        "LINK        "
        f"{name1:<4.4s} "
        f"{resn1:>3.3s} "
        f"{chain1}"
        f"{int(number1):>4d}"
        f"{icode1}"
        f"{'':15s}"
        f"{name2:<4.4s} "
        f"{resn2:>3.3s} "
        f"{chain2}"
        f"{int(number2):>4d}"
        f"{icode2}"
    )


def ssbond(chain1, resi1, chain2, resi2, serial=1) -> str:
    number1, icode1 = _split_resi(resi1)
    number2, icode2 = _split_resi(resi2)
    return (
        "SSBOND"
        f"{serial:>4d} "
        "CYS "
        f"{chain1} "
        f"{int(number1):>4d}"
        f"{icode1}"
        "   CYS "
        f"{chain2} "
        f"{int(number2):>4d}"
        f"{icode2}"
    )


def write_pdb(directory, lines, name="test.pdb") -> Path:
    path = Path(directory) / name
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


def pair(distance, element_a="C", element_b="C", chain_b="B",
         name_a="CA", name_b="CB", resn_a="ALA", resn_b="ALA", **kwargs) -> list[str]:
    """Two atoms a fixed distance apart, on different chains by default."""
    return [
        atom(1, name_a, resn_a, "A", 1, 0.0, 0.0, 0.0, element_a),
        atom(2, name_b, resn_b, chain_b, 1, distance, 0.0, 0.0, element_b, **kwargs),
    ]
