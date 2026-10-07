from __future__ import annotations

import pytest
from records import atom, link, pair, ssbond, write_pdb

from pdbcheck import check_clashes, classify
from pdbcheck.parser import parse_pdb

# --------------------------------------------------------------------------
# the rule table itself
# --------------------------------------------------------------------------


@pytest.mark.parametrize("distance,elements,expected", [
    # below 2.2 A nothing non-bonded can explain the distance
    (1.00, ("C", "C"), "steric"),
    (2.19, ("C", "O"), "steric"),
    (1.90, ("N", "O"), "tight"),
    (2.19, ("O", "O"), "tight"),
    (2.10, ("S", "N"), "tight"),
    # 2.2 to 2.6 A is short for any pair
    (2.20, ("C", "C"), "polar"),
    (2.59, ("N", "O"), "polar"),
    # 2.6 to 3.0 A is a clash only when a carbon takes part
    (2.60, ("C", "C"), "carbon"),
    (2.99, ("C", "O"), "carbon"),
    (2.60, ("N", "O"), None),  # an ordinary hydrogen bond
    (2.90, ("O", "O"), None),
    (2.75, ("N", "N"), None),
    # the thresholds themselves count as healthy
    (3.00, ("C", "C"), None),
    (3.60, ("C", "C"), None),
])
def test_rule_table(distance, elements, expected):
    assert classify(distance, *elements) == expected


def test_carbon_on_either_side_counts(tmp_path):
    assert classify(2.8, "C", "O") == "carbon"
    assert classify(2.8, "O", "C") == "carbon"


# --------------------------------------------------------------------------
# which pairs are looked at
# --------------------------------------------------------------------------


def test_inter_chain_overlap_is_reported(tmp_path):
    contacts = check_clashes(write_pdb(tmp_path, pair(1.0)))
    assert len(contacts) == 1
    assert contacts[0].severity == "steric"
    assert contacts[0].distance == 1.0


def test_healthy_distance_is_not_reported(tmp_path):
    assert check_clashes(write_pdb(tmp_path, pair(2.8, "N", "O"))) == []


def test_same_chain_is_ignored(tmp_path):
    assert check_clashes(write_pdb(tmp_path, pair(1.0, chain_b="A"))) == []


def test_pair_order_is_canonical(tmp_path):
    """The same clash reads the same way no matter the order in the file."""
    path = write_pdb(tmp_path, [
        atom(1, "CB", "ALA", "B", 7, 0.0, 0.0, 0.0, "C"),
        atom(2, "CA", "ALA", "A", 3, 1.0, 0.0, 0.0, "C"),
    ])
    contact = check_clashes(path)[0]
    assert (contact.a.chain, contact.b.chain) == ("A", "B")


def test_results_are_sorted_worst_first(tmp_path):
    path = write_pdb(tmp_path, [
        atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
        atom(2, "CB", "ALA", "B", 1, 2.8, 0.0, 0.0, "C"),
        atom(3, "N", "ALA", "A", 2, 0.0, 50.0, 0.0, "N"),
        atom(4, "O", "ALA", "B", 2, 2.4, 50.0, 0.0, "O"),
        atom(5, "CA", "ALA", "A", 3, 0.0, 100.0, 0.0, "C"),
        atom(6, "CB", "ALA", "B", 3, 1.5, 100.0, 0.0, "C"),
    ])
    assert [c.severity for c in check_clashes(path)] == ["steric", "polar", "carbon"]


# --------------------------------------------------------------------------
# bonds that a distance test alone would report
# --------------------------------------------------------------------------


def test_disulfide_is_not_a_clash(tmp_path):
    path = write_pdb(tmp_path, [
        ssbond("A", 30, "B", 80),
        atom(1, "SG", "CYS", "A", 30, 0.0, 0.0, 0.0, "S"),
        atom(2, "SG", "CYS", "B", 80, 2.05, 0.0, 0.0, "S"),
    ])
    assert parse_pdb(path).ssbonds == frozenset({(("A", "30"), ("B", "80"))})
    assert check_clashes(path) == []


def test_undeclared_disulfide_is_not_a_clash(tmp_path):
    """SSBOND records go missing often enough that SG to SG is excluded outright."""
    path = write_pdb(tmp_path, [
        atom(1, "SG", "CYS", "A", 30, 0.0, 0.0, 0.0, "S"),
        atom(2, "SG", "CYS", "B", 80, 2.05, 0.0, 0.0, "S"),
    ])
    assert check_clashes(path) == []


def test_declared_link_is_not_a_clash(tmp_path):
    path = write_pdb(tmp_path, [
        link("C1", "XYZ", "A", 10, "OG", "SER", "B", 20),
        atom(1, "C1", "XYZ", "A", 10, 0.0, 0.0, 0.0, "C", record="HETATM"),
        atom(2, "OG", "SER", "B", 20, 1.50, 0.0, 0.0, "O"),
    ])
    assert len(parse_pdb(path).links) == 1
    assert check_clashes(path) == []


def test_glycosylation_is_not_a_clash(tmp_path):
    """N-linked sugars are bonded to Asn by a bond LINK rarely declares."""
    path = write_pdb(tmp_path, [
        atom(1, "ND2", "ASN", "A", 50, 0.0, 0.0, 0.0, "N"),
        atom(2, "C1", "NAG", "B", 1, 1.45, 0.0, 0.0, "C", record="HETATM"),
    ])
    assert check_clashes(path) == []


def test_sugar_to_sugar_is_not_a_clash(tmp_path):
    path = write_pdb(tmp_path, [
        atom(1, "O4", "NAG", "A", 1, 0.0, 0.0, 0.0, "O", record="HETATM"),
        atom(2, "C1", "BMA", "B", 2, 1.44, 0.0, 0.0, "C", record="HETATM"),
    ])
    assert check_clashes(path) == []


def test_a_real_overlap_near_a_sugar_still_counts(tmp_path):
    """The glycosidic exemption only reaches as far as a bond length."""
    path = write_pdb(tmp_path, [
        atom(1, "ND2", "ASN", "A", 50, 0.0, 0.0, 0.0, "N"),
        atom(2, "C1", "NAG", "B", 1, 1.95, 0.0, 0.0, "C", record="HETATM"),
    ])
    assert [c.severity for c in check_clashes(path)] == ["steric"]


# --------------------------------------------------------------------------
# atoms that never enter the search
# --------------------------------------------------------------------------


def test_hydrogen_is_dropped(tmp_path):
    path = write_pdb(tmp_path, pair(0.9, element_b="H", name_b="HB"))
    assert [a.element for a in parse_pdb(path).atoms] == ["C"]
    assert check_clashes(path) == []


def test_alternate_conformation_is_dropped(tmp_path):
    path = write_pdb(tmp_path, pair(1.0, altloc="B"))
    assert len(parse_pdb(path).atoms) == 1
    assert check_clashes(path) == []


def test_zero_occupancy_is_dropped(tmp_path):
    path = write_pdb(tmp_path, pair(1.0, occupancy=0.0))
    assert len(parse_pdb(path).atoms) == 1
    assert check_clashes(path) == []


def test_solvent_is_dropped(tmp_path):
    path = write_pdb(tmp_path, [
        atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
        atom(2, "O", "HOH", "B", 1, 1.0, 0.0, 0.0, "O", record="HETATM"),
        atom(3, "ZN", "ZN", "C", 1, 1.2, 0.0, 0.0, "ZN", record="HETATM"),
    ])
    assert [a.resn for a in parse_pdb(path).atoms] == ["ALA"]
    assert check_clashes(path) == []


def test_only_the_first_model_is_read(tmp_path):
    """A second model is another copy of the same chains, not a neighbour."""
    path = write_pdb(tmp_path, [
        "MODEL        1",
        *pair(5.0),
        "ENDMDL",
        "MODEL        2",
        *pair(1.0),
        "ENDMDL",
    ])
    assert len(parse_pdb(path).atoms) == 2
    assert check_clashes(path) == []


def test_insertion_codes_survive(tmp_path):
    """Antibody numbering leans on them, so they must reach the report intact."""
    path = write_pdb(tmp_path, [
        atom(1, "CA", "ALA", "A", "100A", 0.0, 0.0, 0.0, "C"),
        atom(2, "CB", "ALA", "B", "100B", 1.0, 0.0, 0.0, "C"),
    ])
    contact = check_clashes(path)[0]
    assert (contact.a.resi, contact.b.resi) == ("100A", "100B")
