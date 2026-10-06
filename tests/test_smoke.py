from __future__ import annotations

from pathlib import Path

from pdbcheck import __version__
from pdbcheck.cli import build_parser, collect_inputs, default_output


def test_version_is_exposed():
    assert __version__ == "0.1.0"


def test_parser_reads_a_target():
    args = build_parser().parse_args(["/data/lmk/ABCD.pdb"])
    assert args.target == Path("/data/lmk/ABCD.pdb")
    assert args.out is None


def test_directory_scan_is_flat_and_sorted(tmp_path):
    (tmp_path / "b.pdb").touch()
    (tmp_path / "a.PDB").touch()
    (tmp_path / "notes.txt").touch()
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "deep.pdb").touch()

    found = [p.name for p in collect_inputs(tmp_path)]
    assert found == ["a.PDB", "b.pdb"]


def test_csv_lands_beside_the_target(tmp_path):
    f = tmp_path / "x.pdb"
    f.touch()
    assert default_output(f) == tmp_path / "pdb_check.csv"
    assert default_output(tmp_path) == tmp_path / "pdb_check.csv"
