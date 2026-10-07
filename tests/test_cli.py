from __future__ import annotations

import csv

from records import atom, pair, write_pdb

from pdbcheck.cli import main
from pdbcheck.report import COLUMNS


def read_rows(path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_single_file_reports_beside_the_target(tmp_path, capsys):
    target = write_pdb(tmp_path, pair(1.0), name="ABCD.pdb")
    assert main([str(target)]) == 0
    assert (tmp_path / "pdb_check.csv").exists()
    assert "ABCD.pdb analysed" in capsys.readouterr().out


def test_csv_header_is_the_agreed_schema(tmp_path):
    main([str(write_pdb(tmp_path, pair(1.0)))])
    with open(tmp_path / "pdb_check.csv", encoding="utf-8", newline="") as handle:
        assert next(csv.reader(handle)) == COLUMNS


def test_a_clash_row_names_both_atoms(tmp_path):
    target = write_pdb(tmp_path, [
        atom(1, "NZ", "LYS", "A", 45, 0.0, 0.0, 0.0, "N"),
        atom(2, "CZ", "PHE", "B", 112, 1.08, 0.0, 0.0, "C"),
    ])
    main([str(target)])

    row = read_rows(tmp_path / "pdb_check.csv")[0]
    assert row["severity"] == "steric"
    assert row["distance"] == "1.08"
    assert (row["chain1"], row["resn1"], row["resi1"], row["atom1"]) == ("A", "LYS", "45", "NZ")
    assert (row["chain2"], row["resn2"], row["resi2"], row["atom2"]) == ("B", "PHE", "112", "CZ")
    assert row["same_chain"] == "False"


def test_a_clean_file_still_gets_a_row(tmp_path):
    main([str(write_pdb(tmp_path, pair(5.0)))])
    rows = read_rows(tmp_path / "pdb_check.csv")
    assert len(rows) == 1
    assert rows[0]["file"] == "test.pdb"
    assert rows[0]["severity"] == "clean"
    assert rows[0]["chain1"] == ""


def test_directory_scan_covers_every_file(tmp_path, capsys):
    write_pdb(tmp_path, pair(1.0), name="bad.pdb")
    write_pdb(tmp_path, pair(5.0), name="good.PDB")
    (tmp_path / "notes.txt").write_text("ignored", encoding="utf-8")

    assert main([str(tmp_path)]) == 0
    rows = read_rows(tmp_path / "pdb_check.csv")
    assert {row["file"] for row in rows} == {"bad.pdb", "good.PDB"}
    assert "2 PDB files processed" in capsys.readouterr().out


def test_report_is_replaced_not_appended(tmp_path):
    out = tmp_path / "pdb_check.csv"
    out.write_text("stale,content\n1,2\n", encoding="utf-8")
    main([str(write_pdb(tmp_path, pair(1.0)))])
    assert "stale" not in out.read_text(encoding="utf-8")


def test_out_overrides_the_destination(tmp_path):
    target = write_pdb(tmp_path, pair(1.0))
    elsewhere = tmp_path / "reports"
    elsewhere.mkdir()
    out = elsewhere / "report.csv"

    main([str(target), "-o", str(out)])
    assert out.exists()
    assert not (tmp_path / "pdb_check.csv").exists()


def test_quiet_prints_nothing(tmp_path, capsys):
    assert main([str(write_pdb(tmp_path, pair(1.0))), "-q"]) == 0
    assert capsys.readouterr().out == ""


def test_an_unreadable_file_does_not_stop_the_batch(tmp_path):
    write_pdb(tmp_path, pair(1.0), name="good.pdb")
    (tmp_path / "empty.pdb").write_text("", encoding="utf-8")

    assert main([str(tmp_path)]) == 1
    severities = {row["file"]: row["severity"] for row in read_rows(tmp_path / "pdb_check.csv")}
    assert severities["empty.pdb"] == "error"
    assert severities["good.pdb"] == "steric"


def test_missing_path_is_an_argument_error(tmp_path):
    assert main([str(tmp_path / "nope.pdb")]) == 2
