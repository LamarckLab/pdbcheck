"""Command line entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pdbcheck import __version__

OUTPUT_NAME = "pdb_check.csv"
SUFFIXES = (".pdb", ".PDB")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdbcheck",
        description="Screen PDB files for inter-chain steric clashes.",
    )
    p.add_argument("target", type=Path,
                   help="a .pdb file, or a directory containing .pdb files")
    p.add_argument("-o", "--out", type=Path, default=None,
                   help=f"output CSV path (default: {OUTPUT_NAME} beside TARGET)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="print nothing on success")
    p.add_argument("--version", action="version", version=f"pdbcheck {__version__}")
    return p


def collect_inputs(target: Path) -> list[Path]:
    """One file, or every .pdb directly inside a directory. Never recursive."""
    if target.is_file():
        return [target]
    if target.is_dir():
        return sorted(p for p in target.iterdir() if p.is_file() and p.suffix in SUFFIXES)
    raise FileNotFoundError(target)


def default_output(target: Path) -> Path:
    """The CSV always lands in the directory the target lives in."""
    return (target.parent if target.is_file() else target) / OUTPUT_NAME


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        files = collect_inputs(args.target)
    except FileNotFoundError:
        print(f"pdbcheck: no such file or directory: {args.target}", file=sys.stderr)
        return 2

    if not files:
        print(f"pdbcheck: no .pdb files in {args.target}", file=sys.stderr)
        return 2

    out = args.out or default_output(args.target)

    # TODO stage 2: parse, detect clashes, write the CSV
    if not args.quiet:
        print(f"[stub] {len(files)} file(s) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
