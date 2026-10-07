# pdbcheck

Fast inter-chain steric clash screening for PDB structures.

[![CI](https://github.com/LamarckLab/pdbcheck/actions/workflows/ci.yml/badge.svg)](https://github.com/LamarckLab/pdbcheck/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pdbcheck.svg)](https://pypi.org/project/pdbcheck/)
[![Python](https://img.shields.io/pypi/pyversions/pdbcheck.svg)](https://pypi.org/project/pdbcheck/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/LamarckLab/pdbcheck/blob/main/LICENSE)

`pdbcheck` reads one PDB file, or a directory of them, and writes a CSV naming every
pair of atoms from different chains that sit closer together than chemistry allows.

It exists for triaging the output of structure predictors — AlphaFold, Chai-1, Boltz,
RFdiffusion — where the question is which of five hundred models to throw away, not
how a deposited structure scored. For authoritative validation of an experimental
structure, use MolProbity: it places hydrogens, computes a normalised clashscore, and
is what the PDB itself reports. `pdbcheck` deliberately does less, and in exchange
runs on heavy atoms alone, finishes a structure in milliseconds, needs no setup beyond
`pip install`, and emits a table meant for a dataframe rather than a report meant for
a human.

## Installation

```bash
pip install pdbcheck
```

Python 3.9 or newer. `numpy` and `scipy` are the only dependencies, and they are
pinned loosely so that installing into an existing environment leaves whatever is
already there alone.

To get the command without adding anything to an environment:

```bash
pipx install pdbcheck
```

Or straight from source:

```bash
pip install git+https://github.com/LamarckLab/pdbcheck.git
```

## Usage

### Command line

```bash
pdbcheck /data/structures/model.pdb     # one file
pdbcheck /data/structures               # every .pdb in the directory
```

Both forms write `pdb_check.csv` into the directory the target lives in, replacing any
report already sitting there. The directory scan is flat — subdirectories are not
visited — and recognises `.pdb` and `.PDB` only.

| Option | Effect |
| :--- | :--- |
| `-o, --out PATH` | write the CSV somewhere else |
| `-q, --quiet` | print nothing on success |
| `--version` | print the version and exit |

Exit status is `0` on success, `1` when a file could not be read, and `2` when the path
does not exist or holds no PDB files. Clashes do not affect it: a structure riddled with
overlaps is still a successful run.

### Python

```python
from pdbcheck import check_clashes

for contact in check_clashes("model.pdb"):
    print(contact.severity, contact.distance, contact.a.chain, contact.a.resi)
```

`check_clashes` returns a list of `Contact`, worst first. Each carries `severity`,
`distance`, and the two `Atom` records `a` and `b`, which hold `chain`, `resi`, `resn`,
`name`, `element` and `x`/`y`/`z`.

## Output

One row per clashing atom pair, ordered worst first.

| Column | Meaning |
| :--- | :--- |
| `file` | name of the structure the row came from |
| `severity` | `steric`, `tight`, `polar`, `carbon`, `clean` or `error` |
| `distance` | centre-to-centre separation in Angstrom, two decimals |
| `chain1` `resn1` `resi1` `atom1` `elem1` | first atom of the pair |
| `chain2` `resn2` `resi2` `atom2` `elem2` | second atom of the pair |
| `same_chain` | always `False` in this version |

```
file      severity  distance  chain1 resn1 resi1 atom1 elem1  chain2 resn2 resi2 atom2 elem2  same_chain
1IGT.pdb  polar     2.38      D      ASN   314   ND2   N      F      NAG   1     O5    O      False
1IGT.pdb  carbon    2.81      A      GLY   91    O     O      B      TYR   100I  CD2   C      False
```

`resi` keeps the insertion code when the structure carries one, as `100I` above, so
Kabat and IMGT numbered antibody models survive the round trip intact.

A structure with no clashes still produces a row, with `severity` set to `clean` and the
atom columns empty. A structure that could not be read produces one set to `error`. A
batch report therefore accounts for every file it was handed, and a missing row means a
file was never seen rather than that it was fine.

## The criterion

Nitrogen, oxygen and sulfur form hydrogen bonds. Carbon does not. A nitrogen and an
oxygen 2.8 Å apart are a textbook hydrogen bond; a carbon 2.8 Å from anything is a
modelling error. No single distance cut-off can express both, so the thresholds split on
whether a carbon takes part.

| Separation (Å) | Carbon involved | No carbon (N/O/S only) |
| :--- | :--- | :--- |
| `d < 2.2` | **steric** | **tight** |
| `2.2 <= d < 2.6` | polar | polar |
| `2.6 <= d < 3.0` | carbon | healthy |
| `d >= 3.0` | healthy | healthy |

`steric` and `tight` are the verdicts worth acting on. `polar` and `carbon` mark contacts
that are tighter than average but do occur in well refined structures; read them as a
count, not as a list of defects.

The outer thresholds reproduce the 0.4 Å van der Waals overlap criterion MolProbity uses.
Two carbons sum to 3.40 Å, and 3.40 − 0.4 = 3.00. Nitrogen and oxygen sum to 3.07 Å, and
3.07 − 0.4 = 2.67, rounded down to 2.6 here so that the grade errs toward silence. Below
2.2 Å nothing survives: the shortest known low barrier hydrogen bond is near 2.4 Å, so a
non-bonded pair closer than 2.2 Å has no physical account at all.

The split at 2.2 Å separates a certain error from a contact that still deserves an eye.
A carbon that close is overlap and nothing else, so it is graded `steric`; the same
distance between polar atoms is graded `tight`, because an extreme salt bridge or a metal
site can in principle reach it.

## What is excluded

Two filters run before any distance is measured. They are the reason the grades above
mean anything, and most of the work in this package lives in them.

Atoms that never enter the search:

| Dropped | Reason |
| :--- | :--- |
| hydrogens | every threshold above is a heavy-atom threshold |
| alternate conformations other than `A` | two positions of one atom necessarily overlap |
| zero occupancy | modelled, but flagged as not actually present |
| `HOH WAT SO4 EDO GOL PEG CL NA MG NI ZN CA` | crystallisation additives and free ions |
| every model after the first | later models are copies of the same chains, not neighbours |

Pairs that are bonds rather than clashes:

| Excluded | Reason |
| :--- | :--- |
| any pair named in a `LINK` record | a declared covalent link |
| any pair named in an `SSBOND` record | a declared disulfide |
| any `SG`–`SG` pair | SSBOND records go missing often enough to warrant the blanket rule |
| sugar `C1` to `ASN ND2`, `SER OG` or `THR OG1`, below 1.8 Å | glycosylation, which LINK rarely declares |
| sugar to sugar, below 1.8 Å | the glycan's own backbone, likewise undeclared |

## Scope and limitations

### Inter-chain contacts only

Contacts inside a single chain are not reported. The atoms of one chain are joined by
bonds this package does not model — backbone, side chain, the proline ring — and a
distance test blind to them flags every peptide bond in the structure. Doing it properly
means building a connectivity graph and excluding pairs one and two bonds apart. That is
planned, and it is not here yet.

The consequence is worth stating plainly: **a single-chain structure always comes back
`clean`**, because it contains no inter-chain pairs at all. `pdbcheck` has nothing to say
about a monomer in this version.

### Known residual false positives

The bond exclusions are a curated list rather than a connectivity graph, so they exempt
bonded atoms but not their neighbours. Around an N-linked glycan the `C1`–`ND2` bond is
excluded correctly, while `ND2`···`O5` at 2.38 Å and `CG`···`C1` at 2.44 Å — two bonds
away, and held there by the geometry of the bond itself — are reported as `polar`. Six of
the eight `polar` contacts in 1IGT arise this way. The `steric` and `tight` grades are not
affected.

Setting `POLAR` at 2.6 Å also catches unusually short but genuine hydrogen bonds, which
surface as `polar` between 2.5 and 2.6 Å.

### Not implemented

- hydrogen positions, and therefore any MolProbity style clashscore
- crystallographic symmetry mates
- metal coordination geometry, since metals are discarded with the solvent
- mmCIF input
- any modification or repair of the structure

## Validation

Four structures from the PDB, run as one batch:

| PDB | Structure | Chains | Atoms | steric | tight | polar | carbon |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1UBQ | ubiquitin, 1.8 Å | 1 | 602 | 0 | 0 | 0 | 0 |
| 1IGT | intact IgG, 2.8 Å | 6 | 10416 | 0 | 0 | 8 | 8 |
| 3HFM | lysozyme–Fab complex, 3.0 Å | 3 | 4296 | 0 | 0 | 6 | 12 |
| 6M0J | SARS-CoV-2 RBD–ACE2, 2.45 Å | 2 | 6476 | 0 | 0 | 0 | 1 |

No well refined structure reports a `steric` or `tight` contact, which is the property the
thresholds are chosen to have. 1UBQ reports nothing at all because it is a monomer, which
is the limitation above seen from the other side.

Switching the bond exclusions off raises two `steric` and five `tight` contacts on 1IGT
alone — a structure carrying seventeen disulfides and glycans on two chains. That
difference is the measurement that justifies the exclusion list.

The whole batch, 21790 atoms, takes 0.11 s.

## Development

```bash
git clone https://github.com/LamarckLab/pdbcheck.git
cd pdbcheck
pip install -e ".[dev]"
pytest
ruff check .
```

## License

MIT. See [LICENSE](https://github.com/LamarckLab/pdbcheck/blob/main/LICENSE).
