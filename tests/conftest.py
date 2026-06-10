"""Shared pytest fixtures and helpers for the benchmark test suite.

Makes the repo's bin/ directory importable (its scripts are standalone, not
packaged) and provides a small synthetic two-chain PDB builder so the
geometry/constraint tests run on known coordinates without shipping a real
structure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))


def _atom_line(serial, name, resname, chain, resseq, x, y, z,
               occ=1.0, bfac=0.0, element="C"):
    """Build one fixed-column PDB ATOM record (80 chars).

    Only the columns the pipeline parsers read are load-bearing — atom name
    (13-16), chain (22), resSeq (23-26), x/y/z (31-54), and tempFactor
    (61-66, used as pLDDT) — but the whole record is laid out to spec.
    """
    line = [" "] * 80

    def put(text, start):  # start is the 1-based column
        i = start - 1
        line[i:i + len(text)] = list(text)

    put("ATOM", 1)
    put(f"{serial:>5}", 7)
    put(name, 14)            # ≤3-char names start at col 14 → " CA " in 13-16
    put(f"{resname:>3}", 18)
    put(chain, 22)
    put(f"{resseq:>4}", 23)
    put(f"{x:8.3f}", 31)
    put(f"{y:8.3f}", 39)
    put(f"{z:8.3f}", 47)
    put(f"{occ:6.2f}", 55)
    put(f"{bfac:6.2f}", 61)
    put(f"{element:>2}", 77)
    return "".join(line).rstrip()


def make_pdb(atoms):
    """Render a PDB string from a list of atom dicts.

    Each atom dict needs at least ``chain``, ``resseq``, ``x``, ``y``, ``z``;
    ``name`` defaults to ``CA``, ``resname`` to ``ALA``, ``bfac`` to 0.0.
    """
    lines = []
    for i, a in enumerate(atoms, start=1):
        lines.append(_atom_line(
            serial=i,
            name=a.get("name", "CA"),
            resname=a.get("resname", "ALA"),
            chain=a["chain"],
            resseq=a["resseq"],
            x=a["x"], y=a["y"], z=a["z"],
            bfac=a.get("bfac", 0.0),
        ))
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.fixture
def make_pdb_file(tmp_path):
    """Factory: write a synthetic PDB to a temp file and return its path."""
    def _write(atoms, name="synthetic.pdb"):
        p = tmp_path / name
        p.write_text(make_pdb(atoms))
        return p
    return _write


@pytest.fixture
def two_chain_atoms():
    """A deterministic two-chain layout (receptor B, effector C) with known
    inter-chain distances:

        B/1 CA (0,0,0)    C/1 CA (3,0,0)   -> 3.0 Å
        B/2 CA (10,0,0)                     -> B2-C1 = 7.0 Å
        B/3 CA (0,10,0)                     -> B3-C1 = sqrt(109) ≈ 10.44 Å

    So with an 8 Å cutoff: pocket = receptor residues {1, 2}; contact pairs
    (sorted) = [(1,1,3.0), (2,1,7.0)].
    """
    return [
        {"chain": "B", "resseq": 1, "x": 0.0,  "y": 0.0,  "z": 0.0, "bfac": 90.0},
        {"chain": "B", "resseq": 2, "x": 10.0, "y": 0.0,  "z": 0.0, "bfac": 80.0},
        {"chain": "B", "resseq": 3, "x": 0.0,  "y": 10.0, "z": 0.0, "bfac": 70.0},
        {"chain": "C", "resseq": 1, "x": 3.0,  "y": 0.0,  "z": 0.0, "bfac": 60.0},
    ]
