"""Unit tests for bin/_constraint_geometry.py and the two constraint-extractor
CLIs that build on it.

These run on a laptop with no dependencies beyond the standard library — the
geometry is pure Python over a synthetic two-chain PDB with known distances.
"""

import subprocess
import sys
from pathlib import Path

import _constraint_geometry as geom
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin"

pytestmark = pytest.mark.local_unit


# ── shared helpers ────────────────────────────────────────────────────────

def test_read_ca_by_chain_groups_by_chain(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    coords = geom.read_ca_by_chain(pdb)
    assert set(coords) == {"B", "C"}
    assert set(coords["B"]) == {1, 2, 3}
    assert coords["B"][1] == (0.0, 0.0, 0.0)
    assert coords["C"][1] == (3.0, 0.0, 0.0)


def test_read_ca_by_chain_filters_requested_chains(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    coords = geom.read_ca_by_chain(pdb, chains={"B"})
    assert set(coords) == {"B"}


def test_read_ca_by_chain_ignores_non_ca(make_pdb_file):
    atoms = [
        {"chain": "B", "resseq": 1, "x": 0.0, "y": 0.0, "z": 0.0, "name": "CB"},
        {"chain": "B", "resseq": 1, "x": 1.0, "y": 0.0, "z": 0.0, "name": "CA"},
    ]
    coords = geom.read_ca_by_chain(make_pdb_file(atoms))
    # Only the CA atom is read.
    assert coords["B"][1] == (1.0, 0.0, 0.0)


def test_distance():
    assert geom.distance((0, 0, 0), (3, 4, 0)) == pytest.approx(5.0)


def test_pocket_residues_within_cutoff(make_pdb_file, two_chain_atoms):
    coords = geom.read_ca_by_chain(make_pdb_file(two_chain_atoms))
    pocket = geom.pocket_residues(coords["B"], coords["C"], cutoff=8.0)
    assert pocket == [1, 2]  # B3 is 10.44 Å away → excluded


def test_pocket_residues_tight_cutoff(make_pdb_file, two_chain_atoms):
    coords = geom.read_ca_by_chain(make_pdb_file(two_chain_atoms))
    pocket = geom.pocket_residues(coords["B"], coords["C"], cutoff=5.0)
    assert pocket == [1]  # only B1 (3.0 Å); B2 is 7.0 Å


def test_pocket_residues_sorted(make_pdb_file):
    # Receptor residues deliberately out of order; all within cutoff.
    atoms = [
        {"chain": "B", "resseq": 5, "x": 0.0, "y": 0.0, "z": 0.0},
        {"chain": "B", "resseq": 2, "x": 1.0, "y": 0.0, "z": 0.0},
        {"chain": "C", "resseq": 1, "x": 0.5, "y": 0.0, "z": 0.0},
    ]
    coords = geom.read_ca_by_chain(make_pdb_file(atoms))
    assert geom.pocket_residues(coords["B"], coords["C"], cutoff=8.0) == [2, 5]


def test_contact_pairs_sorted_and_capped(make_pdb_file, two_chain_atoms):
    coords = geom.read_ca_by_chain(make_pdb_file(two_chain_atoms))
    pairs = geom.contact_pairs(coords["B"], coords["C"], cutoff=8.0, max_pairs=10)
    assert [(r, e) for r, e, _ in pairs] == [(1, 1), (2, 1)]
    assert pairs[0][2] == pytest.approx(3.0)
    assert pairs[1][2] == pytest.approx(7.0)


def test_contact_pairs_respects_cutoff(make_pdb_file, two_chain_atoms):
    coords = geom.read_ca_by_chain(make_pdb_file(two_chain_atoms))
    pairs = geom.contact_pairs(coords["B"], coords["C"], cutoff=5.0, max_pairs=10)
    assert [(r, e) for r, e, _ in pairs] == [(1, 1)]  # 7.0 Å pair dropped


def test_contact_pairs_max_cap(make_pdb_file, two_chain_atoms):
    coords = geom.read_ca_by_chain(make_pdb_file(two_chain_atoms))
    pairs = geom.contact_pairs(coords["B"], coords["C"], cutoff=8.0, max_pairs=1)
    assert len(pairs) == 1
    assert pairs[0][:2] == (1, 1)  # the closest survives the cap


def test_format_pocket_block():
    block = geom.format_pocket_block("B", "C", [1, 2], 6.0)
    assert block == (
        "  - pocket:\n"
        "      binder: C\n"
        "      contacts: [[B, 1], [B, 2]]\n"
        "      max_distance: 6.0\n"
        "      force: true"
    )


def test_format_contact_block():
    block = geom.format_contact_block("B", "C", 256, 78, 4.4)
    assert block == (
        "  - contact:\n"
        "      token1: [B, 256]\n"
        "      token2: [C, 78]\n"
        "      max_distance: 4.4\n"
        "      force: true"
    )


# ── extractor CLIs end-to-end ─────────────────────────────────────────────

def _run(script, *args, pdb):
    return subprocess.run(
        [sys.executable, str(BIN / script), str(pdb), *map(str, args)],
        capture_output=True, text=True,
    )


def test_boltz1_cli_emits_pocket_only(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    r = _run("extract_constraints_boltz1.py", "B", "C", 8.0, 6.0, pdb=pdb)
    assert r.returncode == 0, r.stderr
    assert "constraints:" in r.stdout
    assert "- pocket:" in r.stdout
    assert "- contact:" not in r.stdout          # Boltz-1 never emits contacts
    assert "max_distance: 6.0" in r.stdout        # pinned for Boltz-1
    assert "contacts: [[B, 1], [B, 2]]" in r.stdout


def test_boltz1_cli_errors_on_bad_chain(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    r = _run("extract_constraints_boltz1.py", "Z", "C", 8.0, 6.0, pdb=pdb)
    assert r.returncode == 1
    assert "No Cα atoms for chain Z" in r.stderr


def test_boltz2_cli_emits_pocket_and_contacts(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    r = _run("extract_constraints_boltz2.py", "B", "C", 8.0, 50, 0.0, 8.0, 8.0, pdb=pdb)
    assert r.returncode == 0, r.stderr
    assert "- pocket:" in r.stdout
    assert r.stdout.count("- contact:") == 2      # two pairs within 8 Å
    # closest pair B1-C1 at 3.0 Å, tolerance 0.0 → max_distance 3.0
    assert "token1: [B, 1]" in r.stdout
    assert "max_distance: 3.0" in r.stdout


def test_boltz2_cli_contact_cap(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    r = _run("extract_constraints_boltz2.py", "B", "C", 8.0, 1, 0.0, 8.0, 8.0, pdb=pdb)
    assert r.returncode == 0, r.stderr
    assert r.stdout.count("- contact:") == 1      # capped at 1
