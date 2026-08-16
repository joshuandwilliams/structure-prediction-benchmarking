"""Unit tests for the structural-RMSD path in compute_metrics.py.

This is the metric the whole benchmark is scored on
(``rmsd_effector_receptor_aligned``, "ra_eff"), so the arithmetic is pinned
here on synthetic structures with known answers rather than trusted.

The geometry is built so each test isolates one behaviour:

* a rigid-body copy of the reference must give ~0 everywhere — a non-zero
  answer would mean the Kabsch fit is not removing rotation/translation;
* displacing ONLY the effector must leave the receptor fit at ~0 while
  ra_eff reports the displacement, which is the exact discrimination the
  benchmark relies on;
* the effector's own independent fit must stay ~0 under that displacement,
  since a rigidly-moved chain still superposes perfectly onto itself. That
  is what separates "wrongly placed" from "wrongly folded".
"""

import math

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit


# ── Geometry helpers ──────────────────────────────────────────────────────

def _rot(ax, ay, az):
    """Rotation matrix from three axis angles (radians)."""
    cx, sx, cy, sy, cz, sz = (math.cos(ax), math.sin(ax), math.cos(ay),
                              math.sin(ay), math.cos(az), math.sin(az))
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return rz @ ry @ rx


def _chain_atoms(chain, resname, coords, start=1):
    return [{"chain": chain, "resseq": start + i, "resname": resname,
             "x": float(x), "y": float(y), "z": float(z), "bfac": 90.0}
            for i, (x, y, z) in enumerate(coords)]


def _two_chain_coords(n_rec=30, n_eff=24):
    """A receptor helix and a separate effector helix, well apart.

    Chains must be long enough to clear the iterative fit's min_pairs=20
    floor, and distinct enough in sequence that greedy chain resolution
    cannot swap them (receptor ALA, effector GLY).
    """
    t = np.arange(n_rec, dtype=float)
    rec = np.column_stack([np.cos(t / 2) * 5, np.sin(t / 2) * 5, t * 1.5])
    u = np.arange(n_eff, dtype=float)
    eff = np.column_stack([np.cos(u / 2) * 4 + 30, np.sin(u / 2) * 4, u * 1.5])
    return rec, eff


@pytest.fixture
def ref_and_coords(make_pdb_file):
    rec, eff = _two_chain_coords()
    ref = make_pdb_file(_chain_atoms("A", "ALA", rec) + _chain_atoms("B", "GLY", eff),
                        name="ref.pdb")
    return ref, rec, eff


# ── Kabsch core ───────────────────────────────────────────────────────────

def test_kabsch_recovers_a_rigid_transform_exactly():
    rec, _ = _two_chain_coords()
    R_true = _rot(0.3, -0.7, 1.1)
    moved = (R_true @ rec.T).T + np.array([12.0, -4.0, 7.5])
    rmsd, R, t, n = cm._kabsch_align(moved, rec)
    assert n == len(rec)
    assert rmsd == pytest.approx(0.0, abs=1e-6)
    # Applying the returned transform must land the moved copy back on rec.
    back = (R @ moved.T).T + t
    assert np.allclose(back, rec, atol=1e-6)


def test_kabsch_reports_known_uniform_displacement():
    rec, _ = _two_chain_coords()
    # A pure translation cannot be fitted away only if it is applied to a
    # subset; applied to all points, Kabsch removes it entirely.
    rmsd, _, _, _ = cm._kabsch_align(rec + np.array([3.0, 0.0, 0.0]), rec)
    assert rmsd == pytest.approx(0.0, abs=1e-6)


def test_kabsch_is_symmetric_in_rmsd():
    rec, eff = _two_chain_coords(n_rec=20, n_eff=20)
    a, _, _, _ = cm._kabsch_align(rec, eff)
    b, _, _, _ = cm._kabsch_align(eff, rec)
    assert a == pytest.approx(b, abs=1e-9)


def test_iterative_kabsch_prunes_a_single_outlier():
    rec, _ = _two_chain_coords()
    moved = rec.copy()
    moved[0] += np.array([25.0, 0.0, 0.0])          # one wildly wrong residue
    whole, _, _, _ = cm._kabsch_align(moved, rec)
    pruned = cm._kabsch_align_iterative(moved, rec, cutoff=2.0, min_pairs=20)
    assert whole > 1.0                               # outlier drags the whole fit
    assert pruned[0] < whole                         # pruning recovers the core


# ── Sequence-alignment pairing ────────────────────────────────────────────

def test_seqalign_pairs_identical_sequences_one_to_one():
    pred, ref = cm._seqalign_pair_indices("ACDEFGHIK", "ACDEFGHIK")
    assert pred == list(range(9))
    assert ref == list(range(9))


def test_seqalign_skips_residues_missing_from_the_reference():
    # The prediction models two extra N-terminal residues the crystal did not.
    # Naive file-order pairing would offset every Ca by two.
    pred, ref = cm._seqalign_pair_indices("MGACDEFGHIK", "ACDEFGHIK")
    assert len(pred) == len(ref) == 9
    assert pred == [2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert ref == list(range(9))


def test_seqalign_returns_empty_on_empty_input():
    assert cm._seqalign_pair_indices("", "ACDE") == ([], [])
    assert cm._seqalign_pair_indices("ACDE", "") == ([], [])


# ── End-to-end structural metrics ─────────────────────────────────────────

def test_rigid_copy_of_reference_scores_zero_everywhere(ref_and_coords, make_pdb_file):
    ref, rec, eff = ref_and_coords
    R = _rot(0.9, 0.2, -0.4)
    shift = np.array([-15.0, 22.0, 3.0])
    pred = make_pdb_file(
        _chain_atoms("A", "ALA", (R @ rec.T).T + shift) +
        _chain_atoms("B", "GLY", (R @ eff.T).T + shift), name="pred_rigid.pdb")

    m = cm.compute_structural_rmsds(str(pred), str(ref), rec_chain="A", eff_chain="B")
    assert m["rmsd_receptor"] == pytest.approx(0.0, abs=1e-4)
    assert m["rmsd_effector_independent"] == pytest.approx(0.0, abs=1e-4)
    assert m["rmsd_effector_receptor_aligned"] == pytest.approx(0.0, abs=1e-4)
    assert m["n_receptor_ca"] == len(rec)
    assert m["n_effector_ca"] == len(eff)


def test_displaced_effector_moves_ra_eff_but_not_the_receptor_fit(
        ref_and_coords, make_pdb_file):
    """The discrimination the benchmark depends on.

    The effector is translated 9 A as a rigid body. The receptor is untouched,
    so its fit stays ~0; the effector still superposes onto itself, so its
    independent fit stays ~0; only ra_eff — the effector under the receptor's
    transform — sees the 9 A. A metric that conflated placement with folding
    would move all three.
    """
    ref, rec, eff = ref_and_coords
    offset = np.array([9.0, 0.0, 0.0])
    pred = make_pdb_file(
        _chain_atoms("A", "ALA", rec) + _chain_atoms("B", "GLY", eff + offset),
        name="pred_shifted.pdb")

    m = cm.compute_structural_rmsds(str(pred), str(ref), rec_chain="A", eff_chain="B")
    assert m["rmsd_receptor"] == pytest.approx(0.0, abs=1e-4)
    assert m["rmsd_effector_independent"] == pytest.approx(0.0, abs=1e-4)
    assert m["rmsd_effector_receptor_aligned"] == pytest.approx(9.0, abs=1e-3)


def test_ra_eff_scales_with_the_size_of_the_displacement(
        ref_and_coords, make_pdb_file):
    ref, rec, eff = ref_and_coords
    seen = []
    for d in (2.0, 5.0, 20.0):
        pred = make_pdb_file(
            _chain_atoms("A", "ALA", rec) +
            _chain_atoms("B", "GLY", eff + np.array([d, 0.0, 0.0])),
            name=f"pred_{d}.pdb")
        m = cm.compute_structural_rmsds(str(pred), str(ref), rec_chain="A", eff_chain="B")
        seen.append(m["rmsd_effector_receptor_aligned"])
    assert seen == sorted(seen)
    assert seen[0] == pytest.approx(2.0, abs=1e-3)
    assert seen[2] == pytest.approx(20.0, abs=1e-3)


def test_chain_ids_are_resolved_by_sequence_not_by_label(
        ref_and_coords, make_pdb_file):
    """Boltz writes A/B regardless of the reference's labels, and AF2M wrote
    B/C on four benchmark targets. Resolution must follow sequence identity,
    so swapping the labels in the prediction must not swap the metrics."""
    ref, rec, eff = ref_and_coords
    pred = make_pdb_file(
        _chain_atoms("B", "ALA", rec) + _chain_atoms("C", "GLY", eff),
        name="pred_relabelled.pdb")

    m = cm.compute_structural_rmsds(str(pred), str(ref), rec_chain="A", eff_chain="B")
    assert m["pred_receptor_chain"] == "B"
    assert m["pred_effector_chain"] == "C"
    assert m["rmsd_receptor"] == pytest.approx(0.0, abs=1e-4)
    assert m["rmsd_effector_receptor_aligned"] == pytest.approx(0.0, abs=1e-4)


def test_missing_prediction_file_does_not_raise(ref_and_coords):
    ref, _, _ = ref_and_coords
    m = cm.compute_structural_rmsds("/nonexistent/pred.pdb", str(ref))
    assert m["rmsd_effector_receptor_aligned"] is None
