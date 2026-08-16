"""Defensive branches in the metric parsers.

Predictor outputs vary between versions and occasionally arrive truncated. The
parsers are written to degrade to a defaulted value rather than abort a
benchmark part-way through, and these pin that behaviour so a refactor cannot
turn a warning into a crash, or worse into a silent zero that reads as a real
measurement.
"""

import json
import pickle

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit

CHAINS = [4, 3]


def _pdb(path, bfactors=(90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0)):
    lines, serial = [], 1
    for chain, n, start in (("A", 4, 0), ("B", 3, 4)):
        for i in range(n):
            lines.append(
                f"ATOM  {serial:>5}  CA  ALA {chain}{i + 1:>4}    "
                f"{i * 3.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}"
                f"{bfactors[start + i]:6.2f}           C")
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _cif(pdb_path, cif_path):
    import gemmi
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()
    cif_path.write_text(st.make_mmcif_document().as_string())
    return cif_path


# ── ipSAE guards ──────────────────────────────────────────────────────────

def test_a_pae_matrix_of_the_wrong_size_is_rescued_when_possible(capsys):
    """Chain lengths can disagree with the emitted matrix. If the receptor
    length still fits, the effector length is inferred rather than discarding
    the matrix."""
    out = cm.compute_ipsae(np.zeros((10, 10)), [4, 3])
    assert out["ipsae_min"] == 1.0
    assert "!=" in capsys.readouterr().out


def test_a_pae_matrix_smaller_than_the_receptor_is_discarded(capsys):
    out = cm.compute_ipsae(np.zeros((2, 2)), [4, 3])
    assert out == {"ipsae_ab": 0.0, "ipsae_ba": 0.0, "ipsae_min": 0.0}


def test_ipae_returns_zero_on_a_size_mismatch():
    assert cm.compute_ipae(np.zeros((5, 5)), [4, 3]) == 0.0


def test_ipae_returns_zero_without_a_matrix():
    assert cm.compute_ipae(None, [4, 3]) == 0.0


# ── pLDDT readers ─────────────────────────────────────────────────────────

def test_plddt_from_pdb_returns_zero_for_a_missing_file():
    assert cm.plddt_from_pdb("/nonexistent/x.pdb") == 0.0


def test_plddt_from_pdb_returns_zero_when_no_ca_atoms_exist(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("END\n")
    assert cm.plddt_from_pdb(p) == 0.0


def test_plddt_from_cif_returns_zero_for_a_missing_file():
    assert cm.plddt_from_cif("/nonexistent/x.cif") == 0.0


def test_plddt_from_cif_returns_zero_for_an_unparseable_file(tmp_path):
    p = tmp_path / "bad.cif"
    p.write_text("this is not mmCIF")
    assert cm.plddt_from_cif(p) == 0.0


# ── RMSD guards ───────────────────────────────────────────────────────────

def test_rmsd_is_nan_when_chains_cannot_be_resolved(tmp_path, capsys):
    """A prediction of an entirely different protein must not be silently
    matched to the reference."""
    ref = _pdb(tmp_path / "ref.pdb")
    other = tmp_path / "other.pdb"
    lines = []
    for i in range(6):
        lines.append(f"ATOM  {i + 1:>5}  CA  TRP C{i + 1:>4}    "
                     f"{i * 3.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           C")
    other.write_text("\n".join(lines) + "\nEND\n")

    m = cm.compute_structural_rmsds(str(other), str(ref))
    assert m["rmsd_effector_receptor_aligned"] is None
    assert "Chain resolution failed" in capsys.readouterr().out


def test_rmsd_is_none_when_the_prediction_has_no_ca(tmp_path, capsys):
    ref = _pdb(tmp_path / "ref.pdb")
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    m = cm.compute_structural_rmsds(str(empty), str(ref))
    assert m["rmsd_receptor"] is None


def test_rmsd_is_none_when_the_reference_has_no_ca(tmp_path, capsys):
    pred = _pdb(tmp_path / "pred.pdb")
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    m = cm.compute_structural_rmsds(str(pred), str(empty))
    assert m["rmsd_receptor"] is None


# ── AF2M pickles ──────────────────────────────────────────────────────────

def test_af2m_reads_metrics_from_the_result_pickles(tmp_path, capsys):
    d = tmp_path / "input"
    d.mkdir()
    _pdb(d / "ranked_0.pdb")
    with open(d / "result_model_1_multimer_v3_pred_0.pkl", "wb") as fh:
        pickle.dump({"ptm": np.float32(0.72), "iptm": np.float32(0.61),
                     "plddt": np.full(7, 85.0),
                     "predicted_aligned_error": np.full((7, 7), 3.0)}, fh)
    (d / "ranking_debug.json").write_text(json.dumps(
        {"order": ["model_1_multimer_v3_pred_0"],
         "iptm+ptm": {"model_1_multimer_v3_pred_0": 0.63}}))

    r = cm.parse_af2m(str(tmp_path), CHAINS)
    assert r[0]["ptm"] == pytest.approx(0.72, abs=1e-4)
    assert r[0]["iptm"] == pytest.approx(0.61, abs=1e-4)
    assert r[0]["avg_plddt"] == pytest.approx(85.0)
    assert r[0]["pae_mean"] == pytest.approx(3.0)


def test_af2m_survives_a_corrupt_pickle(tmp_path, capsys):
    d = tmp_path / "input"
    d.mkdir()
    _pdb(d / "ranked_0.pdb")
    (d / "result_model_1_multimer_v3_pred_0.pkl").write_bytes(b"not a pickle")
    r = cm.parse_af2m(str(tmp_path), CHAINS)
    assert len(r) == 1                      # falls back to B-factors
    assert r[0]["avg_plddt"] == pytest.approx(60.0)


def test_af2m_survives_a_corrupt_ranking_debug(tmp_path, capsys):
    d = tmp_path / "input"
    d.mkdir()
    _pdb(d / "ranked_0.pdb")
    (d / "ranking_debug.json").write_text("{broken")
    assert len(cm.parse_af2m(str(tmp_path), CHAINS)) == 1


# ── AF3 layouts ───────────────────────────────────────────────────────────

def test_af3_reads_the_aggregate_model_alongside_the_seeds(tmp_path, capsys):
    """AF3 writes a top-level aggregate as well as per-seed samples. Both are
    parsed, and the aggregate is tagged so it can be dropped downstream."""
    job = tmp_path / "job"
    job.mkdir()
    _cif(_pdb(tmp_path / "src.pdb"), job / "job_model.cif")
    (job / "job_confidences.json").write_text(json.dumps(
        {"atom_plddts": [90.0] * 7, "pae": np.full((7, 7), 2.0).tolist()}))
    (job / "job_summary_confidences.json").write_text(
        json.dumps({"ptm": 0.8, "iptm": 0.7}))

    r = cm.parse_af3(str(tmp_path), CHAINS)
    assert any("__aggregate" in e["model_name"] for e in r)


def test_af3_handles_a_seed_directory_without_confidences(tmp_path, capsys):
    seed = tmp_path / "job" / "seed-42_sample-0"
    seed.mkdir(parents=True)
    _cif(_pdb(tmp_path / "src.pdb"), seed / "model.cif")
    r = cm.parse_af3(str(tmp_path), CHAINS)
    assert len(r) >= 1
    assert r[0]["ptm"] == 0.0               # defaulted, not invented


def test_af3_confidences_accepts_either_plddt_key(tmp_path):
    """Builds differ between atom_plddts and plddt."""
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"plddt": [80.0] * 7}))
    avg, pae = cm._parse_af3_confidences(str(p), CHAINS)
    assert avg == pytest.approx(80.0)
    assert pae is None


def test_af3_confidences_returns_defaults_for_a_corrupt_file(tmp_path):
    p = tmp_path / "c.json"
    p.write_text("{broken")
    assert cm._parse_af3_confidences(str(p), CHAINS) == (0.0, None)


# ── ColabFold ─────────────────────────────────────────────────────────────

def test_colabfold_falls_back_to_any_scores_file(tmp_path, capsys):
    """Naming varies between ColabFold versions, so an exact match is tried
    first and then any scores file in the same directory."""
    d = tmp_path / "output"
    d.mkdir()
    _pdb(d / "cmplx_unrelaxed_rank_001_model_1_seed_000.pdb")
    (d / "oddly_named_scores.json").write_text(json.dumps(
        {"ptm": 0.6, "iptm": 0.5, "plddt": [70.0] * 7}))
    r = cm.parse_colabfold(str(tmp_path), CHAINS)
    assert r[0]["ptm"] == pytest.approx(0.6)


def test_colabfold_survives_a_corrupt_scores_file(tmp_path, capsys):
    d = tmp_path / "output"
    d.mkdir()
    _pdb(d / "cmplx_unrelaxed_rank_001_model_1_seed_000.pdb")
    (d / "cmplx_scores_rank_001_model_1_seed_000.json").write_text("{broken")
    r = cm.parse_colabfold(str(tmp_path), CHAINS)
    assert len(r) == 1
    assert r[0]["avg_plddt"] == pytest.approx(60.0)


# ── ESMFold2 ──────────────────────────────────────────────────────────────

def test_esmfold2_falls_back_to_cif_b_factors(tmp_path, capsys):
    d = tmp_path / "all_outputs" / "seed42"
    d.mkdir(parents=True)
    _cif(_pdb(tmp_path / "src.pdb"), d / "esmfold2_pred.cif")
    r = cm.parse_esmfold2(str(d.parent), CHAINS)
    assert len(r) == 1
    assert r[0]["avg_plddt"] == pytest.approx(60.0)


def test_esmfold2_returns_nothing_for_an_empty_tree(tmp_path, capsys):
    assert cm.parse_esmfold2(str(tmp_path), CHAINS) == []


# ── Reference detection ───────────────────────────────────────────────────

def test_the_reference_pdb_is_not_parsed_as_a_prediction(tmp_path, capsys):
    """The published reference sits alongside the outputs and would otherwise
    be scored as a perfect prediction of itself."""
    d = tmp_path / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb")
    _pdb(tmp_path / "reference.pdb")
    r = cm.parse_boltz2(str(tmp_path), CHAINS)
    assert [e["model_name"] for e in r] == ["job_model_0"]
