"""Unit tests for the per-model output parsers in compute_metrics.py.

Each predictor writes a different tree, and these parsers are what turn that
into the numbers the benchmark reports. A parser that silently finds nothing
yields an empty result rather than an error, so every test asserts on the
values recovered rather than only that the call succeeded.

Fixtures are minimal but real: a two-chain PDB or mmCIF plus whatever
confidence sidecar that predictor emits.
"""

import json

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit

CHAINS = [4, 3]          # receptor, effector residue counts


# ── Fixture builders ──────────────────────────────────────────────────────

def _pdb(path, bfactors=(90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0)):
    """Two-chain Cα-only PDB, 4 residues in A and 3 in B."""
    lines, serial = [], 1
    for chain, n, start in (("A", 4, 0), ("B", 3, 4)):
        for i in range(n):
            b = bfactors[start + i]
            lines.append(
                f"ATOM  {serial:>5}  CA  ALA {chain}{i + 1:>4}    "
                f"{i * 3.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{b:6.2f}           C")
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _pae(path, n=7, value=5.0):
    np.savez(path, pae=np.full((n, n), value, dtype=float))


def _plddt(path, values=(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)):
    np.savez(path, plddt=np.asarray(values, dtype=float))


# ── Boltz ─────────────────────────────────────────────────────────────────

@pytest.fixture
def boltz_dir(tmp_path):
    d = tmp_path / "all_outputs" / "seed42" / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb")
    (d / "confidence_job_model_0.json").write_text(
        json.dumps({"ptm": 0.81, "iptm": 0.72}))
    _plddt(d / "plddt_job_model_0.npz")
    _pae(d / "pae_job_model_0.npz")
    return tmp_path


def test_boltz2_reads_confidence_plddt_and_pae(boltz_dir):
    r = cm.parse_boltz2(str(boltz_dir), CHAINS)
    assert len(r) == 1
    e = r[0]
    assert e["ptm"] == pytest.approx(0.81)
    assert e["iptm"] == pytest.approx(0.72)
    assert e["avg_plddt"] == pytest.approx(0.6, abs=1e-6)   # mean of the npz
    assert e["pae_mean"] == pytest.approx(5.0)
    assert e["ipsae_min"] > 0                                # PAE was readable


def test_boltz1_shares_the_boltz2_parser(boltz_dir):
    assert cm.parse_boltz1(str(boltz_dir), CHAINS) == cm.parse_boltz2(
        str(boltz_dir), CHAINS)


def test_boltz_falls_back_to_b_factors_when_plddt_npz_is_absent(tmp_path):
    d = tmp_path / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb")
    r = cm.parse_boltz2(str(tmp_path), CHAINS)
    assert r[0]["avg_plddt"] == pytest.approx(60.0)   # mean of the B-factors


def test_boltz_returns_nothing_for_an_empty_directory(tmp_path):
    assert cm.parse_boltz2(str(tmp_path), CHAINS) == []


def test_boltz_survives_a_corrupt_confidence_file(tmp_path):
    d = tmp_path / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb")
    (d / "confidence_job_model_0.json").write_text("{not json")
    r = cm.parse_boltz2(str(tmp_path), CHAINS)
    assert len(r) == 1 and r[0]["ptm"] == 0.0        # defaulted, not crashed


# ── Chai-1 ────────────────────────────────────────────────────────────────

def test_chai1_reads_scores_npz(tmp_path):
    d = tmp_path / "all_outputs" / "seed42"
    d.mkdir(parents=True)
    _pdb(d / "pred.model_idx_0.pdb")
    np.savez(d / "scores.model_idx_0.npz",
             ptm=np.array([0.77]), iptm=np.array([0.66]))
    r = cm.parse_chai1(str(tmp_path), CHAINS)
    assert len(r) == 1
    assert r[0]["ptm"] == pytest.approx(0.77)
    assert r[0]["iptm"] == pytest.approx(0.66)


def test_chai1_has_no_pae_so_ipsae_stays_zero(tmp_path):
    """Chai-1 emits no PAE matrix, so PAE-derived metrics must be 0 rather
    than a fabricated value. The analysis blanks these afterwards."""
    d = tmp_path / "out"
    d.mkdir(parents=True)
    _pdb(d / "pred.model_idx_0.pdb")
    r = cm.parse_chai1(str(tmp_path), CHAINS)
    assert r[0]["ipsae_min"] == 0.0
    assert r[0]["pae_mean"] == 0.0


# ── AlphaFold2-Multimer ───────────────────────────────────────────────────

def test_af2m_reads_ranked_pdbs_and_ranking_debug(tmp_path):
    d = tmp_path / "input"
    d.mkdir(parents=True)
    _pdb(d / "ranked_0.pdb")
    (d / "ranking_debug.json").write_text(json.dumps(
        {"order": ["model_1_multimer_v3_pred_0"],
         "iptm+ptm": {"model_1_multimer_v3_pred_0": 0.83}}))
    r = cm.parse_af2m(str(tmp_path), CHAINS)
    assert len(r) == 1
    assert r[0]["avg_plddt"] == pytest.approx(60.0)   # from B-factors


def test_af2m_returns_nothing_when_no_ranked_pdbs_exist(tmp_path):
    (tmp_path / "input").mkdir()
    assert cm.parse_af2m(str(tmp_path), CHAINS) == []


# ── AlphaFold 3 ───────────────────────────────────────────────────────────

def _cif(pdb_path, cif_path):
    """Convert the fixture PDB to mmCIF, which AF3 and ESMFold2 emit."""
    import gemmi
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()
    cif_path.write_text(st.make_mmcif_document().as_string())
    return cif_path


def _af3_seed(dirpath, plddt=(90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0)):
    dirpath.mkdir(parents=True)
    _cif(_pdb(dirpath / "src.pdb"), dirpath / "model.cif")
    (dirpath / "src.pdb").unlink()
    (dirpath / "confidences.json").write_text(json.dumps({
        "atom_plddts": list(plddt),
        "pae": np.full((7, 7), 4.0).tolist()}))
    (dirpath / "summary_confidences.json").write_text(
        json.dumps({"ptm": 0.79, "iptm": 0.69}))


def test_af3_reads_seed_directories(tmp_path):
    _af3_seed(tmp_path / "job" / "seed-42_sample-0")
    r = cm.parse_af3(str(tmp_path), CHAINS)
    assert len(r) >= 1
    e = r[0]
    assert e["ptm"] == pytest.approx(0.79)
    assert e["iptm"] == pytest.approx(0.69)
    assert e["avg_plddt"] == pytest.approx(60.0)
    assert e["pae_mean"] == pytest.approx(4.0)


def test_af3_returns_nothing_for_an_empty_tree(tmp_path):
    assert cm.parse_af3(str(tmp_path), CHAINS) == []


# ── ColabFold ─────────────────────────────────────────────────────────────

def test_colabfold_reads_its_scores_json(tmp_path):
    d = tmp_path / "output"
    d.mkdir(parents=True)
    _pdb(d / "cmplx_unrelaxed_rank_001_model_1_seed_000.pdb")
    (d / "cmplx_scores_rank_001_model_1_seed_000.json").write_text(json.dumps({
        "ptm": 0.75, "iptm": 0.65,
        "plddt": [90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0],
        "pae": np.full((7, 7), 6.0).tolist()}))
    r = cm.parse_colabfold(str(tmp_path), CHAINS)
    assert len(r) == 1
    assert r[0]["ptm"] == pytest.approx(0.75)
    assert r[0]["iptm"] == pytest.approx(0.65)
    assert r[0]["pae_mean"] == pytest.approx(6.0)


def test_colabfold_returns_nothing_for_an_empty_tree(tmp_path):
    assert cm.parse_colabfold(str(tmp_path), CHAINS) == []


# ── ESMFold2 ──────────────────────────────────────────────────────────────

def test_esmfold2_reads_cif_and_confidences(tmp_path):
    d = tmp_path / "all_outputs" / "seed42"
    d.mkdir(parents=True)
    _cif(_pdb(tmp_path / "src.pdb"), d / "esmfold2_pred.cif")
    (d / "confidences.json").write_text(json.dumps({
        "plddt": [90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0],
        "ptm": 0.7, "iptm": 0.6}))
    r = cm.parse_esmfold2(str(d.parent), CHAINS)
    assert len(r) == 1
    assert r[0]["ptm"] == pytest.approx(0.7)
    assert r[0]["avg_plddt"] == pytest.approx(60.0)


# ── The ranking score ─────────────────────────────────────────────────────

def test_ranking_score_is_the_alphafold_weighting(boltz_dir):
    """0.8*ipTM + 0.2*pTM. Named ranking_score rather than actifpTM, which
    needs a distogram no predictor here saves."""
    e = cm.parse_boltz2(str(boltz_dir), CHAINS)[0]
    assert e["ranking_score"] == pytest.approx(0.8 * 0.72 + 0.2 * 0.81, abs=1e-4)
