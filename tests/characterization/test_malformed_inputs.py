"""Parsers and metrics against malformed sidecar files.

Every branch here fires on input that is present but unusable: a truncated
npz, a PAE matrix of the wrong shape, a chain that aligns to nothing. The
contract is that the run continues with the affected value defaulted or
blanked, never with a fabricated number.
"""

import json
import sys

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit

CHAINS = [4, 3]


def _pdb(path, chains=(("A", 4), ("B", 3)), resname="ALA", b=60.0):
    lines, serial = [], 1
    for chain, n in chains:
        for i in range(n):
            lines.append(
                f"ATOM  {serial:>5}  CA  {resname} {chain}{i + 1:>4}    "
                f"{i * 3.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{b:6.2f}           C")
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _boltz(root, **files):
    d = root / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb")
    for name, writer in files.items():
        writer(d / name)
    return root


# ── Corrupt numpy sidecars ────────────────────────────────────────────────

def test_a_truncated_plddt_npz_falls_back_to_b_factors(tmp_path, capsys):
    root = _boltz(tmp_path,
                  **{"plddt_job_model_0.npz": lambda p: p.write_bytes(b"junk")})
    r = cm.parse_boltz2(str(root), CHAINS)
    assert r[0]["avg_plddt"] == pytest.approx(60.0)
    assert "plddt parse failed" in capsys.readouterr().out


def test_a_truncated_pae_npz_leaves_pae_metrics_at_zero(tmp_path, capsys):
    root = _boltz(tmp_path,
                  **{"pae_job_model_0.npz": lambda p: p.write_bytes(b"junk")})
    r = cm.parse_boltz2(str(root), CHAINS)
    assert r[0]["pae_mean"] == 0.0
    assert "pae parse failed" in capsys.readouterr().out


def test_a_three_dimensional_pae_array_takes_the_first_slice(tmp_path, capsys):
    """Some builds emit a batch dimension."""
    root = _boltz(tmp_path, **{"pae_job_model_0.npz":
                               lambda p: np.savez(p, pae=np.full((1, 7, 7), 3.0))})
    r = cm.parse_boltz2(str(root), CHAINS)
    assert r[0]["pae_mean"] == pytest.approx(3.0)


def test_a_multidimensional_plddt_array_is_flattened(tmp_path, capsys):
    root = _boltz(tmp_path, **{"plddt_job_model_0.npz":
                               lambda p: np.savez(p, plddt=np.full((1, 7), 0.75))})
    r = cm.parse_boltz2(str(root), CHAINS)
    assert r[0]["avg_plddt"] == pytest.approx(0.75)


# ── RMSD guards on unalignable input ──────────────────────────────────────

def test_zero_length_receptor_gives_no_rmsd(tmp_path, capsys):
    """A prediction with an effector but no receptor cannot be superposed."""
    ref = _pdb(tmp_path / "ref.pdb")
    pred = _pdb(tmp_path / "pred.pdb", chains=(("B", 3),))
    m = cm.compute_structural_rmsds(str(pred), str(ref))
    assert m["rmsd_receptor"] is None


def test_a_prediction_missing_the_effector_gives_no_rmsd(tmp_path, capsys):
    ref = _pdb(tmp_path / "ref.pdb")
    pred = _pdb(tmp_path / "pred.pdb", chains=(("A", 4),))
    m = cm.compute_structural_rmsds(str(pred), str(ref))
    assert m["rmsd_effector_receptor_aligned"] is None


# ── ESMFold2 confidences ──────────────────────────────────────────────────

def _esm(root, conf):
    d = root / "seed42"
    d.mkdir(parents=True)
    import gemmi
    st = gemmi.read_structure(str(_pdb(root / "src.pdb")))
    st.setup_entities()
    (d / "esmfold2_pred.cif").write_text(st.make_mmcif_document().as_string())
    (d / "confidences.json").write_text(conf)
    return root


def test_a_corrupt_esmfold2_confidences_file_is_survived(tmp_path, capsys):
    root = _esm(tmp_path, "{broken")
    r = cm.parse_esmfold2(str(root), CHAINS)
    assert len(r) == 1
    assert "confidences parse failed" in capsys.readouterr().out


def test_an_oversized_esmfold2_pae_is_trimmed(tmp_path, capsys):
    """A PAE bigger than the token count is cropped rather than discarded."""
    root = _esm(tmp_path, json.dumps({
        "plddt": [80.0] * 7, "ptm": 0.7,
        "pae": np.full((9, 9), 4.0).tolist()}))
    r = cm.parse_esmfold2(str(root), CHAINS)
    assert r[0]["pae_mean"] == pytest.approx(4.0)


def test_an_undersized_esmfold2_pae_is_discarded(tmp_path, capsys):
    """Too small means the mapping to chains is unknown, so guessing would
    fabricate an interface score."""
    root = _esm(tmp_path, json.dumps({
        "plddt": [80.0] * 7, "ptm": 0.7,
        "pae": np.full((3, 3), 4.0).tolist()}))
    r = cm.parse_esmfold2(str(root), CHAINS)
    assert r[0]["ipsae_min"] == 0.0


def test_the_parser_passes_plddt_through_without_rescaling(tmp_path, capsys):
    """Rescaling to 0-100 is the writer's job, in esmfold2_fold. The parser
    reports what the file holds, and the analysis rescales across models."""
    root = _esm(tmp_path, json.dumps({"plddt": [0.8] * 7, "ptm": 0.7}))
    r = cm.parse_esmfold2(str(root), CHAINS)
    assert r[0]["avg_plddt"] == pytest.approx(0.8)


# ── B-factor reader ───────────────────────────────────────────────────────

def test_unparseable_b_factor_columns_are_skipped(tmp_path):
    """A malformed column must not abort the whole file."""
    p = tmp_path / "odd.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  ABCDE           C\n"
        "ATOM      2  CA  ALA A   2       3.000   0.000   0.000  1.00 70.00           C\n"
        "END\n")
    assert cm.plddt_from_pdb(p) == pytest.approx(70.0)


# ── main() reference handling ─────────────────────────────────────────────

def test_a_reference_that_cannot_be_read_leaves_rmsds_blank(tmp_path, capsys):
    root = _boltz(tmp_path / "pred")
    bad_ref = tmp_path / "bad.pdb"
    bad_ref.write_text("not a pdb at all\n")
    out = tmp_path / "m.csv"
    sys.argv = ["compute_metrics.py", "--model", "boltz2",
                "--prediction-dir", str(root), "--chain-lengths", "4", "3",
                "--output-csv", str(out), "--reference-pdb", str(bad_ref)]
    cm.main()
    import csv
    rows = list(csv.DictReader(open(out)))
    assert rows[0]["rmsd_receptor"] == ""
