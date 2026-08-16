"""End-to-end tests of the compute_metrics CLI.

This is what Nextflow invokes for every predictor, so its output CSV is the
sole input to everything downstream. These drive main() rather than the
parsers directly, covering the argument handling, the reference-free
representative selection, and the RMSD path when a reference is supplied.
"""

import csv
import json
import sys

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit

CHAINS = ["4", "3"]


def _pdb(path, coords_offset=0.0, bfactors=(90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0)):
    """Two-chain Cα PDB, receptor A (4 residues) and effector B (3)."""
    lines, serial = [], 1
    for chain, n, start in (("A", 4, 0), ("B", 3, 4)):
        for i in range(n):
            x = i * 3.0 + (coords_offset if chain == "B" else 0.0)
            lines.append(
                f"ATOM  {serial:>5}  CA  ALA {chain}{i + 1:>4}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{bfactors[start + i]:6.2f}"
                f"           C")
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _boltz_tree(root, n_models=3):
    """A Boltz-style output tree with n_models predictions of varying pLDDT."""
    d = root / "all_outputs" / "seed42" / "predictions" / "job"
    d.mkdir(parents=True)
    for i in range(n_models):
        _pdb(d / f"job_model_{i}.pdb")
        (d / f"confidence_job_model_{i}.json").write_text(
            json.dumps({"ptm": 0.5 + i * 0.1, "iptm": 0.4 + i * 0.1}))
        np.savez(d / f"plddt_job_model_{i}.npz",
                 plddt=np.full(7, 0.5 + i * 0.1))
        np.savez(d / f"pae_job_model_{i}.npz", pae=np.full((7, 7), 5.0))
    return root


def _run(tmp_path, pred_dir, out_csv, *extra):
    sys.argv = ["compute_metrics.py", "--model", "boltz2",
                "--prediction-dir", str(pred_dir),
                "--chain-lengths", *CHAINS,
                "--output-csv", str(out_csv), *map(str, extra)]
    cm.main()
    with open(out_csv, newline="") as fh:
        return list(csv.DictReader(fh))


# ── Basic CLI ─────────────────────────────────────────────────────────────

def test_one_row_is_written_per_prediction(tmp_path, capsys):
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv")
    assert len(rows) == 3
    assert {r["model_name"] for r in rows} == {f"job_model_{i}" for i in range(3)}


def test_confidence_columns_are_populated(tmp_path, capsys):
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv")
    r = rows[0]
    for col in ("avg_plddt", "ptm", "iptm", "ipsae_min", "ranking_score",
                "pae_mean", "ipae"):
        assert r[col] not in ("", None)


def test_rmsd_columns_are_blank_without_a_reference(tmp_path, capsys):
    """No reference means no structural comparison, and the columns must be
    empty rather than zero, which would read as a perfect prediction."""
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv")
    assert rows[0]["rmsd_effector_receptor_aligned"] == ""


def test_an_empty_prediction_dir_writes_no_rows(tmp_path, capsys):
    empty = tmp_path / "empty"
    empty.mkdir()
    rows = _run(tmp_path, empty, tmp_path / "m.csv")
    assert rows == []


# ── Representative selection ──────────────────────────────────────────────

def test_the_published_model_is_the_most_confident(tmp_path, capsys):
    """Selection must be reference-free. Picking the lowest RMSD would need
    the answer the benchmark is trying to measure."""
    best = tmp_path / "best"
    _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv",
         "--best-model-dir", best)
    copied = list(best.iterdir())
    assert len(copied) == 1
    assert "boltz2_best" in copied[0].name
    assert "highest avg_plddt" in capsys.readouterr().out


def test_selection_falls_back_when_no_plddt_is_usable(tmp_path, capsys):
    d = tmp_path / "pred" / "predictions" / "job"
    d.mkdir(parents=True)
    _pdb(d / "job_model_0.pdb", bfactors=(0.0,) * 7)
    best = tmp_path / "best"
    _run(tmp_path, tmp_path / "pred", tmp_path / "m.csv", "--best-model-dir", best)
    assert "no usable pLDDT" in capsys.readouterr().out


# ── With a reference ──────────────────────────────────────────────────────

def test_rmsds_are_computed_against_a_reference(tmp_path, capsys):
    ref = _pdb(tmp_path / "ref.pdb")
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv",
                "--reference-pdb", ref)
    r = rows[0]
    assert float(r["rmsd_receptor"]) == pytest.approx(0.0, abs=1e-4)
    assert float(r["rmsd_effector_receptor_aligned"]) == pytest.approx(0.0, abs=1e-4)
    assert int(r["n_receptor_ca"]) == 4
    assert int(r["n_effector_ca"]) == 3


def test_a_displaced_effector_shows_up_only_in_ra_eff(tmp_path, capsys):
    ref = _pdb(tmp_path / "ref.pdb")
    pred = tmp_path / "pred" / "predictions" / "job"
    pred.mkdir(parents=True)
    _pdb(pred / "job_model_0.pdb", coords_offset=7.0)

    rows = _run(tmp_path, tmp_path / "pred", tmp_path / "m.csv",
                "--reference-pdb", ref)
    r = rows[0]
    assert float(r["rmsd_receptor"]) == pytest.approx(0.0, abs=1e-4)
    assert float(r["rmsd_effector_independent"]) == pytest.approx(0.0, abs=1e-4)
    assert float(r["rmsd_effector_receptor_aligned"]) == pytest.approx(7.0, abs=1e-3)


def test_a_missing_reference_file_is_reported_not_fatal(tmp_path, capsys):
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv",
                "--reference-pdb", tmp_path / "nope.pdb")
    assert len(rows) == 3
    assert "not found" in capsys.readouterr().out


def test_custom_chain_ids_are_honoured(tmp_path, capsys):
    ref = _pdb(tmp_path / "ref.pdb")
    rows = _run(tmp_path, _boltz_tree(tmp_path / "pred"), tmp_path / "m.csv",
                "--reference-pdb", ref,
                "--receptor-chain", "A", "--effector-chain", "B")
    assert rows[0]["pred_receptor_chain"] == "A"
    assert rows[0]["pred_effector_chain"] == "B"


# ── Argument validation ───────────────────────────────────────────────────

def test_an_unknown_model_is_rejected(tmp_path):
    sys.argv = ["compute_metrics.py", "--model", "nonesuch",
                "--prediction-dir", str(tmp_path),
                "--chain-lengths", *CHAINS,
                "--output-csv", str(tmp_path / "m.csv")]
    with pytest.raises(SystemExit):
        cm.main()


def test_every_parser_is_reachable_from_the_cli():
    """PARSERS is what --model validates against, so a parser missing here
    cannot be run by Nextflow at all."""
    assert set(cm.PARSERS) >= {"boltz1", "boltz2", "chai1", "af2m", "af3",
                               "colabfold", "esmfold2"}
