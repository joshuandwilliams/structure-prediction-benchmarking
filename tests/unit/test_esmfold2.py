"""Unit tests for the ESMFold2 integration:
  - bin/esmfold2_fold.py::build_confidences (the pure, GPU-free part of the runner)
  - bin/compute_metrics.py::parse_esmfold2 (reads cif + confidences.json)

Both run on a laptop — no esm/torch/GPU needed (the heavy imports in
esmfold2_fold.py are deferred into main()).
"""

import json

import compute_metrics as cm
import esmfold2_fold as ef
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit


# ── build_confidences ─────────────────────────────────────────────────────

def test_build_confidences_rescales_plddt_0_1_to_0_100():
    conf = ef.build_confidences([0.0, 0.5, 1.0])
    assert conf["plddt"] == [0.0, 50.0, 100.0]


def test_build_confidences_leaves_0_100_plddt_alone():
    conf = ef.build_confidences([10.0, 80.0, 95.0])
    assert conf["plddt"] == [10.0, 80.0, 95.0]


def test_build_confidences_passes_through_ptm_iptm():
    conf = ef.build_confidences([0.9], ptm=0.42, iptm=0.31)
    assert conf["ptm"] == pytest.approx(0.42)
    assert conf["iptm"] == pytest.approx(0.31)


def test_build_confidences_omits_absent_scores():
    conf = ef.build_confidences([0.9])
    assert "ptm" not in conf and "iptm" not in conf and "pae" not in conf


def test_build_confidences_includes_2d_pae_only():
    conf = ef.build_confidences([0.9], pae=[[0.0, 1.0], [1.0, 0.0]])
    assert conf["pae"] == [[0.0, 1.0], [1.0, 0.0]]
    # a 1-D "pae" is not a valid matrix → dropped
    assert "pae" not in ef.build_confidences([0.9], pae=[1.0, 2.0, 3.0])


# ── parse_esmfold2 ────────────────────────────────────────────────────────

_MINIMAL_CIF = """\
data_pred
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.B_iso_or_equiv
ATOM CA A 90.0
ATOM CA B 80.0
"""


def _write_seed(seed_dir, plddt, ptm=None, iptm=None, pae=None):
    seed_dir.mkdir(parents=True, exist_ok=True)
    (seed_dir / "esmfold2_pred.cif").write_text(_MINIMAL_CIF)
    conf = {"plddt": plddt}
    if ptm is not None:
        conf["ptm"] = ptm
    if iptm is not None:
        conf["iptm"] = iptm
    if pae is not None:
        conf["pae"] = pae
    (seed_dir / "confidences.json").write_text(json.dumps(conf))


def test_parse_reads_confidences(tmp_path):
    _write_seed(tmp_path / "output_seed42", plddt=[80.0, 90.0], ptm=0.6, iptm=0.4)
    results = cm.parse_esmfold2(str(tmp_path), [1, 1])
    assert len(results) == 1
    r = results[0]
    assert r["avg_plddt"] == 85.0
    assert r["ptm"] == 0.6
    assert r["iptm"] == 0.4
    # actifptm = 0.8*iptm + 0.2*ptm
    assert r["actifptm"] == pytest.approx(0.8 * 0.4 + 0.2 * 0.6, abs=1e-4)
    assert r["pdb_path"].endswith("esmfold2_pred.cif")


def test_parse_collects_all_seeds(tmp_path):
    for s in (42, 123, 456):
        _write_seed(tmp_path / f"output_seed{s}", plddt=[70.0], ptm=0.5)
    results = cm.parse_esmfold2(str(tmp_path), [1, 1])
    assert len(results) == 3


def test_parse_pae_drives_ipsae(tmp_path):
    # Perfect 2x2-chain interface (all-zero inter-chain PAE) → ipsae_min == 1.0
    pae = np.zeros((4, 4)).tolist()
    _write_seed(tmp_path / "output_seed42", plddt=[90.0], ptm=0.9, iptm=0.9, pae=pae)
    r = cm.parse_esmfold2(str(tmp_path), [2, 2])[0]
    assert r["ipsae_min"] == 1.0
    assert r["ipae"] == 0.0


def test_parse_missing_confidences_falls_back_to_cif_plddt(tmp_path):
    # No confidences.json → avg_plddt read from the cif B-factor column.
    seed = tmp_path / "output_seed42"
    seed.mkdir(parents=True)
    (seed / "esmfold2_pred.cif").write_text(_MINIMAL_CIF)
    r = cm.parse_esmfold2(str(tmp_path), [1, 1])[0]
    assert r["avg_plddt"] == 85.0   # mean of 90 and 80 from the cif
    assert r["ipsae_min"] == 0.0    # no PAE → interface metrics fall back to 0


def test_parse_empty_dir_returns_nothing(tmp_path):
    assert cm.parse_esmfold2(str(tmp_path), [1, 1]) == []


def test_esmfold2_registered_in_parsers():
    assert cm.PARSERS["esmfold2"] is cm.parse_esmfold2
