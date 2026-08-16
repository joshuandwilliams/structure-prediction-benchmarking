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

# A real gemmi-written mmCIF with two Cα atoms at B-factor 90 and 80.
# Written by gemmi rather than hand-rolled, so it parses the way an
# actual ESMFold2 output does.
_MINIMAL_CIF = """data_pred
_entry.id pred

_cell.entry_id pred
_cell.length_a 1
_cell.length_b 1
_cell.length_c 1
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90

_symmetry.entry_id pred
_symmetry.space_group_name_H-M ''

loop_
_entity.id
_entity.type
ALA! non-polymer
GLY! non-polymer




loop_
_chem_comp.id
_chem_comp.type
ALA .
GLY .

loop_
_struct_asym.id
_struct_asym.entity_id
Ax1 ALA!
Bx1 GLY!



loop_
_atom_type.symbol
C


loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C CA . ALA Ax1 ALA! . ? 0 0 0 1 90 ? 1 A 1
HETATM 2 C CA . GLY Bx1 GLY! . ? 10 0 0 1 80 ? 1 B 1
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
    # ranking_score = 0.8*iptm + 0.2*ptm
    assert r["ranking_score"] == pytest.approx(0.8 * 0.4 + 0.2 * 0.6, abs=1e-4)
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


# ── HPC tier: the fold itself ─────────────────────────────────────────────

@pytest.mark.hpc
def test_a_real_fold_writes_a_structure_and_confidences(tmp_path):
    """The model-loading path needs a 40 GB GPU and the esm package, so it
    cannot run on a laptop. Everything above tests the pure-Python parts
    around it."""
    pytest.importorskip("esm", reason="ESMFold2 fold needs the esm package")
    pytest.importorskip("torch", reason="ESMFold2 fold needs torch")
    import json
    import subprocess
    import sys
    from pathlib import Path

    repo = Path(__file__).resolve().parents[2]
    spec = tmp_path / "input.json"
    spec.write_text(json.dumps({
        "receptor": {"id": "A", "sequence": "MKFLVAAALLLGAVSA"},
        "effector": {"id": "B", "sequence": "GTALPPWWQDFAERLK"}}))

    out = tmp_path / "seed42"
    r = subprocess.run(
        [sys.executable, str(repo / "bin" / "esmfold2_fold.py"),
         "--input-json", str(spec), "--out-dir", str(out), "--seed", "42"],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (out / "esmfold2_pred.cif").is_file()

    conf = json.loads((out / "confidences.json").read_text())
    assert len(conf["plddt"]) > 0
    assert max(conf["plddt"]) > 1.5     # rescaled to 0-100, not left as 0-1
