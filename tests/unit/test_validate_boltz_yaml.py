"""Unit tests for bin/validate_boltz_yaml.py.

The script imports PyYAML, so the whole module skips cleanly when PyYAML is
not installed (it is part of the .[test] extra). Exercised via the CLI.
"""

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("yaml", reason="validate_boltz_yaml.py requires PyYAML")

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "validate_boltz_yaml.py"

pytestmark = pytest.mark.local_unit


def _run(yaml_text, model, tmp_path):
    p = tmp_path / "input.yaml"
    p.write_text(yaml_text)
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(p), "--model", model],
        capture_output=True, text=True,
    )


POCKET_6 = """\
version: 1
sequences:
  - protein: {id: A, sequence: MKFL}
  - protein: {id: B, sequence: GTAL}
constraints:
  - pocket:
      binder: B
      contacts: [[A, 1]]
      max_distance: 6.0
      force: true
"""


def test_boltz1_valid_pocket(tmp_path):
    r = _run(POCKET_6, "boltz1", tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr


def test_boltz1_rejects_wrong_max_distance(tmp_path):
    bad = POCKET_6.replace("max_distance: 6.0", "max_distance: 8.0")
    r = _run(bad, "boltz1", tmp_path)
    assert r.returncode == 1
    assert "requires exactly 6.0" in r.stdout


def test_boltz1_rejects_contact_constraints(tmp_path):
    with_contact = POCKET_6 + """\
  - contact:
      token1: [A, 1]
      token2: [B, 1]
      max_distance: 6.0
      force: true
"""
    r = _run(with_contact, "boltz1", tmp_path)
    assert r.returncode == 1
    assert "Contact constraints present" in r.stdout


BOLTZ2_CONTACT = """\
version: 1
sequences:
  - protein: {id: A, sequence: MKFL}
  - protein: {id: B, sequence: GTAL}
constraints:
  - contact:
      token1: [A, 1]
      token2: [B, 2]
      max_distance: 4.4
      force: true
"""


def test_boltz2_valid_contact(tmp_path):
    r = _run(BOLTZ2_CONTACT, "boltz2", tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr


def test_boltz2_rejects_malformed_token(tmp_path):
    bad = BOLTZ2_CONTACT.replace("token1: [A, 1]", "token1: A")
    r = _run(bad, "boltz2", tmp_path)
    assert r.returncode == 1
    assert "malformed token1" in r.stdout
