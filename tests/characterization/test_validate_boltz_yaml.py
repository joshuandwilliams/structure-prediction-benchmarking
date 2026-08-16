"""Unit tests for bin/validate_boltz_yaml.py.

Skips cleanly when PyYAML is absent. Run in-process rather than through the
CLI so coverage records the module.
"""

import sys
from types import SimpleNamespace

import pytest

pytest.importorskip("yaml", reason="validate_boltz_yaml.py requires PyYAML")

import validate_boltz_yaml as vby  # noqa: E402

pytestmark = pytest.mark.local_unit


def _run(yaml_text, model, tmp_path, capsys=None):
    """Call main in-process, returning (returncode, stdout)."""
    p = tmp_path / "input.yaml"
    p.write_text(yaml_text)
    sys.argv = ["validate_boltz_yaml.py", str(p), "--model", model]
    code = 0
    try:
        vby.main()
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else 1
    out = capsys.readouterr() if capsys else SimpleNamespace(out="", err="")
    return SimpleNamespace(returncode=code, stdout=out.out, stderr=out.err)


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


def test_boltz1_valid_pocket(tmp_path, capsys):
    r = _run(POCKET_6, "boltz1", tmp_path, capsys)
    assert r.returncode == 0, r.stdout + r.stderr


def test_boltz1_rejects_wrong_max_distance(tmp_path, capsys):
    bad = POCKET_6.replace("max_distance: 6.0", "max_distance: 8.0")
    r = _run(bad, "boltz1", tmp_path, capsys)
    assert r.returncode == 1
    assert "requires exactly 6.0" in r.stdout


def test_boltz1_rejects_contact_constraints(tmp_path, capsys):
    with_contact = POCKET_6 + """\
  - contact:
      token1: [A, 1]
      token2: [B, 1]
      max_distance: 6.0
      force: true
"""
    r = _run(with_contact, "boltz1", tmp_path, capsys)
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


def test_boltz2_valid_contact(tmp_path, capsys):
    r = _run(BOLTZ2_CONTACT, "boltz2", tmp_path, capsys)
    assert r.returncode == 0, r.stdout + r.stderr


def test_boltz2_rejects_malformed_token(tmp_path, capsys):
    bad = BOLTZ2_CONTACT.replace("token1: [A, 1]", "token1: A")
    r = _run(bad, "boltz2", tmp_path, capsys)
    assert r.returncode == 1
    assert "malformed token1" in r.stdout
