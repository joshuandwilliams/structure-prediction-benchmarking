"""Make the repo's bin/ directory importable from experiment analysis scripts.

The pipeline's Python modules live under bin/ as standalone scripts that
Nextflow invokes. They are not packaged, so bin/ is not on sys.path by
default and modules like compute_metrics cannot be imported directly from
analysis scripts elsewhere in the repo.

Analysis scripts under experiments/ should reuse those modules rather than
copying their logic (the ipSAE / RMSD / parsing code in compute_metrics.py
is the single source of truth for what the metrics mean).

Import this module once at the top of an analysis script:

    import _path_setup  # noqa: F401  — adds bin/ to sys.path

After that, `from compute_metrics import compute_ipsae` and similar imports
from bin/ work normally.

REPO_ROOT (a pathlib.Path) is also exposed for callers that need the repo
root for other purposes (locating params files, results trees, etc.).
"""

from __future__ import annotations

import sys
from pathlib import Path


def _discover_repo_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "main.nf").is_file() and (candidate / "bin").is_dir():
            return candidate
    raise RuntimeError(
        "Could not locate repo root: no ancestor of "
        f"{start} contains both 'main.nf' and 'bin/'. "
        "experiments/_path_setup.py expects to live inside the "
        "structure-prediction-benchmarking repo."
    )


REPO_ROOT: Path = _discover_repo_root(Path(__file__).resolve().parent)

_BIN_DIR = str(REPO_ROOT / "bin")
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
