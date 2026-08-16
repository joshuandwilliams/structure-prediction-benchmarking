"""The argv contract of bin/run_boltz_seeds.sh.

All nine Boltz processes route through this script, and it is shell rather than
Python, so the coverage gate over bin/ says nothing about it. A regression here
is expensive: `boltz predict` rejects the whole command and Nextflow records an
ignored error, so the run completes with that variant silently missing.

singularity is stubbed with a script that echoes the argv it was handed, which
is enough to pin how arguments reach boltz without a GPU or a container.
"""

import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.local_unit

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "run_boltz_seeds.sh"
SEEDS = ["42", "123", "456", "789", "1024"]

STUB = """#!/bin/bash
# Drop everything up to and including the `boltz` token, then record the rest.
while [ $# -gt 0 ]; do
    case "$1" in boltz) shift; break ;; *) shift ;; esac
done
echo "BOLTZ $*" >> "${ARGV_LOG}"

# Write a PDB where boltz was told to, so the script's own "did this seed
# produce anything" guard sees real output.
out=""
prev=""
for a in "$@"; do
    [ "${prev}" = "--out_dir" ] && out="${a}"
    prev="${a}"
done
[ -n "${out}" ] && { mkdir -p "${out}"; : > "${out}/model.pdb"; }
exit 0
"""


@pytest.fixture
def run(tmp_path):
    """Run the script with singularity and the aggregator stubbed out."""
    bindir = tmp_path / "stub_bin"
    bindir.mkdir()
    log = tmp_path / "argv.log"

    sing = bindir / "singularity"
    sing.write_text(STUB)
    sing.chmod(0o755)

    # The script hands off to aggregate_seed_outputs.sh and then lists
    # all_outputs/. The stub only has to create that directory, which keeps
    # this test about argv rather than about aggregation.
    repo = tmp_path / "repo" / "bin"
    repo.mkdir(parents=True)
    agg = repo / "aggregate_seed_outputs.sh"
    agg.write_text('#!/bin/bash\nmkdir -p all_outputs\n: > all_outputs/m.pdb\n')
    agg.chmod(0o755)

    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "input.fasta").write_text(">A\nMKF\n")

    def _run(*extra):
        env = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}",
                   ARGV_LOG=str(log))
        proc = subprocess.run(
            ["bash", str(SCRIPT), "input.fasta", "/fake.img",
             str(tmp_path / "repo"), "Boltz-test", *extra],
            cwd=workdir, capture_output=True, text=True, env=env)
        lines = log.read_text().splitlines() if log.exists() else []
        return proc, lines

    return _run


def test_extra_flags_reach_boltz_without_positional_leakage(run):
    """The regression that shipped. `local seed="$1"` does not consume the
    positional parameters, so "$@" still began with the seed and the output
    directory and boltz rejected them as extra arguments."""
    proc, lines = run("--model", "boltz1")
    assert proc.returncode == 0, proc.stderr
    assert lines, "singularity was never invoked"
    for line in lines:
        args = line.split()
        assert "--model" in args and "boltz1" in args
        # Nothing between the input and the first flag, and no bare seed.
        assert args[:4] == ["BOLTZ", "predict", "input.fasta", "--out_dir"]
        for seed in SEEDS:
            assert seed not in args[:5], f"seed leaked into argv: {line}"


def test_every_seed_runs_once_into_its_own_output_dir(run):
    proc, lines = run("--model", "boltz1")
    assert proc.returncode == 0, proc.stderr
    assert len(lines) == len(SEEDS)

    seeds_used, out_dirs = [], []
    for line in lines:
        args = line.split()
        seeds_used.append(args[args.index("--seed") + 1])
        out_dirs.append(args[args.index("--out_dir") + 1])

    assert seeds_used == SEEDS
    # The first seed writes to output/, the rest to output_seed<N>/, which is
    # what aggregate_seed_outputs.sh expects to find.
    assert out_dirs[0] == "output"
    assert out_dirs[1:] == [f"output_seed{s}" for s in SEEDS[1:]]


def test_it_runs_with_no_extra_flags(run):
    """Boltz-2 passes no --model, so the no-extra-args path has to work too."""
    proc, lines = run()
    assert proc.returncode == 0, proc.stderr
    assert len(lines) == len(SEEDS)
    for line in lines:
        assert "--model" not in line.split()


def test_a_seed_producing_no_pdb_fails_loudly(run, tmp_path):
    """boltz predict exits 0 even when it silently skips an input, so the
    script checks for output itself. Without that a benchmark reports success
    for a variant that produced nothing."""
    sing = tmp_path / "stub_bin" / "singularity"
    sing.write_text("#!/bin/bash\nexit 0\n")   # writes no PDB
    sing.chmod(0o755)

    proc, _ = run("--model", "boltz1")
    assert proc.returncode != 0
    assert "produced no PDB" in proc.stderr


def test_the_seed_list_matches_the_documented_sampling(run):
    """25 predictions per variant is 5 seeds by 5 diffusion samples, and both
    numbers are quoted in the methods."""
    proc, lines = run()
    assert proc.returncode == 0, proc.stderr
    for line in lines:
        args = line.split()
        assert args[args.index("--diffusion_samples") + 1] == "5"
    assert len(lines) == 5
