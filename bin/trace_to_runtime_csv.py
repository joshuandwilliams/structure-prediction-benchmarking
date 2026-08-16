#!/usr/bin/env python3
"""Filter a Nextflow trace.txt down to the predictor processes.

Emits predictor_runtime_stats.csv with wall-clock and memory per predictor,
dropping preprocessing, metrics and aggregation rows.

standalone_elapsed_s folds the shared COLABFOLD_SEARCH step into the variants
that cannot run without it, and equals elapsed_s for everything else. Comparing
raw elapsed_s would give those variants a free MSA search.

Usage:
    trace_to_runtime_csv.py <trace.txt> <output.csv>
"""

import csv
import re
import sys

PREDICTORS = {
    "BOLTZ1", "BOLTZ1_MSA", "BOLTZ1_CONSTRAINED",
    "BOLTZ2", "BOLTZ2_MSA", "BOLTZ2_CONSTRAINED",
    "CHAI1", "AF2M", "AF3", "AF3_NOMSA",
    "COLABFOLD", "COLABFOLD_NOMSA", "ESMFOLD2",
    "BOLTZ2_TEMPLATE", "BOLTZ2_MSA_TEMPLATE", "BOLTZ2_CONSTRAINED_TEMPLATE",
}

NEEDS_MSA = {"BOLTZ1_MSA", "BOLTZ2_MSA", "COLABFOLD", "BOLTZ2_MSA_TEMPLATE"}

_UNITS = {"B": 1 / 1024**3, "KB": 1 / 1024**2, "MB": 1 / 1024, "GB": 1.0, "TB": 1024.0}
_DURATION = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0}

FIELDS = ["model", "status", "exit_code", "queue", "elapsed_hms", "elapsed_s",
          "standalone_elapsed_hms", "standalone_elapsed_s", "pct_cpu",
          "rss_gb", "vmem_gb", "peak_rss_gb", "peak_vmem_gb"]


def parse_mem_gb(text):
    """Nextflow memory string such as '6.5 GB' to GB, or None."""
    if not text or text in {"-", "0"}:
        return None
    m = re.fullmatch(r"([\d.]+)\s*([KMGT]?B)", text.strip())
    return float(m.group(1)) * _UNITS[m.group(2)] if m else None


def parse_duration_s(text):
    """Nextflow duration to whole seconds, or None.

    Handles both forms the observer emits: plain integer milliseconds, and
    human-readable spans like '1h 23m 45s'.
    """
    if not text or text == "-":
        return None
    text = text.strip()
    if text.isdigit():
        return int(int(text) / 1000)
    total, found = 0.0, False
    for value, unit in re.findall(r"(\d+(?:\.\d+)?)(ms|s|m|h)", text):
        total += float(value) * _DURATION[unit]
        found = True
    return int(total) if found else None


def hms(seconds):
    if seconds is None:
        return ""
    return f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def rows_from_trace(lines):
    """Predictor rows from the trace lines, MSA time already folded in."""
    header = lines[0].rstrip("\n").split("\t")
    col = {name: i for i, name in enumerate(header)}

    def field(fields, name):
        i = col.get(name, -1)
        return fields[i] if 0 <= i < len(fields) else ""

    parsed = [ln.rstrip("\n").split("\t") for ln in lines[1:] if ln.strip()]

    msa_s = max((parse_duration_s(field(f, "realtime")) or 0
                 for f in parsed if field(f, "process") == "COLABFOLD_SEARCH"),
                default=0)

    out = []
    for f in parsed:
        name = field(f, "process")
        if name not in PREDICTORS:
            continue
        elapsed = parse_duration_s(field(f, "realtime"))
        standalone = None if elapsed is None else elapsed + (msa_s if name in NEEDS_MSA else 0)
        mem = {k: parse_mem_gb(field(f, k)) for k in ("rss", "vmem", "peak_rss", "peak_vmem")}
        out.append({
            "model": name.lower(),
            "status": field(f, "status"),
            "exit_code": field(f, "exit"),
            "queue": field(f, "queue"),
            "elapsed_hms": hms(elapsed),
            "elapsed_s": "" if elapsed is None else elapsed,
            "standalone_elapsed_hms": hms(standalone),
            "standalone_elapsed_s": "" if standalone is None else standalone,
            "pct_cpu": field(f, "%cpu").replace("%", ""),
            **{f"{k}_gb": "" if v is None else f"{v:.2f}" for k, v in mem.items()},
        })
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: trace_to_runtime_csv.py <trace.txt> <output.csv>")
    trace, output = sys.argv[1], sys.argv[2]

    try:
        with open(trace) as fh:
            lines = fh.readlines()
    except OSError as e:
        sys.exit(f"cannot read {trace}: {e}")
    if len(lines) < 2:
        sys.exit(f"{trace} has no data rows")

    rows = rows_from_trace(lines)
    if not rows:
        sys.exit(f"no predictor rows found in {trace}")

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} predictor rows to {output}")


if __name__ == "__main__":
    main()
