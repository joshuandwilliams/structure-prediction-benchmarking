#!/bin/bash
# Run boltz predict across the five benchmark seeds and aggregate the outputs.
#
# Shared by all nine Boltz processes, which differ only in the input file they
# build and whether they pass --model boltz1.
#
# boltz predict exits 0 even when it silently skips an input, so each seed is
# guarded on having produced at least one PDB.
#
#   run_boltz_seeds.sh <input> <container> <pipeline_dir> <label> [extra flags...]

set -euo pipefail

INPUT="${1:?input yaml or fasta}"
CONTAINER="${2:?singularity image}"
PIPELINE_DIR="${3:?repo root}"
LABEL="${4:?label for log lines}"
shift 4

SEEDS=(42 123 456 789 1024)

run_seed() {
    local seed="$1" out_dir="$2"

    singularity exec --nv \
        --bind "${PWD}:${PWD}" \
        "${CONTAINER}" boltz predict \
            "${INPUT}" \
            --out_dir "${out_dir}" \
            --recycling_steps 20 \
            --diffusion_samples 5 \
            --sampling_steps 20 \
            --seed "${seed}" \
            --num_workers 0 \
            --output_format pdb \
            --write_full_pae \
            --use_potentials \
            --no_kernels \
            --override \
            "$@"

    local n_pdb
    n_pdb=$(find "${out_dir}" -name "*.pdb" 2>/dev/null | wc -l)
    if [ "${n_pdb}" -eq 0 ]; then
        echo "ERROR: seed ${seed} produced no PDB files, boltz skipped the input." >&2
        exit 1
    fi
    echo "  seed ${seed}: ${n_pdb} PDB(s) OK"
}

for i in "${!SEEDS[@]}"; do
    seed="${SEEDS[$i]}"
    echo "Running ${LABEL} with seed ${seed}..."
    if [ "$i" -eq 0 ]; then
        run_seed "${seed}" output "$@"
    else
        run_seed "${seed}" "output_seed${seed}" "$@"
    fi
done

bash "${PIPELINE_DIR}/bin/aggregate_seed_outputs.sh" pdb npz json

echo "Aggregated PDB files:"
find all_outputs -name "*.pdb" | sort
