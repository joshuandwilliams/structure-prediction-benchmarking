/*
 * ESMFold2 (ESMC-6B backbone), single-sequence and MSA-free, natively
 * multi-chain. One structure per seed rather than the 25 the others produce,
 * so any per-prediction comparison must divide by 5 rather than 25.
 *
 * ESMC-6B is ~25 GB in fp32, so this needs a 40 GB GPU and ~96 GB RAM.
 */

/* Fold the receptor+effector complex across 5 seeds, one mmCIF per seed. */
process ESMFOLD2 {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/esmfold2", mode: 'copy',
        pattern: '{esmfold2_input.json,all_outputs/**}'

    input:
    val receptor_seq
    val effector_seq
    val rec_chain
    val eff_chain

    output:
    path "all_outputs",        emit: prediction_dir
    path "esmfold2_input.json"

    script:
    def num_loops          = 20
    def num_sampling_steps = 100
    """
    set -euo pipefail

    cat > esmfold2_input.json << JSONEOF
{
  "receptor": {"id": "${rec_chain}", "sequence": "${receptor_seq}"},
  "effector": {"id": "${eff_chain}", "sequence": "${effector_seq}"}
}
JSONEOF

    echo "Input spec:"
    cat esmfold2_input.json
    echo ""

    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.esmfold2_container} \\
                python ${projectDir}/bin/esmfold2_fold.py \\
                    --input-json \${PWD}/esmfold2_input.json \\
                    --out-dir "\${out_dir}" \\
                    --seed "\${seed}" \\
                    --num-loops ${num_loops} \\
                    --num-sampling-steps ${num_sampling_steps}

        local n_cif
        n_cif=\$(find "\${out_dir}" -name "*.cif" 2>/dev/null | wc -l)
        if [ "\${n_cif}" -eq 0 ]; then
            echo "ERROR: Seed \${seed} produced no CIF files." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_cif} CIF OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running ESMFold2 with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    bash ${projectDir}/bin/aggregate_seed_outputs.sh cif json

    echo "Aggregated ESMFold2 CIF files:"
    find all_outputs -name "*.cif" | sort
    """
}
