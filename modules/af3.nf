/*
 * =============================================================================
 * AlphaFold 3 module
 * =============================================================================
 * Three processes:
 *
 *   AF3_SETUP_DB  — idempotently build the combined AF3 database directory
 *                   at $HOME/af3_db via symlinks to the shared reference
 *                   data location.  Runs exactly once; AF3 and AF3_NOMSA
 *                   both consume its output channel.
 *
 *   AF3           — full MSA + template search (runs the data pipeline)
 *   AF3_NOMSA     — --norun_data_pipeline with empty MSA/template fields
 *
 * Both predictor variants take JSON input with 5 modelSeeds
 * (42, 123, 456, 789, 1024).  The NBI AF3 installation (Dec 2024 build)
 * does not expose --num_diffusion_samples / --num_seeds / --num_recycles
 * CLI flags; sampling is controlled only via the modelSeeds JSON field.
 * Default is 5 diffusion samples per seed and 10 recycles, giving
 * 5 × 5 = 25 output structures.
 *
 * The AF3 data directory lives at $HOME/af3_db because AF3's run_alphafold.py
 * requires all databases to sit in a single --db_dir.  The shared
 * reference-data location cannot itself serve as --db_dir because the
 * BFD/MGnify files from db-v2.3.2 need to appear alongside the v3.0.0 files,
 * so we build a flat union via symlinks.
 *
 * As with AF2M, no env-snapshot/restore is needed because COMPUTE_METRICS is
 * a separate Nextflow process with a clean environment — any PATH pollution
 * from `source package` stays confined to the AF3 work directory.
 */


/*
 * AF3_SETUP_DB
 * ------------
 * Create/refresh the $HOME/af3_db symlink farm exactly once.  Emits a
 * value channel containing the dir path, which AF3 and AF3_NOMSA both
 * consume.  Using `-resume` will skip this process cleanly because the
 * output channel content (a single path string) is deterministic.
 */
process AF3_SETUP_DB {
    tag "${params.project_name}"
    label 'cpu'

    output:
    path "af3_db_ready.flag", emit: ready_flag
    val  "${System.getenv('HOME')}/af3_db", emit: db_dir

    script:
    """
    set -euo pipefail

    AF3_DATA_DIR="\${HOME}/af3_db"
    mkdir -p "\${AF3_DATA_DIR}"

    # ─── Link v3.0.0 databases ────────────────────────────────────────────
    for f in ${params.af3_db_v3}/*; do
        ln -sfn "\$f" "\${AF3_DATA_DIR}/\$(basename "\$f")"
    done

    # ─── Link BFD/MGnify from v2.3.2 (AF3 still needs these) ─────────────
    ln -sfn ${params.af2_data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta \\
        "\${AF3_DATA_DIR}/bfd-first_non_consensus_sequences.fasta"
    ln -sfn ${params.af2_data_dir}/mgnify/mgy_clusters_2022_05.fa \\
        "\${AF3_DATA_DIR}/mgy_clusters_2022_05.fa"

    echo "AF3 database dir ready: \${AF3_DATA_DIR}"
    ls -la "\${AF3_DATA_DIR}" | head -20

    # Touch a sentinel file that Nextflow can pass as a path dependency.
    touch af3_db_ready.flag
    """
}


/*
 * AF3
 * ---
 * AlphaFold 3 with full data pipeline (MSA + template search).
 */
process AF3 {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/af3", mode: 'copy',
        pattern: '{input.json,output/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  af3_db_dir
    path db_ready_flag

    output:
    path "output",      emit: prediction_dir
    path "input.json"

    script:
    """
    set -euo pipefail

    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2

    mkdir -p output

    # ─── Load AF3 environment ─────────────────────────────────────────────
    source package ${params.af3_package_id}

    # ─── Write AF3 JSON input ─────────────────────────────────────────────
    cat > input.json << EOF
{
  "name": "benchmark_test",
  "modelSeeds": [42, 123, 456, 789, 1024],
  "dialect": "alphafold3",
  "version": 1,
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "${receptor_seq}"
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "${effector_seq}"
      }
    }
  ]
}
EOF

    echo "Running AlphaFold 3 prediction..."

    run_alphafold.py \\
        --json_path=\${PWD}/input.json \\
        --model_dir="${params.af3_model_dir}" \\
        --db_dir="${af3_db_dir}" \\
        --output_dir=\${PWD}/output

    echo "AF3 output files:"
    find output -type f | head -20 || true
    """
}


/*
 * AF3_NOMSA
 * ---------
 * AlphaFold 3 in single-sequence mode (--norun_data_pipeline + empty
 * unpairedMsa / pairedMsa / templates fields).  Runs inference only on
 * GPU; no MSA or template search is performed.
 */
process AF3_NOMSA {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/af3_nomsa", mode: 'copy',
        pattern: '{input.json,output/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  af3_db_dir
    path db_ready_flag

    output:
    path "output",      emit: prediction_dir
    path "input.json"

    script:
    """
    set -euo pipefail

    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2

    mkdir -p output

    source package ${params.af3_package_id}

    # ─── AF3 JSON with empty MSA + no templates ───────────────────────────
    cat > input.json << EOF
{
  "name": "benchmark_test",
  "modelSeeds": [42, 123, 456, 789, 1024],
  "dialect": "alphafold3",
  "version": 1,
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "${receptor_seq}",
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "${effector_seq}",
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
      }
    }
  ]
}
EOF

    echo "Running AlphaFold 3 prediction (no MSA, inference only)..."

    run_alphafold.py \\
        --json_path=\${PWD}/input.json \\
        --model_dir="${params.af3_model_dir}" \\
        --db_dir="${af3_db_dir}" \\
        --output_dir=\${PWD}/output \\
        --norun_data_pipeline

    echo "AF3 (no MSA) output files:"
    find output -type f | head -20 || true
    """
}
