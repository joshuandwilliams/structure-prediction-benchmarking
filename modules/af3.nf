/*
 * =============================================================================
 * AlphaFold 3 module
 * =============================================================================
 * Three processes:
 *
 *   AF3_SETUP_DB  — idempotently build the combined AF3 database directory
 *                   at params.af3_db_dir via symlinks to the shared reference
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
 * The AF3 data directory lives at params.af3_db_dir because AF3's run_alphafold.py
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
 * Create/refresh the params.af3_db_dir symlink farm exactly once.  Emits a
 * value channel containing the dir path, which AF3 and AF3_NOMSA both
 * consume.  Using `-resume` will skip this process cleanly because the
 * output channel content (a single path string) is deterministic.
 */
process AF3_SETUP_DB {
    tag "${params.project_name}"
    label 'cpu'

    output:
    path "af3_db_ready.flag", emit: ready_flag
    val  "${params.af3_db_dir}", emit: db_dir

    script:
    """
    set -euo pipefail

    AF3_DATA_DIR="${params.af3_db_dir}"

    # ─── Use the existing symlink farm if it's already populated ─────────
    # The farm is maintained at params.af3_db_dir.  If the critical BFD link
    # is already there we trust it; otherwise (re)build the farm from the
    # shared reference data.  ln -sfn is idempotent and never deletes links
    # the user added by hand, so this is safe to run against a curated farm.
    if [ -e "\${AF3_DATA_DIR}/bfd-first_non_consensus_sequences.fasta" ]; then
        echo "Using existing AF3 database farm at \${AF3_DATA_DIR}"
    else
        echo "Building AF3 database farm at \${AF3_DATA_DIR}..."
        mkdir -p "\${AF3_DATA_DIR}"
        for f in ${params.af3_db_v3}/*; do
            ln -sfn "\$f" "\${AF3_DATA_DIR}/\$(basename "\$f")"
        done
        ln -sfn ${params.af2_data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta \\
            "\${AF3_DATA_DIR}/bfd-first_non_consensus_sequences.fasta"
        ln -sfn ${params.af2_data_dir}/mgnify/mgy_clusters_2022_05.fa \\
            "\${AF3_DATA_DIR}/mgy_clusters_2022_05.fa"
    fi

    # ─── Verify the critical databases resolve (fail early and clearly) ──
    # Catches a dangling symlink here instead of deep inside run_alphafold.py.
    for req in bfd-first_non_consensus_sequences.fasta mgy_clusters_2022_05.fa mmcif_files; do
        if [ ! -e "\${AF3_DATA_DIR}/\${req}" ]; then
            echo "ERROR: \${AF3_DATA_DIR}/\${req} is missing or a dangling symlink." >&2
            echo "       Check params.af3_db_dir and the shared reference-data paths." >&2
            exit 1
        fi
    done

    echo "AF3 database dir ready: \${AF3_DATA_DIR}"
    ls -la "\${AF3_DATA_DIR}" | head -20

    # Touch a sentinel file that Nextflow can pass as a path dependency.
    touch af3_db_ready.flag
    """
}


/*
 * AF3
 * ---
 * AlphaFold 3 with the full data pipeline (MSA + PDB template search) for both
 * chains — vanilla AF3.  No custom effector template is injected here: AF3
 * rejects a chain that combines an auto-built MSA with a custom template
 * (ValueError: "...set only partially"), so effector-template steering is
 * applied only in AF3_NOMSA (and ColabFold), not in this with-MSA condition.
 *
 * Note: AF2M cannot inject custom templates via its CLI either; its templates
 * are limited to structures deposited before 2020-05-14.
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

    # ─── Generate AF3 JSON (vanilla AF3: full pipeline, no template) ──────
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/af3_input.py \\
            "${receptor_seq}" "${effector_seq}" full \\
            --output input.json

    # ─── Load AF3 environment and run ────────────────────────────────────
    source package ${params.af3_package_id}

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
 * AlphaFold 3 with --norun_data_pipeline: empty MSA, no PDB template search.
 * The effector structural template is still injected if available — giving
 * inference-only folding guided by the known effector structure but no MSA.
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
    path effector_template_cif

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

    # ─── Generate AF3 JSON (no MSA; effector template if available) ───────
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/af3_input.py \\
            "${receptor_seq}" "${effector_seq}" nomsa \\
            --template ${effector_template_cif} \\
            --output input.json

    source package ${params.af3_package_id}

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
