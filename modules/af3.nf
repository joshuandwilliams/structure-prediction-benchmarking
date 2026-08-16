/*
 * AlphaFold 3, with the full data pipeline and with --norun_data_pipeline.
 *
 * AF3_SETUP_DB builds a flat symlink union of the v3.0.0 databases plus BFD and
 * MGnify from v2.3.2, because run_alphafold.py needs them all under one
 * --db_dir. It runs once and both predictor processes consume its output.
 *
 * Both take 5 modelSeeds. The NBI build exposes no sampling flags, so it runs
 * its defaults of 5 diffusion samples per seed and 10 recycles, giving 25
 * structures. Ten recycles is half what Boltz, Chai-1 and ColabFold use here.
 */

/*
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

    # Use the existing symlink farm if it's already populated
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

    # Verify the critical databases resolve (fail early and clearly)
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
 * AlphaFold 3 with the full data pipeline (MSA search) for both chains.
 * Templates are disabled by pinning --max_template_date to 1900-01-01, which
 * returns zero hits.  Left at AF3's own default of 2021-09-30, the template
 * search can retrieve 13 of the 18 Tier-1 benchmark complexes themselves —
 * the answer the benchmark is scoring.  Every variant in the model
 * comparison is template-free, so this keeps the inputs uniform across
 * predictors.
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

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/af3_input.py \\
            "${receptor_seq}" "${effector_seq}" full \\
            --output input.json

    # Sourced by absolute path. `source` resolves its argument against PATH,
    # and /nbi/software/production/bin is only on PATH in a login shell, so a
    # bare `source package` fails whenever the run was submitted from a
    # non-interactive shell.
    source /nbi/software/production/bin/package ${params.af3_package_id}

    echo "Running AlphaFold 3 prediction..."

    run_alphafold.py \\
        --json_path=\${PWD}/input.json \\
        --model_dir="${params.af3_model_dir}" \\
        --db_dir="${af3_db_dir}" \\
        --output_dir=\${PWD}/output \\
        --max_template_date=1900-01-01

    echo "AF3 output files:"
    find output -type f | head -20 || true
    """
}

/*
 * AlphaFold 3 with --norun_data_pipeline: empty MSA, no PDB template search,
 * no injected template.  Pure inference-only folding from sequence alone.
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

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/af3_input.py \\
            "${receptor_seq}" "${effector_seq}" nomsa \\
            --output input.json

    # Sourced by absolute path. `source` resolves its argument against PATH,
    # and /nbi/software/production/bin is only on PATH in a login shell, so a
    # bare `source package` fails whenever the run was submitted from a
    # non-interactive shell.
    source /nbi/software/production/bin/package ${params.af3_package_id}

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
