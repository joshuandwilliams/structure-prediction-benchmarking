/*
 * Sequence extraction from the reference complex, and the effector structural
 * template used only by the Run-2 *_template variants.
 */

/*
 * Parse chain sequences from the reference complex PDB and emit a
 * sequences.json file that the workflow parses with JsonSlurper. Also
 * publishes the reference PDB to the results root so downstream
 * metric/aggregation steps can find it in a predictable location.
 */
process EXTRACT_SEQUENCES {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}", mode: 'copy',
        pattern: '{sequences.json,reference.pdb}'

    input:
    path ref_pdb
    val  receptor_chain
    val  effector_chain

    output:
    path "sequences.json",   emit: sequences_json
    path "reference.pdb",    emit: reference_pdb

    script:
    """
    # Keep an unambiguous copy of the reference PDB alongside the sequences.
    cp ${ref_pdb} reference.pdb

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_sequences.py \\
            reference.pdb \\
            --output sequences.json \\
            --env sequences.env \\
            --chains ${receptor_chain} ${effector_chain}
    """
}

/*
 * Extract the effector chain from the reference PDB as an mmCIF structural
 * template, relabelled to chain A. Used ONLY by the *_template Boltz-2
 * variants (Run 2).  Every other variant in the benchmark is template-free,
 * so that the model-vs-model comparison holds inputs constant.  The
 * *_template variants exist to measure what supplying the known effector
 * fold is worth, as a within-model comparison against the matching template-
 * free Boltz-2 variant.
 */
process EXTRACT_EFFECTOR_TEMPLATE {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}", mode: 'copy',
        pattern: 'effector_template.*'

    input:
    path ref_pdb
    val  effector_chain

    output:
    path "effector_template.cif", emit: template_cif

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_effector_template.py \\
            ${ref_pdb} \\
            ${effector_chain}
    """
}
