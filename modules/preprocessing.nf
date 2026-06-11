/*
 * =============================================================================
 * Preprocessing module
 * =============================================================================
 * Runs extract_sequences.py inside the Boltz1_Boltz2_Chai1.img container to
 * pull chain sequences out of the reference PDB.  Downstream processes read
 * chain_A / chain_B sequences as value channels parsed from the produced
 * sequences.json.
 *
 * The --chains flag pins which PDB chains become CHAIN_A/CHAIN_B aliases.
 * For 7B1I that is "B C", for 6G10 it is "A B", etc.
 */


/*
 * EXTRACT_SEQUENCES
 * -----------------
 * Parse chain sequences from the reference complex PDB and emit a
 * sequences.json file that the workflow parses with JsonSlurper.
 *
 * Also publishes the reference PDB to the results root so downstream
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
 * EXTRACT_EFFECTOR_TEMPLATE
 * -------------------------
 * Extract the effector chain from the reference PDB and write it as both
 * PDB and mmCIF for use as a structural template in AF3 and ColabFold.
 * Only runs in PDB input mode; FASTA mode passes a no-template sentinel.
 *
 * The extracted chain is relabelled to chain A (single-chain convention).
 * AF3 uses the mmCIF form embedded in the JSON template block; ColabFold
 * uses the PDB form via --custom-template-path.
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
    path "effector_template.pdb", emit: template_pdb
    path "effector_template.cif", emit: template_cif

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_effector_template.py \\
            ${ref_pdb} \\
            ${effector_chain}
    """
}
