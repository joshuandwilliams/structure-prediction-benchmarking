/*
 * =============================================================================
 * Preprocessing module
 * =============================================================================
 * Runs extract_sequences.py inside the benchmark_models.img container to
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
