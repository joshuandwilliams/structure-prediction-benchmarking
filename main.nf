#!/usr/bin/env nextflow

/*
 * Benchmark protein-complex structure predictors against a reference complex.
 *
 * Run 1 is the 13 template-free variants, the model comparison. Run 2 is the
 * three Boltz-2 variants that receive the reference effector as a structural
 * template, which pair with their template-free twins to measure what the known
 * effector fold is worth. Select which to run via the models list in params.
 *
 * Each variant writes to ${outdir}/<model>/. The ColabFold MSA is computed once
 * and shared. GPU seed loops are serial to respect the small jic-gpu queue.
 */

nextflow.enable.dsl = 2

// Parameter defaults

// Project
params.project_name     = "structure_prediction_benchmark"
params.outdir           = "${launchDir}/${params.project_name}_results"

// A two-chain reference complex plus the chain IDs pinning receptor and
// effector. Required, since every run scores RMSDs against it.
params.reference_pdb    = null
params.receptor_chain   = null
params.effector_chain   = null

// Model selection
// A list of models to run; null means 'run them all'.
params.models           = null

// Infrastructure paths
// Containers
params.benchmark_container = "/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/Boltz1_Boltz2_Chai1.img"
params.colabfold_container = "/hpc-home/jowillia/singularity/ColabFold/colabfold.img"
params.esmfold2_container  = "/hpc-home/jowillia/singularity/ESMFold2/esmfold2.img"

// Databases
params.colabfold_db     = "/nbi/Reference-Data/AlphaFold/colabfold_databases"
params.af2_data_dir     = "/nbi/Reference-Data/AlphaFold/db-v2.3.2"
params.af3_db_v3        = "/nbi/Reference-Data/AlphaFold/db-v3.0.0"

// AF3 model weights directory
params.af3_model_dir    = "/hpc-home/jowillia/singularity/AlphaFold3"

// AF3 database directory: a single --db_dir containing symlinks to all the
// reference databases AF3 needs (v3.0.0 dbs + BFD/MGnify from v2.3.2).
// AF3_SETUP_DB uses this farm if it's already populated, else builds it.
params.af3_db_dir       = "/hpc-home/jowillia/singularity/AlphaFold3/af3_db"

// HPC source package UUIDs
params.af2_package_id   = "be036146-04b2-42a1-b7f8-26af883f2281"
params.af3_package_id   = "e8edb411-7374-4342-b9f1-408da41fc197"

// GPU concurrency (matches max_gpu_parallel in nextflow.config)
params.max_gpu_parallel = 1

// Wrapped in a function because newer Nextflow rejects top-level statements.
def validate_params() {
    if (!params.reference_pdb) {
        error "Please provide --reference_pdb (a two-chain reference complex)."
    }
    if (!params.receptor_chain || !params.effector_chain) {
        error "Please provide --receptor_chain and --effector_chain (PDB chain IDs)."
    }
    if (params.receptor_chain == params.effector_chain) {
        error "--receptor_chain and --effector_chain must differ " +
              "(both were '${params.receptor_chain}')."
    }
}

// Model selection

// MODEL_TO_PARSER names the output format, so variants sharing a format
// (boltz1_msa parses as boltz1) still keep their own tag in the output CSV.
def ALL_MODELS = [
    'boltz1',              // parser: boltz1
    'boltz1_msa',          // parser: boltz1
    'boltz1_constrained',  // parser: boltz1
    'boltz2',              // parser: boltz2
    'boltz2_msa',          // parser: boltz2
    'boltz2_constrained',  // parser: boltz2
    'chai1',               // parser: chai1
    'af2m',                // parser: af2m
    'af3',                 // parser: af3
    'af3_nomsa',           // parser: af3
    'colabfold',           // parser: colabfold
    'colabfold_nomsa',     // parser: colabfold
    'esmfold2',            // parser: esmfold2
    // Run 2, not part of the model comparison
    'boltz2_template',              // parser: boltz2
    'boltz2_msa_template',          // parser: boltz2
    'boltz2_constrained_template',  // parser: boltz2
]

def MODEL_TO_PARSER = [
    'boltz1'             : 'boltz1',
    'boltz1_msa'         : 'boltz1',
    'boltz1_constrained' : 'boltz1',
    'boltz2'             : 'boltz2',
    'boltz2_msa'         : 'boltz2',
    'boltz2_constrained' : 'boltz2',
    'chai1'              : 'chai1',
    'af2m'               : 'af2m',
    'af3'                : 'af3',
    'af3_nomsa'          : 'af3',
    'colabfold'          : 'colabfold',
    'colabfold_nomsa'    : 'colabfold',
    'esmfold2'           : 'esmfold2',
    'boltz2_template'             : 'boltz2',
    'boltz2_msa_template'         : 'boltz2',
    'boltz2_constrained_template' : 'boltz2',
]

// ALL_MODELS is passed in rather than read from the enclosing scope: a
// top-level `def` is local to the script body and is not visible here.
def selected_models(all_models) {
    if (params.models == null) {
        return all_models
    }
    // May arrive as a YAML list or a comma/space-separated CLI string.
    def raw = params.models
    def requested = (raw instanceof List)
        ? raw.collect { it.toString().trim() }
        : raw.toString().replaceAll(',', ' ').split()*.trim()
    requested = requested.findAll { it }

    def unknown = requested.findAll { !(it in all_models) }
    if (unknown) {
        error "Unknown model(s) in --models: ${unknown.join(', ')}\nValid: ${all_models.join(', ')}"
    }
    return requested
}

def SELECTED_MODELS = selected_models(ALL_MODELS)

def needs_shared_msa      = SELECTED_MODELS.any { it in ['boltz1_msa', 'boltz2_msa', 'colabfold',
                                                         'boltz2_msa_template'] }
def needs_af3_db          = SELECTED_MODELS.any { it in ['af3', 'af3_nomsa'] }

// The effector structural template is used ONLY by the Run-2 *_template variants.
// Every variant in the model-vs-model comparison is template-free.
def needs_eff_template    = SELECTED_MODELS.any { it in ['boltz2_template', 'boltz2_msa_template',
                                                         'boltz2_constrained_template'] }

log.info """
=============================================================
Structure Prediction Benchmark — v0.1.0
=============================================================
Project:         ${params.project_name}
Reference PDB:   ${params.reference_pdb}
Receptor chain:  ${params.receptor_chain}
Effector chain:  ${params.effector_chain}
Output:          ${params.outdir}
Models:          ${SELECTED_MODELS.join(', ')}
Shared MSA:      ${needs_shared_msa ? 'yes' : 'no'}
AF3 DB setup:    ${needs_af3_db ? 'yes' : 'no'}
=============================================================
""".stripIndent()

// Include modules

include { EXTRACT_SEQUENCES            } from './modules/preprocessing'
include { EXTRACT_EFFECTOR_TEMPLATE    } from './modules/preprocessing'
include { COLABFOLD_SEARCH             } from './modules/msa'

include { BOLTZ1                  } from './modules/boltz1'
include { BOLTZ1_MSA              } from './modules/boltz1'
include { BOLTZ1_CONSTRAINED      } from './modules/boltz1'

include { BOLTZ2                  } from './modules/boltz2'
include { BOLTZ2_MSA              } from './modules/boltz2'
include { BOLTZ2_CONSTRAINED      } from './modules/boltz2'
include { BOLTZ2_TEMPLATE             } from './modules/boltz2'
include { BOLTZ2_MSA_TEMPLATE         } from './modules/boltz2'
include { BOLTZ2_CONSTRAINED_TEMPLATE } from './modules/boltz2'

include { CHAI1                   } from './modules/chai1'
include { AF2M                    } from './modules/af2m'

include { AF3_SETUP_DB            } from './modules/af3'
include { AF3                     } from './modules/af3'
include { AF3_NOMSA               } from './modules/af3'

include { COLABFOLD               } from './modules/colabfold'
include { COLABFOLD_NOMSA         } from './modules/colabfold'

include { ESMFOLD2                } from './modules/esmfold2'

// Metrics: one aliased instance per model so publishDir subdirs stay distinct
include { COMPUTE_METRICS as METRICS_BOLTZ1              } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ1_MSA          } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ1_CONSTRAINED  } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2              } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2_MSA          } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2_CONSTRAINED  } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_CHAI1               } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_AF2M                } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_AF3                 } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_AF3_NOMSA           } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_COLABFOLD           } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_COLABFOLD_NOMSA     } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_ESMFOLD2            } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2_TEMPLATE             } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2_MSA_TEMPLATE         } from './modules/metrics'
include { COMPUTE_METRICS as METRICS_BOLTZ2_CONSTRAINED_TEMPLATE } from './modules/metrics'

include { AGGREGATE_RESULTS       } from './modules/aggregate'

// Helper: build the input tuple for a COMPUTE_METRICS alias
// Each alias consumes a queue channel of exactly one tuple:
//   [model_tag, parser_tag, prediction_dir,
//    reference_pdb, chain_a_len, chain_b_len, rec_chain, eff_chain]
//
// prediction_dir is a queue channel emitting a single directory path from
// the upstream predictor.  reference_pdb, chain lengths and chain IDs are
// value channels so they fan out correctly across all metric calls.
//
// NOTE: parser_tag is passed in explicitly rather than looked up from
// MODEL_TO_PARSER here, because script-level `def` bindings are not
// visible inside top-level function definitions in Groovy.  The workflow
// body (below) is a closure and CAN see MODEL_TO_PARSER, so the lookup
// happens at the call site.

def build_metric_input(model_tag, parser_tag, pred_dir_ch, ref_pdb_ch, ch_a_len_ch, ch_b_len_ch, rec_ch, eff_ch) {
    return pred_dir_ch.combine(ref_pdb_ch).combine(ch_a_len_ch).combine(ch_b_len_ch).map { tuple ->
        def (pdir, rpdb, alen, blen) = tuple
        return [model_tag, parser_tag, pdir, rpdb, alen, blen, rec_ch, eff_ch]
    }
}

// Workflow

workflow {

    validate_params()

    // Emits sequences.json and republishes the reference to a predictable
    // path so every COMPUTE_METRICS alias can find it.

    ref_pdb_file = Channel.fromPath(params.reference_pdb, checkIfExists: true)

    EXTRACT_SEQUENCES(
        ref_pdb_file,
        params.receptor_chain,
        params.effector_chain
    )
    def sequences_json_ch    = EXTRACT_SEQUENCES.out.sequences_json
    def reference_pdb_src_ch = EXTRACT_SEQUENCES.out.reference_pdb

    // One combined map avoids consuming the queue channel twice. .first()
    // promotes it to a value channel so every metric alias can read it.
    seqs_info_ch = sequences_json_ch.map { json_file ->
        def data = new groovy.json.JsonSlurper().parseText(json_file.text)
        // extract_sequences.py populates `chains` with id/sequence/length
        // entries — build a map keyed by chain ID so we can look up by
        // receptor_chain and effector_chain.
        def by_id = [:]
        data.chains.each { c -> by_id[c.id] = c }
        def rec = by_id[params.receptor_chain]
        def eff = by_id[params.effector_chain]
        if (!rec || !eff) {
            error "Chain ${params.receptor_chain}/${params.effector_chain} not found in sequences.json: available=${data.chain_ids}"
        }
        return [
            receptor_seq: rec.sequence,
            effector_seq: eff.sequence,
            chain_a_len : rec.length,
            chain_b_len : eff.length,
        ]
    }.first()

    receptor_seq_ch = seqs_info_ch.map { it.receptor_seq }
    effector_seq_ch = seqs_info_ch.map { it.effector_seq }
    chain_a_len_ch  = seqs_info_ch.map { it.chain_a_len }
    chain_b_len_ch  = seqs_info_ch.map { it.chain_b_len }

    // Reference PDB (real or placeholder) as a reusable value channel for
    // every metric alias.
    ref_pdb_ch = reference_pdb_src_ch.first()

    // Stage 0b, conditional: effector structural template, Run 2 only
    //
    // Extracted only when a *_template variant is selected.  No variant in the
    // model-vs-model comparison receives it, so nothing derived from the
    // reference structure reaches those predictions.

    def effector_template_cif_ch = Channel.empty()
    if (needs_eff_template) {
        EXTRACT_EFFECTOR_TEMPLATE(ref_pdb_ch, params.effector_chain)
        effector_template_cif_ch = EXTRACT_EFFECTOR_TEMPLATE.out.template_cif.first()
    }

    // Stage 1 (conditional): shared ColabFold MSA

    if (needs_shared_msa) {
        COLABFOLD_SEARCH(receptor_seq_ch, effector_seq_ch)
        shared_msa_dir_ch     = COLABFOLD_SEARCH.out.msa_dir.first()
        shared_a3m_a_ch       = COLABFOLD_SEARCH.out.chain_a_a3m.first()
        shared_a3m_b_ch       = COLABFOLD_SEARCH.out.chain_b_a3m.first()
        shared_complex_a3m_ch = COLABFOLD_SEARCH.out.complex_a3m.first()
    }

    // Stage 2 (conditional): AF3 database symlink farm

    if (needs_af3_db) {
        AF3_SETUP_DB()
        af3_db_dir_ch    = AF3_SETUP_DB.out.db_dir.first()
        af3_db_flag_ch   = AF3_SETUP_DB.out.ready_flag.first()
    }

    // Stage 3: run each selected model and its metrics alias
    //
    // Collect per-model [tagged_metrics, tagged_best_dir] channels into
    // two lists, then AGGREGATE_RESULTS concatenates them at the end.

    all_tagged_metrics = Channel.empty()
    all_tagged_best    = Channel.empty()

    // BOLTZ1
    if ('boltz1' in SELECTED_MODELS) {
        BOLTZ1(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain
        )
        METRICS_BOLTZ1(
            build_metric_input('boltz1', MODEL_TO_PARSER['boltz1'],
                BOLTZ1.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ1.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ1.out.tagged_best_dir)
    }

    // BOLTZ1_MSA
    if ('boltz1_msa' in SELECTED_MODELS) {
        BOLTZ1_MSA(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            shared_a3m_a_ch, shared_a3m_b_ch
        )
        METRICS_BOLTZ1_MSA(
            build_metric_input('boltz1_msa', MODEL_TO_PARSER['boltz1_msa'],
                BOLTZ1_MSA.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ1_MSA.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ1_MSA.out.tagged_best_dir)
    }

    // BOLTZ1_CONSTRAINED
    if ('boltz1_constrained' in SELECTED_MODELS) {
        BOLTZ1_CONSTRAINED(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            ref_pdb_ch
        )
        METRICS_BOLTZ1_CONSTRAINED(
            build_metric_input('boltz1_constrained', MODEL_TO_PARSER['boltz1_constrained'],
                BOLTZ1_CONSTRAINED.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ1_CONSTRAINED.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ1_CONSTRAINED.out.tagged_best_dir)
    }

    // BOLTZ2
    if ('boltz2' in SELECTED_MODELS) {
        BOLTZ2(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain
        )
        METRICS_BOLTZ2(
            build_metric_input('boltz2', MODEL_TO_PARSER['boltz2'],
                BOLTZ2.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2.out.tagged_best_dir)
    }

    // BOLTZ2_MSA
    if ('boltz2_msa' in SELECTED_MODELS) {
        BOLTZ2_MSA(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            shared_a3m_a_ch, shared_a3m_b_ch
        )
        METRICS_BOLTZ2_MSA(
            build_metric_input('boltz2_msa', MODEL_TO_PARSER['boltz2_msa'],
                BOLTZ2_MSA.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2_MSA.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2_MSA.out.tagged_best_dir)
    }

    // BOLTZ2_CONSTRAINED
    if ('boltz2_constrained' in SELECTED_MODELS) {
        BOLTZ2_CONSTRAINED(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            ref_pdb_ch
        )
        METRICS_BOLTZ2_CONSTRAINED(
            build_metric_input('boltz2_constrained', MODEL_TO_PARSER['boltz2_constrained'],
                BOLTZ2_CONSTRAINED.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2_CONSTRAINED.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2_CONSTRAINED.out.tagged_best_dir)
    }

    // CHAI1
    if ('chai1' in SELECTED_MODELS) {
        CHAI1(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain
        )
        METRICS_CHAI1(
            build_metric_input('chai1', MODEL_TO_PARSER['chai1'],
                CHAI1.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_CHAI1.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_CHAI1.out.tagged_best_dir)
    }

    // AF2M
    if ('af2m' in SELECTED_MODELS) {
        AF2M(receptor_seq_ch, effector_seq_ch)
        METRICS_AF2M(
            build_metric_input('af2m', MODEL_TO_PARSER['af2m'],
                AF2M.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_AF2M.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_AF2M.out.tagged_best_dir)
    }

    // AF3
    if ('af3' in SELECTED_MODELS) {
        AF3(
            receptor_seq_ch, effector_seq_ch,
            af3_db_dir_ch, af3_db_flag_ch
        )
        METRICS_AF3(
            build_metric_input('af3', MODEL_TO_PARSER['af3'],
                AF3.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_AF3.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_AF3.out.tagged_best_dir)
    }

    // AF3_NOMSA
    if ('af3_nomsa' in SELECTED_MODELS) {
        AF3_NOMSA(
            receptor_seq_ch, effector_seq_ch,
            af3_db_dir_ch, af3_db_flag_ch
        )
        METRICS_AF3_NOMSA(
            build_metric_input('af3_nomsa', MODEL_TO_PARSER['af3_nomsa'],
                AF3_NOMSA.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_AF3_NOMSA.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_AF3_NOMSA.out.tagged_best_dir)
    }

    // COLABFOLD
    if ('colabfold' in SELECTED_MODELS) {
        COLABFOLD(shared_complex_a3m_ch)
        METRICS_COLABFOLD(
            build_metric_input('colabfold', MODEL_TO_PARSER['colabfold'],
                COLABFOLD.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_COLABFOLD.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_COLABFOLD.out.tagged_best_dir)
    }

    // COLABFOLD_NOMSA
    if ('colabfold_nomsa' in SELECTED_MODELS) {
        COLABFOLD_NOMSA(receptor_seq_ch, effector_seq_ch)
        METRICS_COLABFOLD_NOMSA(
            build_metric_input('colabfold_nomsa', MODEL_TO_PARSER['colabfold_nomsa'],
                COLABFOLD_NOMSA.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_COLABFOLD_NOMSA.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_COLABFOLD_NOMSA.out.tagged_best_dir)
    }

    // ESMFOLD2
    // Single-sequence diffusion complex predictor; no MSA / no shared deps.
    if ('esmfold2' in SELECTED_MODELS) {
        ESMFOLD2(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain
        )
        METRICS_ESMFOLD2(
            build_metric_input('esmfold2', MODEL_TO_PARSER['esmfold2'],
                ESMFOLD2.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_ESMFOLD2.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_ESMFOLD2.out.tagged_best_dir)
    }

    // Run 2. Not part of the model comparison. Each pairs with its
    // template-free twin above, differing only in the supplied effector fold.

    // BOLTZ2_TEMPLATE
    if ('boltz2_template' in SELECTED_MODELS) {
        BOLTZ2_TEMPLATE(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            effector_template_cif_ch
        )
        METRICS_BOLTZ2_TEMPLATE(
            build_metric_input('boltz2_template', MODEL_TO_PARSER['boltz2_template'],
                BOLTZ2_TEMPLATE.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2_TEMPLATE.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2_TEMPLATE.out.tagged_best_dir)
    }

    // BOLTZ2_MSA_TEMPLATE
    if ('boltz2_msa_template' in SELECTED_MODELS) {
        BOLTZ2_MSA_TEMPLATE(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            shared_a3m_a_ch, shared_a3m_b_ch,
            effector_template_cif_ch
        )
        METRICS_BOLTZ2_MSA_TEMPLATE(
            build_metric_input('boltz2_msa_template', MODEL_TO_PARSER['boltz2_msa_template'],
                BOLTZ2_MSA_TEMPLATE.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2_MSA_TEMPLATE.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2_MSA_TEMPLATE.out.tagged_best_dir)
    }

    // BOLTZ2_CONSTRAINED_TEMPLATE
    if ('boltz2_constrained_template' in SELECTED_MODELS) {
        BOLTZ2_CONSTRAINED_TEMPLATE(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            ref_pdb_ch, effector_template_cif_ch
        )
        METRICS_BOLTZ2_CONSTRAINED_TEMPLATE(
            build_metric_input('boltz2_constrained_template',
                MODEL_TO_PARSER['boltz2_constrained_template'],
                BOLTZ2_CONSTRAINED_TEMPLATE.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_BOLTZ2_CONSTRAINED_TEMPLATE.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_BOLTZ2_CONSTRAINED_TEMPLATE.out.tagged_best_dir)
    }

    // Stage 4: Aggregate

    AGGREGATE_RESULTS(
        all_tagged_metrics.collect(),
        all_tagged_best.collect()
    )
}

// On completion

workflow.onComplete {
    // Runs in the hook rather than a process so trace.txt is guaranteed
    // complete. The parsing lives in Python so it can be unit-tested.
    def trace = file("${params.outdir}/trace.txt")
    def stats = file("${params.outdir}/predictor_runtime_stats.csv")

    if (!trace.exists()) {
        log.warn "trace.txt not found at ${trace}, skipping predictor_runtime_stats.csv"
    } else {
        def proc = ["python3", "${projectDir}/bin/trace_to_runtime_csv.py",
                    trace.toString(), stats.toString()].execute()
        proc.waitFor()
        if (proc.exitValue() != 0) {
            log.warn "trace_to_runtime_csv.py failed: ${proc.err.text.trim()}"
        } else {
            log.info proc.text.trim()
        }
    }

    log.info """
    =============================================================
    Benchmark Complete
    =============================================================
    Project:    ${params.project_name}
    Output:     ${params.outdir}
    Duration:   ${workflow.duration}
    Success:    ${workflow.success}

    Key outputs:
      ${params.outdir}/all_metrics.csv
      ${params.outdir}/all_metrics_ranked_by_effector_rmsd.csv
      ${params.outdir}/predictor_runtime_stats.csv
      ${params.outdir}/best_models/
      ${params.outdir}/trace.txt
      ${params.outdir}/pipeline_report.html
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
