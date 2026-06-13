#!/usr/bin/env nextflow

/*
 * =============================================================================
 * Structure Prediction Benchmark Pipeline (v0.1.0)
 * =============================================================================
 *
 * Benchmarks multiple protein complex structure predictors against a
 * reference PDB:
 *
 *   Boltz-1, Boltz-1 + MSA, Boltz-1 + pocket constraint
 *   Boltz-2, Boltz-2 + MSA, Boltz-2 + pocket/contact constraints
 *   Chai-1
 *   AlphaFold2-Multimer (with MSA)
 *   AlphaFold 3, AlphaFold 3 (no MSA)
 *   ColabFold, ColabFold (no MSA)
 *   ESMFold2 (single-sequence; diffusion complex predictor)
 *
 * Each model writes predictions + metrics to ${outdir}/<model>/.  The
 * shared ColabFold MSA is computed once and reused by every model that
 * needs it (boltz1_msa, boltz2_msa, colabfold).  All seed loops inside
 * each model are serial to respect the small jic-gpu queue at NBI/JIC.
 *
 * Resource accounting is handled natively by Nextflow's trace / report /
 * timeline outputs — no custom collect_slurm_stats.sh step is needed.
 *
 * Author: Josh Williams
 * Institute: John Innes Centre / The Sainsbury Laboratory
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults
// ---------------------------------------------------------------------------

// ── Project ─────────────────────────────────────────────────────────────
params.project_name     = "structure_prediction_benchmark"
params.outdir           = "${launchDir}/${params.project_name}_results"

// ── Inputs ──────────────────────────────────────────────────────────────
// Two mutually-exclusive input modes:
//
//   PDB mode (original):  --reference_pdb + --receptor_chain + --effector_chain
//                         Provides an experimental reference complex; the
//                         pipeline computes structural RMSD/DockQ against it
//                         in addition to confidence-based metrics.
//
//   FASTA mode (new):     --input_fasta  (a 2-entry protein FASTA: receptor
//                         then effector).  No experimental reference structure
//                         is required; structural metrics are skipped and
//                         only confidence metrics (pLDDT/ipTM/pTM/...) are
//                         emitted.  --receptor_chain / --effector_chain are
//                         optional in this mode and default to A / B — they
//                         become the synthetic chain IDs assigned to the two
//                         FASTA entries when handed to the predictors.
//
params.reference_pdb    = null      // PDB mode: two-chain reference complex
params.input_fasta      = null      // FASTA mode: two-entry FASTA file
params.receptor_chain   = null      // PDB mode: receptor chain ID (required)
                                    // FASTA mode: synthetic chain ID (default 'A')
params.effector_chain   = null      // PDB mode: effector chain ID (required)
                                    // FASTA mode: synthetic chain ID (default 'B')

// ── Model selection ─────────────────────────────────────────────────────
// A list of models to run; null means 'run them all'.
params.models           = null

// ── Infrastructure paths ────────────────────────────────────────────────
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

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------
//
// Exactly one of --reference_pdb / --input_fasta must be provided.  In PDB
// mode the chain IDs are required (they pin which chains in the reference
// complex are receptor vs effector).  In FASTA mode chain IDs default to
// A/B and act as synthetic labels for the two FASTA entries.

if (params.reference_pdb && params.input_fasta) {
    error "Provide either --reference_pdb (PDB mode) OR --input_fasta " +
          "(FASTA mode), not both."
}
if (!params.reference_pdb && !params.input_fasta) {
    error "Please provide either --reference_pdb (PDB mode) or " +
          "--input_fasta (FASTA mode)."
}

def INPUT_MODE = params.reference_pdb ? 'pdb' : 'fasta'

if (INPUT_MODE == 'pdb') {
    if (!params.receptor_chain || !params.effector_chain) {
        error "PDB mode requires --receptor_chain and --effector_chain " +
              "(PDB chain IDs)."
    }
} else {
    // FASTA mode: chain IDs are optional synthetic labels.
    if (!params.receptor_chain) { params.receptor_chain = 'A' }
    if (!params.effector_chain) { params.effector_chain = 'B' }
    if (params.receptor_chain == params.effector_chain) {
        error "FASTA mode: --receptor_chain and --effector_chain must differ " +
              "(both were '${params.receptor_chain}')."
    }
}

// ---------------------------------------------------------------------------
// Model selection
// ---------------------------------------------------------------------------

// The full catalogue of models this pipeline can run.  Parser tag names the
// underlying output format so compute_metrics.py picks the right parser;
// variants that share a parser (e.g. boltz1_msa → boltz1) keep their own
// model_tag in the output CSV.
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
    'esmfold2',            // parser: esmfold2 (single-sequence diffusion complex predictor)
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
]

def SELECTED_MODELS
if (params.models == null) {
    SELECTED_MODELS = ALL_MODELS
} else {
    // params.models may arrive as a List (from YAML) or a whitespace/comma
    // -separated string (from the CLI).  Normalise to a List<String>.
    def raw = params.models
    def requested
    if (raw instanceof List) {
        requested = raw.collect { it.toString().trim() }
    } else {
        requested = raw.toString().replaceAll(',', ' ').split()*.trim()
    }
    requested = requested.findAll { it }  // drop empty strings

    def unknown = requested.findAll { !(it in ALL_MODELS) }
    if (unknown) {
        error "Unknown model(s) in --models: ${unknown.join(', ')}\nValid: ${ALL_MODELS.join(', ')}"
    }
    SELECTED_MODELS = requested
}

def needs_shared_msa      = SELECTED_MODELS.any { it in ['boltz1_msa', 'boltz2_msa', 'colabfold'] }
def needs_af3_db          = SELECTED_MODELS.any { it in ['af3', 'af3_nomsa'] }
def needs_eff_template    = INPUT_MODE == 'pdb' &&
                            SELECTED_MODELS.any { it in ['af3_nomsa', 'colabfold', 'colabfold_nomsa',
                                                          'boltz2', 'boltz2_msa'] }

log.info """
=============================================================
Structure Prediction Benchmark — v0.1.0
=============================================================
Project:         ${params.project_name}
Input mode:      ${INPUT_MODE}
${ INPUT_MODE == 'pdb' \
    ? "Reference PDB:   ${params.reference_pdb}" \
    : "Input FASTA:     ${params.input_fasta}" }
Receptor chain:  ${params.receptor_chain}${ INPUT_MODE == 'fasta' ? '  (synthetic label for FASTA entry 1)' : '' }
Effector chain:  ${params.effector_chain}${ INPUT_MODE == 'fasta' ? '  (synthetic label for FASTA entry 2)' : '' }
Output:          ${params.outdir}
Models:          ${SELECTED_MODELS.join(', ')}
Shared MSA:      ${needs_shared_msa ? 'yes' : 'no'}
AF3 DB setup:    ${needs_af3_db ? 'yes' : 'no'}
${ INPUT_MODE == 'fasta' ? 'NOTE: FASTA mode — structural RMSD/DockQ skipped, confidence metrics only.' : '' }
=============================================================
""".stripIndent()

// ---------------------------------------------------------------------------
// Include modules
// ---------------------------------------------------------------------------

include { EXTRACT_SEQUENCES            } from './modules/preprocessing'
include { EXTRACT_EFFECTOR_TEMPLATE    } from './modules/preprocessing'
include { EXTRACT_SEQUENCES_FROM_FASTA } from './modules/preprocessing_fasta'
include { COLABFOLD_SEARCH             } from './modules/msa'

include { BOLTZ1                  } from './modules/boltz1'
include { BOLTZ1_MSA              } from './modules/boltz1'
include { BOLTZ1_CONSTRAINED      } from './modules/boltz1'

include { BOLTZ2                  } from './modules/boltz2'
include { BOLTZ2_MSA              } from './modules/boltz2'
include { BOLTZ2_CONSTRAINED      } from './modules/boltz2'

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

include { AGGREGATE_RESULTS       } from './modules/aggregate'

// ---------------------------------------------------------------------------
// Helper: build the input tuple for a COMPUTE_METRICS alias
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // =====================================================================
    // Stage 0: extract sequences from reference PDB OR from input FASTA
    // =====================================================================
    //
    // Both preprocessing processes emit identically-shaped outputs:
    //   .sequences_json — JSON with chains[].id/sequence/length
    //   .reference_pdb  — either the real reference (PDB mode) or a tiny
    //                     placeholder containing only REMARK lines (FASTA
    //                     mode).  compute_metrics.py treats the placeholder
    //                     as "no reference available" and skips structural
    //                     RMSD/DockQ while still emitting confidence metrics.
    //
    // The downstream stages consume `sequences_json_ch` and
    // `reference_pdb_src_ch` regardless of which branch ran, so no other
    // part of the workflow needs to know about the input mode.

    def sequences_json_ch
    def reference_pdb_src_ch

    if (INPUT_MODE == 'pdb') {
        ref_pdb_file = Channel.fromPath(params.reference_pdb, checkIfExists: true)

        EXTRACT_SEQUENCES(
            ref_pdb_file,
            params.receptor_chain,
            params.effector_chain
        )
        sequences_json_ch    = EXTRACT_SEQUENCES.out.sequences_json
        reference_pdb_src_ch = EXTRACT_SEQUENCES.out.reference_pdb
    } else {
        input_fasta_file = Channel.fromPath(params.input_fasta, checkIfExists: true)

        EXTRACT_SEQUENCES_FROM_FASTA(
            input_fasta_file,
            params.receptor_chain,
            params.effector_chain
        )
        sequences_json_ch    = EXTRACT_SEQUENCES_FROM_FASTA.out.sequences_json
        reference_pdb_src_ch = EXTRACT_SEQUENCES_FROM_FASTA.out.reference_pdb
    }

    // Parse sequences.json into value channels.  One combined map step
    // avoids consuming the queue channel twice; .first() promotes the
    // result into a value channel so downstream processes (each metric
    // alias) can read it repeatedly.
    seqs_info_ch = sequences_json_ch.map { json_file ->
        def data = new groovy.json.JsonSlurper().parseText(json_file.text)
        // Both extract_sequences.py (PDB) and extract_sequences_from_fasta.py
        // (FASTA) populate `chains` with id/sequence/length entries — build a
        // map keyed by chain ID so we can look up by receptor_chain and
        // effector_chain regardless of input mode.
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
    // every metric alias.  In FASTA mode this is a comment-only stub and
    // compute_metrics.py will fall through to confidence-only metrics.
    ref_pdb_ch = reference_pdb_src_ch.first()

    // =====================================================================
    // Stage 0b (conditional): extract effector structural template
    // =====================================================================
    //
    // PDB mode only, and only when AF3 or ColabFold variants are selected.
    // AF3 uses the mmCIF form (injected into the JSON template block).
    // ColabFold uses the PDB form (--custom-template-path).
    // In FASTA mode or when no template-capable model is selected, the
    // sentinel file (empty, zero bytes) is passed instead — processes check
    // with `[ -s ]` and skip template injection when the file is empty.
    //
    // AF2M note: run_alphafold.py has no custom-template-path flag; AF2M
    // templates come from its built-in database search (cutoff 2020-05-14).

    def no_template_sentinel = Channel.value(file("${projectDir}/bin/no_template.sentinel"))
    def effector_template_pdb_ch = no_template_sentinel
    def effector_template_cif_ch = no_template_sentinel

    if (needs_eff_template) {
        EXTRACT_EFFECTOR_TEMPLATE(ref_pdb_ch, params.effector_chain)
        effector_template_pdb_ch = EXTRACT_EFFECTOR_TEMPLATE.out.template_pdb.first()
        effector_template_cif_ch = EXTRACT_EFFECTOR_TEMPLATE.out.template_cif.first()
    }

    // =====================================================================
    // Stage 1 (conditional): shared ColabFold MSA
    // =====================================================================

    if (needs_shared_msa) {
        COLABFOLD_SEARCH(receptor_seq_ch, effector_seq_ch)
        shared_msa_dir_ch     = COLABFOLD_SEARCH.out.msa_dir.first()
        shared_a3m_a_ch       = COLABFOLD_SEARCH.out.chain_a_a3m.first()
        shared_a3m_b_ch       = COLABFOLD_SEARCH.out.chain_b_a3m.first()
        shared_complex_a3m_ch = COLABFOLD_SEARCH.out.complex_a3m.first()
    }

    // =====================================================================
    // Stage 2 (conditional): AF3 database symlink farm
    // =====================================================================

    if (needs_af3_db) {
        AF3_SETUP_DB()
        af3_db_dir_ch    = AF3_SETUP_DB.out.db_dir.first()
        af3_db_flag_ch   = AF3_SETUP_DB.out.ready_flag.first()
    }

    // =====================================================================
    // Stage 3: run each selected model and its metrics alias
    // =====================================================================
    //
    // Collect per-model [tagged_metrics, tagged_best_dir] channels into
    // two lists, then AGGREGATE_RESULTS concatenates them at the end.

    all_tagged_metrics = Channel.empty()
    all_tagged_best    = Channel.empty()

    // ── BOLTZ1 ──────────────────────────────────────────────────────────
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

    // ── BOLTZ1_MSA ──────────────────────────────────────────────────────
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

    // ── BOLTZ1_CONSTRAINED ──────────────────────────────────────────────
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

    // ── BOLTZ2 ──────────────────────────────────────────────────────────
    if ('boltz2' in SELECTED_MODELS) {
        BOLTZ2(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            effector_template_cif_ch
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

    // ── BOLTZ2_MSA ──────────────────────────────────────────────────────
    if ('boltz2_msa' in SELECTED_MODELS) {
        BOLTZ2_MSA(
            receptor_seq_ch, effector_seq_ch,
            params.receptor_chain, params.effector_chain,
            shared_a3m_a_ch, shared_a3m_b_ch,
            effector_template_cif_ch
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

    // ── BOLTZ2_CONSTRAINED ──────────────────────────────────────────────
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

    // ── CHAI1 ───────────────────────────────────────────────────────────
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

    // ── AF2M ────────────────────────────────────────────────────────────
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

    // ── AF3 ─────────────────────────────────────────────────────────────
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

    // ── AF3_NOMSA ───────────────────────────────────────────────────────
    if ('af3_nomsa' in SELECTED_MODELS) {
        AF3_NOMSA(
            receptor_seq_ch, effector_seq_ch,
            af3_db_dir_ch, af3_db_flag_ch,
            effector_template_cif_ch
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

    // ── COLABFOLD ───────────────────────────────────────────────────────
    if ('colabfold' in SELECTED_MODELS) {
        COLABFOLD(shared_complex_a3m_ch, effector_template_pdb_ch)
        METRICS_COLABFOLD(
            build_metric_input('colabfold', MODEL_TO_PARSER['colabfold'],
                COLABFOLD.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_COLABFOLD.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_COLABFOLD.out.tagged_best_dir)
    }

    // ── COLABFOLD_NOMSA ─────────────────────────────────────────────────
    if ('colabfold_nomsa' in SELECTED_MODELS) {
        COLABFOLD_NOMSA(receptor_seq_ch, effector_seq_ch, effector_template_pdb_ch)
        METRICS_COLABFOLD_NOMSA(
            build_metric_input('colabfold_nomsa', MODEL_TO_PARSER['colabfold_nomsa'],
                COLABFOLD_NOMSA.out.prediction_dir,
                ref_pdb_ch, chain_a_len_ch, chain_b_len_ch,
                params.receptor_chain, params.effector_chain)
        )
        all_tagged_metrics = all_tagged_metrics.mix(METRICS_COLABFOLD_NOMSA.out.tagged_metrics)
        all_tagged_best    = all_tagged_best.mix(METRICS_COLABFOLD_NOMSA.out.tagged_best_dir)
    }

    // ── ESMFOLD2 ────────────────────────────────────────────────────────
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

    // =====================================================================
    // Stage 4: Aggregate
    // =====================================================================

    AGGREGATE_RESULTS(
        all_tagged_metrics.collect(),
        all_tagged_best.collect()
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    // ── Write predictor_runtime_stats.csv ──────────────────────────────
    // Parse trace.txt and emit a filtered CSV with only the predictor
    // processes — no preprocessing, no metrics, no aggregation.  Produces
    // the same shape of information collect_slurm_stats.sh used to emit
    // in the bash pipeline, derived purely from Nextflow's observer data.
    //
    // Runs in the onComplete hook rather than a Nextflow process so it
    // has guaranteed read access to trace.txt after every task finishes.
    def PREDICTOR_NAMES = [
        'BOLTZ1', 'BOLTZ1_MSA', 'BOLTZ1_CONSTRAINED',
        'BOLTZ2', 'BOLTZ2_MSA', 'BOLTZ2_CONSTRAINED',
        'CHAI1', 'AF2M',
        'AF3', 'AF3_NOMSA',
        'COLABFOLD', 'COLABFOLD_NOMSA',
        'ESMFOLD2',
    ] as Set

    // Models whose standalone wall-clock time must include the shared
    // COLABFOLD_SEARCH MSA step — they cannot run without it in practice.
    def MSA_MODELS = ['BOLTZ1_MSA', 'BOLTZ2_MSA', 'COLABFOLD'] as Set

    // Parse Nextflow-formatted memory strings like "6.5 GB" / "120 MB"
    // into GB floats.  Returns null on unrecognised / missing values.
    def parse_mem_gb = { String s ->
        if (!s || s == '-' || s == '0') return null
        def m = (s.trim() =~ /^([\d.]+)\s*([KMGT]?B)$/)
        if (!m) return null
        def val  = m[0][1] as double
        def unit = m[0][2]
        switch (unit) {
            case 'B':  return val / (1024d * 1024d * 1024d)
            case 'KB': return val / (1024d * 1024d)
            case 'MB': return val / 1024d
            case 'GB': return val
            case 'TB': return val * 1024d
        }
        return null
    }

    // Parse Nextflow duration strings like "1h 23m 45s", "45.3s", "2m 10s".
    // Returns (elapsed_s as long, hms as String).  Nextflow emits the
    // realtime column as milliseconds when the DSL2 observer writes it,
    // but falls back to hms strings in some paths; handle both.
    def parse_duration = { String s ->
        if (!s || s == '-') return [null, null]
        s = s.trim()
        // Plain integer milliseconds
        if (s ==~ /^\d+$/) {
            long ms = s as long
            long secs = (ms / 1000d) as long
            long h = secs / 3600
            long m = (secs % 3600) / 60
            long sec = secs % 60
            return [secs, String.format('%02d:%02d:%02d', h, m, sec)]
        }
        // Human-readable: match any combination of  Nh Nm N(.N)s  Nms
        long total_ms = 0
        def patterns = [
            (~/(\d+)ms/)          : 1L,
            (~/(\d+(?:\.\d+)?)s/) : 1000L,
            (~/(\d+)m(?!s)/)      : 60_000L,
            (~/(\d+)h/)           : 3_600_000L,
        ]
        boolean any = false
        patterns.each { pat, mult ->
            def mm = (s =~ pat)
            while (mm.find()) {
                any = true
                total_ms += ((mm.group(1) as double) * mult) as long
            }
        }
        if (!any) return [null, null]
        long secs = (total_ms / 1000d) as long
        long h = secs / 3600
        long m = (secs % 3600) / 60
        long sec = secs % 60
        return [secs, String.format('%02d:%02d:%02d', h, m, sec)]
    }

    // Format seconds into HH:MM:SS without going through parse_duration.
    def fmt_hms = { long secs ->
        long h   = secs / 3600
        long m   = (secs % 3600) / 60
        long sec = secs % 60
        String.format('%02d:%02d:%02d', h, m, sec)
    }

    def trace_path = file("${params.outdir}/trace.txt")
    def out_path   = file("${params.outdir}/predictor_runtime_stats.csv")

    if (!trace_path.exists()) {
        log.warn "trace.txt not found at ${trace_path} — skipping predictor_runtime_stats.csv"
    } else {
        def lines = trace_path.readLines()
        if (lines.size() < 2) {
            log.warn "trace.txt has no data rows — skipping predictor_runtime_stats.csv"
        } else {
            def header = lines[0].split('\t') as List
            def col = { String name -> header.indexOf(name) }

            def i_process  = col('process')
            def i_status   = col('status')
            def i_exit     = col('exit')
            def i_realtime = col('realtime')
            def i_cpus     = col('cpus')      // may not be present
            def i_rss      = col('rss')
            def i_vmem     = col('vmem')
            def i_peak_rss = col('peak_rss')
            def i_peak_vm  = col('peak_vmem')
            def i_pcpu     = col('%cpu')
            def i_queue    = col('queue')

            // First pass: collect COLABFOLD_SEARCH wall-clock time so we can
            // add it to MSA-dependent predictor runtimes (standalone_elapsed_s).
            // If COLABFOLD_SEARCH ran multiple times, take the max.
            long msa_elapsed_s = 0L
            lines.drop(1).each { String msa_line ->
                def msa_fields = msa_line.split('\t', -1) as List
                if (i_process < 0 || i_process >= msa_fields.size()) return
                if (msa_fields[i_process] != 'COLABFOLD_SEARCH') return
                def rt = i_realtime >= 0 && i_realtime < msa_fields.size() ? msa_fields[i_realtime] : ''
                def (s, _hms) = parse_duration(rt)
                if (s != null && s > msa_elapsed_s) msa_elapsed_s = s
            }

            def rows_out = []
            rows_out << [
                'model', 'status', 'exit_code', 'queue',
                'elapsed_hms', 'elapsed_s',
                'standalone_elapsed_hms', 'standalone_elapsed_s',
                'pct_cpu',
                'rss_gb', 'vmem_gb', 'peak_rss_gb', 'peak_vmem_gb',
            ].join(',')

            lines.drop(1).each { String line ->
                def fields = line.split('\t', -1) as List
                if (i_process < 0 || i_process >= fields.size()) return
                def proc_name = fields[i_process]
                if (!(proc_name in PREDICTOR_NAMES)) return

                def status = i_status >= 0 && i_status < fields.size() ? fields[i_status] : ''
                def exit_c = i_exit   >= 0 && i_exit   < fields.size() ? fields[i_exit]   : ''
                def queue  = i_queue  >= 0 && i_queue  < fields.size() ? fields[i_queue]  : ''

                def realtime_raw = i_realtime >= 0 && i_realtime < fields.size() ? fields[i_realtime] : ''
                def (elapsed_s, elapsed_hms) = parse_duration(realtime_raw)

                def pct_cpu = i_pcpu >= 0 && i_pcpu < fields.size() ? fields[i_pcpu] : ''

                def rss_gb      = parse_mem_gb(i_rss      >= 0 && i_rss      < fields.size() ? fields[i_rss]      : '')
                def vmem_gb     = parse_mem_gb(i_vmem     >= 0 && i_vmem     < fields.size() ? fields[i_vmem]     : '')
                def peak_rss_gb = parse_mem_gb(i_peak_rss >= 0 && i_peak_rss < fields.size() ? fields[i_peak_rss] : '')
                def peak_vm_gb  = parse_mem_gb(i_peak_vm  >= 0 && i_peak_vm  < fields.size() ? fields[i_peak_vm]  : '')

                def fmt = { Double v -> v == null ? '' : String.format('%.2f', v) }

                // Standalone time = GPU time + MSA time (for MSA-dependent models)
                def standalone_s   = (elapsed_s != null)
                    ? (elapsed_s + ((proc_name in MSA_MODELS) ? msa_elapsed_s : 0L))
                    : null
                def standalone_hms = standalone_s != null ? fmt_hms(standalone_s) : ''

                rows_out << [
                    proc_name.toLowerCase(),
                    status,
                    exit_c,
                    queue,
                    elapsed_hms ?: '',
                    elapsed_s   ?: '',
                    standalone_hms,
                    standalone_s != null ? standalone_s.toString() : '',
                    pct_cpu?.replace('%', '') ?: '',
                    fmt(rss_gb),
                    fmt(vmem_gb),
                    fmt(peak_rss_gb),
                    fmt(peak_vm_gb),
                ].join(',')
            }

            if (rows_out.size() > 1) {
                out_path.text = rows_out.join('\n') + '\n'
                log.info "Wrote ${rows_out.size() - 1} predictor rows to ${out_path}"
            } else {
                log.warn "No predictor rows found in trace.txt — predictor_runtime_stats.csv not written"
            }
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
      ${params.outdir}/all_metrics.csv                           (prediction quality)
      ${params.outdir}/all_metrics_ranked_by_effector_rmsd.csv   (sorted by eff RMSD)
      ${params.outdir}/predictor_runtime_stats.csv               (SLURM resource usage)
      ${params.outdir}/best_models/
      ${params.outdir}/trace.txt                                 (full Nextflow trace)
      ${params.outdir}/pipeline_report.html                      (HTML dashboard)
      ${params.outdir}/timeline.html                             (Gantt chart)
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
