#!/usr/bin/env nextflow

// Load workflows
// Stage
include { stage } from './workflows/stage.nf'

// Run pipeline
include { run_pipeline } from './workflows/run_pipeline.nf'

// Nextflow's strict syntax parser no longer supports selecting an entry workflow via
// the `-entry` CLI flag - dispatch on a param instead. Run with `--workflow stage` or
// `--workflow run_pipeline` (see README.md).
workflow {
    if (params.workflow == "stage") {
        stage()
    } else if (params.workflow == "run_pipeline") {
        run_pipeline()
    } else {
        error "Unknown or missing --workflow '${params.workflow}'. Must be one of: stage, run_pipeline"
    }
}
