#!/usr/bin/env nextflow
include { validateParameters } from 'plugin/nf-schema'

// Validate parameters plugin
validateParameters()

// Load workflows
// Stage
include { stage } from './workflows/stage.nf'

// Run pipeline
include { run_pipeline } from './workflows/run_pipeline.nf'

// Run subcell
include { run_subcell } from './workflows/run_subcell.nf'

