#!/usr/bin/env nextflow

include { stage_as_plate as qc_stage_as_plate } from '../processes/intensity.nf'
include { qc_decon_samples; qc_render_report } from '../processes/qc.nf'

// Builds the self-contained QC report (params.qc_run). Tab 1 (general QC) is always
// produced from measure_intensity's output (guaranteed present - qc_run requires
// rn_cache_images, see checks.nf); tabs 2-6 each independently degrade to "not
// available" when their upstream stage didn't run, rather than failing the report.
//
// Registration/flatfield PNGs and scaling_index.tsv are read directly off the
// filesystem at their well-known publish locations (there's no clean Nextflow channel
// output for "all the PNGs a process wrote as a storeDir side effect") - the
// *_barrier values below exist purely to force qc_render_report to wait until those
// upstream stages have actually finished writing before it reads that directory.
workflow qc_report {

    take:
        measurements                 // tuple(well, [object/image_features.parquet]) - finalize_images.out.measurements
        manifest_file                // file(params.rn_manifest)
        manifest_registration        // channel of RegistrationRecord maps, or Channel.empty()
        manifest_registration_file   // file(params.rn_manifest_registration), or NO_REGISTRATION_MANIFEST
        has_registration              // bool: params.rn_manifest_registration != null
        blacklist_file                // file(params.rn_blacklist), or NO_BLACKLIST
        well_channel                  // tuple(plate, Well) - setup.out.well_channel
        flatfield_out                 // flatfield_estimation.out.flatfield_out (barrier only), or Channel.empty() if !ff_run
        decon_out                     // tuple(well, manifest, path) from deconvolute, or Channel.empty() if !dc_run
        registration_out              // tuple(well, RegistrationRecord, path) from register, or Channel.empty() if !has_registration
        scaling_index_file            // estimate_scaling_factors.out.scaling_index_file (real file, or NO_SCALING_INDEX)

    main:
        //------------------------------------------------------------------------
        // Tab 1/2/5 data source: aggregate measure_intensity's per-well parquet
        // output into per-plate dirs (same pattern finalize_images.nf uses for
        // sc_autoscale, aliased since a process can only be invoked once per
        // workflow scope)
        //------------------------------------------------------------------------
        plate_measurements = qc_stage_as_plate(
            measurements
                .map { well, files -> tuple(well.plate, files) }
                .groupTuple(by: 0)
                .map { plate, files -> tuple(plate, files.flatten()) }
        )
        measurements_dir = plate_measurements.map { plate, dir -> dir }.collect()

        //------------------------------------------------------------------------
        // Tab 2: registration sample images
        //------------------------------------------------------------------------
        if (has_registration) {
            registration_images_dir = file(params.rn_publish_dir + "/registration")
            registration_barrier = registration_out.collect()
        } else {
            registration_images_dir = file("NO_REGISTRATION")
            registration_barrier = Channel.value(true)
        }

        //------------------------------------------------------------------------
        // Tab 3: flatfields
        //------------------------------------------------------------------------
        show_flatfield = params.ff_run
        if (show_flatfield) {
            flatfields_dir = file(params.rn_publish_dir + "/flatfields")
            flatfield_barrier = flatfield_out.collect()
        } else {
            flatfields_dir = file("NO_FLATFIELDS")
            flatfield_barrier = Channel.value(true)
        }

        ff_params_json = groovy.json.JsonOutput.toJson([
            ff_mode: params.ff_mode,
            ff_degree: params.ff_degree,
            ff_nimg: params.ff_nimg,
            ff_merge_n: params.ff_merge_n,
            ff_global_flatfield: params.ff_global_flatfield,
            ff_pseudoreplicates: params.ff_pseudoreplicates,
            ff_no_tune: params.ff_no_tune,
            ff_threshold: params.ff_threshold,
            ff_autosegment: params.ff_autosegment,
            ff_use_ridge: params.ff_use_ridge,
        ])

        //------------------------------------------------------------------------
        // Tab 4: decon before/after samples
        //
        // Deterministic (stride-based, no RNG) sample of qc_n_sample_decon wells -
        // keeps -resume stable and avoids needing a fixed seed across Groovy/Python.
        // With registration, samples are drawn from the reference plate's wells only
        // (using the first registration group when several exist - a QC convenience
        // sample, not exhaustive coverage). Without registration, samples are drawn
        // across all plates.
        //------------------------------------------------------------------------
        show_decon = params.dc_run
        if (show_decon) {
            if (has_registration) {
                reg_first = manifest_registration.first()
                sample_candidates = well_channel
                    .combine(reg_first.map { r -> tuple(r.ref_plate, r.qry_plates) })
                    .filter { plate, well, ref_plate, qry_plates -> plate == ref_plate }
                    .map { plate, well, ref_plate, qry_plates -> tuple(well.well, ref_plate, qry_plates) }
                    .unique { it[0] }
                    .toSortedList { a, b -> a[0] <=> b[0] }
            } else {
                sample_candidates = well_channel
                    .map { plate, well -> tuple(well.key, plate, well.well) }
                    .unique { it[0] }
                    .map { key, plate, well_id -> tuple(well_id, plate, []) }
                    .toSortedList { a, b -> (a[1] + a[0]) <=> (b[1] + b[0]) }
            }

            decon_sample_tuples = sample_candidates.flatMap { candidates ->
                def n = Math.min(params.qc_n_sample_decon as int, candidates.size())
                if (n == 0) return []
                def stride = Math.max(1, (candidates.size() / n) as int)
                (0..<n).collect { i -> candidates[Math.min(i * stride, candidates.size() - 1)] }
            }

            decon_barrier = decon_out.collect()
            decon_images = qc_decon_samples(
                decon_sample_tuples,
                file(params.rn_image_dir),
                file(params.rn_decon_dir),
                registration_images_dir,
                decon_barrier,
            ).images.collect()
        } else {
            // Empty list + stageAs glob in qc_render_report means the "decon_samples"
            // dir simply won't exist - same convention already used elsewhere for a
            // process input that's a no-op in this branch (see e.g. run_pipeline.nf's
            // estimate_scaling_factors(cellpose_out.last(), Channel.value([]), ...)).
            decon_images = Channel.value([])
        }

        //------------------------------------------------------------------------
        // Tab 6: scaling factors - scaling_index_file is already a real channel
        // dependency (no barrier needed)
        //------------------------------------------------------------------------
        show_scaling = params.sc_autoscale

        //------------------------------------------------------------------------
        // Render
        //------------------------------------------------------------------------
        qc_render_report(
            measurements_dir,
            manifest_file,
            blacklist_file,
            manifest_registration_file,
            registration_barrier,
            registration_images_dir,
            show_flatfield,
            flatfield_barrier,
            flatfields_dir,
            params.ff_global_flatfield,
            ff_params_json,
            show_decon,
            decon_images,
            show_scaling,
            scaling_index_file,
        )

    emit:
        qc_report_html = qc_render_report.out
}
