#!/usr/bin/env nextflow

include { calculate_scaling_factors } from '../processes/scaling.nf'


workflow estimate_scaling_factors {

    take:
        images_ready
        plate_dirs
        sc_autoscale
        sc_manualscale
        sc_scale_slope
        sc_scale_bias
        sc_control_list
        sc_channel_map
    main:

        // Get scaling string
        if (sc_manualscale != null) {
            scaling_channel = Channel.value(file(sc_manualscale))
        } else {
            scaling_channel = Channel.value(file("NO_SCALE"))
        }

        // Optional slope for the sigmoid curve
        if (sc_scale_slope != null) {
            slope_channel = Channel.value(file(sc_scale_slope))
        } else {
            slope_channel = Channel.value(file("NO_SLOPE"))
        }

        // Optional bias for the sigmoid curve
        if (sc_scale_bias != null) {
            bias_channel = Channel.value(file(sc_scale_bias))
        } else {
            bias_channel = Channel.value(file("NO_BIAS"))
        }

        // scaling_index.tsv (the full per-plate/channel table behind scaling_factors.txt -
        // used by the QC report's Tab 6) only exists on the autoscale path; manual scaling
        // never produces it.
        scaling_index_channel = Channel.value(file("NO_SCALING_INDEX"))

        // Autoscale derives factors (and sigmoid bias/slope) from the measured intensities.
        // sc_channel_map/sc_control_list are both optional - when absent, calculate_scaling_factors
        // (via processes/scaling.nf) falls back to dynamic-range-only scaling from sc_autoscale_q1/q2.
        if (sc_autoscale) {
            calc_out = calculate_scaling_factors(
                images_ready,
                plate_dirs,
                Channel.value(file(sc_channel_map ?: "NO_CHANNEL_MAP")),
                Channel.value(file(sc_control_list ?: "NO_CONTROLS"))
            )

            scaling_channel = calc_out.scaling_factors.first()
            scaling_index_channel = calc_out.scaling_index.first()

            // Autoscale fits its own sigmoid bias/slope from the control images - use those
            // unless the user explicitly supplied manual overrides
            if (sc_scale_slope == null && sc_scale_bias == null) {
                slope_channel = calc_out.scaling_slope.first()
                bias_channel = calc_out.scaling_bias.first()
            }
        }
    emit:
        scaling_file=scaling_channel
        slope_file=slope_channel
        bias_file=bias_channel
        scaling_index_file=scaling_index_channel
}
