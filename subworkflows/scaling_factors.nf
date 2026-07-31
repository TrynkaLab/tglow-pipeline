#!/usr/bin/env nextflow

include { calculate_scaling_factors } from '../processes/scaling.nf'


workflow estimate_scaling_factors {

    take:
        images_ready
        plate_dirs
        rn_autoscale
        rn_manualscale
        rn_scale_slope
        rn_scale_bias
        rn_control_list
        rn_channel_map
    main:

        // Get scaling string
        if (rn_manualscale != null) {
            scaling_channel = Channel.value(file(rn_manualscale))
        } else {
            scaling_channel = Channel.value(file("NO_SCALE"))
        }

        // Optional slope for the sigmoid curve
        if (rn_scale_slope != null) {
            slope_channel = Channel.value(file(rn_scale_slope))
        } else {
            slope_channel = Channel.value(file("NO_SLOPE"))
        }

        // Optional bias for the sigmoid curve
        if (rn_scale_bias != null) {
            bias_channel = Channel.value(file(rn_scale_bias))
        } else {
            bias_channel = Channel.value(file("NO_BIAS"))
        }

        // Autoscale derives factors (and sigmoid bias/slope) from the measured intensities -
        // requires rn_control_list and rn_channel_map (enforced in checks.nf)
        if (rn_autoscale) {
            calc_out = calculate_scaling_factors(
                images_ready,
                plate_dirs,
                Channel.value(file(rn_channel_map)),
                Channel.value(file(rn_control_list))
            )

            scaling_channel = calc_out.scaling_factors.first()

            // Autoscale fits its own sigmoid bias/slope from the control images - use those
            // unless the user explicitly supplied manual overrides
            if (rn_scale_slope == null && rn_scale_bias == null) {
                slope_channel = calc_out.scaling_slope.first()
                bias_channel = calc_out.scaling_bias.first()
            }
        }
    emit:
        scaling_file=scaling_channel
        slope_file=slope_channel
        bias_file=bias_channel
}
