include { validateManifestContent; validateBlacklistContent; validateRegistrationManifestContent; validateControlListContent; formatValidationErrors } from './validateInputs.nf'

// Sanity check the run parameters, erroring out early on invalid combinations
def checkParamsMain(params) {


    // Basic IO checks for input files
    checkParamsBase(params)

    // Checks for scaling parameters
    checkParamsScaling(params)

    // Check cellprofiler parameters
    checkParamsCpr(params)

    // Check flatfield estimation parameters
    checkParamsFlatfield(params)

    // Check deconvolution parameters
    checkParamsDecon(params)

    // Check registration parameters
    checkParamsRegistration(params)

    // Check cellcrops
    if (params.rn_make_cellcrops && !params.rn_cache_images) {
        error("rn_cache_images must be true when rn_make_cellcrops is true")
    }
    

    // Check cellpose is run
    if ((!params.cp_run && params.rn_cache_images) || (!params.cp_run && params.cpr_run)) {
        error("rn_cache_images & cpr_run can only be true when cp_run is also true")
    }
}


def checkParamsBase(params) {
        //------------------------------------------------------------
        // Defaults & sanity checks
        //------------------------------------------------------------
        if (params.rn_publish_dir == null) {
            error "rn_publish_dir file parameter is required: --rn_publish_dir"
        }
        
        if (params.rn_manifest == null) {
            error "rn_manifest file parameter is required: --rn_manifest"
        }
        
        if (!file(params.rn_manifest).exists()) {
            error("Specified manifest file does not exist or is not readable")
        }

        def manifestErrors = validateManifestContent(params.rn_manifest)
        if (manifestErrors) {
            error(formatValidationErrors("rn_manifest", manifestErrors))
        }

        // Check blacklist
        if (params.rn_blacklist != null) {
            if (!file(params.rn_blacklist).exists()) {
                error("Specified blacklist file does not exist or is not readable")
            }

            def blacklistErrors = validateBlacklistContent(params.rn_blacklist)
            if (blacklistErrors) {
                error(formatValidationErrors("rn_blacklist", blacklistErrors))
            }
        }

        // Check registration manifest
        if (params.rn_manifest_registration != null) {
            if (!file(params.rn_manifest_registration).exists()) {
                error("Specified registration manifest file does not exist or is not readable")
            }

            def registrationErrors = validateRegistrationManifestContent(params.rn_manifest_registration)
            if (registrationErrors) {
                error(formatValidationErrors("rn_manifest_registration", registrationErrors))
            }
        }

        // Check well level manifest(s) - comma separated list of paths
        if (params.rn_manifest_well != null) {
            params.rn_manifest_well.split(',').each { path ->
                if (!file(path).exists()) {
                    error("Specified rn_manifest_well file does not exist or is not readable: ${path}")
                }
            }
        }

        if (params.rn_cache_images && !(params.rn_max_project || params.rn_hybrid)) {
            log.warn "Caching images in 3d mode can take up a lot of space, are you sure this is what you want?"
        }

        // rn_hybrid ignores rn_max_project - warn so it's not mistaken for both being active
        if (params.rn_hybrid && params.rn_max_project) {
            log.warn "Both rn_hybrid and rn_max_project are true - rn_hybrid ignores rn_max_project, so rn_max_project has no effect"
        }

}


def checkParamsScaling(params) {
    if (!params.sc_autoscale && params.sc_control_list != null) {
        error "Provided --sc_control_list but --sc_autoscale false. Either drop --sc_control_list or set --sc_autoscale true"
    }
    
    if ((params.sc_manualscale != null) && params.sc_autoscale) {
        log.warn "Both sc_autoscale and sc_manualscale are provided, sc_manualscale will be overridden"
    }
        
    // Check control list
    if (params.sc_control_list != null) {
        if (!file(params.sc_control_list).exists()) {
            error("Specified control list file does not exist or is not readable")
        }

        def controlListErrors = validateControlListContent(params.sc_control_list)
        if (controlListErrors) {
            error(formatValidationErrors("sc_control_list", controlListErrors))
        }
    }

    // Check channel map
    if (params.sc_channel_map != null) {
        if (!file(params.sc_channel_map).exists()) {
            error("Specified channel map file does not exist or is not readable")
        }
    }

    // sc_control_list and sc_channel_map are both optional: with neither set, autoscale
    // falls back to dynamic-range-only scaling (sc_autoscale_q1/q2, no plate offsets or
    // sigmoid fitting). But plate-offset correction needs a feature to average per channel,
    // which only comes from sc_channel_map's rep_features_offset column, so sc_control_list
    // requires sc_channel_map to also be set.
    if (params.sc_autoscale) {
        if (params.sc_control_list != null && params.sc_channel_map == null) {
            error("sc_control_list requires sc_channel_map to be set - calculate_scaling_factors needs it to know which measured feature to average per channel for the plate offset")
        }

        if (!params.rn_cache_images) {
            error("Autoscaling is set to true, but image caching is disabled. Please set params.rn_cache_images=true")
        }
    }


    // Check manual scale
    if (params.sc_manualscale != null) {
        if (!file(params.sc_manualscale).exists()) {
            error("Specified manual scale file does not exist or is not readable")
        }
    }

    // sc_scale_slope and sc_scale_bias must be supplied together, or neither -
    // rescale_images.py errors out on this deep in the run, so catch it here instead
    if ((params.sc_scale_slope != null) != (params.sc_scale_bias != null)) {
        error("Must supply both --sc_scale_slope and --sc_scale_bias, or neither")
    }

    // Check scale slope/bias files
    if (params.sc_scale_slope != null) {
        if (!file(params.sc_scale_slope).exists()) {
            error("Specified sc_scale_slope file does not exist or is not readable")
        }
    }

    if (params.sc_scale_bias != null) {
        if (!file(params.sc_scale_bias).exists()) {
            error("Specified sc_scale_bias file does not exist or is not readable")
        }
    }

    // sc_scale_in_finalize applies scaling directly in finalize, before any unscaled images
    // exist to measure/estimate factors from - so it can only work with pre-supplied manual
    // factors, never with autoscale (which needs those unscaled images first).
    if (params.sc_scale_in_finalize) {
        if (params.sc_manualscale == null) {
            error("sc_scale_in_finalize requires sc_manualscale to be set - autoscale factors cannot be estimated without first producing unscaled finalized images, which sc_scale_in_finalize skips")
        }
        if (params.sc_autoscale) {
            error("sc_scale_in_finalize cannot be combined with sc_autoscale - autoscale factors require measuring the unscaled finalized images first, which sc_scale_in_finalize skips producing")
        }
    }

    // Autoscale factors are estimated from the unscaled finalized images, which only exist
    // when images are cached
    if (params.sc_autoscale && !params.rn_cache_images) {
        error("sc_autoscale requires rn_cache_images to be true - autoscale factors are estimated from the cached, finalized unscaled images")
    }

}


def checkParamsCpr(params) {
    if (params.cpr_run) {
        if (params.rn_max_project || params.rn_hybrid) {
            if (params.cpr_pipeline_2d == null) {
                error("Running in hybrid or 2d mode with cpr_run=true, but cpr_pipeline_2d is not set")
            }
            if (!file(params.cpr_pipeline_2d).exists()) {
                error("Specified 2d cellprofiler pipeline is not accessible")
            }
            if (params.cpr_pipeline_3d != null) {
                log.warn "cpr_pipeline_3d is set but running in hybrid/2d mode - cpr_pipeline_3d will be ignored"
            }
        } else {
            if (params.cpr_pipeline_3d == null) {
                error("Running in 3d mode with cpr_run=true, but cpr_pipeline_3d is not set")
            }
            if (!file(params.cpr_pipeline_3d).exists()) {
                error("Specified 3d cellprofiler pipeline is not accessible")
            }
            if (params.cpr_pipeline_2d != null) {
                log.warn "cpr_pipeline_2d is set but running in 3d mode - cpr_pipeline_2d will be ignored"
            }
        }

        // cpr_plugins is passed straight to --plugins-directory with no validation otherwise.
        // Warn rather than error: it defaults to a path that doesn't exist in a fresh checkout
        // (no plugins installed), which is a valid "no plugins" state, not a misconfiguration.
        if (params.cpr_plugins != null && !file(params.cpr_plugins).exists()) {
            log.warn "cpr_plugins directory does not exist, plugins will not be loaded: ${params.cpr_plugins}"
        }
    }

}


def checkParamsFlatfield(params) {
    if (params.ff_run) {
        def valid_modes = ["BASICPY", "POLY", "PE"]
        if (!(params.ff_mode in valid_modes)) {
            error("ff_mode must be one of ${valid_modes} but got: ${params.ff_mode}")
        }
    }
}


def checkParamsDecon(params) {
    if (params.dc_run) {
        def valid_modes = ["clij2", "clij2_nc", "rlf_cpu", "rlf_gpu"]
        if (!(params.dc_mode in valid_modes)) {
            error("dc_mode must be one of ${valid_modes} but got: ${params.dc_mode}")
        }

        // Decon's storeDir is rn_decon_dir - if it matches rn_image_dir, decon output
        // would land on top of the raw images
        if (params.rn_image_dir == params.rn_decon_dir) {
            error("rn_image_dir and rn_decon_dir must not be the same path when dc_run is true")
        }
    }
}


def checkParamsRegistration(params) {
    if (params.rn_manifest_registration != null) {
        // STACKREG is disabled due to dependency issues with pystackreg (see nextflow.config)
        if (params.rg_mode != "CROSS") {
            error("rg_mode '${params.rg_mode}' is not supported - only 'CROSS' is currently enabled (STACKREG is disabled due to dependency issues)")
        }
    }
}