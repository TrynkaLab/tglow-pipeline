class ManifestRecord {
    String plate
    String index_xml
    def channels
    def ff_channels
    def cp_nucl_channel
    def cp_cell_channel
    def dc_channels
    def dc_psfs
    def mask_channels

    ManifestRecord(Map args) {
        this.plate = args.plate
        
        this.index_xml = args.index_xml

        this.channels = (args.channels == null || args.channels == "none") ? "none" : args.channels.split(',').collect { it as Integer - 1 }

        this.ff_channels = (args.ff_channels == null || args.ff_channels == "none") ? "none" : args.ff_channels.split(',').collect { it as Integer - 1 }

        this.cp_nucl_channel = (args.cp_nucl_channel == null || args.cp_nucl_channel == "none") ? "none" : (args.cp_nucl_channel as Integer) - 1
        
        this.cp_cell_channel =  (args.cp_cell_channel == null || args.cp_cell_channel == "none") ? "none" : (args.cp_cell_channel as Integer) - 1

        // dc_psfs is a comma-separated list of <1-indexed channel>=<path> pairs, e.g.
        // "1=psf1.tif,4=psf4.tif". dc_channels is derived from its keys rather than
        // being its own column, so the two can never fall out of sync.
        this.dc_psfs = (args.dc_psfs == null || args.dc_psfs == "none") ? "none" : args.dc_psfs.split(',').collectEntries { pair ->
            def (channel, path) = pair.split('=', 2)
            [(channel as Integer - 1): path]
        }

        this.dc_channels = (this.dc_psfs == "none") ? "none" : this.dc_psfs.keySet().toList()

        this.mask_channels = (args.mask_channels == null || args.mask_channels == "none") ? "none" : args.mask_channels.split(',').collect { it as Integer - 1 }
    }
    
    String toString() {
        return "ManifestRecord(plate: ${plate}, ...)"
    }

    // @Override
    // String toString() {
    //     return(
    //     """
    //     ManifestRecord(
    //         plate: ${plate},
    //         index_xml: ${index_xml},
    //         channels: ${channels},
    //         ff_channels: ${ff_channels},
    //         cp_nucl_channel: ${cp_nucl_channel},
    //         cp_cell_channel: ${cp_cell_channel},
    //         dc_channels: ${dc_channels},
    //         dc_psfs: ${dc_psfs},
    //         mask_channels: ${mask_channels},
    //         scale_factors: ${scale_factors}
    //     )
    //     """)
    // }
}