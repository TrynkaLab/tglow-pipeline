class ManifestRecord {
    String plate
    String index_xml
    def channels
    def bp_channels
    String cp_nucl_channel
    String cp_cell_channel
    def dc_channels
    def dc_psfs
    def mask_channels

    ManifestRecord(Map args) {
        this.plate = args.plate
        
        this.index_xml = args.index_xml

        this.channels = (args.channels == null || args.channels == "none") ? "none" : args.channels.split(',').collect { it as Integer - 1 }

        this.bp_channels = (args.bp_channels == null || args.bp_channels == "none") ? "none" : args.bp_channels.split(',').collect { it as Integer - 1 }

        this.cp_nucl_channel = args.cp_nucl_channel
        
        this.cp_cell_channel = args.cp_cell_channel

        this.dc_channels = (args.dc_channels == null || args.dc_channels == "none") ? "none" : args.dc_channels.split(',').collect { it as Integer - 1 }

        this.dc_psfs = (args.dc_psfs == null) ? "none" : args.dc_psfs.split(',')

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
    //         bp_channels: ${bp_channels},
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