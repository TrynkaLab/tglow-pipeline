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

        this.dc_channels = (args.dc_channels == null || args.dc_channels == "none") ? "none" : args.dc_channels.split(',').collect { it as Integer - 1 }

        this.dc_psfs = (args.dc_psfs == null || args.dc_psfs == "none") ? "none" : args.dc_psfs.split(',')

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