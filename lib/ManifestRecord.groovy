import nextflow.io.ValueObject

@ValueObject
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

    static ManifestRecord fromRow(Map row) {
        // dc_psfs is a comma-separated list of <1-indexed channel>=<path> pairs, e.g.
        // "1=psf1.tif,4=psf4.tif". dc_channels is derived from its keys rather than
        // being its own column, so the two can never fall out of sync.
        def dc_psfs = (row.dc_psfs == null || row.dc_psfs == "none") ? "none" : row.dc_psfs.split(',').collectEntries { pair ->
            def (channel, path) = pair.split('=', 2)
            [(channel as Integer - 1): path]
        }

        new ManifestRecord(
            plate:           row.plate,
            index_xml:       row.index_xml,
            channels:        (row.channels == null || row.channels == "none") ? "none" : row.channels.split(',').collect { it as Integer - 1 },
            ff_channels:     (row.ff_channels == null || row.ff_channels == "none") ? "none" : row.ff_channels.split(',').collect { it as Integer - 1 },
            cp_nucl_channel: (row.cp_nucl_channel == null || row.cp_nucl_channel == "none") ? "none" : (row.cp_nucl_channel as Integer) - 1,
            cp_cell_channel: (row.cp_cell_channel == null || row.cp_cell_channel == "none") ? "none" : (row.cp_cell_channel as Integer) - 1,
            dc_psfs:         dc_psfs,
            dc_channels:     (dc_psfs == "none") ? "none" : dc_psfs.keySet().toList(),
            mask_channels:   (row.mask_channels == null || row.mask_channels == "none") ? "none" : row.mask_channels.split(',').collect { it as Integer - 1 }
        )
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