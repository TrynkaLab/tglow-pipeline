import nextflow.io.ValueObject

@ValueObject
class RegistrationRecord {

    String ref_plate
    def ref_channel
    def qry_plates
    def qry_channels

    static RegistrationRecord fromVars(String ref_plate, Integer ref_channel, List qry_plates, List qry_channels) {
        new RegistrationRecord(
            ref_plate:    ref_plate,
            ref_channel:  ref_channel,
            qry_plates:   qry_plates,
            qry_channels: qry_channels
        )
    }

    static RegistrationRecord fromRow(Map row) {
        new RegistrationRecord(
            ref_plate:    row.ref_plate,
            ref_channel:  (row.ref_channel as Integer) - 1,
            qry_plates:   row.qry_plates.split(',').toList(),
            qry_channels: row.qry_channels.split(',').collect { (it as Integer) - 1 }
        )
    }

    String toString() {
        return "RegistrationRecord(ref_plate=${ref_plate}, ref_channel=${ref_channel}, qry_plates=${qry_plates}, qry_channels=${qry_channels})"
    }    

}