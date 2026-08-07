class RegistrationRecord {

    static Map fromVars(String ref_plate, Integer ref_channel, List qry_plates, List qry_channels) {
        [
            ref_plate:    ref_plate,
            ref_channel:  ref_channel,
            qry_plates:   qry_plates,
            qry_channels: qry_channels
        ]
    }

    static Map fromRow(Map row) {
        [
            ref_plate:    row.ref_plate,
            ref_channel:  (row.ref_channel as Integer) - 1,
            qry_plates:   row.qry_plates.split(',').toList(),
            qry_channels: row.qry_channels.split(',').collect { (it as Integer) - 1 }
        ]
    }

}
