

class RegistrationRecord {
    
    String ref_plate
    String ref_channel
    def qry_plates
    def qry_channels

    RegistrationRecord(Map args) {
        
        this.ref_plate = args.ref_plate
        this.ref_channel = args.ref_channel
        this.qry_plates = args.qry_plates.split(',')
        this.qry_channels = args.qry_channels.split(',').collect { it as Integer - 1 }
    }  

    String toString() {
        return "RegistrationRecord(ref_plate=${ref_plate}, ref_channel=${ref_channel}, qry_plates=${qry_plates}, qry_channels=${qry_channels})"
    }    
}