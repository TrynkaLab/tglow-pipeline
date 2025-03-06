
// Convert a 1 indexed integer string to zero indexed integer
def convertChannelType(String input) {
    if (input.isInteger()) {
        return (input.toInteger() - 1)
    } else {
        return -9
    }        
}


// Parse the manfiest and extract what things should get flatfields
def parseManifestFlatfields(String rn_manifest) {
    //def csvFile = new File(params.rn_manifest)
    def csvFile = new File(rn_manifest)
    def csvData = csvFile.readLines()

    def plate_channel = []
    def plate_channel_global = []

    // Skip header
    // Substract one from the channel is python is 0 indexed
    // Returns a list of tuples, one for each plate channel with (plate, channel, pe_index)
    for (int i = 1; i < csvData.size(); i++) {
        def curLine = csvData[i].split('\t')
        
        if (curLine[3] != "none") {
            def curChannels = curLine[3].split(',')
            for (channel in curChannels) {
                plate_channel << tuple(curLine[0], channel.toInteger()-1, curLine[1])
                
                // Save the first line as the global channel
                if (i == 1) {
                    plate_channel_global << tuple(curLine[0], channel.toInteger()-1, curLine[1])
                }
            }
        }
    }
    return [plate_channel, plate_channel_global]
}


// Read blacklist as list of plate:well
def readBlacklist(String path){
    
    blacklist=[]
    
    new File(path).splitEachLine("\t") {fields ->
        blacklist.add(fields[0] + ":" + fields[1])
    }
    
    log.info("Blacklist consists of " + blacklist.size() + " items")
            
    return blacklist
}


// Parse the registration manifest and return a list of lists
// one item per qry-ref pair
def parseRegManfiestAsQryRef(String rn_manifest_registration) {
    qry_ref=[]
    new File(rn_manifest_registration).splitEachLine("\t") {fields ->
        ref = fields[0]  
        
        // Add the ref plate
        qry_ref.add([ref, ref])
        
        // Add the qry plates      
        qry_items = fields[2].split(",")
        for (qry in qry_items){
            qry_ref.add([qry, ref])
        }
    }
    
    return qry_ref         
}


def parseRegManfiestAsQry(String rn_manifest_registration) {
    qry_ref=[]
    new File(rn_manifest_registration).splitEachLine("\t") {fields ->
        ref = fields[0]  
        qry_ref.add(ref)
        // Add the ref plate
        //qry_ref.add([ref, ref])
        
        // Add the qry plates      
        //qry_items = fields[2].split(",")
        //for (qry in qry_items){
        //    qry_ref.add([qry, ref])
       // }
    }
    
    return qry_ref         
}


