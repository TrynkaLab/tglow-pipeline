class Well {
    String plate
    String row
    String col
    String well
    String key
    String relpath
    
    Well(Map args) {
        this.plate = args.plate
        this.well = args.well
        this.row = args.row
        this.col = args.col

        if (!plate) throw new IllegalArgumentException("plate must be provided")
        if (!well) throw new IllegalArgumentException("well must be provided")
        if (!row) this.row = well.replaceAll(/\d+/,'')
        if (!col) this.col = well.replaceAll(/[A-Z]+/,'').replaceAll(/^0+/, '')
        this.key = plate + ":" + well
        this.relpath = plate + "/" + row + "/" + col
    }
    
    Well(String plate, String well) {
        this.plate = plate
        this.well = well
        this.row = well.replaceAll(/\d+/,'')
        this.col = well.replaceAll(/[A-Z]+/,'').replaceAll(/^0+/, '')
        this.key = plate + ":" + well
    }
    
    Well(String plate, String row, String col, String well) {
        this.plate = plate
        this.row = row
        this.col = col
        this.well = well
        this.key = plate + ":" + well
    }
    
    String toString() {
        return "Well(key=$key, plate=$plate, well=$well, row=$row, col=$col)"
    }
}