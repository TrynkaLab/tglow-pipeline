import nextflow.io.ValueObject

@ValueObject
class Well {
    String plate
    String row
    String col
    String well
    String key
    String relpath

    static Well fromMap(Map args) {
        def plate = args.plate
        def well = args.well
        if (!plate) throw new IllegalArgumentException("plate must be provided")
        if (!well) throw new IllegalArgumentException("well must be provided")
        def row = args.row ?: well.replaceAll(/\d+/,'')
        def col = args.col ?: well.replaceAll(/[A-Z]+/,'').replaceAll(/^0+/, '')

        new Well(
            plate:   plate,
            well:    well,
            row:     row,
            col:     col,
            key:     plate + ":" + well,
            relpath: plate + "/" + row + "/" + col
        )
    }

    static Well fromPlateWell(String plate, String well) {
        def row = well.replaceAll(/\d+/,'')
        def col = well.replaceAll(/[A-Z]+/,'').replaceAll(/^0+/, '')

        new Well(
            plate: plate,
            well:  well,
            row:   row,
            col:   col,
            key:   plate + ":" + well
        )
    }

    static Well fromVars(String plate, String row, String col, String well) {
        new Well(
            plate: plate,
            row:   row,
            col:   col,
            well:  well,
            key:   plate + ":" + well
        )
    }
    
    String toString() {
        return "Well(key=$key, plate=$plate, well=$well, row=$row, col=$col)"
    }
}