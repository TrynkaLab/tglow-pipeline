#!/usr/bin/env nextflow

// Content validation for the user-supplied TSV input files (main manifest,
// blacklist, registration manifest, control list). lib/checks.nf only checks
// existence - this file checks headers and per-row value shape, so a
// malformed file fails fast with one aggregated error instead of silently
// falling through to a "none" sentinel deep in a downstream script.

// The only accepted "disabled/not set" marker for optional fields.
def SENTINEL() { return "none" }

def SINGLE_INT_RE() { return ~/^[0-9]+$/ }
def INT_LIST_RE() { return ~/^[0-9]+(,[0-9]+)*$/ }
// Long well format: 1-2 letter row + zero-padded 2-3 digit column, e.g. A01 (not A1)
def WELL_RE() { return ~/^[A-Z]{1,2}[0-9]{2,3}$/ }

// Values people commonly mistype in place of the "none" sentinel - these are
// deliberately rejected everywhere rather than silently tolerated
def BLANKISH() { return ["", "na", "n/a", SENTINEL(), "null", "nan", "false", "0"] }

def isBlankish(String value) {
    return value == null || BLANKISH().contains(value.trim().toLowerCase())
}

// value must be exactly the sentinel, or satisfy pattern - used for optional columns
def checkOptional(String value, def pattern, String colName, String hint) {
    if (value == SENTINEL()) return null
    if (isBlankish(value)) {
        return "${colName} '${value}' is invalid - use ${hint} or '${SENTINEL()}'"
    }
    if (pattern != null && !(value ==~ pattern)) {
        return "${colName} '${value}' is invalid - use ${hint} or '${SENTINEL()}'"
    }
    return null
}

// value must satisfy pattern (or just be non-blank if pattern is null) - used for required columns
def checkRequired(String value, def pattern, String colName, String hint) {
    if (isBlankish(value)) {
        return "required field '${colName}' is missing/blank"
    }
    if (pattern != null && !(value ==~ pattern)) {
        return "${colName} '${value}' is invalid - use ${hint}"
    }
    return null
}

// value must be a comma-separated list of <positive integer>=<non-blank path> pairs,
// with no channel number repeated - used for dc_psfs
def checkPsfMap(String value, String colName) {
    def hint = "comma-separated <channel>=<path> pairs (e.g. 1=psf1.tif,4=psf4.tif)"
    def tokens = value.split(',')

    def malformed = tokens.find { token ->
        def parts = token.split('=', 2)
        parts.size() != 2 || !(parts[0] ==~ SINGLE_INT_RE()) || parts[1].trim().isEmpty()
    }
    if (malformed != null) {
        return "${colName} '${value}' is invalid - use ${hint} or '${SENTINEL()}'"
    }

    def channels = tokens.collect { it.split('=', 2)[0] }
    def duplicates = channels.findAll { channels.count(it) > 1 }.unique()
    if (duplicates) {
        return "${colName} '${value}' has duplicate channel(s) ${duplicates} - each channel may only appear once"
    }

    return null
}


// Main manifest (rn_manifest)
def validateManifestContent(String path) {
    def MANIFEST_COLUMNS = ["plate", "index_xml", "channels", "ff_channels", "cp_nucl_channel",
                            "cp_cell_channel", "dc_psfs", "mask_channels"] as Set

    def errors = []
    def lines = new File(path).readLines()
    if (lines.isEmpty()) {
        return ["file is empty"]
    }

    def header = lines[0].split('\t', -1)*.trim()
    def missing = MANIFEST_COLUMNS - (header as Set)
    def extra = (header as Set) - MANIFEST_COLUMNS
    if (missing || extra) {
        def msg = "header does not match expected columns ${MANIFEST_COLUMNS}."
        if (missing) msg += " Missing: ${missing}."
        if (extra) msg += " Unexpected: ${extra}."
        return [msg]
    }

    def colIndex = header.withIndex().collectEntries { name, i -> [(name): i] }

    lines.drop(1).eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        def fields = line.split('\t', -1)*.trim()
        def rowNum = i + 1
        def plate = fields[colIndex.plate]
        def context = isBlankish(plate) ? "row ${rowNum}" : "row ${rowNum} (plate=${plate})"

        def rowErrors = []
        rowErrors << checkRequired(plate, null, "plate", "a plate name")
        rowErrors << checkRequired(fields[colIndex.channels], INT_LIST_RE(), "channels", "a comma-separated list of integers")
        rowErrors << checkOptional(fields[colIndex.index_xml], null, "index_xml", "a file path")
        rowErrors << checkOptional(fields[colIndex.ff_channels], INT_LIST_RE(), "ff_channels", "a comma-separated list of integers")
        rowErrors << checkOptional(fields[colIndex.cp_nucl_channel], SINGLE_INT_RE(), "cp_nucl_channel", "a single integer")
        rowErrors << checkOptional(fields[colIndex.cp_cell_channel], SINGLE_INT_RE(), "cp_cell_channel", "a single integer")
        rowErrors << checkOptional(fields[colIndex.mask_channels], INT_LIST_RE(), "mask_channels", "a comma-separated list of integers")

        def dcPsfs = fields[colIndex.dc_psfs]
        if (dcPsfs != SENTINEL()) {
            def psfError = isBlankish(dcPsfs) ?
                "dc_psfs '${dcPsfs}' is invalid - use comma-separated <channel>=<path> pairs (e.g. 1=psf1.tif,4=psf4.tif) or '${SENTINEL()}'" :
                checkPsfMap(dcPsfs, "dc_psfs")
            rowErrors << psfError
        }

        rowErrors.findAll { it != null }.each { errors << "${context}: ${it}" }
    }

    return errors
}


// Blacklist (rn_blacklist) - no header, <plate> <well>
def validateBlacklistContent(String path) {
    def errors = []
    new File(path).readLines().eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        def fields = line.split('\t', -1)*.trim()
        def rowNum = i + 1
        if (fields.size() != 2) {
            errors << "row ${rowNum}: expected 2 tab-separated columns (plate, well), got ${fields.size()}"
            return
        }
        def (plate, well) = fields
        def err = checkRequired(plate, null, "plate", "a plate name")
        if (err) errors << "row ${rowNum}: ${err}"
        err = checkRequired(well, WELL_RE(), "well", "long format, e.g. A01")
        if (err) errors << "row ${rowNum} (plate=${plate}): ${err}"
    }
    return errors
}


// Registration manifest (rn_manifest_registration)
def validateRegistrationManifestContent(String path) {
    def REGISTRATION_COLUMNS = ["reference_plate", "reference_channel", "query_plates", "query_channels"]

    def errors = []
    def lines = new File(path).readLines()
    if (lines.isEmpty()) {
        return ["file is empty"]
    }

    def header = lines[0].split('\t', -1)*.trim()
    if (header != REGISTRATION_COLUMNS) {
        return ["header must be exactly ${REGISTRATION_COLUMNS.join('\t')}, got ${header.join('\t')}".toString()]
    }

    lines.drop(1).eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        def fields = line.split('\t', -1)*.trim()
        def rowNum = i + 1
        def (refPlate, refChannel, qryPlates, qryChannels) = fields
        def context = isBlankish(refPlate) ? "row ${rowNum}" : "row ${rowNum} (reference_plate=${refPlate})"

        def rowErrors = []
        rowErrors << checkRequired(refPlate, null, "reference_plate", "a plate name")
        rowErrors << checkRequired(refChannel, SINGLE_INT_RE(), "reference_channel", "a single integer")
        rowErrors << checkRequired(qryPlates, null, "query_plates", "a plate name or comma-separated list of plate names")
        rowErrors << checkRequired(qryChannels, INT_LIST_RE(), "query_channels", "an integer or comma-separated list of integers")

        if (!isBlankish(qryPlates) && !isBlankish(qryChannels) && (qryChannels ==~ INT_LIST_RE())) {
            def plateCount = qryPlates.split(',').size()
            def channelCount = qryChannels.split(',').size()
            if (plateCount != channelCount) {
                rowErrors << "query_plates has ${plateCount} entries but query_channels has ${channelCount} - they must pair up one-to-one".toString()
            }
        }

        rowErrors.findAll { it != null }.each { errors << "${context}: ${it}" }
    }

    return errors
}


// Control list (sc_control_list)
def validateControlListContent(String path) {
    def CONTROL_LIST_COLUMNS = ["plate", "well", "control_type"]

    def errors = []
    def lines = new File(path).readLines()
    if (lines.isEmpty()) {
        return ["file is empty"]
    }

    def header = lines[0].split('\t', -1)*.trim()
    if (header != CONTROL_LIST_COLUMNS) {
        return ["header must be exactly ${CONTROL_LIST_COLUMNS.join('\t')}, got ${header.join('\t')}".toString()]
    }

    lines.drop(1).eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        def fields = line.split('\t', -1)*.trim()
        def rowNum = i + 1
        def (plate, well, controlType) = fields
        def context = isBlankish(plate) ? "row ${rowNum}" : "row ${rowNum} (plate=${plate})"

        def rowErrors = []
        rowErrors << checkRequired(plate, null, "plate", "a plate name")
        rowErrors << checkRequired(well, WELL_RE(), "well", "long format, e.g. A01")
        rowErrors << checkRequired(controlType, null, "control_type", "a control type label")

        rowErrors.findAll { it != null }.each { errors << "${context}: ${it}" }
    }

    return errors
}


def formatValidationErrors(String label, List errors) {
    if (!errors) return null
    def lines = errors.collect { "  ${it}" }.join('\n')
    return "${label} has ${errors.size()} problem(s):\n${lines}"
}
