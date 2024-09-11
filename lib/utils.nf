
// Convert a 1 indexed integer string to zero indexed integer
def convertChannelType(String input) {
    if (input.isInteger()) {
        return (input.toInteger() - 1)
    } else {
        return -9
    }        
}
