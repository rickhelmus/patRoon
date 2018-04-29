checkHasNames <- function(x, n, type = "unique") checkmate::checkNames(names(x), must.include = n, type = type)
assertHasNames <- checkmate::makeAssertionFunction(checkHasNames)

checkAnalysisInfo <- function(anaInfo)
{
    ret <- checkmate::checkDataFrame(anaInfo, types = "character", any.missing = FALSE, min.rows = 1)
    if (!isTRUE(ret))
        return(ret)
    
    return(checkHasNames(anaInfo, c("path", "analysis", "group", "ref")))
}
assertAnalysisInfo <- checkmate::makeAssertionFunction(checkAnalysisInfo)
