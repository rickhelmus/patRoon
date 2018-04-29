checkHasNames <- function(x, n, type = "unique") checkmate::checkNames(names(x), must.include = n, type = type)
assertHasNames <- checkmate::makeAssertionFunction(checkHasNames)

assertAnalysisInfo <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertDataFrame(x, types = "character", any.missing = FALSE, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, c("path", "analysis", "group", "ref"), .var.name = .var.name, add = add)
}
