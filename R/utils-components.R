#' @include main.R
NULL

# from https://stackoverflow.com/a/17531678
orderComponentsNames <- function(n) order(nchar(n), n)

mergeComponents <- function(compList, compNames, nameColumn)
{
    retCInfo <- copy(componentInfo(compList[[1]]))
    retCInfo[, (nameColumn) := compNames[[1]]]
    retCTable <- componentTable(compList[[1]])
    
    for (mi in seq(2, length(compList)))
    {
        if (length(compList[[mi]]) == 0)
            next
        rci <- copy(componentInfo(compList[[mi]]))
        rci[, (nameColumn) := compNames[[mi]]]
        retCInfo <- rbind(retCInfo, rci, fill = TRUE)
        retCTable <- c(retCTable, componentTable(compList[[mi]]))
    }
    
    retCInfo[, name := paste0(name, "-", get(nameColumn))][]
    names(retCTable) <- retCInfo[["name"]]
    
    return(list(components = retCTable, componentInfo = retCInfo))
}
