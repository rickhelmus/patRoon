#' @include main.R
NULL

# from https://stackoverflow.com/a/17531678
orderComponentsNames <- function(n) order(nchar(n), n)

mergeComponents <- function(compList, compNames, nameColumn)
{
    notEmpty <- lengths(compList) > 0
    compList <- compList[notEmpty]; compNames <- compNames[notEmpty]
    
    if (length(compList) == 0)
        return(list(components = list(), componentInfo = data.table()))
    
    retCInfo <- copy(componentInfo(compList[[1]]))
    retCInfo[, (nameColumn) := compNames[[1]]]
    retCTable <- componentTable(compList[[1]])

    if (length(compList) > 1)
    {
        for (mi in seq(2, length(compList)))
        {
            if (length(compList[[mi]]) == 0)
                next
            rci <- copy(componentInfo(compList[[mi]]))
            rci[, (nameColumn) := compNames[[mi]]]
            retCInfo <- rbind(retCInfo, rci, fill = TRUE)
            retCTable <- c(retCTable, componentTable(compList[[mi]]))
        }
    }               
    
    retCInfo[, name := paste0(name, "-", get(nameColumn))][]
    names(retCTable) <- retCInfo[["name"]]
    
    return(list(components = retCTable, componentInfo = retCInfo))
}

calculateComponentIntensities <- function(comps, fGroups)
{
    getGroupInt <- function(grp)
    {
        ints <- fGroups[[grp]]
        return(mean(ints[ints != 0]))
    }
    return(lapply(comps, function(cmp)
    {
        cmp <- copy(cmp)
        cmp[, intensity := sapply(group, getGroupInt)]
        cmp[, intensity_rel := intensity / max(intensity)]
        return(cmp[])
    }))
}
