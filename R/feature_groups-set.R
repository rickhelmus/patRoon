# UNDONE: combine featureGroupsComparison interface (perhaps this could replace it?)
# UNDONE: design common set base class for this, compounds etc?

#' @include main.R
#' @include feature_groups-comparison.R
NULL

neutralizeFGroups <- function(fGroups, adduct)
{
    adductMZ <- adductMZDelta(adduct)
    fGroups@groupInfo$mzs <- fGroups@groupInfo$mzs - adductMZ
    fGroups@features@features <- lapply(fGroups@features@features, function(fTab)
    {
        fTab <- copy(fTab)
        fTab[, mz := mz - adductMZ]
        fTab[, mzmin := mz - adductMZ]
        fTab[, mzmax := mz - adductMZ]
        return(fTab)
    })
    return(fGroups)
}

featureGroupsSet <- setClass("featureGroupsSet",
                             slots = c(adducts = "list", setObjects = "list",
                                       neutralizedFGroups = "list"),
                             contains = "featureGroups")

setMethod("initialize", "featureGroupsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("c", "featureGroups", function(x, ..., adducts, groupAlgo, groupArgs = list(rtalign = FALSE))
{
    fGroupsList <- list(x, ...)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(fGroupsList, types = "featureGroups", any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assert(checkmate::checkCharacter(adducts, any.missing = FALSE, min.len = 1,
                                                max.len = length(fGroupsList)),
                      checkmate::checkList(adducts, types = c("adduct", "character"), any.missing = FALSE,
                                           min.len = 1, max.len = length(fGroupsList)),
                      .var.name = "adducts")
    checkmate::assertChoice(groupAlgo, c("xcms", "openms"), add = ac)
    checkmate::assertList(groupArgs, any.missing = FALSE, names = "unique", add = ac)
    checkmate::reportAssertions(ac)
    
    n <- getArgNames(..., def = sapply(fGroupsList, algorithm))
    names(fGroupsList) <- make.unique(n)
    
    adductNamed <- checkmate::testNames(names(adducts))
    if (adductNamed)
        checkmate::assertNames(names(adducts), type = "unique", must.include = names(fGroupsList))
    adducts <- lapply(adducts, checkAndToAdduct, .var.name = "adducts")
    adducts <- rep(adducts, length.out = length(fGroupsList))
    
    if (!adductNamed)
        names(adducts) <- names(fGroupsList)
    else
        adducts <- adducts[names(fGroupsList)] # synchronize order
    
    # neutralize featureGroups
    neutralizedFGroups <- mapply(fGroupsList, adducts, FUN = neutralizeFGroups,
                                 SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    # convert feature groups to features
    featsFromGroups <- convertFeatureGroupsToFeatures(neutralizedFGroups)
    
    if (groupAlgo == "xcms")
        groupArgs <- c(list(exportedData = FALSE), groupArgs)
    compGroups <- do.call(groupFeatures, c(list(featsFromGroups, groupAlgo), groupArgs))
    
    anaInfo <- do.call(rbind, c(mapply(fGroupsList, names(fGroupsList), FUN = function(fG, n)
    {
        fG@analysisInfo$set <- n
        return(fG@analysisInfo)
    }, SIMPLIFY = FALSE), list(stringsAsFactors = FALSE, make.row.names = FALSE)))
    
    return(featureGroupsSet(adducts = adducts, setObjects = fGroupsList, neutralizedFGroups = neutralizedFGroups,
                            groups = groups(compGroups), analysisInfo = anaInfo, groupInfo = groupInfo(compGroups),
                            features = getFeatures(compGroups), ftindex = groupFeatIndex(compGroups)))
})
