#' @include utils-IPO.R
#' @include doe-optimizer.R
#' @include features.R
#' @include feature_groups.R
#' @include main.R
NULL

featureGroupsOptimizer <- setRefClass("featureGroupsOptimizer", contains = c("DoEOptimizer", "VIRTUAL"),
                                      fields = list(features = "features"))

featureGroupsOptimizer$methods(

    defaultParamRanges = function(...) getDefFGroupsOptParamRanges(algorithm),

    calculateResponse = function(params, task, keepObject)
    {
        # UNDONE: error handling necessary as is done in IPO?
        # retCorFailed <- if (!is.null(params[["rtalign"]]) && params$rtalign) 1.1 else 1
        retCorFailed <- 1

        fg <- do.call(groupFeatures, c(list(features, algorithm, verbose = FALSE), params))

        ret <- utilsIPO$getRGTVValues(getXCMSSet(fg, verbose = FALSE, exportedData = TRUE), task, retCorFailed)

        if (keepObject)
            ret <- list(response = ret, object = fg)

        return(ret)
    },

    getResponseScores = function(response)
    {
        GS <- response$GS
        RCS <- response$RCS

        # give penalty when retcor failed
        RCS_penalty <- 1 / response$retcor_done
        RCS <- RCS / RCS_penalty

        # normalize
        norm_GS <- (GS - min(GS)) / (max(GS) - min(GS))
        norm_RCS <- (RCS - min(RCS)) / (max(RCS) - min(RCS))
        norm_GS[is.na(norm_GS)] <- 0
        norm_RCS[is.na(norm_RCS)] <- 0

        return(norm_GS + norm_RCS)
    },

    getFinalScore = function(oldr, newr)
    {
        rs <- getResponseScores(rbind(oldr[, -"experiment"], newr[names(newr) != "object"]))
        return(rs[length(rs)])
    },

    resultIncreased = function(history)
    {
        index = length(history)
        if (index < 2)
            return(TRUE)

        prevFR <- history[[index-1]]$finalResult$response
        curFR <- history[[index]]$finalResult$response

        if (curFR$bad_groups == 0)
        {
            curFR$bad_groups = 1
            curFR$good_groups = curFR$good_groups + 1
        }

        if (prevFR$bad_groups == 0)
        {
            prevFR$bad_groups = 1
            prevFR$good_groups = prevFR$good_groups + 1
        }

        if ((curFR$good_groups^2/curFR$bad_groups <= prevFR$good_groups^2/prevFR$bad_groups) ||
            (curFR$RCS <= prevFR$RCS))
            return(FALSE)

        return(TRUE)
    }
)


#' @rdname feature-optimization
#' @param features A \code{\link{features}} object with the features that should
#'   be used to optimize grouping.
#' @export
optimizeFeatureGrouping <- function(features, algorithm, ..., templateParams = list(),
                                    paramRanges = list(), maxIterations = 50, maxModelDeviation = 0.1,
                                    parallel = TRUE)
{
    params <- list(...)

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(features, "features", add = ac)
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "kpic2"), add = ac)
    assertOptimArgs(params, templateParams, paramRanges, maxIterations, maxModelDeviation, parallel, ac)
    checkmate::reportAssertions(ac)

    go <- switch(algorithm,
                 openms = featureGroupsOptimizerOpenMS,
                 xcms = featureGroupsOptimizerXCMS,
                 xcms3 = featureGroupsOptimizerXCMS3,
                 kpic2 = featureGroupsOptimizerKPIC2)

    go <- go$new(features = features, algorithm = algorithm, parallel = parallel)
    result <- go$optimize(params, templateParams, paramRanges, maxIterations, maxModelDeviation)

    return(optimizationResult(algorithm = algorithm, paramSets = result$paramSets,
                              bestParamSet = result$bestParamSet))
}

#' @rdname feature-optimization
#' @export
generateFGroupsOptPSet <- function(algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "kpic2"))

    f <- switch(algorithm,
                openms = generateFGroupsOptPSetOpenMS,
                xcms = generateFGroupsOptPSetXCMS,
                xcms3 = generateFGroupsOptPSetXCMS3,
                kpic2 = generateFGroupsOptPSetKPIC2)

    defs <- f(...)
    return(modifyList(defs, list(...)))
}

#' @rdname feature-optimization
#' @export
getDefFGroupsOptParamRanges <- function(algorithm)
{
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "kpic2"))

    if (algorithm == "openms")
        return(getDefFGroupsOptParamRangesOpenMS())
    if (algorithm == "xcms")
        return(getDefFGroupsOptParamRangesXCMS())
    if (algorithm == "xcms3")
        return(getDefFGroupsOptParamRangesXCMS3())
    # if (algorithm == "kpic2")
    return(getDefFGroupsOptParamRangesKPIC2())
}
