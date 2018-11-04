#' @include utils-IPO.R
#' @include doe-optimizer.R
#' @include features.R
#' @include feature_groups.R
#' @include main.R
NULL

featureGroupsOptimizer <- setRefClass("featureGroupsOptimizer", contains = c("DoEOptimizer", "VIRTUAL"),
                                      fields = list(features = "features"))

featureGroupsOptimizer$methods(

    # dummy methods to be potentially overrided
    convertOptToCallParams = function(params) params,

    calculateResponse = function(params, task, keepObject)
    {
        # UNDONE: do we want to keep caching this?

        # UNDONE: error handling necessary as is done in IPO?
        # retCorFailed <- if (!is.null(params[["rtalign"]]) && params$rtalign) 1.1 else 1
        retCorFailed <- 1

        params <- convertOptToCallParams(params)
        fg <- do.call(groupFeatures, c(list(features, algorithm, verbose = FALSE), params))

        ret <- utilsIPO$getRGTVValues(getXCMSSet(fg, TRUE, verbose = FALSE), task, retCorFailed)

        if (keepObject)
            ret <- list(response = ret, object <- fg)

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


optimizeFeatureGrouping <- function(features, algorithm, ..., templateParams = list(),
                                    paramRanges = list(), maxIterations = 50, maxModelDeviation = 0.1)
{
    params <- list(...)

    ac <- checkmate::makeAssertCollection()
    checkmate::assert_class(features, "features", add = ac)
    checkmate::assertChoice(algorithm, c("openms", "xcms"), add = ac)
    assertOptimArgs(params, templateParams, paramRanges, maxIterations, maxModelDeviation, ac)
    checkmate::reportAssertions(ac)

    go <- switch(algorithm,
                 openms = featureGroupsOptimizerOpenMS,
                 xcms = featureGroupsOptimizerXCMS)

    go <- go$new(features = features, algorithm = algorithm)
    result <- go$optimize(params, templateParams, paramRanges, maxIterations, maxModelDeviation)

    return(optimizationResult(algorithm = algorithm, paramSets = result$paramSets,
                              bestParamSet = result$bestParamSet))
}
