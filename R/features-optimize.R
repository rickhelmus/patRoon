#' @include utils-IPO.R
#' @include doe-optimizer.R
#' @include features.R
#' @include main.R
NULL

featuresOptimizer <- setRefClass("featuresOptimizer", contains = c("DoEOptimizer", "VIRTUAL"),
                                 fields = list(anaInfo = "data.frame", isoIdent = "character",
                                               checkPeakShape = "character", CAMERAOpts = "list"))

featuresOptimizer$methods(

    defaultParamRanges = function(params) getDefFeaturesOptParamRanges(algorithm, params[["method"]]),

    # Adapted from IPO: add OpenMS isotope detection
    calcPPS = function(feat)
    {
        fTable <- featureTable(feat)

        ret <- list(featureCount = length(feat), nonRP = 0, RP = 0, PPS = 0)

        if (length(feat) == 0)
            return(ret)

        doOpenMS <- isoIdent == "OpenMS"
        if (!doOpenMS) # no need to find isotopes with OpenMS algo
        {
            xset <- getXCMSSet(feat, verbose = FALSE, exportedData = TRUE)
            peak_source <- utilsIPO$peaks_IPO(xset)[, c("mz", "rt", "sample", "into", "mzmin",
                                                        "mzmax", "rtmin", "rtmax"), drop = FALSE]
            if(isoIdent == "IPO")
                iso_mat <- utilsIPO$findIsotopes.IPO(xset, checkPeakShape)
            else
                iso_mat <- do.call(utilsIPO$findIsotopes.CAMERA, c(list(xset), CAMERAOpts))
        }

        isotope_abundance = 0.01108

        #calculating low intensity peaks
        for (anai in seq_along(analyses(feat)))
        {
            if (doOpenMS)
            {
                # intensities/masses of features without assigned isotopes
                intensities <- fTable[[anai]][isocount == 1, intensity]
                masses <- fTable[[anai]][isocount == 1, mz]
            }
            else
            {
                non_isos_peaks <- peak_source

                if (nrow(iso_mat) > 0)
                    non_isos_peaks <- peak_source[-unique(c(iso_mat)), , drop = FALSE]

                speaks <- non_isos_peaks[non_isos_peaks[,"sample"]==anai, , drop = FALSE]
                intensities <- speaks[,"into"]
                na_int <- is.na(intensities)
                intensities <- intensities[!na_int]

                if (length(intensities) > 0)
                {
                    masses <- speaks[!na_int, "mz"]
                    #floor((masses-2*CH3)/CH2) + 2
                }
            }

            if (length(intensities) > 0)
            {
                tmp <- intensities[order(intensities)]
                int_cutoff <- mean(tmp[1:max(round((length(tmp)/33),0),1)])

                maximum_carbon <- utilsIPO$calcMaximumCarbon(masses)
                carbon_probability <- maximum_carbon * isotope_abundance

                iso_int <- intensities * carbon_probability

                not_loq_peaks <- sum(iso_int > int_cutoff)
                ret$nonRP <- ret$nonRP + not_loq_peaks
            }
        }

        if (doOpenMS)
        {
            # isocount represent the number of isotopes collapsed in a feature. When
            # it's one, no other isotopes are present and hence its a NP. The
            # remaining (except the last, hence - 1) are assumed to be RP.
            ret$RP <- sum(unlist(lapply(fTable, function(ft) ft[isocount > 1, isocount - 1])))
        }
        else
            ret$RP <- length(unique(c(iso_mat)))

        if (ret[3] == 0)
            ret$PPS <- (ret$RP+1)^2/(ret$nonRP+1)
        else
            ret$PPS <- ret$RP^2/ret$nonRP

        return(ret)
    },

    getResponseScores = function(response) response$PPS,
    getFinalScore = function(oldr, newr) newr$PPS,

    calculateResponse = function(params, task, keepObject)
    {
        feat <- do.call(findFeatures, c(list(anaInfo, algorithm, verbose = FALSE), params))
        ret <- calcPPS(feat)

        if (keepObject)
            ret <- list(response = ret, object = feat)

        return(ret)
    },

    resultIncreased = function(history)
    {
        index <- length(history)
        if (history[[index]]$finalResult$response$PPS == 0 & index == 1)
        {
            warning("No isotopes have been detected!")
            return(FALSE)
        }
        if (index < 2)
            return(TRUE)
        if (history[[index-1]]$finalResult$response$PPS >= history[[index]]$finalResult$response$PPS)
            return(FALSE)
        return(TRUE)
    }
)


#' @param isoIdent Sets the algorithm used to identify isotopes. Valid values
#'   are: \code{"IPO"}, \code{"CAMERA"} and \code{"OpenMS"}. The latter can only
#'   be used when OpenMS is used to find features, and is highly recommended in
#'   this situation.
#' @param checkPeakShape Additional peak shape checking of isotopes. Only used
#'   if \code{isoIdent="IPO"}. Valid values: \code{"none"},
#'   \code{"borderIntensity"}, \code{"sinusCurve"} or \code{"normalDistr"}.
#' @param CAMERAOpts A \code{list} with additional arguments passed to
#'   \code{\link[CAMERA:findIsotopes-methods]{CAMERA::findIsotopes}} when \code{isoIdent="CAMERA"}.
#'
#' @rdname feature-optimization
#' @export
optimizeFeatureFinding <- function(anaInfo, algorithm, ..., templateParams = list(),
                                   paramRanges = list(),
                                   isoIdent = if (algorithm == "openms") "OpenMS" else "IPO",
                                   checkPeakShape = "none", CAMERAOpts = list(), maxIterations = 50, maxModelDeviation = 0.1)
{
    params <- list(...)

    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, add = ac)
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "envipick", "kpic2"), add = ac)
    assertOptimArgs(params, templateParams, paramRanges, maxIterations, maxModelDeviation, ac)
    checkmate::assertChoice(isoIdent, c("IPO", "CAMERA", "OpenMS"), add = ac)
    checkmate::assertChoice(checkPeakShape, c("none", "borderIntensity", "sinusCurve", "normalDistr"), add = ac)
    checkmate::assertList(CAMERAOpts, add = ac)
    checkmate::reportAssertions(ac)

    if (algorithm != "openms" && isoIdent == "OpenMS")
        stop("OpenMS isotope identification can only be used when OpenMS is used to find features.")

    fo <- switch(algorithm,
                 openms = featuresOptimizerOpenMS,
                 xcms = featuresOptimizerXCMS,
                 xcms3 = featuresOptimizerXCMS3,
                 envipick = featuresOptimizerEnviPick,
                 kpic2 = featuresOptimizerKPIC2)

    fo <- fo$new(anaInfo = anaInfo, algorithm = algorithm, isoIdent = isoIdent,
                 checkPeakShape = checkPeakShape, CAMERAOpts = CAMERAOpts)
    result <- fo$optimize(params, templateParams, paramRanges, maxIterations, maxModelDeviation)

    return(optimizationResult(algorithm = algorithm, paramSets = result$paramSets,
                              bestParamSet = result$bestParamSet))
}

#' @param method Method used by XCMS to find features (only if \code{algorithm="xcms"}).
#' @rdname feature-optimization
#' @export
generateFeatureOptPSet <- function(algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "envipick", "kpic2"))

    f <- switch(algorithm,
                openms = generateFeatureOptPSetOpenMS,
                xcms = generateFeatureOptPSetXCMS,
                xcms3 = generateFeatureOptPSetXCMS3,
                envipick = generateFeatureOptPSetEnviPick,
                kpic2 = generateFeatureOptPSetKPIC2)

    defs <- f(...)
    return(modifyList(defs, list(...)))
}

#' @rdname feature-optimization
#' @export
getDefFeaturesOptParamRanges <- function(algorithm, method = "centWave")
{
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "envipick", "kpic2"))

    if (algorithm == "openms")
        return(getDefFeaturesOptParamRangesOpenMS())
    else if (algorithm == "xcms")
        return(getDefFeaturesOptParamRangesXCMS(method))
    else if (algorithm == "xcms3")
        return(getDefFeaturesOptParamRangesXCMS3(method))
    else if (algorithm == "envipick")
        return(getDefFeaturesOptParamRangesEnviPick())
    return(getDefFeaturesOptParamRangesKPIC2())
}
