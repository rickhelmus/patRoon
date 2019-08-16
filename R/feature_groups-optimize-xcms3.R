#' @include main.R
#' @include feature_groups-optimize-xcms.R
NULL

# UNDONE: do the subset/subsetAdjust/minFraction parameters currently work for ret correction?

featureGroupsOptimizerXCMS3 <- setRefClass("featureGroupsOptimizerXCMS3", contains = "featureGroupsOptimizerXCMS")

featureGroupsOptimizerXCMS3$methods(
    
    convertOptToCallParams = function(params)
    {
        params <- callSuper(params)
        gMethod <- params$groupArgs$method; rtMethod <- params$retcorArgs$method
        ret <- list()
        
        if (!is.null(gMethod))
        {
            gParams <- params$groupArgs[names(params$groupArgs) != "method"]
            gParams$sampleGroups <- analysisInfo(features)$group
            if (gMethod == "density")
            {
                # update to newstyle naming
                gParams <- setListNamesIfPresent(gParams,
                                                 c("minfrac", "minsamp", "mzwid", "max"),
                                                 c("minFraction", "minSamples", "binSize", "maxFeatures"))
                ret$groupParam <- do.call(xcms::PeakDensityParam, gParams)
            }
            else if (gMethod == "nearest")
            {
                gParams <- setListNamesIfPresent(gParams,
                                                 c("mzCheck", "rtCheck"),
                                                 c("absMz", "absRt"))
                ret$groupParam <- do.call(xcms::NearestPeaksParam, gParams)
            }
            else if (gMethod == "mzClust") # UNDONE: not really supported in XCMS base class
            {
                gParams <- setListNamesIfPresent(gParams,
                                                 c("mzppm", "mzabs", "minsamp", "minfrac"),
                                                 c("ppm", "absMz", "minSamples", "minFraction"))
                ret$groupParam <- do.call(xcms::MzClustParam, gParams)
            }
        }
        
        if (!is.null(rtMethod))
        {
            rtParams <- params$retcorArgs[names(params$retcorArgs) != "method"]
            if (rtMethod == "obiwarp")
            {
                rtParams <- setListNamesIfPresent(rtParams,
                                                  c("profStep", "center", "distFunc"),
                                                  c("binSize", "centerSample", "distFun"))
                ret$retAlignParam <- do.call(xcms::ObiwarpParam, rtParams)
            }
            else if (rtMethod == "peakgroups")
            {
                # UNDONE: old style missing param doesn't translate
                rtParams <- setListNamesIfPresent(rtParams,
                                                  c("extra"),
                                                  c("extraPeaks"))
                ret$retAlignParam <- do.call(xcms::PeakGroupsParam, rtParams)
            }
        }
        
        return(ret)
    },
    
    checkInitialParamsDontUse = function(params)
    {
        ret <- list()
        
        # convert to old style
        if (isXCMSClass(params$groupParam, "PeakDensityParam"))
            ret$groupArgs <- list(minfrac = xcms::minFraction(params$groupParam),
                                  minsamp = xcms::minSamples(params$groupParam),
                                  bw = xcms::bw(params$groupParam),
                                  mzwid = xcms::binSize(params$groupParam),
                                  max = xcms::maxFeatures(params$groupParam))
        else if (isXCMSClass(params$groupParam, "NearestPeaksParam"))
        {
            ret$groupArgs <- list(mzVsRTbalance = xcms::mzVsRtBalance(params$groupParam),
                                  mzCheck = xcms::absMz(params$groupParam),
                                  rtCheck = xcms::absRt(params$groupParam),
                                  kNN = xcms::kNN(params$groupParam))
        }
        else if (isXCMSClass(params$groupParam, "MzClustParam"))
        {
            ret$groupArgs <- list(mzppm = xcms::ppm(params$groupParam),
                                  mzabs = xcms::absMz(params$groupParam),
                                  minsamp = xcms::minSamples(params$groupParam),
                                  minfrac = xcms::minFraction(params$groupParam))
        }
        
        if (isXCMSClass(params$retAlignParam, "ObiwarpParam"))
            ret$retcorArgs <- list(profStep = xcms::binSize(params$retAlignParam),
                                   center = xcms::centerSample(params$retAlignParam),
                                   response = xcms::response(params$retAlignParam),
                                   distFunc = xcms::distFun(params$retAlignParam),
                                   gapInit = xcms::gapInit(params$retAlignParam),
                                   gapExtend = xcms::gapExtend(params$retAlignParam),
                                   factorDiag = xcms::factorDiag(params$retAlignParam),
                                   factorGap = xcms::factorGap(params$retAlignParam),
                                   localAlignment = xcms::localAlignment(params$retAlignParam),
                                   initPenalty = xcms::initPenalty(params$retAlignParam))
        else if (isXCMSClass(params$retAlignParam, "peakgroups"))
            ret$retcorArgs <- list(extra = xcms::extraPeaks(params$retAlignParam),
                                   smooth = xcms::smooth(params$retAlignParam),
                                   span = xcms::span(params$retAlignParam),
                                   family = xcms::family(params$retAlignParam))
        
        return(callSuper(ret))
    }
)

generateFGroupsOptPSetXCMS3 <- generateFGroupsOptPSetXCMS
getDefFGroupsOptParamRangesXCMS3 <- getDefFGroupsOptParamRangesXCMS
