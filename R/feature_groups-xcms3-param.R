# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeatureGroupsXCMS3ParamDefs <- paramConfigDefsFact(list(
    # NOTE: PeakDensityParam() and LamaParama() cannot be constructed completely from defaults. For these, we set a
    # dummy value during construction and remove it after list conversion. The actual values are filled in prior to the
    # actual feature grouping.
    
    rtalign = list(
        default = TRUE,
        description = "If TRUE then retention times are aligned",
        type = "flag"
    ),
    loadRawData = list(
        default = TRUE,
        description = "If TRUE then raw data is loaded into memory",
        type = "flag"
    ),
    groupAlgo = list(
        default = "PeakDensity",
        description = "Algorithm to use for feature grouping (after alignment)",
        type = "choice",
        typeCheckArgs = list(choices = c("PeakDensity", "NearestPeaks", "MzClust"))
    ),
    alignAlgo = list(
        default = "Obiwarp",
        description = "Algorithm to use for retention time alignment",
        type = "choice",
        typeCheckArgs = list(choices = c("Obiwarp", "PeakGroups", "Lama"))
    ),
    preGroupAlgo = list(
        default = "PeakDensity",
        description = "Algorithm to use for feature grouping prior to alignment (only with peak groups alignment)",
        type = "choice",
        typeCheckArgs = list(choices = c("PeakDensity", "NearestPeaks", "MzClust"))
    ),
    groupPeakDensityParam = list(
        default = removeListEntries(as(xcms::PeakDensityParam(sampleGroups = "dummy"), "list"), "sampleGroups"),
        description = "PeakDensity XCMS parameters (use xcms::PeakDensityParam())",
        type = "list"
    ),
    groupNearestPeaksParam = list(
        default = as(xcms::NearestPeaksParam(), "list"),
        description = "NearestPeaks XCMS parameters (use xcms::NearestPeaksParam())",
        type = "list"
    ),
    groupMzClustParam = list(
        default = as(xcms::MzClustParam(), "list"),
        description = "MzClust XCMS parameters (use xcms::MzClustParam())",
        type = "list"
    ),
    alignObiwarpParam = list(
        default = as(xcms::ObiwarpParam(), "list"),
        description = "Obiwarp XCMS parameters (use xcms::ObiwarpParam())",
        type = "list"
    ),
    alignPeakGroupsParam = list(
        default = as(xcms::PeakGroupsParam(), "list"),
        description = "PeakGroups XCMS parameters (use xcms::PeakGroupsParam())",
        type = "list"
    ),
    alignLamaParama = list(
        default = removeListEntries(as(xcms::LamaParama(lamas = data.frame(mz = 1, rt = 1)), "list"), "lamas"),
        description = "Lama XCMS parameters (use xcms::LamaParama())",
        type = "list"
    ),
    preGroupPeakDensityParam = list(
        default = removeListEntries(as(xcms::PeakDensityParam(sampleGroups = "dummy"), "list"), "sampleGroups"),
        description = "PeakDensity XCMS parameters for pre-alignment grouping (use xcms::PeakDensityParam())",
        type = "list"
    ),
    preGroupNearestPeaksParam = list(
        default = as(xcms::NearestPeaksParam(), "list"),
        description = "NearestPeaks XCMS parameters for pre-alignment grouping (use xcms::NearestPeaksParam())",
        type = "list"
    ),
    preGroupMzClustParam = list(
        default = as(xcms::MzClustParam(), "list"),
        description = "MzClust XCMS parameters for pre-alignment grouping (use xcms::MzClustParam())",
        type = "list"
    ),    
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeatureGroupsXCMS3Param <- setClass("FeatureGroupsXCMS3Param", slots = c(lamas = "ANY"), contains = "param")
setMethod("initialize", "FeatureGroupsXCMS3Param", function(.Object, ..., lamas = NULL)
{
    ret <- callNextMethod(.Object, name = "FeatureGroupsXCMS3Param", baseName = "FeatureGroupsXCMS3Param",
                          description = "Parameters for XCMS3 feature grouping", version = "1.0",
                          definitions = getFeatureGroupsXCMS3ParamDefs(), ...)
    ret@lamas <- lamas # NOTE: assign here, as ... are passed to data slot in parent constructor
    return(ret)
})

setValidity("FeatureGroupsXCMS3Param", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    if (parsFilled$alignAlgo == "Lama") # UNDONE: needs XCMS4 interface
        return("Lama alignment algorithm is not yet supported.")
    return(TRUE)
})
