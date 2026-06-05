# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

# UNDONE: some XCMS objects cannot be default constructed, set to NULL? Construct with given parameter? How to handle this in the GUI?
# --> this is only for xcms::PeakDensityParam() (needs repl groups) and Lama (needs ref mz/rt values) --> set density automatically and not support lama for now?
getFeatureGroupsXCMS3ParamDefs <- paramConfigDefsFact(list(
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
        default = "density",
        description = "Algorithm to use for feature grouping (after alignment)",
        type = "choice",
        choices = c("PeakDensity", "NearestPeaks", "MzClust")
    ),
    alignAlgo = list(
        default = "obiwarp",
        description = "Algorithm to use for retention time alignment",
        type = "choice",
        choices = c("Obiwarp", "PeakGroups", "Lama")
    ),
    preGroupAlgo = list(
        default = "density",
        description = "Algorithm to use for feature grouping prior to alignment (only with peak groups alignment)",
        type = "choice",
        choices = c("PeakDensity", "NearestPeaks", "MzClust")
    ),
    groupPeakDensityParam = list(
        default = as(xcms::PeakDensityParam(), "list"),
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
    alignLamaParam = list(
        default = as(xcms::LamaParama(), "list"),
        description = "Lama XCMS parameters (use xcms::LamaParama())",
        type = "list"
    ),
    preGroupPeakDensityParam = list(
        default = as(xcms::PeakDensityParam(), "list"),
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

FeatureGroupsXCMS3Param <- setClass("FeatureGroupsXCMS3Param", contains = "param")
setMethod("initialize", "FeatureGroupsXCMS3Param", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeatureGroupsXCMS3Param", baseName = "FeatureGroupsXCMS3Param",
                   description = "Parameters for XCMS3 feature grouping", version = "1.0",
                   definitions = getFeatureGroupsXCMS3ParamDefs(), ...)
})