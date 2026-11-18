# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesXCMS3ParamDefs <- paramConfigDefsFact(list(
    # BUG/UNDONE: the as() calls can sometimes fail
    algo = list(
        advanced = FALSE,
        default = "CentWave",
        description = "Algorithm to use for feature detection",
        type = "choice",
        typeCheckArgs = list(choices = c("CentWave", "CentWavePredIso", "MatchedFilter", "Massifquant", "MSW"))
    ),
    CentWaveParam = list(
        advanced = FALSE,
        default = as(xcms::CentWaveParam(), "list"),
        description = "CentWave XCMS parameters (use xcms::CentWaveParam())",
        type = "list"
    ),
    CentWavePredIsoParam = list(
        advanced = FALSE,
        default = as(xcms::CentWavePredIsoParam(), "list"),
        description = "CentWavePredIso XCMS parameters (use xcms::CentWavePredIsoParam())",
        type = "list"
    ),
    MatchedFilterParam = list(
        advanced = FALSE,
        default = as(xcms::MatchedFilterParam(), "list"),
        description = "MatchedFilter XCMS parameters (use xcms::MatchedFilterParam())",
        type = "list"
    ),
    MassifQuantParam = list(
        advanced = FALSE,
        default = as(xcms::MassifquantParam(), "list"),
        description = "MassifQuant XCMS parameters (use xcms::MassifquantParam())",
        type = "list"
    ),
    MSWParam = list(
        advanced = FALSE,
        default = as(xcms::MSWParam(), "list"),
        description = "MSW XCMS parameters (use xcms::MSWParam())",
        type = "list"
    ),
    extraOpts = list(
        advanced = TRUE,
        default = NULL,
        description = "List with extra options passed to xcms::findChromPeaks()",
        type = "list",
        typeCheckArgs = list(names = "unique", null.ok = TRUE)
    ),
    verbose = list(
        advanced = TRUE,
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeaturesXCMS3Param <- setClass("FeaturesXCMS3Param", contains = "param")
setMethod("initialize", "FeaturesXCMS3Param", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesXCMS3Param", baseName = "xcms3", subDir = "features",
                   description = "Parameters for XCMS3 feature detection", version = "1.0",
                   definitions = getFeaturesXCMS3ParamDefs(), ...)
})
