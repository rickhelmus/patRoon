# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getCompoundsLibraryParamDefs <- paramConfigDefsFact(list(
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for annotation similarity calculations",
        type = "specSimParams"
    ),
    minSim = list(
        default = 0.75,
        description = "The minimum spectral similarity for candidate records",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minAnnSim = list(
        default = 0.75,
        description = "The minimum spectral similarity of a record for it to be used to find annotations",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    absMzDev = list(
        default = defaultLim("mz", "narrow"),
        description = "The maximum absolute m/z deviation between the feature group and library record m/z values",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    checkIons = list(
        default = "adduct",
        description = "Exclude library records with different adduct or MS ionization polarity",
        type = "choice",
        typeCheckArgs = list(choices = c("adduct", "polarity", "none"))
    ),
    spectrumType = list(
        default = "MS2",
        description = "The spectrum type(s) to use from the library (NULL to use all)",
        type = "character",
        typeCheckArgs = list(min.len = 1, min.chars = 1, null.ok = TRUE)
    ),
    specSimParamsLib = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for pre-treatment of library spectra",
        type = "specSimParams"
    ),
    minIMSSpecSim = list(
        default = 0,
        description = "Minimum spectral similarity to copy annotations from the IMS precursor",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    setThreshold = list(
        default = 0,
        description = "Minimum abundance for a candidate among all sets",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setThresholdAnn = list(
        default = 0,
        description = "As setThreshold, but only taking sets into account in which features of the feature group are present",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setAvgSpecificScores = list(
        default = FALSE,
        description = "If TRUE, average set specific scores (e.g. MS/MS match)",
        type = "flag"
    )
))

#' @export
CompoundsLibraryParam <- setClass("CompoundsLibraryParam", contains = "param")

setMethod("initialize", "CompoundsLibraryParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "CompoundsLibraryParam", baseName = "CompoundsLibraryParam",
                   description = "Parameters for MS library compound generation", version = "1.0",
                   definitions = getCompoundsLibraryParamDefs(), ...)
})
