# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getComponentsRAMClustRParamDefs <- paramConfigDefsFact(list(
    st = list(
        default = NULL,
        description = "sigma t - time similarity decay value (when NULL computed as the half median for all chrom peak widths, passed to RAMClustR::ramclustR)",
        type = "number",
        lower = 0,
        null.ok = TRUE
    ),
    sr = list(
        default = NULL,
        description = "sigma r - correlational similarity decay value (NULL = auto, passed to RAMClustR::ramclustR)",
        type = "number",
        lower = 0,
        null.ok = TRUE
    ),
    maxt = list(
        default = 12,
        description = "maximum time difference to calculate retention similarity for - all values beyond this are assigned similarity of zero (passed to RAMClustR::ramclustR)",
        type = "number",
        lower = 0
    ),
    hmax = list(
        default = 0.3,
        description = "precut the tree at this height. See ?cutreeDynamicTree (passed to RAMClustR::ramclustR)",
        type = "number",
        lower = 0
    ),
    normalize = list(
        default = "TIC",
        description = "normalization of feature intensities (passed to RAMClustR::ramclustR)",
        type = "string"
    ),
    absMzDev = list(
        default = defaultLim("mz", "narrow"),
        description = "Absolute m/z deviation (pased as mzabs.error to RAMClustR::do.findmain)",
        type = "number",
        lower = 0
    ),
    relMzDev = list(
        default = defaultLim("mz", "narrow_rel"),
        description = "Relative m/z deviation (passed as ppm.error to RAMClustR::do.findmain)",
        type = "number",
        lower = 0
    ),
    minSize = list(
        default = 2,
        description = "Minimum number of feature groups to form a component",
        type = "count",
        positive = TRUE
    ),
    relMinReplicates = list(
        default = 0.5,
        description = "The features of a feature group in a component must be present in at least this relative fraction of replicates",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    RCExperimentVals = list(
        default = list(design = list(platform = "LC-MS"), instrument = list(MSlevs = 1)),
        description = "Experiment / instrument values used to construct ExpDes for RAMClustR",
        type = "list",
        names = "unique",
        any.missing = FALSE
    ),
    extraOptsRC = list(
        default = NULL,
        description = "Extra options passed to RAMClustR::ramclustR as a named character list",
        type = "list",
        types = "character",
        names = "unique",
        null.ok = TRUE
    ),
    extraOptsFM = list(
        default = NULL,
        description = "Extra options passed to RAMClustR::do.findmain as a named character list",
        type = "list",
        types = "character",
        names = "unique",
        null.ok = TRUE
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

#' @export
ComponentsRAMClustRParam <- setClass("ComponentsRAMClustRParam", contains = "param")

setMethod("initialize", "ComponentsRAMClustRParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsRAMClustRParam", baseName = "ComponentsRAMClustRParam",
                   description = "Parameters for RAMClustR component generation", version = "1.0",
                   definitions = getComponentsRAMClustRParamDefs(), ...)
})
