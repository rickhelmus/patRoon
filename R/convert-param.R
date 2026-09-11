# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getConvertMSFilesPWizParamDefs <- paramConfigDefsFact(list(
    centroidVendor = list(
        default = TRUE,
        description = "Use vendor centroiding algorithms",
        type = "flag"
    ),
    minIntensity = list(
        default = 0,
        description = "Minimum mass peak intensity",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    filters = list(
        default = NULL,
        description = "MSConvert filters",
        type = "character",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra command line options",
        type = "character",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    ),
    PWizBatchSize = list(
        default = 1,
        description = "Number of analyses per MSConvert call",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    )
))

#' @export
ConvertMSFilesPWizParam <- setClass("ConvertMSFilesPWizParam", contains = "param")
setMethod("initialize", "ConvertMSFilesPWizParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ConvertMSFilesPWizParam", baseName = "ConvertMSFilesPWizParam",
                   description = "Parameters for ProteoWizard/MSConvert MS file conversion",
                   version = "1.0", definitions = getConvertMSFilesPWizParamDefs(), ...)
})

getConvertMSFilesOpenMSParamDefs <- paramConfigDefsFact(list(
    extraOpts = list(
        default = NULL,
        description = "Extra command line options",
        type = "character",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    )
))

#' @export
ConvertMSFilesOpenMSParam <- setClass("ConvertMSFilesOpenMSParam", contains = "param")
setMethod("initialize", "ConvertMSFilesOpenMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ConvertMSFilesOpenMSParam", baseName = "ConvertMSFilesOpenMSParam",
                   description = "Parameters for OpenMS MS file conversion", version = "1.0",
                   definitions = getConvertMSFilesOpenMSParamDefs(), ...)
})

getConvertMSFilesBrukerParamDefs <- paramConfigDefsFact(list())

#' @export
ConvertMSFilesBrukerParam <- setClass("ConvertMSFilesBrukerParam", contains = "param")
setMethod("initialize", "ConvertMSFilesBrukerParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ConvertMSFilesBrukerParam", baseName = "ConvertMSFilesBrukerParam",
                   description = "Parameters for Bruker MS file conversion", version = "1.0",
                   definitions = getConvertMSFilesBrukerParamDefs(), ...)
})

getConvertMSFilesIMSCollapseParamDefs <- paramConfigDefsFact(list(
    mzRange = list(
        default = NULL,
        description = "m/z range",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    mobilityRange = list(
        default = NULL,
        description = "Ion mobility range",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    smoothWindow = list(
        default = 0,
        description = "Smoothing window",
        type = "count"
    ),
    halfWindow = list(
        default = 2,
        description = "Centroiding half window",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    maxGap = list(
        default = 0.005,
        description = "Maximum centroiding gap",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    clusterMethod = list(
        default = "distance_mean",
        description = "MS/MS clustering method",
        type = "choice",
        typeCheckArgs = list(choices = c("bin", "distance_mean", "distance_point", "hclust"))
    ),
    mzWindow = list(
        default = defaultLim("mz", "medium"),
        description = "MS/MS clustering m/z window",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minIntensityIMS = list(
        default = 0,
        description = "Minimum IMS peak intensity",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    includeMSMS = list(
        default = FALSE,
        description = "Include MS/MS spectra in the export",
        type = "flag"
    )
))

#' @export
ConvertMSFilesIMSCollapseParam <- setClass("ConvertMSFilesIMSCollapseParam", contains = "param")
setMethod("initialize", "ConvertMSFilesIMSCollapseParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ConvertMSFilesIMSCollapseParam", baseName = "ConvertMSFilesIMSCollapseParam",
                   description = "Parameters for IMS-collapsing MS file conversion", version = "1.0",
                   definitions = getConvertMSFilesIMSCollapseParamDefs(), ...)
})

getConvertMSFilesTIMSCONVERTParamDefs <- paramConfigDefsFact(list(
    centroidRaw = list(
        default = FALSE,
        description = "Use raw mode for centroided data",
        type = "flag"
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra command line options",
        type = "character",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    ),
    virtualenv = list(
        default = "patRoon-TIMSCONVERT",
        description = "TIMSCONVERT Python virtual environment",
        type = "string",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    )
))

#' @export
ConvertMSFilesTIMSCONVERTParam <- setClass("ConvertMSFilesTIMSCONVERTParam", contains = "param")
setMethod("initialize", "ConvertMSFilesTIMSCONVERTParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ConvertMSFilesTIMSCONVERTParam", baseName = "ConvertMSFilesTIMSCONVERTParam",
                   description = "Parameters for TIMSCONVERT MS file conversion", version = "1.0",
                   definitions = getConvertMSFilesTIMSCONVERTParamDefs(), ...)
})
