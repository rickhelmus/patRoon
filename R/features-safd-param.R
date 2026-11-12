# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesSAFDParamDefs <- paramConfigDefsFact(list(
    prefCentroid = list(
        advanced = TRUE,
        default = FALSE,
        description = "Set to TRUE to prefer centroided data over other MS data specified in analysisInfo",
        type = "flag"
    ),
    mzRange = list(
        default = c(0, 400),
        description = "The m/z window to be imported",
        type = "numeric",
        typeCheckArgs = list(lower = 0, finite = TRUE, any.missing = FALSE, len = 2)
    ),
    maxNumbIter = list(
        advanced = TRUE,
        default = 1000L,
        description = "Maximum number of iterations (max_numb_iter parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    maxTPeakW = list(
        advanced = TRUE,
        default = 300L,
        description = "Maximum number of scans within a feature (max_t_peak_w parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    resolution = list(
        advanced = TRUE,
        default = 30000L,
        description = "Expected MS resolution (res parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    minMSW = list(
        advanced = TRUE,
        default = 0.02,
        description = "Minimum MS peak width (min_ms_w parameter of safd_s3D)",
        type = "number"
    ),
    RThreshold = list(
        advanced = TRUE,
        default = 0.75,
        description = "R2 threshold for Gaussian fit (r_thresh parameter of safd_s3D)",
        type = "number"
    ),
    minInt = list(
        advanced = TRUE,
        default = 2000L,
        description = "Minimum feature intensity (min_int parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    sigIncThreshold = list(
        advanced = TRUE,
        default = 5L,
        description = "Parameter for detecting overlapping features in time domain (sig_inc_thresh parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    S2N = list(
        advanced = TRUE,
        default = 2L,
        description = "Signal to noise threshold (S2N parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    minPeakWS = list(
        advanced = TRUE,
        default = 3L,
        description = "Minimum peak chrom. peak width (min_peak_w_s parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    centroidMethod = list(
        advanced = TRUE,
        default = "RFM",
        description = "Centroid method (method parameter of safd_s3d_cent)",
        type = "choice",
        typeCheckArgs = list(choices = c("RFM", "BG", "CT"))
    ),
    centroidDM = list(
        advanced = TRUE,
        default = 0.005,
        description = "mdm parameter of safd_s3d_cent",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    verbose = list(
        advanced = TRUE,
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeaturesSAFDParam <- setClass("FeaturesSAFDParam", contains = "param")
setMethod("initialize", "FeaturesSAFDParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesSAFDParam", baseName = "safd", subDir = "features",
                   description = "Parameters for SAFD feature detection", version = "1.0",
                   definitions = getFeaturesSAFDParamDefs(), ...)
})
