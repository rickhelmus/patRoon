# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesSAFDParamDefs <- paramConfigDefsFact(list(
    prefCentroid = list(
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
        default = 1000L,
        description = "Maximum number of iterations (max_numb_iter parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    maxTPeakW = list(
        default = 300L,
        description = "Maximum number of scans within a feature (max_t_peak_w parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    resolution = list(
        default = 30000L,
        description = "Expected MS resolution (res parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    minMSW = list(
        default = 0.02,
        description = "Minimum MS peak width (min_ms_w parameter of safd_s3D)",
        type = "number"
    ),
    RThreshold = list(
        default = 0.75,
        description = "R2 threshold for Gaussian fit (r_thresh parameter of safd_s3D)",
        type = "number"
    ),
    minInt = list(
        default = 2000L,
        description = "Minimum feature intensity (min_int parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    sigIncThreshold = list(
        default = 5L,
        description = "Parameter for detecting overlapping features in time domain (sig_inc_thresh parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    S2N = list(
        default = 2L,
        description = "Signal to noise threshold (S2N parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    minPeakWS = list(
        default = 3L,
        description = "Minimum peak chrom. peak width (min_peak_w_s parameter of safd_s3D)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    centroidMethod = list(
        default = "RFM",
        description = "Centroid method (method parameter of safd_s3d_cent)",
        type = "choice",
        typeCheckArgs = list(choices = c("RFM", "BG", "CT"))
    ),
    centroidDM = list(
        default = 0.005,
        description = "mdm parameter of safd_s3d_cent",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeaturesSAFDParam <- setClass("FeaturesSAFDParam", contains = "param")
setMethod("initialize", "FeaturesSAFDParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesSAFDParam", baseName = "FeaturesSAFDParam",
                   description = "Parameters for SAFD feature detection", version = "1.0",
                   definitions = getFeaturesSAFDParamDefs(), ...)
})