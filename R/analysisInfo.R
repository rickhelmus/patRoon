# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' AnalysisInfo data.frame methods
#'
#' Various parsing and plotting functions for the analysisInfo data.frame.
#'
#' @param obj An \code{analysisInfo} data.frame object as obtained by \link[patRoon]{generateAnalysisInfo} function.
#' @param MSLevel Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to obtain traces.
#' @param retentionRange Range of retention time (in seconds), m/z, respectively. Should be a numeric vector with length
#'   of two containing the min/max values. The maximum can be Inf to specify no maximum range. Set to NULL to skip this
#'   step.
#' @param \dots Further arguments passed to \code{\link[graphics]{plot}}.
#'
#' @author Ricardo Cunha (\email{cunha@@iuta.de}) and Rick Helmus (\email{r.helmus@@uva.nl})
#'
#' @name analysisinfo-dataframe
NULL

#' @describeIn analysisinfo-dataframe Obtain the total ion chromatogram/s (TICs) of the analyses.
#' @export
setMethod("getTICs", "data.frame", function(obj, retentionRange = NULL, MSLevel = 1)
{
    ac <- checkmate::makeAssertCollection()
    obj <- assertAndPrepareAnaInfo(obj, add = ac)
    assertRange(retentionRange, null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    checkmate::reportAssertions(ac)
    
    return(rbindlist(applyMSData(obj, showProgress = FALSE, func = function(ana, path, backend)
    {
        openMSReadBackend(backend, path)
        md <- getMSMetadata(backend, MSLevel)
        if (!is.null(retentionRange))
            md <- md[numGTE(md$time, retentionRange[1]) & numLTE(md$time, retentionRange[2]), ]
        return(data.table(analysis = ana, replicate = obj$replicate[match(ana, obj$analysis)], ret = md$time,
                          MSLevel = MSLevel, intensity = md$TIC))
    })))
})

#' @describeIn analysisinfo-dataframe Obtain the base peak chromatogram/s (BPCs) of the analyses.
#' @export
setMethod("getBPCs", "data.frame", function(obj, retentionRange = NULL, MSLevel = 1)
{
    ac <- checkmate::makeAssertCollection()
    obj <- assertAndPrepareAnaInfo(obj, add = ac)
    assertRange(retentionRange, null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    checkmate::reportAssertions(ac)
    
    return(rbindlist(applyMSData(obj, showProgress = FALSE, func = function(ana, path, backend)
    {
        openMSReadBackend(backend, path)
        md <- getMSMetadata(backend, MSLevel)
        if (!is.null(retentionRange))
            md <- md[numGTE(md$time, retentionRange[1]) & numLTE(md$time, retentionRange[2]), ]
        return(data.table(analysis = ana, replicate = obj$replicate[match(ana, obj$analysis)], ret = md$time,
                          MSLevel = MSLevel, intensity = md$BPC))
    })))
})

#' @describeIn analysisinfo-dataframe Plots the TICs of the analyses.
#' @param retMin Plot retention time in minutes (instead of seconds).
#' @param title Character string used for title of the plot. If \code{NULL} a title will be automatically generated.
#' @param groupBy Sets the automatic colour selection: "none" for a single colour or "analyses"/"rGroups" for a
#'   distinct colour per analysis or analysis replicate.
#' @param showLegend Plot a legend if TRUE.
#' @template plot-lim

#' @export
setMethod("plotTICs", "data.frame", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL, 
                                             groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    doPlotHeaders(as.data.table(obj), what = "tic", retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim,
                  ylim, ...)
})


#' @describeIn analysisinfo-dataframe Plots the BPCs of the analyses.
#' @export
setMethod("plotBPCs", "data.frame", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL,
                                             groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    doPlotHeaders(as.data.table(obj), what = "bpc", retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim,
                  ylim, ...)
})

#' @rdname analysisinfo-dataframe
#' @export
setMethod("plotTICs", "data.table", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL, 
                                             groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    doPlotHeaders(obj, what = "tic", retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
})

#' @rdname analysisinfo-dataframe
#' @describeIn analysisinfo-dataframe Plots the BPCs of the analyses.
#' @export
setMethod("plotBPCs", "data.table", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL,
                                             groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    doPlotHeaders(obj, what = "bpc", retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
})
