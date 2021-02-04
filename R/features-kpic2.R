#' @include features.R
NULL

makeKPIC2PeakInfo <- function(ft)
{
    cols <- intersect(names(ft), c("ret", "retmin", "retmax", "mz", "mzmin", "mzmax", "intensity", "area", "sn"))
    ret <- ft[, cols, with = FALSE]
    setnames(ret,
             c("ret", "retmin", "retmax", "intensity"),
             c("rt", "rtmin", "rtmax", "maxo"))
    if (!is.null(ret[["sn"]]))
        setnames(ret, "sn", "snr")
    else
        ret[, snr := NA_real_]
    ret[, mzrsd := NA_real_]
    setcolorder(ret, c("rt", "rtmin", "rtmax", "mz", "mzmin", "mzmax", "mzrsd", "maxo", "area", "snr"))
    ret <- as.matrix(ret)
    return(ret)
}

updatePICSet <- function(old, new)
{
    new@picsList <- Map(new@picsList, featureTable(old), featureTable(new), f = function(pics, oft, nft)
    {
        removed <- oft[!ID %in% nft$ID]$ID
        if (length(removed) > 0)
        {
            pics$pics <- pics$pics[-removed]
            pics$peaks <- pics$peaks[-removed]
            pics$peakinfo <- pics$peakinfo[-removed, ]
        }

        # update peakinfo
        pics$peakinfo <- makeKPIC2PeakInfo(nft)
        
        return(pics)
    })
    
    # ensure IDs are along rows in case features were removed
    new@features <- lapply(new@features, function(ft)
    {
        if (nrow(ft) > 0 && nrow(ft) < last(ft$ID))
            set(ft, j = "ID", value = seq_len(nrow(ft)))
        return(ft)
    })
    
    return(new)
}

#' @rdname features-class
#' @export
featuresKPIC2 <- setClass("featuresKPIC2", slots = list(picsList = "ANY"), contains = "features")

setMethod("initialize", "featuresKPIC2",
          function(.Object, ...) callNextMethod(.Object, algorithm = "kpic2", ...))

#' @rdname features-class
#' @export
setReplaceMethod("featureTable", "featuresKPIC2", function(obj, value)
{
    ret <- callNextMethod()
    ret <- updatePICSet(obj, ret)
    return(ret)
})

#' @rdname features-class
#' @export
setMethod("[", c("featuresKPIC2", "ANY", "missing", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., drop = drop)
    x@picsList <- x@picsList[names(x@features)]
    return(x)
})

#' @rdname features-class
#' @export
setReplaceMethod("[", c("featuresKPIC2", "ANY", "missing"), function(x, i, j, value)
{
    ret <- callNextMethod()
    ret <- updatePICSet(x, ret)
    return(ret)
})

#' @rdname features-class
#' @export
setReplaceMethod("[[", c("featuresKPIC2", "ANY", "missing"), function(x, i, j, value)
{
    ret <- callNextMethod()
    ret <- updatePICSet(x, ret)
    return(ret)
})

#' @rdname features-class
#' @export
setReplaceMethod("$", "featuresKPIC2", function(x, name, value)
{
    ret <- callNextMethod()
    ret <- updatePICSet(x, ret)
    return(ret)
})


#' @rdname feature-finding
#' @export
findfeaturesKPIC2 <- function(analysisInfo, kmeans, level = 1000, ..., parallel = TRUE, verbose = TRUE)
{
    # UNDONE: docs
    #   - mention that filter() doesn't alter KPIC object, but IDs can be used to retrace
    #       - or make function that gives synchronized object?
    
    checkPackage("KPIC", "https://github.com/hcji/KPIC2")
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    aapply(checkmate::assertFlag, . ~ kmeans + parallel + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    baseHash <- makeHash(kmeans, level, list(...))
    
    if (verbose)
        printf("Finding features with KPIC2 for %d analyses ...\n", nrow(analysisInfo))

    doKP <- function(ana, path)
    {
        inFile <- getMzMLOrMzXMLAnalysisPath(ana, path)
        hash <- makeHash(baseHash, makeFileHash(inFile))
        cachef <- loadCacheData("featuresKPIC2", hash)
        if (!is.null(cachef))
            return(cachef)
        
        raw <- KPIC::LoadData(inFile)
        pics <- do.call(if (kmeans) KPIC::getPIC.kmeans else KPIC::getPIC,
                        c(list(raw = raw, level = level), ...), verbose)
        pics <- KPIC::PICsplit(pics) # UNDONE: make optional?
        pics <- KPIC::getPeaks(pics)
        
        patRoon:::doProgress()
        
        return(pics)
    }
    
    if (parallel)
        allPics <- withProg(nrow(analysisInfo), future.apply::future_Map(doKP, analysisInfo$analysis,
                                                                          analysisInfo$path))
    else
        allPics <- withProg(nrow(analysisInfo), Map(doKP, analysisInfo$analysis, analysisInfo$path))
    
    ret <- importfeaturesKPIC2(allPics, analysisInfo)

    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(ret@features)
    }

    return(ret)
}

#' @details \code{importfeaturesKPIC2} converts features obtained with with the
#'   \pkg{KPIC2} package to a new \code{\link{features}} object.
#'
#' @param picsList A \code{list} with a \code{pics} objects obtained with
#'   \code{\link{getPIC}} or \code{\link{getPIC.kmeans}} for each analysis.
#'
#' @rdname feature-finding
#' @export
importfeaturesKPIC2 <- function(picsList, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(picsList, "list", min.len = 1, add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(picsList) != nrow(analysisInfo))
        stop("Length of picsList should be equal to number analyses")

    feat <- setNames(lapply(picsList, function(pics)
    {
        ret <- as.data.table(pics$peakinfo)
        setnames(ret, c("rt", "rtmin", "rtmax", "maxo", "snr"),
                 c("ret", "retmin", "retmax", "intensity", "sn"))
        ret[, ID := seq_len(.N)][]
        return(ret)
    }), analysisInfo$analysis)
    
    return(featuresKPIC2(picsList = picsList, features = feat, analysisInfo = analysisInfo))
}

#' @export
setMethod("getPICSet", "featuresKPIC2", function(obj, ...)
{
    return(unname(obj@picsList))
})

#' @export
setMethod("getPICSet", "features", function(obj, exportedData = TRUE)
{
    checkmate::assertFlag(exportedData)
    
    anaInfo <- analysisInfo(obj)
    fTable <- featureTable(obj)
    EICs <- if (exportedData) getEICsForFeatures(obj) else NULL
    return(lapply(seq_along(fTable), function(anai)
    {
        ret <- list()
        if (exportedData)
        {
            ret$path = getMzMLOrMzXMLAnalysisPath(anaInfo$analysis[anai], anaInfo$path[anai])
            ret$scantime <- loadSpectra(ret$path, verbose = FALSE)$header$retentionTime
            ret$pics <- Map(EICs[[anai]], fTable[[anai]]$mz, f = function(eic, mz)
            {
                setDT(eic)
                setnames(eic, "intensity", "int")
                eic[, mz := mz] # UNDONE? Could add actual m/z for each scan...
                eic[, scan := sapply(time, function(t) which.min(abs(t - ret$scantime)))]
                return(as.matrix(eic[, c("scan", "int", "mz"), with = FALSE]))
            })
            ret$peaks <- Map(ret$pics, fTable[[anai]]$intensity, f = function(pic, int)
            {
                # UNDONE: some dummy values here
                return(list(peakIndex = which.min(abs(pic[, "int"] - int)),
                            snr = NA_real_,
                            signals = int,
                            peakScale = 10))
            })
        }
        else
        {
            ret$path <- anaInfo$analysis[anai]
            ret$scantime <- integer()
            ret$pics <- ret$peaks <- list()
        }
        
        ret$peakinfo <- makeKPIC2PeakInfo(fTable[[anai]])

        return(ret)
    }))
})
