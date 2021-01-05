#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresKPIC2 <- setClass("featuresKPIC2", slots = list(picsList = "ANY"), contains = "features")

setMethod("initialize", "featuresKPIC2",
          function(.Object, ...) callNextMethod(.Object, algorithm = "kpic2", ...))

#' @rdname features-class
#' @export
setMethod("[", c("featuresKPIC2", "ANY", "missing", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., drop = drop)
    x@picsList <- x@picsList[names(x@features)]
    return(x)
})


#' @rdname feature-finding
#' @export
findfeaturesKPIC2 <- function(analysisInfo, kmeans, level = 1000, ..., verbose = TRUE)
{
    # UNDONE: docs
    #   - mention that filter() doesn't alter KPIC object, but IDs can be used to retrace
    #       - or make function that gives synchronized object?
    
    checkPackage("KPIC", "https://github.com/hcji/KPIC2")
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    aapply(checkmate::assertFlag, . ~ kmeans + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    baseHash <- makeHash(kmeans, level, list(...))
    
    if (verbose)
        printf("Finding features with KPIC2 for %d analyses ...\n", nrow(analysisInfo))

    prog <- progressr::progressor(steps = nrow(analysisInfo), enable = verbose)

    allPics <- future.apply::future_Map(analysisInfo$analysis, analysisInfo$path, f = function(ana, path)
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
        
        prog()
        
        return(pics)
    })
    
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
        setnames(ret, c("rt", "rtmin", "rtmax", "maxo"),
                 c("ret", "retmin", "retmax", "intensity"))
        ret[, ID := seq_len(.N)][]
        return(ret)
    }), analysisInfo$analysis)
    
    return(featuresKPIC2(picsList = picsList, features = feat, analysisInfo = analysisInfo))
}
