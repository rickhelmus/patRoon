#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresEnviPick <- setClass("featuresEnviPick", contains = "features")

setMethod("initialize", "featuresEnviPick",
          function(.Object, ...) callNextMethod(.Object, algorithm = "envipick", ...))


#' @details \code{findFeaturesEnviPick} uses the
#'   \code{\link[enviPick]{enviPickwrap}}. function from the \pkg{enviPick} R
#'   package to extract features.
#'
#' @note \code{findFeaturesEnviPick} Requires analysis files to be in the
#'   \code{mzXML} format.
#'
#' @rdname feature-finding
#' @export
findFeaturesEnviPick <- function(analysisInfo, ..., verbose = TRUE)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzXML")

    anaCount <- nrow(analysisInfo)
    ret <- featuresEnviPick(analysisInfo = analysisInfo)

    if (verbose)
    {
        printf("Finding features with enviPick for %d analyses ...\n", anaCount)
        prog <- openProgBar(0, anaCount)
    }

    fts <- list()
    for (i in seq_len(anaCount))
    {
        fp <- getMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i])
        hash <- makeHash(makeFileHash(fp), list(...))
        f <- loadCacheData("featuresEnviPick", hash)

        if (is.null(f))
        {
            invisible(utils::capture.output(ep <- enviPick::enviPickwrap(fp, ...)))
            f <- importEnviPickPeakList(ep$Peaklist)
            saveCacheData("featuresEnviPick", f, hash)
        }

        fts[[analysisInfo$analysis[i]]] <- f

        if (verbose)
            setTxtProgressBar(prog, i)
    }

    ret@features <- fts

    if (verbose)
    {
        close(prog)
        printf("Done!\n")
        printFeatStats(fts)
    }

    return(ret)
}

#' @details \code{importFeaturesEnviMass} imports features from a project
#'   generated by the \pkg{enviMass} package. NOTE: this functionality has only
#'   been tested with older versions of \pkg{enviMass}.
#'
#' @param enviProjPath The path of the enviMass project.
#' @rdname feature-finding
#' @export
importFeaturesEnviMass <- function(analysisInfo, enviProjPath)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzXML", add = ac)
    checkmate::assertDirectoryExists(enviProjPath, "r", add = ac)
    checkmate::assertDirectoryExists(file.path(enviProjPath, "peaklist"), "r", .var.name = "enviProjPath", add = ac)
    checkmate::reportAssertions(ac)

    cat("Importing features from enviMass...\n")

    ret <- featuresEnviPick(analysisInfo = analysisInfo)

    fts <- list()
    scount <- nrow(analysisInfo)
    prog <- openProgBar(0, scount)

    for (i in seq_len(nrow(analysisInfo)))
    {
        load(file.path(enviProjPath, "peaklist", analysisInfo$analysis[i])) # load into 'peaklist'
        fts[[analysisInfo$analysis[i]]] <- importEnviPickPeakList(as.data.frame(peaklist))
        setTxtProgressBar(prog, i)
    }

    ret@features <- fts

    setTxtProgressBar(prog, scount)
    close(prog)

    cat("Done!\n")

    return(ret)
}

importEnviPickPeakList <- function(peaklist)
{
    # peaklist is a single number (zero) when no results
    if (length(peaklist) == 1 || nrow(peaklist) == 0)
        return(data.table(ID = character(), ret = numeric(), mz = numeric(), intensity = numeric(),
                          area = numeric(), retmin = numeric(), retmax = numeric(), mzmin = numeric(),
                          mzmax = numeric()))

    ft <- as.data.table(peaklist)

    setnames(ft, c("m/z", "max_int", "sum_int", "RT", "minRT", "maxRT", "peak_ID"),
             c("mz", "intensity", "area", "ret", "retmin", "retmax", "ID"))

    # Estimate mzrange from variance
    s <- sqrt(ft$`var_m/z`)
    ft[, mzmin := mz - 2*s]
    ft[, mzmax := mz + 2*s]

    return(ft[, c("ID", "ret", "mz", "intensity", "area", "retmin", "retmax", "mzmin", "mzmax")])
}
