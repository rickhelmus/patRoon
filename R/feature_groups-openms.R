#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname feature-grouping
#' @export
featureGroupsOpenMS <- setClass("featureGroupsOpenMS", contains = "featureGroups")

#' @details \code{groupFeaturesOpenMS} uses the OpenMS tools for grouping of
#'   features (see \url{http://www.openms.de}). Retention times may be aligned
#'   by the
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_MapAlignerPoseClustering.html}{MapAlignerPoseClustering}
#'    TOPP tool. Grouping is achieved by either the
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureLinkerUnlabeled.html}{FeatureLinkerUnlabeled}
#'    or
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureLinkerUnlabeledQT.html}{FeatureLinkerUnlabeledQT}
#'    TOPP tools.
#'
#' @param QT If enabled, use \code{FeatureLinkerUnlabeledQT} instead of
#'   \code{FeatureLinkerUnlabeled} for feature grouping.
#' @param maxAlignRT,maxAlignMZ Used for retention alignment. Maximum retention
#'   time or m/z difference (seconds/Dalton) for feature pairing. Sets
#'   \code{-algorithm:pairfinder:distance_RT:max_difference} and
#'   \code{-algorithm:pairfinder:distance_MZ:max_difference} otpions,
#'   respectively.
#' @param maxGroupRT,maxGroupMZ as \code{maxAlignRT} and \code{maxAlignMZ}, but
#'   for grouping of features. Sets \code{-algorithm:distance_RT:max_difference}
#'   and \code{-algorithm:distance_MZ:max_difference} options, respectively.
#'
#' @references \insertRef{Rst2016}{patRoon}
#'
#' @rdname feature-grouping
#' @export
groupFeaturesOpenMS <- function(feat, rtalign = TRUE, QT = FALSE, maxAlignRT = 30, maxAlignMZ = 0.005, maxGroupRT = 12,
                                maxGroupMZ = 0.005)
{
    # UNDONE: allow extra options for aligning/grouping?

    hash <- makeHash(feat, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ)
    cachefg <- loadCacheData("featureGroupsOpenMS", hash)
    if (!is.null(cachefg))
        return(cachefg)

    cat("Grouping features with OpenMS...\n===========\n")

    cfile <- tempfile("cons", fileext = ".consensusXML")
    generateConsensusXML(feat, cfile, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ)
    fgimp <- importConsensusXML(feat, cfile)

    ret <- featureGroupsOpenMS(groups = fgimp$groups, groupInfo=fgimp$gInfo, analysisInfo = analysisInfo(feat),
                               features=feat, ftindex=fgimp$ftindex)

    saveCacheData("featureGroupsOpenMS", ret, hash)

    cat("\n===========\nDone!\n")
    return(ret)
}

setMethod("generateConsensusXML", "features", function(feat, out, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ)
{
    sGroup <- analysisInfo(feat)
    fts <- featureTable(feat)

    edtaFiles <- list()
    featFiles <- list()
    for (datafile in sGroup$analysis)
    {
        featFiles[[datafile]] <- tempfile(datafile, fileext = ".featureXML")
        writeFeatureXML(fts[[datafile]], featFiles[[datafile]])
    }

    if (rtalign)
    {
        settings <- c("-algorithm:max_num_peaks_considered", -1,
                      "-algorithm:superimposer:mz_pair_max_distance", maxAlignMZ,
                      "-algorithm:superimposer:num_used_points", 10000,
                      "-algorithm:pairfinder:distance_RT:max_difference", maxAlignRT,
                      "-algorithm:pairfinder:distance_MZ:max_difference", maxAlignMZ,
                      "-algorithm:pairfinder:distance_MZ:unit", "Da")
        executeCommand(getCommandWithOptPath("MapAlignerPoseClustering", "OpenMS"),
                       c(settings, "-in", featFiles, "-out", featFiles))
    }

    settings <- c("-algorithm:distance_RT:max_difference", maxGroupRT,
                  "-algorithm:distance_MZ:max_difference", maxGroupMZ,
                  "-algorithm:distance_MZ:unit", "Da")
    executeCommand(getCommandWithOptPath(if (QT) "FeatureLinkerUnlabeledQT" else "FeatureLinkerUnlabeled", "OpenMS"),
                   c(settings, "-in", featFiles, "-out", out))
})

# generating XML via package is too slow...http://r.789695.n4.nabble.com/Creating-XML-document-extremely-slow-td4376088.html
# generate by simply writing text to file instead
writeFeatureXML <- function(ft, out)
{
    con <- file(out, "w")

    writeLine <- function(level, txt, ...)
    {
        if (level > 0)
            cat(rep("    ", level), file=con)
        cat(paste0(sprintf(txt, ...), "\n"), file=con)
    }

    writeLine(0, '<?xml version="1.0" encoding="ISO-8859-1"?>')

    writeLine(0, '<featureMap version="1.9" id="fm" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_9.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">')

    writeLine(1, '<featureList count="%d">', nrow(ft))

    for (i in seq_len(nrow(ft)))
    {
        writeLine(2, '<feature id="f_%d">', i)

        writeLine(3, '<position dim="0">%f</position>', ft[i, ret])
        writeLine(3, '<position dim="1">%f</position>', ft[i, mz])
        writeLine(3, '<intensity>%f</intensity>', ft[i, area])
        writeLine(3, '<quality dim="0">0</quality>')
        writeLine(3, '<quality dim="1">0</quality>')
        writeLine(3, '<overallquality>0</overallquality>')
        writeLine(3, '<charge>0</charge>')

        writeLine(2, '</feature>')
    }

    writeLine(1, "</featureList>")

    writeLine(0, "</featureMap>")
    close(con)
}

importConsensusXML <- function(feat, cfile)
{
    cat("Importing consensus XML...")

    sGroup <- analysisInfo(feat)
    fTable <- featureTable(feat)

    doc <- XML::xmlTreeParse(cfile)
    docrt <- XML::xmlRoot(doc)

    ccount <- XML::xmlSize(docrt[["consensusElementList"]])
    groups <- data.table(matrix(0, nrow = nrow(sGroup), ncol = ccount))

    ftindex <- copy(groups)

    gInfoDT <- data.table(rts=numeric(ccount), mzs=numeric(ccount))
    gCount <- 0L

    XML::xmlSApply(docrt[["consensusElementList"]], function(xn)
    {
        gCount <<- gCount + 1L

        set(gInfoDT, gCount, "rts", as.numeric(XML::xmlAttrs(xn[["centroid"]])[["rt"]]))
        set(gInfoDT, gCount, "mzs", as.numeric(XML::xmlAttrs(xn[["centroid"]])[["mz"]]))

        XML::xmlSApply(xn[["groupedElementList"]], function(xsn)
        {
            a <- XML::xmlAttrs(xsn)
            fid <- as.integer(a[["map"]]) + 1L
            ftid <- as.integer(a[["id"]])
            set(ftindex, fid, gCount, ftid)
            set(groups, fid, gCount, fTable[[fid]][["intensity"]][ftid])
        })
    })

    gInfo <- as.data.frame(gInfoDT, stringsAsFactors = FALSE)
    gNames <- sapply(seq_len(nrow(gInfo)), function(grpi) makeFGroupName(grpi, gInfo$rts[grpi], gInfo$mzs[grpi]))
    rownames(gInfo) <- gNames
    setnames(groups, gNames)
    setnames(ftindex, gNames)

    cat("Done!\n")

    return(list(groups = groups, gInfo = gInfo, ftindex = ftindex))
}

# test with xml2 package: nicer code but featureXML code with xml2 is slower, so don't use for now
importConsensusXML2 <- function(feat, cfile)
{
    fTable <- featureTable(feat)

    cat("Importing consensus XML...")

    xml <- xml2::read_xml(cfile)
    consPath <- xml_find_all(xml, "/consensusXML/consensusElementList/consensusElement")

    gInfo <- data.frame(rts <- xml_double(xml_find_all(consPath, "centroid/@rt")),
                        mzs <- xml_double(xml_find_all(consPath, "centroid/@mz")),
                        stringsAsFactors = FALSE)

    elementPath <- xml_find_all(consPath, "groupedElementList")
    consSizes <- xml_length(elementPath)
    consTable <- data.table(maps = xml_integer(xml_find_all(elementPath, "element/@map")),
                            ftInds = xml_integer(xml_find_all(elementPath, "element/@id")),
                            consensus = unlist(sapply(seq_along(consSizes), function(e) rep(e, consSizes[e]))))

    gCount <- length(consSizes)
    groups <- data.table(matrix(0, nrow = nrow(analysisInfo(feat)), ncol = gCount))
    ftindex <- data.table(matrix(0L, nrow = nrow(analysisInfo(feat)), ncol = gCount))

    for (grpi in seq_len(gCount))
    {
        cons <- consTable[consensus == grpi]
        anaInds <- cons[["maps"]] + 1L
        ftInds <- cons[["ftInds"]]
        set(ftindex, anaInds, grpi, ftInds)
        set(groups, anaInds, grpi,
            sapply(seq_len(nrow(cons)), function(i) fTable[[anaInds[i]]][["intensity"]][ftInds[i]]))
    }

    gNames <- sapply(seq_len(nrow(gInfo)), function(grpi) makeFGroupName(grpi, gInfo$rts[grpi], gInfo$mzs[grpi]))
    rownames(gInfo) <- gNames
    setnames(groups, gNames)
    setnames(ftindex, gNames)

    cat("Done!\n")

    return(list(groups = groups, gInfo = gInfo, ftindex = ftindex))






    sGroup <- analysisInfo(feat)
    fTable <- featureTable(feat)

    doc <- XML::xmlTreeParse(cfile)
    docrt <- XML::xmlRoot(doc)

    ccount <- XML::xmlSize(docrt[["consensusElementList"]])
    groups <- data.table(matrix(0, nrow = nrow(sGroup), ncol = ccount))
    # setnames(groups, XML::xmlSApply(docrt[["consensusElementList"]], function(xn) XML::xmlAttrs(xn)[["id"]]))

    ftindex <- copy(groups)

    gInfo <- data.frame(rts=numeric(ccount), mzs=numeric(ccount), row.names = names(groups), stringsAsFactors = F)
    gCount <- 0

    XML::xmlSApply(docrt[["consensusElementList"]], function(xn)
    {
        # cid <- XML::xmlAttrs(xn)[["id"]]

        gCount <<- gCount + 1

        gInfo[gCount, "rts"] <<- as.numeric(XML::xmlAttrs(xn[["centroid"]])[["rt"]])
        gInfo[gCount, "mzs"] <<- as.numeric(XML::xmlAttrs(xn[["centroid"]])[["mz"]])

        XML::xmlSApply(xn[["groupedElementList"]], function(xsn)
        {
            a <- XML::xmlAttrs(xsn)
            fid <- as.numeric(a[["map"]]) + 1
            ftid <- as.numeric(a[["id"]])
            ftindex[[fid, gCount]] <<- ftid
            groups[[fid, gCount]] <<- fTable[[fid]][["intensity"]][ftid]
        })
    })

    gNames <- sapply(seq_len(nrow(gInfo)), function(grpi) makeFGroupName(grpi, gInfo$rts[grpi], gInfo$mzs[grpi]))
    rownames(gInfo) <- gNames
    setnames(groups, gNames)
    setnames(ftindex, gNames)

    cat("Done!\n")

    return(list(groups = groups, gInfo = gInfo, ftindex = ftindex))
}
