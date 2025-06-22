#' @include main.R
NULL

doInitBrukerLib <- function()
{
    libp <- getOption("patRoon.path.BrukerTIMS", "")
    return(nzchar(libp) && initBrukerLibrary(libp))
}

MSFileExtensions <- function()
{
    list(thermo = "raw",
         bruker = "d",
         bruker_ims = "d",
         agilent = "d",
         agilent_ims = "d",
         ab = "wiff",
         waters = "raw",
         mzML = "mzML",
         mzXML = "mzXML")
}

getMSReadBackends <- function() c("opentims", "mzr", "mstoolkit", "streamcraft")

getMSDataFileHash <- function(path)
{
    # NOTE: for hashing limit length as this function may be called frequently. The path is also included to make it
    # hopefully sufficiently unique.
    return(makeHash(path, makeFileHash(path, length = 8192)))
}

verifyFileForFormat <- function(path, format)
{
    isDir <- file.info(path, extra_cols = FALSE)$isdir
    if (isDir)
    {
        if (format == "agilent" && file.exists(file.path(path, "AcqData")))
            return(TRUE)
        if (format == "agilent_ims" && file.exists(file.path(path, "AcqData")))
            return(TRUE) # UNDONE: more checks?
        if (format == "bruker" && file.exists(file.path(path, "analysis.baf")))
            return(TRUE)
        if (format == "bruker_ims" && file.exists(file.path(path, "analysis.tdf")))
            return(TRUE)
        if (format == "waters")
            return(TRUE) # UNDONE: also more checks?
    }
    return(!isDir && !format %in% c("agilent", "agilent_ims", "bruker", "bruker_ims", "waters"))
}

filterMSFileDirs <- function(files, from)
{
    if (length(files) == 0)
        return(files)
    
    allFromExts <- MSFileExtensions()[from]
    keep <- sapply(files, function(file)
    {
        fExt <- tolower(tools::file_ext(file))
        return(any(mapply(allFromExts, names(allFromExts), FUN = function(e, f)
        {
            tolower(e) == fExt && verifyFileForFormat(file, f)
        })))
    })
    
    return(files[keep])
}

listMSFiles <- function(dirs, from)
{
    dirs <- normalizePath(unique(dirs), mustWork = FALSE, winslash = "/")
    
    allExts <- MSFileExtensions()
    allExts <- unique(unlist(allExts[from]))
    
    files <- list.files(dirs, full.names = TRUE, pattern = paste0(paste0(".+\\.", allExts, collapse = "|"), "$"),
                        ignore.case = TRUE)
    
    # try to get to back to the original directory order
    ord <- match(dirname(files), dirs)
    files <- files[order(ord)]
    
    return(filterMSFileDirs(files, from))
}

getPathsFromAnaInfo <- function(anaInfo, type) return(anaInfo[[paste0("path_", type)]])
getAnaInfoPathCols <- function(anaInfo) intersect(paste0("path_", getMSFileTypes()), names(anaInfo))

getMSFilesFromAnaInfo <- function(anaInfo, types, formats, mustExist = TRUE)
{
    if (!is.list(formats))
    {
        # formats can be a vector if only one type was specified (common scenario)
        stopifnot(length(types) == 1)
        formats <- list(formats)
    }
    
    missing <- list()
    for (i in seq_along(types))
    {
        # go through allowed formats one by one: otherwise the resulting file types could be mixed, which might be unwanted(?)
        for (form in formats[[i]])
        {
            aip <- getPathsFromAnaInfo(anaInfo, types[i])
            if (is.null(aip))
            {
                missing[[types[i]]][[form]] <- seq_along(anaInfo$analysis)
                next
            }
            msFilePaths <- listMSFiles(aip, form)
            msFilesNoExt <- tools::file_path_sans_ext(basename(msFilePaths))
            found <- anaInfo$analysis %in% msFilesNoExt
            if (all(found)) # success!
                return(msFilePaths[match(anaInfo$analysis, msFilesNoExt)])
            missing[[types[i]]][[form]] <- which(!found)
        }
    }

    if (!mustExist)
        return(NULL)

    missingAll <- Reduce(intersect, lapply(missing, Reduce, f = intersect))
    stop(sprintf("The following analyses are not found with any valid type (%s) and/or data format (%s): %s",
                 paste0(types, collapse = ", "), paste0(unique(unlist(formats)), collapse = ", "),
                 paste0(anaInfo$analysis[missingAll], collapse = ", ")),
         call. = FALSE)
}

getMSFileHashesFromAvailBackend <- function(anaInfo, types = getMSFileTypes(), formats = names(MSFileExtensions()),
                                            needIMS = FALSE)
{
    backends <- getOption("patRoon.MS.backends", character())
    
    for (bn in backends)
    {
        if (!backendAvailable(bn))
            next
        filePaths <- maybeGetMSFiles(bn, anaInfo, types, formats, needIMS)
        if (!is.null(filePaths))
        {
            if (bn == "opentims")
                filePaths <- file.path(filePaths, "analysis.tdf_bin")
            return(setNames(sapply(filePaths, getMSDataFileHash), anaInfo$analysis))
        }
    }
    
    stop("Failed to load a correct MS read backend. Please ensure patRoon.MS.backends is configured properly. See ?patRoon",
         call. = FALSE)
}

# shortcut for common case
getCentroidedMSFilesFromAnaInfo <- function(anaInfo, formats = c("mzML", "mzXML"), mustExist = TRUE)
{
    msf <- getMSFilesFromAnaInfo(anaInfo, "centroid", formats, mustExist)
    if (is.null(msf))
        return(NULL)
    return(setNames(msf, anaInfo$analysis))
}

doGetEICsForAna <- function(...)
{
    # NOTE: getEICList() return lists, which are converted to data.frames and is a lot faster than returning
    # data.frames directly.
    EICs <- getEICList(...)
    EICs[] <- lapply(EICs, setDF)
    return(EICs)
}

doGetEICs <- function(anaInfo, EICInfoList, mzExpIMSWindow = 0, minIntensityIMS = 0, mode = "simple",
                      minEICIntensity = 0, minEICAdjTime = 0, minEICAdjPoints = 0, minEICAdjIntensity = 0,
                      pad = FALSE, doCache = TRUE, cacheDB = NULL)
{
    if (length(EICInfoList) == 0)
        return(list())
    
    anaInfo <- anaInfo[analysis %in% names(EICInfoList)]
    
    needIMS <- !is.null(EICInfoList[[1]][["mobmin"]])
    
    anaHashes <- baseHash <- NULL
    if (doCache)
    {
        if (is.null(cacheDB))
            cacheDB <- openCacheDBScope()
        anaHashes <- getMSFileHashesFromAvailBackend(anaInfo, needIMS = needIMS)
        baseHash <- makeHash(mzExpIMSWindow, minIntensityIMS, mode, minEICIntensity, minEICAdjTime, minEICAdjPoints,
                             minEICAdjIntensity)
    }
    
    allEICs <- applyMSData(anaInfo, EICInfoList, showProgress = TRUE, needIMS = needIMS,
                           func = function(ana, path, backend, EICInfo)
    {
        EICInfo <- copy(EICInfo)
        for (col in c("mobmin", "mobmax"))
        {
            if (is.null(EICInfo[[col]]))
                set(EICInfo, j = col, value = 0)
            else
                setnafill(EICInfo, fill = 0, cols = col)
        }
        
        EICs <- isCached <- anaHash <- anaEICHash <- NULL
        if (doCache)
        {
            anaEICHash <- makeHash(anaHashes[[ana]], baseHash)
            # NOTE: subset columns here, so any additional columns from e.g. feature tables are not considered
            hashes <- EICInfo[, makeHash(anaEICHash, .SD), by = seq_len(nrow(EICInfo)),
                              .SDcols = c("retmin", "retmax", "mzmin", "mzmax", "mobmin", "mobmax")][[2]]
            
            EICs <- unname(loadCacheData(category = "EICs", hashes, dbArg = cacheDB, simplify = FALSE, fixDTs = FALSE))
            ax <- loadCacheData(category = "EICAllTimes", anaEICHash, dbArg = cacheDB)
            isCached <- !is.null(ax) & !sapply(EICs, is.null)
            if (!is.null(ax))
                attr(EICs, "allXValues") <- ax
        }
        else
        {
            EICs <- vector("list", nrow(EICInfo))
            isCached <- rep(FALSE, nrow(EICInfo))
        }

        if (any(!isCached))
        {

            ToDo <- EICInfo[isCached == FALSE]
            openMSReadBackend(backend, path)
            
            newEICs <- doGetEICsForAna(backend, ToDo$mzmin, ToDo$mzmax, ToDo$retmin, ToDo$retmax, ToDo$mobmin,
                                       ToDo$mobmax, mzExpIMSWindow, minIntensityIMS, mode, minEICIntensity,
                                       minEICAdjTime, minEICAdjPoints, minEICAdjIntensity)
            EICs[!isCached] <- newEICs
            attr(EICs, "allXValues") <- attr(newEICs, "allXValues")
            
            if (doCache)
            {
                saveCacheDataList("EICs", EICs[!isCached], hashes[!isCached], cacheDB)
                saveCacheData("EICAllTimes", attr(EICs, "allXValues"), anaEICHash, cacheDB)
            }
        }
        if (pad)
        {
            EICs[] <- Map(EICs, EICInfo$retmin, EICInfo$retmax, f = function(eic, rmin, rmax)
            {
                px <- padEIX(attr(EICs, "allXValues"), rmin, rmax, eic$time, eic$intensity)
                data.frame(time = px[[1]], intensity = px[[2]])
            })
        }
        
        doProgress()
        return(EICs)
    })
    
    return(allEICs)
}

doGetEIMs <- function(anaInfo, EIMInfoList, IMSWindow, clusterMethod, minIntensity, mzExpIMSWindow = 0,
                      compress = TRUE, cacheDB = NULL)
{
    if (length(EIMInfoList) == 0)
        return(list())
    
    anaInfo <- anaInfo[analysis %in% names(EIMInfoList)]
    
    anaHashes <- NULL
    if (is.null(cacheDB))
        cacheDB <- openCacheDBScope()
    anaHashes <- getMSFileHashesFromAvailBackend(anaInfo, needIMS = TRUE)
    
    allEIMs <- applyMSData(anaInfo, EIMInfoList, needIMS = TRUE, func = function(ana, path, backend, EIMInfo)
    {
        EIMInfo <- copy(EIMInfo)
        for (col in c("mobmin", "mobmax"))
        {
            if (is.null(EIMInfo[[col]]))
                set(EIMInfo, j = col, value = 0)
            else
                setnafill(EIMInfo, fill = 0, cols = col)
        }
        
        # NOTE: subset columns here, so any additional columns from e.g. feature tables are not considered
        hashes <- EIMInfo[, makeHash(anaHashes[[ana]], IMSWindow, clusterMethod, minIntensity, mzExpIMSWindow,
                                     compress, .SD),
                          by = seq_len(nrow(EIMInfo)), .SDcols = c("retmin", "retmax", "mzmin", "mzmax", "mobmin",
                                                                   "mobmax")][[2]]
        
        EIMs <- loadCacheData(category = "EIMs", hashes, dbArg = cacheDB, simplify = FALSE)
        isCached <- !sapply(EIMs, is.null)
        if (all(isCached))
        {
            # everything is in the cache
            doProgress()
            return(unname(EIMs))
        }

        ToDo <- EIMInfo[isCached == FALSE]
        
        openMSReadBackend(backend, path)
        
        # NOTE: getEIMList() return lists, which are converted to data.frames and is a lot faster than returning
        # data.frames directly.
        newEIMs <- getEIMList(backend, ToDo$mzmin, ToDo$mzmax, ToDo$retmin, ToDo$retmax, ToDo$mobmin, ToDo$mobmax,
                              clusterMethod, IMSWindow, minIntensity, mzExpIMSWindow, compress)
        newEIMs <- lapply(newEIMs, setDF)
        EIMs[!isCached] <- newEIMs
        
        saveCacheDataList("EIMs", EIMs[!isCached], hashes[!isCached], cacheDB)
        
        doProgress()
        return(EIMs)
    })
    
    return(allEIMs)
}

prepareAgilentIMSCalib <- function(calibrant, massGas)
{
    if (is.list(calibrant))
    {
        if (is.null(calibrant[["massGas"]]))
            calibrant$massGas <- massGas # UNDONE
    }
    else
    {
        if (checkmate::testDirectory(calibrant))
        {
            calibrant <- file.path(calibrant, "AcqData", "OverrideImsCal.xml")
            if (!file.exists(calibrant))
                stop(sprintf("Cannot load '%s': file does not exists", calibrant), call. = FALSE)
        }
        calibrant <- loadAgilentIMSCalibration(calibrant)
    }
    return(calibrant)
}
