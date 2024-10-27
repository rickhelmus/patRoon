#' @include main.R
NULL

MSFileExtensions <- function()
{
    list(thermo = "raw",
         bruker = "d",
         bruker_ims = "d",
         agilent = "d",
         agilent_ims = "d",
         ab = "wiff",
         waters = "raw",
         mzXML = "mzXML",
         mzML = "mzML")
}

getMSFileTypes <- function() c("centroid", "profile", "raw", "ims")

getMSReadBackends <- function() c("opentims", "streamcraft", "mstoolkit", "mzr")

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

getAllMSFilesFromAnaInfo <- function(anaInfo, types, formats)
{
    if (!is.list(formats))
    {
        # formats can be a vector if only one type was specified (common scenario)
        stopifnot(length(types) == 1)
        formats <- list(formats)
    }
    
    ret <- Map(types, formats, f = function(type, forms)
    {
        lapply(forms, function(f)
        {
            aip <- getPathsFromAnaInfo(anaInfo, type)
            if (is.null(aip))
                next
            msFilePaths <- listMSFiles(aip, f)
            msFilesNoExt <- tools::file_path_sans_ext(basename(msFilePaths))
            found <- anaInfo$analysis %in% msFilesNoExt
            return(msFilePaths[match(anaInfo$analysis[found], msFilesNoExt[found])])
        })
    })
    
    return(unique(unlist(ret)))
}

getMSFileHashesFromAvailBackend <- function(anaInfo, types = getMSFileTypes(), formats = names(MSFileExtensions()))
{
    backends <- getOption("patRoon.MSBackends", character())
    
    for (bn in backends)
    {
        if (!backendAvailable(bn))
            next
        filePaths <- maybeGetMSFiles(bn, anaInfo, types, formats)
        if (!is.null(filePaths))
        {
            if (bn == "opentims")
                filePaths <- file.path(filePaths, "analysis.tdf_bin")
            return(setNames(sapply(filePaths, getMSDataFileHash), anaInfo$analysis))
        }
    }
    
    stop("Failed to load a correct MS read backend. Please ensure patRoon.MSBackends is configured properly. See ?patRoon",
         call. = FALSE)
}

# shortcut for common case
getCentroidedMSFilesFromAnaInfo <- function(anaInfo, formats = c("mzML", "mzXML"), mustExist = TRUE)
{
    getMSFilesFromAnaInfo(anaInfo, "centroid", formats, mustExist)
}

doGetEICs <- function(anaInfo, EICInfoList, minIntensityIMS = 0, compress = TRUE, withBP = FALSE, cacheDB = NULL)
{
    # HACK: for now we _don't_ cache EICs if compres==FALSE: the resulting data is very large and takes a long time to
    # be stored. Hence, the caller should cache the final results.
    doCache <- compress

    anaHashes <- NULL
    if (doCache)
    {
        if (is.null(cacheDB))
            cacheDB <- openCacheDBScope()
        anaHashes <- getMSFileHashesFromAvailBackend(anaInfo)
    }
    
    allEICs <- applyMSData(anaInfo, EICInfoList, func = function(ana, path, backend, EICInfo)
    {
        EICInfo <- copy(EICInfo)
        for (col in c("mobmin", "mobmax"))
        {
            if (is.null(EICInfo[[col]]))
                set(EICInfo, j = col, value = 0)
            else
                setnafill(EICInfo, fill = 0, cols = col)
        }
        
        EICs <- vector("list", nrow(EICInfo))
        
        if (doCache)
        {
            # NOTE: subset columns here, so any additional columns from e.g. feature tables are not considered
            hashes <- EICInfo[, makeHash(anaHashes[[ana]], compress, withBP, .SD), by = seq_len(nrow(EICInfo)),
                              .SDcols = c("retmin", "retmax", "mzmin", "mzmax", "mobmin", "mobmax")][[2]]
            
            cachedData <- loadCacheData(category = "EICs", hashes, dbArg = cacheDB, simplify = FALSE)
            if (!is.null(cachedData) && length(cachedData) == nrow(EICInfo))
            {
                doProgress()
                return(unname(cachedData)) # everything is in the cache
            }
            
            cachedInds <- if (!is.null(cachedData)) match(names(cachedData), hashes) else integer()
            isCached <- if (!is.null(cachedData)) hashes %chin% names(cachedData) else rep(FALSE, nrow(EICInfo))
            # NOTE: cachedData is 'subset' below to make sure any duplicate hashes are properly assigned
            EICs[isCached] <- cachedData[match(hashes, names(cachedData), nomatch = 0)]
        }
        else
            isCached <- rep(FALSE, length(EICs))
        
        ToDo <- EICInfo[isCached == FALSE]
        openMSReadBackend(backend, path)
        
        # NOTE: getEICList() return lists, which are converted to data.frames and is a lot faster than returning
        # data.frames directly.
        newEICs <- getEICList(backend, ToDo$mzmin, ToDo$mzmax, ToDo$retmin, ToDo$retmax, ToDo$mobmin, ToDo$mobmax,
                              minIntensityIMS, compress, withBP)
        newEICs <- lapply(newEICs, setDF)
        EICs[!isCached] <- newEICs

        if (doCache)
            saveCacheDataList("EICs", EICs[!isCached], hashes[!isCached], cacheDB)
        
        doProgress()
        return(EICs)
    })
    
    return(allEICs)
}
