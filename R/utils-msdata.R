#' @include main.R
NULL

MSFileExtensions <- function()
{
    list(thermo = "raw",
         bruker = "d",
         bruker_tims = "d",
         agilent = "d",
         ab = "wiff",
         waters = "raw",
         mzXML = "mzXML",
         mzML = "mzML")
}

getMSFileTypes <- function() c("centroid", "profile", "raw", "ims")

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
        if (format == "bruker" && file.exists(file.path(path, "analysis.baf")))
            return(TRUE)
        if (format == "bruker_tims" && file.exists(file.path(path, "analysis.tdf")))
            return(TRUE)
        if (format == "waters")
            return(TRUE) # UNDONE: also more checks?
    }
    return(!isDir && !format %in% c("agilent", "bruker", "bruker_tims", "waters"))
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

getMSFilesFromAnaInfo <- function(anaInfo, type, formats, mustExist = TRUE)
{
    missing <- list()
    
    # go through allowed formats one by one: otherwise the resulting file types could be mixed, which might be unwanted(?)
    for (form in formats)
    {
        aip <- getPathsFromAnaInfo(anaInfo, type)
        if (is.null(aip))
            next
        msFilePaths <- listMSFiles(aip, form)
        msFilesNoExt <- tools::file_path_sans_ext(basename(msFilePaths))
        found <- anaInfo$analysis %in% msFilesNoExt
        if (all(found)) # success!
            return(msFilePaths[match(anaInfo$analysis, msFilesNoExt)])
        missing[[form]] <- which(!found)
    }

    if (!mustExist)
        return(NULL)
    
    missingAll <- Reduce(intersect, missing)
    stop(sprintf("The following analyses with type %s are not found with a correct data format (valid: %s): %s",
                 type, paste0(formats, collapse = ", "), paste0(anaInfo$analysis[missingAll], collapse = ", ")),
         call. = FALSE)
}

# shortcut for common case
getCentroidedMSFilesFromAnaInfo <- function(anaInfo, formats = c("mzML", "mzXML"), mustExist = TRUE)
{
    getMSFilesFromAnaInfo(anaInfo, "centroid", formats, mustExist)
}

doGetEICs <- function(anaInfo, EICInfoList, cacheDB = NULL)
{
    # UNDONE: handle mobilities
    
    if (is.null(cacheDB))
        cacheDB <- openCacheDBScope()
    
    EICs <- applyMSData(anaInfo, EICInfoList, func = function(ana, path, backend, EICInfo)
    {
        anaHash <- makeHash(getMSDataFileHash(path))
        
        # NOTE: subset columns here, so any additional columns from e.g. feature tables are not considered
        hashes <- EICInfo[, makeHash(anaHash, .SD), by = seq_len(nrow(EICInfo)),
                          .SDcols = c("retmin", "retmax", "mzmin", "mzmax")][[2]]
        
        cachedData <- loadCacheData(category = "EICs", hashes, dbArg = cacheDB, simplify = FALSE)
        if (!is.null(cachedData) && length(cachedData) == nrow(EICInfo))
        {
            doProgress()
            return(unname(cachedData)) # everything is in the cache
        }
        
        EICs <- vector("list", nrow(EICInfo))
        cachedInds <- if (!is.null(cachedData)) match(names(cachedData), hashes) else integer()
        isCached <- if (!is.null(cachedData)) hashes %chin% names(cachedData) else rep(FALSE, nrow(EICInfo))
        # NOTE: cachedData is 'subset' below to make sure any duplicate hashes are properly assigned
        EICs[isCached] <- cachedData[match(hashes, names(cachedData), nomatch = 0)]
        
        ToDo <- EICInfo[isCached == FALSE]
        openMSReadBackend(backend, path)
        EICs[!isCached] <- getEICList(backend, ToDo$mzmin, ToDo$mzmax, ToDo$retmin, ToDo$retmax, rep(0, nrow(ToDo)),
                                      rep(0, nrow(ToDo)), TRUE)

        saveCacheDataList("EICs", EICs[!isCached], hashes[!isCached], cacheDB)
        
        doProgress()
        return(EICs)
    })
    
    return(EICs)
}
