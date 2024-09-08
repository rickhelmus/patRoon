#' @include main.R
NULL

MSFileExtensions <- function()
{
    list(thermo = "raw",
         bruker = c("d", "yep", "baf", "fid"),
         agilent = "d",
         ab = "wiff",
         waters = "raw",
         mzXML = "mzXML",
         mzML = "mzML")
}

getMSFileTypes <- function() c("centroid", "profile", "raw", "ims")

MSFileFormatIsDir <- function(format, ext)
{
    # UNDONE: is agilent .d also a directory?
    return((format == "bruker" && ext == "d") || (format == "waters" && ext == "raw"))
}

filterMSFileDirs <- function(files, from)
{
    if (length(files) == 0)
        return(files)
    
    allFromExts <- MSFileExtensions()[from]
    keep <- sapply(files, function(file)
    {
        fExt <- tools::file_ext(file)
        
        fromExts <- pruneList(lapply(allFromExts, function(f) f[tolower(f) %in% tolower(fExt)]), checkEmptyElements = TRUE)
        if (length(fromExts) == 0)
            return(FALSE)
        
        fromCheck <- names(fromExts)
        shouldBeDir <- mapply(fromCheck, fromExts, SIMPLIFY = TRUE,
                              FUN = function(format, exts) sapply(exts, MSFileFormatIsDir, format = format))
        
        if (!allSame(shouldBeDir))
            return(TRUE) # can be either
        
        isDir <- file.info(file, extra_cols = FALSE)$isdir
        if (all(shouldBeDir))
            return(isDir)
        return(!isDir)
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

getMSFilesFromAnaInfo <- function(anaInfo, type, formats)
{
    missing <- list()
    
    # go through allowed formats one by one: otherwise the resulting file types could be mixed, which might be unwanted(?)
    for (form in formats)
    {
        msFilePaths <- listMSFiles(getPathsFromAnaInfo(anaInfo, type), form)
        msFilesNoExt <- tools::file_path_sans_ext(basename(msFilePaths))
        found <- anaInfo$analysis %in% msFilesNoExt
        if (all(found)) # success!
            return(msFilePaths[match(anaInfo$analysis, msFilesNoExt)])
        missing[[form]] <- which(!found)
    }

    missingAll <- Reduce(intersect, missing)
    stop(sprintf("The following analyses with type %s are not found with a correct data format (valid: %s): %s",
                 type, paste0(formats, collapse = ", "), paste0(anaInfo$analysis[missingAll], collapse = ", ")),
         call. = FALSE)
}

# shortcut for common case
getCentroidedMSFilesFromAnaInfo <- function(anaInfo, formats = c("mzML", "mzXML")) getMSFilesFromAnaInfo(anaInfo, "centroid", formats)
