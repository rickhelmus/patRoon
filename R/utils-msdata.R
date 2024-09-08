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
