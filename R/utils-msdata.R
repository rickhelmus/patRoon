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
