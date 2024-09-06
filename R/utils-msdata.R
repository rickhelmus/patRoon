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

getMSFilesFromAnaInfo <- function(anaInfo, type, formats, mustExist = FALSE)
{
    msFilePaths <- listMSFiles(getPathsFromAnaInfo(anaInfo, type), formats)
    msFilesNoExt <- tools::file_path_sans_ext(basename(msFilePaths))
    found <- anaInfo$analysis %in% msFilesNoExt
    
    if (mustExist && any(!found))
    {
        stop(sprintf("The following analyses with type %s are not found with a correct data format (valid: %s): %s",
                     type, paste0(formats, collapse = ", "), paste0(files[!found], collapse = ", ")), call. = FALSE)
    }
    
    return(msFilePaths[match(anaInfo$analysis, msFilesNoExt, nomatch = 0)])
}
