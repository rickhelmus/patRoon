#' @include main.R
NULL

# Fix from R DescTools
#' Internal fix for \pkg{RDCOMClient}, ignore.
#' @param ref,className ignore
#' @keywords internal
#' @export
createCOMReference <- function(ref, className)
{
    RDCOMClient::createCOMReference(ref, className)
}

#' Internal fix for \pkg{RDCOMClient}, ignore.
#' @param ... ignore
#' @keywords internal
#' @export
COMStop <- function(...)
{
    RDCOMClient::COMStop(...)
}

#' @details \code{generateAnalysisInfo} is an utility function that
#'   automatically generates an analysis information object. It will collect all
#'   datafiles from given file paths and convert the filenames into valid
#'   analysis names (\emph{i.e.} without extensions such as \file{.d} and
#'   \file{.mzML}). Duplicate analyses, which may appear when datafiles with
#'   different file extension (\file{.d}, \file{.mzXML} and/or \file{.mzML}) are
#'   present, will be automatically removed.
#'
#' @param paths A character vector containing one or more file paths that should
#'   be used for finding the analyses.
#' @param groups,blanks An (optional) character vector containing replicate
#'   groups and references, respectively (will be recycled). If \code{groups} is
#'   an empty character string (\code{""}) the analysis name will be set as
#'   replicate group.
#' @param concs An optional numeric vector containing concentration values for
#'   each analysis. Can be \code{NA} if unknown. If the length of \code{concs}
#'   is less than the number of analyses the remainders will be set to
#'   \code{NA}. Set to \code{NULL} to not include concentration data.
#' @param formats A character vector of analyses file types. Valid values are:
#'   \code{Bruker}, \code{mzXML} and \code{mzML}.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(paths, groups = "", blanks = "", concs = NULL,
                                 formats = MSFileFormats())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(groups, min.len = 1, add = ac)
    checkmate::assertCharacter(blanks, min.len = 1, add = ac)
    checkmate::assertNumeric(concs, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertSubset(formats, MSFileFormats(), empty.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    files <- listMSFiles(paths, formats)
    
    if (length(files) == 0)
    {
        warning(sprintf("No analyses found in %s!", paste(paths, collapse = ", ")))
        return(NULL)
    }
    
    ret <- data.frame(path = dirname(files), analysis = simplifyAnalysisNames(files), stringsAsFactors = FALSE)
    ret <- ret[!duplicated(ret[, c("path", "analysis")]), ]
    
    # set after duplicate removal
    groups <- rep(groups, length.out = nrow(ret))
    ret$group <- ifelse(!nzchar(groups), ret$analysis, groups)
    ret$blank <- blanks
    
    if (!is.null(concs))
    {
        if (length(concs) >= nrow(ret))
            ret$conc <- concs[seq_len(nrow(ret))]
        else
        {
            ret$conc <- NA
            ret$conc[seq_along(concs)] <- concs
        }
    }
    
    return(ret)
}

# nocov start

# UNDONE: delete some day...
#' @details \code{generateAnalysisInfoFromEnviMass} loads analysis information
#'   from an \pkg{enviMass} project. Note: this funtionality has only been
#'   tested with older versions of \pkg{enviMass}.
#'
#' @param path The path of the enviMass project.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfoFromEnviMass <- function(path)
{
    checkmate::assertDirectoryExists(path, access = "r")
    
    enviSInfo <- fread(file.path(path, "dataframes", "measurements"))[Type %in% c("sample", "blank")]
    
    enviSInfo[, DT := as.POSIXct(paste(Date, Time, sep = " "))]
    enviSInfo[, group := tag3]
    
    blanks <- enviSInfo[enviSInfo$Type == "blank", ]
    if (nrow(blanks) > 0)
    {
        blanks[, blank := paste0("blank", match(DT, unique(DT)))] # add unique date+time identifier
        enviSInfo[, blank := sapply(DT, function(dt)
        {
            bls <- blanks[DT <= dt]
            if (nrow(bls) > 0)
                return(bls[which.max(bls$DT), blank])
            else
                return("")
        })]
        
        enviSInfo[Type == "blank", group := blank]
    }
    else
        enviSInfo[, blank := ""]
    
    enviSInfo[group == "FALSE", group := ""]
    
    ret <- data.frame(path = file.path(path, "files"), analysis = enviSInfo$ID, group = enviSInfo$group,
                      blank = enviSInfo$blank, stringsAsFactors = FALSE)
    
    return(ret)
}
# nocov end


#' Temporarily changes package options
#'
#' This function is inspired by
#' \code{\link[withr:with_options]{withr::with_options}}: it can be used to
#' execute some code where package options are temporarily changed. This
#' function uses a shortened syntax, especially when changing options for
#' \code{patRoon}.
#'
#' @param \dots Named arguments with options to change.
#' @param code The code to be executed.
#' @param prefix A \code{character} that will be used to prefix given option
#'   names.
#'
#' @examples \dontrun{
#' # Set max parallel processes to five while performing formula calculations
#' withOpt(MP.maxProcs = 5, {
#'     formulas <- generateFormulas(fGroups, "genform", ...)
#' })
#' }
#'
#' @export
withOpt <- withr::with_(function(..., prefix = "patRoon.")
{
    checkmate::assertString(prefix, null.ok = TRUE)
    opts <- list(...)
    checkmate::assertNamed(opts, "unique", .var.name = "...")
    
    if (!is.null(prefix))
        names(opts) <- paste0(prefix, names(opts))
    
    options(opts)
}, function(old) options(old), new = FALSE)

#' Prints all the package options of \code{patRoon} and their currently set values.
#' @export
printPackageOpts <- function() dumpPkgOpts(function(s) cat(s, "\n", sep = ""))

#' Verifies if all dependencies are installed properly and instructs the user if
#' this is not the case.
#' @export
verifyDependencies <- function()
{
    # UNDONE: for now just check one command-line tool of a software package
    # UNDONE: skip GenForm for now? Should be present as embedded binary.
    
    OK <- TRUE
    check <- function(name, path, opt, isDir = FALSE)
    {
        pleaseSet <- sprintf("Please set the '%s' option.", opt)
        printf("Checking %s... ", name)
        if (is.null(path) || !nzchar(path))
        {
            cat("not found or configured.", pleaseSet, "\n")
            OK <<- FALSE
        }
        else if (isDir)
        {
            if (!dir.exists(path))
            {
                cat("configured directory path does not exist!", pleaseSet, "\n")
                OK <<- FALSE
            }
            else
                printf("found! (directory '%s')\n", path)
        }
        else
        {
            # NOTE: dirname point to current path if getCommandWithOptPath() found it in PATH
            dn <- dirname(path)
            # if ((nzchar(dn) && dn != "." && !file.exists(path)) || !nzchar(Sys.which(path)))
            if (nzchar(dn) && dn != "." && !file.exists(path))
            {
                cat("configured path does not exist!", pleaseSet, "\n")
                OK <<- FALSE
            }
            else if (nzchar(dn) && dn != ".")
                printf("found! (in '%s')\n", dn)
            else
                cat("found!\n")
        }
    }
    
    check("ProteoWizard", findPWizPath(), "patRoon.path.pwiz", isDir = TRUE)
    check("OpenMS", getCommandWithOptPath("FeatureFinderMetabo", "OpenMS", verify = FALSE), "patRoon.path.OpenMS")
    check("pngquant", getCommandWithOptPath("pngquant", "pngquant", verify = FALSE), "patRoon.path.pngquant")
    check("SIRIUS", getCommandWithOptPath(getSiriusBin(), "SIRIUS", verify = FALSE), "patRoon.path.SIRIUS")
    check("MetFrag CL", getOption("patRoon.path.MetFragCL"), "patRoon.path.MetFragCL")
    check("MetFrag CompTox Database", getOption("patRoon.path.MetFragCompTox"), "patRoon.path.MetFragCompTox")
    check("MetFrag PubChemLite Database", getOption("patRoon.path.MetFragPubChemLite"), "patRoon.path.MetFragPubChemLite")
    check("OpenBabel", getCommandWithOptPath("obabel", "obabel", verify = FALSE), "patRoon.path.obabel")
    check("BioTransformer", getOption("patRoon.path.BioTransformer"), "patRoon.path.BioTransformer")
    
    if (!OK)
        cat("\nSome dependencies were not found. Please make sure that their file locations are configured properly.",
            "For instance, run the following to set the location of MetFragCL:",
            sprintf("options(patRoon.path.MetFragCL = \"C:/MetFrag2.4.5-CL.jar\")"),
            "\nPlease see ?patRoon for more information on how to configure patRoon options.",
            sep = "\n")
    
    invisible(NULL)
}

#' Returns chromatographic peak quality and score names for features and/or feature groups.
#'
#' @param feat If \code{TRUE} then names specific to features are returned.
#' @param group If \code{TRUE} then names specific to groups are returned.
#' @param scores If \code{TRUE} the score names are returned, otherwise the quality names.
#' @param totScore If \code{TRUE} (and \code{scores=TRUE}) then the name of the total score is included.
#'
#' @export
featureQualityNames <- function(feat = TRUE, group = TRUE, scores = FALSE, totScore = TRUE)
{
    ret <- character()
    if (feat)
        ret <- names(featureQualities())
    if (group)
        ret <- c(ret, names(featureGroupQualities()))
    if (scores)
    {
        ret <- paste0(ret, "Score")
        if (totScore)
            ret <- c(ret, "totalScore")
    }
    return(ret)
}
