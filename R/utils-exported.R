#' @include main.R
NULL

# nocov start

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

# nocov end

#' @details \code{generateAnalysisInfo} is an utility function that automatically generates a \code{data.frame} with
#'   analysis information. It scans the directories specified from the \code{paths} argument for analyses, and uses this
#'   to automatically fill in the \code{analysis} and \code{path} columns. Furthermore, this function also correctly
#'   handles analyses which are available in multiple formats.
#'
#' @param paths A character vector containing one or more file paths that should be used for finding the analyses.
#' @param groups,blanks An (optional) character vector containing replicate groups and blanks, respectively (will be
#'   recycled). If \code{groups} is an empty character string (\code{""}) the analysis name will be set as replicate
#'   group.
#' @param concs An optional numeric vector containing concentration values for each analysis. Can be \code{NA} if
#'   unknown. If the length of \code{concs} is less than the number of analyses the remainders will be set to \code{NA}.
#'   Set to \code{NULL} to not include concentration data.
#' @param norm_concs An optional numeric vector containing concentrations used for \emph{feature normalization} (see the
#'   \verb{Feature intensity normalization} section in the \link[=featureGroups-class]{featureGroups documentation}).
#'   \code{NA} values are allowed for analyses that should not be normalized (\emph{e.g.} because no IS is present). If
#'   the length of \code{norm_concs} is less than the number of analyses the remainders will be set to \code{NA}. Set to
#'   \code{NULL} to not include normalization concentration data.
#' @param formats A character vector of analyses file types to consider. Analyses not present in these formats will be
#'   ignored. For valid values see \code{\link{MSFileFormats}}.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(paths, groups = "", blanks = "", concs = NULL, norm_concs = NULL,
                                 formats = MSFileFormats())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(groups, min.len = 1, add = ac)
    checkmate::assertCharacter(blanks, min.len = 1, add = ac)
    aapply(checkmate::assertNumeric, . ~ concs + norm_concs, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
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
    
    getConcs <- function(x)
    {
        if (length(x) >= nrow(ret))
            return(x[seq_len(nrow(ret))])
        ret <- rep(NA_real_, nrow(ret))
        ret[seq_along(x)] <- x
        return(ret)
    }
    
    if (!is.null(concs))
        ret$conc <- getConcs(concs)
    if (!is.null(norm_concs))
        ret$norm_conc <- getConcs(norm_concs)
    
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
    check("OpenMS", getExtDepPath("openms", "FeatureFinderMetabo", verify = FALSE), "patRoon.path.OpenMS")
    check("pngquant", getExtDepPath("pngquant", verify = FALSE), "patRoon.path.pngquant")
    check("SIRIUS", getExtDepPath("sirius", verify = FALSE), "patRoon.path.SIRIUS")
    check("MetFrag CL", getExtDepPath("metfragcl"), "patRoon.path.MetFragCL")
    check("MetFrag CompTox Database", getExtDepPath("metfragct"), "patRoon.path.MetFragCompTox")
    check("MetFrag PubChemLite Database", getExtDepPath("metfragpcl"), "patRoon.path.MetFragPubChemLite")
    check("OpenBabel", getExtDepPath("openbabel"), "patRoon.path.obabel")
    check("BioTransformer", getExtDepPath("biotransformer"), "patRoon.path.BioTransformer")
    
    if (!OK)
        cat("\nSome dependencies were not found. Please make sure that their file locations are configured properly.",
            "For instance, run the following to set the location of MetFragCL:",
            sprintf("options(patRoon.path.MetFragCL = \"C:/MetFragCommandLine-2.4.8.jar\")"),
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
        ret <- c("ApexBoundaryRatio", "FWHM2Base", "Jaggedness", "Modality", "Symmetry", "GaussianSimilarity",
                 "Sharpness", "TPASR", "ZigZag")
    if (group)
        ret <- c(ret, "ElutionShift", "RetentionTimeCorrelation")
    if (scores)
    {
        ret <- paste0(ret, "Score")
        if (totScore)
            ret <- c(ret, "totalScore")
    }
    return(ret)
}

#' Fold change calculation
#'
#' @details Fold change calculation can be used to easily identify significant changes between replicate groups. The
#'   calculation process is configured through a paramater list, which can be constructed with the \code{getFCParams}
#'   function. The parameter list has the following entries: \itemize{
#'
#'   \item \code{rGroups} the name of the two replicate groups to compare (taken from the \code{rGroups} argument to
#'   \code{getFCParams}).
#'
#'   \item \code{thresholdFC}: the threshold log FC for a feature group to be classified as increasing/decreasing.
#'
#'   \item \code{thresholdPV}: the threshold log P for a feature group to be significantly different.
#'
#'   \item \code{zeroMethod},\code{zeroValue}: how to handle zero values when calculating the FC: \code{add} adds an
#'   offset to zero values, \code{"fixed"} sets zero values to a fixed number and \code{"omit"} removes zero data. The
#'   number that is added/set by the former two options is defined by \code{zeroValue}.
#'
#'   \item \code{PVTestFunc}: a function that is used to calculate P values (usually using \code{\link{t.test}}).
#'
#'   \item \code{PVAdjFunc}: a function that is used to adjust P values (usually using \code{\link{p.adjust}})
#'
#'   }
#'   
#' @param rGroups A \code{character} vector with the names of the two replicate groups to be compared.
#' @param \dots Optional named arguments that override defaults.
#' 
#' @author The code to calculate and plot Fold change data was created by Bas van de Velde.
#'
#' @seealso \code{\link{featureGroups-class}} and \link{feature-plotting}
#'
#' @export
getFCParams <- function(rGroups, ...)
{
    checkmate::assertCharacter(rGroups, min.chars = 1, len = 2, any.missing = FALSE)
    
    def <- list(rGroups = rGroups,
                thresholdFC = 0.25,
                thresholdPV = 0.05,
                zeroMethod = "add",
                zeroValue = 0.01,
                PVTestFunc = function(x, y) t.test(x, y, paired = TRUE)$p.value,
                PVAdjFunc = function(pv) p.adjust(pv, "BH"))
    return(modifyList(def, list(...)))
}

#' Obtains extracted ion chromatograms (EICs)
#'
#' This function generates one or more EIC(s) for given retention time and \emph{m/z} ranges.
#'
#' @param file The file path to the sample analysis data file (\file{.mzXML} or \file{.mzML}).
#' @param ranges A \code{data.frame} with \code{numeric} columns \code{"retmin"}, \code{"retmin"}, \code{"mzmin"},
#'   \code{"mzmax"} with the lower/upper ranges of the retention time and \emph{m/z}.
#'
#' @return A \code{list} with EIC data for each of the rows in \code{ranges}. 
#'
#' @export
getEICs <- function(file, ranges)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    checkmate::assertDataFrame(ranges, types = "numeric", any.missing = FALSE, add = ac)
    assertHasNames(ranges, c("mzmin", "mzmax", "retmin", "retmax"))
    checkmate::reportAssertions(ac)
    
    return(doGetEICs(file, as.data.table(ranges)))
}

#' Extracted Ion Chromatogram parameters
#'
#' Parameters for creation of extracted ion chromatograms.
#'
#' To configure the creation of extracted ion chromatograms (EICs) several parameters exist:
#'
#' \itemize{
#'
#' \item \code{rtWindow} Retention time (in seconds) that will be subtracted/added to respectively the minimum and
#' maximum retention time of the feature. Thus, setting this value to \samp{>0} will 'zoom out' on the retention time
#' axis.
#' 
#' \item \code{topMost} Only create EICs for this number of top most intense features. If \code{NULL} then EICs are
#' created for all features.
#' 
#' \item \code{topMostByRGroup} If set to \code{TRUE} and \code{topMost} is set: only create EICs for the top most
#' features in each replicate group. For instance, when \code{topMost=1} and \code{topMostByRGroup=TRUE}, then EICs will
#' be plotted for the most intense feature of each replicate group.
#' 
#' \item \code{onlyPresent} If \code{TRUE} then EICs are created only for analyses in which a feature was detected. If
#' \code{onlyPresent=FALSE} then EICs are generated for \strong{all} analyses. The latter is handy to evaluate if a peak
#' was 'missed' during feature detection or removed during \emph{e.g.} filtering.
#' 
#' }
#' 
#' if \code{onlyPresent=FALSE} then the following parameters are also relevant: \itemize{
#'
#' \item \code{mzExpWindow} To create EICs for analyses in which no feature was found, the \emph{m/z} value is derived
#' from the min/max values of all features in the feature group. The value of \code{mzExpWindow} further expands this
#' window.
#' 
#' \item \code{setsAdductPos},\code{setsAdductNeg} \setsWF In sets workflows the adduct must be known to calculate the
#' ionized \emph{m/z}. If a feature is completely absent in a particular set then it follows no adduct annotations are
#' available and the value of \code{setsAdductPos} (positive ionization data) or \code{setsAdductNeg} (negative
#' ionization data) will be used instead.
#'
#' }
#'
#' These parameters are passed as a named \code{list} as the \code{EICParams} argument to functions that use EICs. The
#' \code{getDefEICParams} function can be used to generate such parameter list with defaults.
#'
#' @param \dots optional named arguments that override defaults.
#'
#' @name EICParams
#' @export
getDefEICParams <- function(...)
{
    def <- list(
        rtWindow = 30,
        topMost = NULL,
        topMostByRGroup = FALSE,
        onlyPresent = TRUE,
        mzExpWindow = 0.001,
        setsAdductPos = "[M+H]+",
        setsAdductNeg = "[M-H]-"
    )
    
    return(modifyList(def, list(...), keep.null = TRUE))
}

#' Obtains a SIRIUS refresh token
#'
#' This function is used to obtain a \command{SIRIUS} refresh token with your login details, which allows
#' \code{\link{generateCompoundsSIRIUS}} to automatically log in to use \emph{e.g.} \command{CSI:FingerID}.
#'
#' \command{SIRIUS} version \samp{5} requires the user to be logged in when using web services such as
#' \command{CSI:FingerID}. To allow secure log in by external software tools such as \pkg{patRoon}, a \emph{referesh
#' token} is needed. This function uses your user name (email) and asks for a password (using \code{\link{getPass}}) to
#' obtain such a token more easily.
#'
#' More details for \emph{e.g.} account creation can be found on
#' \url{https://boecker-lab.github.io/docs.sirius.github.io/account-and-license/}.
#'
#' @param user A \code{character} string with the user name/email.
#'
#' @return A \code{character} string with the (secret) refresh token.
#'
#' @references \insertRef{Dhrkop2019}{patRoon}
#'
#' @export
getSIRIUSToken <- function(user)
{
    checkmate::assertString(user, min.chars = 1)
    
    cmd <- getCommandWithOptPath(getSiriusBin(), "SIRIUS")
    pw <- getPass::getPass("SIRIUS password:", noblank = TRUE)
    
    gotRefHeader <- FALSE
    ret <- NULL
    # NOTE: processx::run() is used as it allows correctly setting the environment, which doesn't seem to work very well
    # with base::system2()
    runv <- processx::run(cmd, c("login", "--user-env=SIRUSER", "--password-env=SIRPW", "--request-token-only"),
                          env = c("current", SIRUSER = user, SIRPW = pw), stdout_line_callback = function(out, px)
                          {
                              if (gotRefHeader && is.null(ret))
                                  ret <<- out
                              else
                                  gotRefHeader <<- grepl("^(\\#)+ Refresh token (\\#)+$", out)
                              
                          })
    
    if (is.null(ret))
    {
        cat(runv$stderr)
        stop("Failed to retrieve SIRIUS token! See error output above for details.", call. = FALSE)
    }
    
    return(ret)
}
