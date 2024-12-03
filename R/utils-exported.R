# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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

#' @export
getMSFileFormats <- function(fileType = NULL)
{
    assertMSFileType(fileType, null.ok = TRUE)
    ret <- names(MSFileExtensions())
    if (!is.null(fileType))
    {
        if (fileType == "raw")
            ret <- setdiff(ret, c("mzXML", "mzML"))
        else if (fileType == "ims")
            ret <- "mzML"
        else
            ret <- c("mzXML", "mzML")
    }
    return(ret)
}

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
#'   ignored. For valid values see \code{\link{getMSFileFormats}}.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(paths, groups = "", blanks = "", concs = NULL, norm_concs = NULL,
                                 formats = getMSFileFormats())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(groups, min.len = 1, add = ac)
    checkmate::assertCharacter(blanks, min.len = 1, add = ac)
    aapply(checkmate::assertNumeric, . ~ concs + norm_concs, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertSubset(formats, getMSFileFormats(), empty.ok = FALSE, add = ac)
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
            # NOTE: dirname point to current path if getExtDepPath() found it in PATH
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
    check("MetFrag CL", getExtDepPath("metfragcl", verify = FALSE), "patRoon.path.MetFragCL")
    check("MetFrag CompTox Database", getExtDepPath("metfragct", verify = FALSE), "patRoon.path.MetFragCompTox")
    check("MetFrag PubChemLite Database", getExtDepPath("metfragpcl", verify = FALSE), "patRoon.path.MetFragPubChemLite")
    check("OpenBabel", getExtDepPath("openbabel", verify = FALSE), "patRoon.path.obabel")
    check("BioTransformer", getExtDepPath("biotransformer", verify = FALSE), "patRoon.path.BioTransformer")
    
    if (!OK)
    {
        cat("\nSome dependencies were not found. Please make sure that their file locations are configured properly.",
            "Most dependencies can be easily installed and configured by installing the patRoonExt package.",
            "Otherwise, options() can be set to configure file locations. For instance, run the following to set the location of MetFragCL:",
            sprintf("options(patRoon.path.MetFragCL = \"C:/MetFragCommandLine-2.4.8.jar\")"),
            "\nPlease see the installation chapter in the handbook and ?patRoon for more details.",
            sep = "\n")
    }
    
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

#' @export
availableBackends <- function(anaInfo = NULL)
{
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, null.ok = TRUE)
    
    allBackends <- getMSReadBackends()
    
    unselected <- setdiff(allBackends, getOption("patRoon.MSBackends", character()))
    notCompiled <- allBackends[!sapply(allBackends, backendAvailable)]
    noAnas <- if (is.null(anaInfo))
        character()
    else
        allBackends[sapply(allBackends, function(b) is.null(maybeGetMSFiles(b, anaInfo, getMSFileTypes(), names(MSFileExtensions()))))]

    checkAvail <- function(b)
    {
        stat <- character()
        if (b %in% unselected)
            stat <- c(stat, "not in patRoon.MSBackends")
        if (b %in% notCompiled)
            stat <- c(stat, "not compiled during installation or unavailable on your system")
        if (b %in% noAnas)
            stat <- c(stat, "no suitable analyses found")
        if (length(stat) == 0)
            return("yes")
        return(sprintf("no (%s)", paste0(stat, collapse = ", ")))
    }
    
    for (b in allBackends)
        printf("Backend '%s': %s\n", b, checkAvail(b))
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
    
    # UNDONE: set file type
    dummyAI <- data.table(analysis = basename(tools::file_path_sans_ext(file)), path_centroid = dirname(file))
    return(doGetEICs(dummyAI, list(as.data.table(ranges)))[[1]])
}

#' @export
getBGMSMSPeaks <- function(anaInfo, rGroups = NULL, MSLevel = 2, retentionRange = NULL, mobilityRange = NULL,
                           minBPIntensity = 5000,
                           avgSpectraParams = getDefAvgPListParams(minAbundance = 0.1, topMost = 25),
                           avgAnalysesParams = getDefAvgPListParams(minAbundance = 0.8, topMost = 25))
{
    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, add = ac)
    checkmate::assertChoice(rGroups, unique(anaInfo$group), null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    assertRange(retentionRange, null.ok = TRUE, add = ac)
    checkmate::assertNumber(minBPIntensity, na.ok = FALSE, add = ac)
    aapply(assertAvgPListParams, . ~ avgSpectraParams + avgAnalysesParams, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(MSLevel, retentionRange, minBPIntensity, avgSpectraParams, avgAnalysesParams)
    
    if (!is.null(rGroups))
        anaInfo <- anaInfo[group %in% rGroups]
    if (is.null(retentionRange))
        retentionRange <- c(0, 0)
    if (is.null(mobilityRange))
        mobilityRange <- c(0, 0)
    
    printf("Averaging the spectra for each of the %d analyses\n", nrow(anaInfo))
    blSpecs <- applyMSData(anaInfo, function(ana, path, backend)
    {
        hash <- makeHash(baseHash, getMSDataFileHash(path))
        cached <- loadCacheData("avgBGMSMS", hash, cacheDB)
        if (!is.null(cached))
            return(cached)
        
        openMSReadBackend(backend, path)
        
        avgsp <- getMSPeakLists(backend, retentionRange[1], retentionRange[2], 0,
                                withPrecursor = FALSE, retainPrecursor = FALSE,
                                MSLevel = MSLevel, method = avgSpectraParams$method,
                                mzWindow = avgSpectraParams$clusterMzWindow,
                                startMobs = mobilityRange[1], endMobs = mobilityRange[2],
                                minAbundance = avgSpectraParams$minAbundance,
                                topMost = avgSpectraParams$topMost,
                                minIntensityIMS = avgSpectraParams$minIntensityIMS,
                                minIntensityPre = avgSpectraParams$minIntensityPre,
                                minIntensityPost = avgSpectraParams$minIntensityPost,
                                minBPIntensity = minBPIntensity)[[1]]
        setDT(avgsp)
        saveCacheData("avgBGMSMS", avgsp, hash, cacheDB)
        doProgress()
        return(avgsp)
    })
    blSpecs <- pruneList(blSpecs, checkZeroRows = TRUE)
    
    printf("Averaging analyses averaged spectra... ")
    ret <- averageSpectraList(list(blSpecs), avgAnalysesParams$clusterMzWindow, avgAnalysesParams$topMost,
                              avgAnalysesParams$minIntensityPre, avgAnalysesParams$minIntensityPost,
                              avgAnalysesParams$minAbundance, avgAnalysesParams$method, FALSE, FALSE, FALSE, FALSE)[[1]]
    ret[, precursor := NULL]
    printf("Done!\n")
    
    return(ret[])
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
    def <- getDefEIXParams()
    def$window <- 30
    
    mod <- list(...)
    if (!is.null(mod[["rtWindow"]]))
    {
        warning("The \"rtWindow\" parameter is deprecated and was renamed to \"window\"", call. = FALSE)
        mod$window <- mod$rtWindow
        mod$rtWindow <- NULL
    }
    
    return(modifyList(def, mod, keep.null = TRUE))
}

#' @name EICParams
#' @export
getDefEIMParams <- function(...)
{
    def <- getDefEIXParams()
    def <- modifyList(def, list(
        window = 0.2,
        maxRTWindow = 2,
        IMSWindow = 0.01,
        clusterMethod = "distance",
        minIntensity = 0
    ))
    return(modifyList(def, list(...), keep.null = TRUE))
}

#' @export
getDefPeakParams <- function(type, algorithm, ...)
{
    checkmate::assertChoice(type, c("chrom", "ims"))
    checkmate::assertChoice(algorithm, c("openms", "xcms3", "envipick", "dietrich"))
    
    def <- NULL
    # UNDONE: put all this in eg an internal YAML?
    if (type == "chrom")
    {
        def <- switch(algorithm,
               openms = list(
                   minPeakWidth = -1,
                   backgroundSubtraction = "none",
                   SGolayFrameLength = 15,
                   SGolayPolyOrder = 3,
                   useGauss = TRUE,
                   SN = 1.0,
                   SNWinLen = 1000,
                   SNBinCount = 30,
                   method = "corrected",
                   integrationType = "intensity_sum",
                   baselineType = "base_to_base",
                   fitEMG = FALSE,
                   extraOpts = list()
               ),
               xcms3 = list(
                   peakwidth = c(20, 50),
                   snthresh = 10,
                   prefilter = c(3, 100),
                   integrate = 1,
                   fitgauss = FALSE,
                   noise = 0
               ),
               envipick = list(
                   minpeak = 4,
                   drtsmall = 20,
                   drtfill = 10,
                   drttotal = 200,
                   recurs = 4,
                   weight = 2,
                   SB = 3,
                   SN = 2,
                   minint = 1E4,
                   maxint = 1E7,
                   ended = 2
               ),
               dietrich = list(
                   minIntensity = 1,
                   SN = 3,
                   peakWidth = c(5, 60),
                   RTRange = c(0, Inf),
                   maxPeaksPerSignal = 10
               )
        )
        def$forcePeakRange <- c(0, 0)
        def$relMinIntensity <- 0
    }
    else # IMS
    {
        # UNDONE!!
        def <- getDefPeakParams(algorithm, type = "chrom", ...)
        def$forcePeakRange = c(0.01, 0.1)
        def$relMinIntensity <- 0.25
    }
    
    return(modifyList(def, c(list(...), algorithm = algorithm)))
}

#' Parameters to aggregate concentrations/toxicity values assigned to feature groups
#'
#' Parameters that are used by method functions such \code{\link[=as.data.table,featureGroups-method]{as.data.table}} to
#' aggregate \link[=pred-quant]{predicted concentrations} or \link[=pred-tox]{toxicities}.
#'
#' Multiple concentration or toxicity values may be assigned to a single feature group. To ease the interpretation and
#' data handling, several functions aggregate these values prior their use. Aggregation occurs by the following data:
#'
#' \itemize{
#'
#' \item The candidate (\emph{i.e.} suspect or annotation candidate). This is mainly relevant for sets workflows, where
#' calculations among sets may yield different results for the same candidate.
#'
#' \item The prediction type, \emph{e.g.} all values that were obtained from suspect or compound annotation data.
#'
#' \item The feature group.
#'
#' }
#'
#' The aggregation of all data first occurs by the same candidate/type/feature group, then the same type/feature group
#' and finally for each feature group. This ensures that \emph{e.g.} large numbers of data points for a prediction type
#' do not bias results.
#'
#' The \code{candidateFunc}, \code{typeFunc} and \code{groupFunc} parameters specify the function that should be used to
#' aggregate data. Commonly, functions such \code{\link{mean}}, \code{\link{min}} or \code{\link{max}} can be used here.
#' Note that the function does not need to handle \code{NA} values, as these are removed in advance.
#'
#' The \code{preferType} parameters specifies the \emph{preferred} prediction type. Any values from other prediction
#' types will be ignored unless the preferred type is not available for a feature group. Valid values are
#' \code{"suspect"} (the default), \code{"compound"} (results from compound annotation by \acronym{SMILES}),
#' \code{"SIRIUS_FP"} (results from formula/compound annotation with \command{SIRIUS+CSI:FingerID}) or \code{"none"}.
#'
#' These parameters should be stored inside a \code{list}. The \code{getDefPredAggrParams} function can be used to
#' generate such parameter list with defaults.
#'
#' @param all The default aggregation function for all types, \emph{e.g.} \code{mean}.
#' @param \dots optional named arguments that override defaults.
#'
#' @name pred-aggr-params
#' @export
getDefPredAggrParams <- function(all = mean, ...)
{
    def <- list(
        candidateFunc = all,
        typeFunc = all,
        groupFunc = all,
        preferType = "suspect"
    )
    
    return(modifyList(def, list(...)))
}

#' @export
getCCSParams <- function(method, ..., calibrant = NULL)
{
    ret <- list(
        method = method,
        defaultCharge = 1,
        temperature = 305,
        massGas = 28.0134,
        MasonSchampConstant = 18500,
        calibrant = calibrant
    )
    
    return(modifyList(ret, list(...)))
}

#' @rdname pred-quant
#' @export
getQuantCalibFromScreening <- function(fGroups, concs, areas = FALSE, average = FALSE, IMS = "maybe")
{
    # UNDONE: mention that duplicate suspects (name.x) are ignored
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkClass(fGroups, "featureGroupsScreening"),
        checkmate::checkClass(fGroups, "featureGroupsScreeningSet"),
        .var.name = "fGroups", add = ac)
    checkmate::assertDataFrame(concs, min.cols = 2, min.rows = 1, add = ac)
    checkmate::assertNames(names(concs), must.include = "name", add = ac)
    checkmate::assertCharacter(concs$name, any.missing = FALSE, min.chars = 1, unique = TRUE, add = ac)
    concRGs <- setdiff(names(concs), "name")
    concRGs <- intersect(concRGs, replicateGroups(fGroups))
    if (length(concRGs) == 0)
        stop("No concentration columns for (relevant) replicate groups found.", call. = FALSE)
    for (col in concRGs)
        checkmate::assertNumeric(concs[[col]], finite = TRUE, .var.name = paste0("concs[['", col, "']]"), add = ac)
    aapply(checkmate::assertFlag, . ~ areas + average, fixed = list(add = ac))
    assertIMSArg(IMS, add = ac)
    checkmate::reportAssertions(ac)

    fGroups <- selectIMSFilter(fGroups, IMS, verbose = FALSE, warn = FALSE)
    
    anaInfo <- analysisInfo(fGroups)
    
    if (!is.data.table(concs))
        concs <- as.data.table(concs)
    
    tab <- as.data.table(fGroups, areas = areas, average = average, collapseSuspects = NULL)
    if (any(!concs$name %chin% tab$susp_name))
    {
        warning("Ignoring suspects not present in screening results: ",
                paste0(setdiff(concs$name, tab$susp_name), collapse = ", "), call. = FALSE)
        concs <- concs[name %chin% tab$susp_name]
        if (nrow(concs) == 0)
            stop("None of the suspects present in the screening results, aborting...", call. = FALSE)
    }
    
    if (is.null(tab[["susp_SMILES"]]))
        stop("Screening results lack SMILES data.", call. = FALSE)
    
    ret <- data.table(name = concs$name, SMILES = tab[match(concs$name, susp_name)]$susp_SMILES,
                      group = tab[match(concs$name, susp_name)]$group)
    ret[, rt := groupInfo(fGroups)$ret[match(ret$group, groupInfo(fGroups)$group)]]
    
    if (anyNA(ret$SMILES))
    {
        warning("Ignoring suspects without SMILES: ", paste0(ret[is.na(SMILES)]$name, collapse = ", "), call. = FALSE)
        concs <- concs[!is.na(SMILES)]
        if (nrow(concs) == 0)
            stop("No suspects with SMILES present, aborting...", call. = FALSE)
    }
    
    intCols <- if (average)
        concRGs
    else
        unlist(lapply(concRGs, function(rg) anaInfo[group == rg]$analysis))
    intCols <- getADTIntCols(intCols)
    ret <- merge(ret, tab[, c("susp_name", intCols), with = FALSE], by.x = "name", by.y = "susp_name")
    
    vname <- if (average) "rGroup" else "analysis"
    ret <- melt(ret, measure.vars = intCols, variable.name = vname, value.name = "intensity")
    ret[, (vname) := stripADTIntSuffix(get(vname))]
    
    if (!average)
        ret[, rGroup := anaInfo$group[match(analysis, anaInfo$analysis)]]
    
    mconcs <- melt(concs, measure.vars = concRGs, variable.name = "rGroup", value.name = "conc") 
    
    ret <- merge(ret, mconcs[, c("name", "rGroup", "conc"), with = FALSE], by = c("name", "rGroup"))
    
    ret <- removeDTColumnsIfPresent(ret, c("analysis", "rGroup", "group"))
    setcolorder(ret, c("name", "SMILES", "rt", "conc", "intensity"))
    
    return(ret[])
}

#' Obtain default rules for metabolic logic
#'
#' This function returns a \code{data.frame} with the default rules for metabolic logic, which can be used by
#' \code{\link{generateTPsLogic}} and \code{\link{genFormulaTPLibrary}}.
#'
#' @return A \code{data.frame} with columns describing each transformation rule.
#'
#' @section Source: The table is based on the work done by Schollee \emph{et al.} (see references).
#'
#' @references \insertRef{Scholle2015}{patRoon}
#' @export
TPLogicTransformations <- function()
{
    return(patRoon:::TPsLogicTransformations) # stored inside R/sysdata.rda
}

#' @export
convertMobilityToCCS <- function(mobility, mz, CCSParams, charge = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertMobilityConversionArgs(mobility, mz, CCSParams, charge, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(charge))
        charge <- rep(CCSParams$defaultCharge, length(mobility))
    charge <- abs(charge)
    
    u = ((mz * CCSParams$massGas) / (mz + CCSParams$massGas))
    
    if (CCSParams$method == "bruker")
    {
        if (!doInitBrukerLib())
            stop("Cannot convert mobility data with Bruker TIMS library: failed to initialize library", call. = FALSE)
        return(getBrukerCCS(mobility, charge, mz))
    }
    else if (CCSParams$method %in% c("mason-schamp_k", "mason-schamp_1/k"))
    {
        if (CCSParams$method == "mason-schamp_k")
            mobility <- 1 / mobility
        return(CCSParams$MasonSchampConst * (charge / (sqrt(u * CCSParams$temperature))) * mobility)
    }
    else if (CCSParams$method == "agilent")
    {
        calibrant <- prepareAgilentIMSCalib(CCSParams$calibrant, CCSParams$massGas)
        return((mobility - calibrant$TFix) * charge / calibrant$beta / sqrt(mz / (mz + calibrant$massGas)))
    }
    # UNDONE: support waters later, procedure seems unclear how it affects different instruments and need test data
    # else # waters
    # {
    #     coeff <- t0 <- exponent <- EDCDelay <- NA_real_ # UNDONE
    #     mobCorr <- mobility - (0.001 * EDCDelay * sqrt(mz))
    #     omagec <- coeff * (mobCorr + t0)^exponent
    #     return((omagec * charges) / sqrt(u))
    # }
}

#' @export
convertCCSToMobility <- function(ccs, mz, CCSParams, charge = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertMobilityConversionArgs(ccs, mz, CCSParams, charge, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(charge))
        charge <- rep(CCSParams$defaultCharge, length(mobility))
    charge <- abs(charge)
    
    u = ((mz * CCSParams$massGas) / (mz + CCSParams$massGas))
    
    if (CCSParams$method == "bruker")
    {
        if (!doInitBrukerLib())
            stop("Cannot convert mobility data with Bruker TIMS library: failed to initialize library", call. = FALSE)
        return(getBrukerMob(ccs, charge, mz))
    }
    else if (CCSParams$method %in% c("mason-schamp_k", "mason-schamp_1/k"))
    {
        mob <- ccs / CCSParams$MasonSchampConst / (charge / (sqrt(u * CCSParams$temperature)))
        if (CCSParams$method == "mason-schamp_k")
            mob <- 1 / mob
        return(mob)
    }
    else if (CCSParams$method == "agilent")
    {
        calibrant <- prepareAgilentIMSCalib(CCSParams$calibrant, CCSParams$massGas)
        return((ccs * sqrt(mz / (mz + calibrant$massGas)) * calibrant$beta / charge) + calibrant$TFix)
    }
    # UNDONE: support waters later, procedure seems unclear how it affects different instruments and need test data
    # else # waters
    # {
    # }
}
