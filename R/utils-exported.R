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
#' @param groups,blanks An (optional) character vector containing replicates and blanks, respectively (will be
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
generateAnalysisInfo <- function(paths, replicates = "", blanks = "", concs = NULL, norm_concs = NULL,
                                 formats = getMSFileFormats())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(replicates, min.len = 1, add = ac)
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
    replicates <- rep(replicates, length.out = nrow(ret))
    ret$replicate <- ifelse(!nzchar(replicates), ret$analysis, replicates)
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
        
        enviSInfo[Type == "blank", replicate := blank]
    }
    else
        enviSInfo[, blank := ""]
    
    enviSInfo[replicate == "FALSE", replicate := ""]
    
    ret <- data.frame(path = file.path(path, "files"), analysis = enviSInfo$ID, replicate = enviSInfo$replicate,
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
#' @details Fold change calculation can be used to easily identify significant changes between replicates. The
#'   calculation process is configured through a paramater list, which can be constructed with the \code{getFCParams}
#'   function. The parameter list has the following entries: \itemize{
#'
#'   \item \code{replicates} the name of the two replicates to compare (taken from the \code{replicates} argument to
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
#' @param replicates A \code{character} vector with the names of the two replicates to be compared.
#' @param \dots Optional named arguments that override defaults.
#' 
#' @author The code to calculate and plot Fold change data was created by Bas van de Velde.
#'
#' @seealso \code{\link{featureGroups-class}} and \link{feature-plotting}
#'
#' @export
getFCParams <- function(replicates, ...)
{
    checkmate::assertCharacter(replicates, min.chars = 1, len = 2, any.missing = FALSE)
    
    def <- list(replicates = replicates,
                thresholdFC = 0.25,
                thresholdPV = 0.05,
                zeroMethod = "add",
                zeroValue = 0.01,
                PVTestFunc = function(x, y) t.test(x, y, paired = TRUE)$p.value,
                PVAdjFunc = function(pv) p.adjust(pv, "BH"))
    return(modifyList(def, list(...)))
}

#' @export
availableBackends <- function(anaInfo = NULL, verbose = TRUE)
{
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, null.ok = TRUE)
    checkmate::assertFlag(verbose)
    
    allBackends <- getMSReadBackends()
    
    unselected <- setdiff(allBackends, getOption("patRoon.MS.backends", character()))
    notCompiled <- allBackends[!sapply(allBackends, backendAvailable)]
    noAnas <- if (is.null(anaInfo))
        character()
    else
        allBackends[sapply(allBackends, function(b) is.null(maybeGetMSFiles(b, anaInfo, getMSFileTypes(), names(MSFileExtensions()))))]

    checkAvail <- function(b)
    {
        stat <- character()
        if (b %in% unselected)
            stat <- c(stat, "not in patRoon.MS.backends")
        if (b %in% notCompiled)
            stat <- c(stat, "not compiled during installation or unavailable on your system")
        if (b %in% noAnas)
            stat <- c(stat, "no suitable analyses found")
        if (length(stat) == 0)
            return("yes")
        return(sprintf("no (%s)", paste0(stat, collapse = ", ")))
    }

    check <- sapply(allBackends, checkAvail)
        
    if (verbose)
    {
        for (b in allBackends)
            printf("Backend '%s': %s\n", b, check[[b]])
    }
    
    return(invisible(names(check[sapply(check, identical, "yes")])))
}

#' Obtains extracted ion chromatograms (EICs)
#'
#' This function generates one or more EIC(s) for given retention time, \emph{m/z} and optionally mobility ranges.
#'
#' @template analysisInfo-arg
#'
#' @param ranges A \code{list} with for each analysis a \code{data.frame} with \code{numeric} columns \code{"retmin"},
#'   \code{"retmax"}, \code{"mzmin"}, \code{"mzmax"} with the lower/upper ranges of the retention time and \emph{m/z}.
#'   Furthermore, columns \code{"mobmin"} and \code{"mobmax"} can be added for IMS data.
#'
#' @return A \code{list} with for each analysis a \code{list} with EIC data for each of the rows in \code{ranges}.
#'
#' @export
getEICs <- function(analysisInfo, ranges, output = "fill", minIntensityIMS = 25)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assert(
        checkmate::checkList(ranges, len = nrow(analysisInfo)),
        checkmate::checkDataFrame(ranges, types = "numeric", any.missing = FALSE),
        .var.name = "ranges", add = ac
    )
    checkmate::assertChoice(output, c("fill", "pad", "raw"), add = ac)
    checkmate::assertNumber(minIntensityIMS, lower = 0, finite = TRUE, na.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)

    if (checkmate::testDataFrame(ranges))
        ranges <- rep(list(as.data.table(ranges)), nrow(analysisInfo))
    
    if (!checkmate::testNamed(ranges))
        names(ranges) <- analysisInfo$analysis
    else
        checkmate::assertSetEqual(names(ranges), analysisInfo$analysis)
    
    for (r in ranges)
    {
        checkmate::assertDataFrame(r, types = "numeric", any.missing = FALSE)
        assertHasNames(r, c("mzmin", "mzmax", "retmin", "retmax"))
    }
    ret <- doGetEICs(analysisInfo, ranges, withBP = TRUE, minIntensityIMS = minIntensityIMS)
    if (output != "raw")
    {
        ret <- lapply(ret, function(anaEICs)
        {
            at <- attr(anaEICs, "allXValues")
            if (is.null(at))
                return(anaEICs) # no EICs
        
            if (output == "fill")
                return(lapply(anaEICs, function(eic) data.frame(time = at,
                                                                intensity = doFillEIXIntensities(at, eic$time, eic$intensity))))
            
            # output == "pad"
            return(lapply(anaEICs, function(eic)
            {
                eic <- setDF(padEIX(at, eic$time, eic$intensity))
                setnames(eic, "xvalue", "time")
            }))
        })
    }
    return(ret)
}

#' @export
getBGMSMSPeaks <- function(anaInfo, replicates = NULL, MSLevel = 2, retentionRange = NULL, mobilityRange = NULL,
                           minBPIntensity = 5000,
                           avgSpectraParams = getDefAvgPListParams(minAbundanceRel = 0.1, topMost = 25),
                           avgAnalysesParams = getDefAvgPListParams(minAbundanceAbs = 0.8, topMost = 25))
{
    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, add = ac)
    checkmate::assertChoice(replicates, unique(anaInfo$replicate), null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    assertRange(retentionRange, null.ok = TRUE, add = ac)
    checkmate::assertNumber(minBPIntensity, na.ok = FALSE, add = ac)
    aapply(assertAvgPListParams, . ~ avgSpectraParams + avgAnalysesParams, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(MSLevel, retentionRange, minBPIntensity, avgSpectraParams, avgAnalysesParams)
    
    if (!is.null(replicates))
        anaInfo <- anaInfo[replicate %in% replicates]
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
                                minAbundanceRel = avgSpectraParams$minAbundanceRel,
                                minAbundanceAbs = avgSpectraParams$minAbundanceAbs,
                                minAbundanceIMSRel = avgSpectraParams$minAbundanceIMSRel,
                                minAbundanceIMSAbs = avgSpectraParams$minAbundanceIMSAbs,
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
                              avgAnalysesParams$minAbundanceRel, avgAnalysesParams$minAbundanceAbs,
                              avgAnalysesParams$method, FALSE, FALSE, FALSE, FALSE)[[1]]
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
#' \item \code{topMostByReplicate} If set to \code{TRUE} and \code{topMost} is set: only create EICs for the top most
#' features in each replicate. For instance, when \code{topMost=1} and \code{topMostByReplicate=TRUE}, then EICs will
#' be plotted for the most intense feature of each replicate.
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
        clusterMethod = "distance"
    ))
    return(modifyList(def, list(...), keep.null = TRUE))
}

#' @export
getDefPeakParams <- function(type, algorithm, ...)
{
    checkmate::assertChoice(type, c("chrom", "bruker_ims", "agilent_ims"))
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
                   gaussWidth = 30,
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
                   peakwidth = c(10, 30),
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
    else if (type == "bruker_ims")
    {
        def <- getDefPeakParams(algorithm, type = "chrom", ...)
        def <- modifyList(def, switch(
            algorithm,
            openms = list(gaussWidth = 0.02),
            xcms3 = list(peakwidth = c(0.01, 0.5), prefilter = c(3, 10), firstBaselineCheck = FALSE),
            envipick = list(drtsmall = 0.2, drtfill = 0.02, drttotal = 1, minint = 10),
            dietrich = list(peakWidth = c(0.02, 0.5))
        ))
        def$forcePeakRange = c(0.01, 0.1)
        def$relMinIntensity <- 0.25
    }
    else if (type == "agilent_ims")
    {
        def <- getDefPeakParams(algorithm, type = "chrom", ...)
        def <- modifyList(def, switch(
            algorithm,
            openms = list(gaussWidth = 0.05),
            xcms3 = list(peakwidth = c(0.1, 5), prefilter = c(3, 10), firstBaselineCheck = FALSE),
            envipick = list(drtsmall = 2, drtfill = 0.4, drttotal = 1, minint = 10),
            dietrich = list(peakWidth = c(0.2, 5))
        ))
        def$forcePeakRange = c(0.2, 1)
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

#' @export
getIMSRangeParams <- function(param, lower, upper, mzRelative = FALSE)
{
    return(list(param = param, lower = lower, upper = upper, mzRelative = mzRelative))
}

#' @export
getIMSMatchParams <- function(param, ...)
{
    ret <- list(
        param = param,
        relative = TRUE,
        minMatches = 0
    )
    
    ret <- modifyList(ret, list(...))
    
    if (is.null(ret[["window"]]))
    {
        # UNDONE: tweak
        ret$window <- if (ret$relative)
            0.05
        else if (ret$param == "mobility")
            0.03
        else # CCS
            10
    }
    
    return(ret)
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
    concReps <- setdiff(names(concs), "name")
    concReps <- intersect(concReps, replicates(fGroups))
    if (length(concReps) == 0)
        stop("No concentration columns for (relevant) replicates found.", call. = FALSE)
    for (col in concReps)
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
        concReps
    else
        unlist(lapply(concReps, function(rg) anaInfo[group == rg]$analysis))
    intCols <- getADTIntCols(intCols)
    ret <- merge(ret, tab[, c("susp_name", intCols), with = FALSE], by.x = "name", by.y = "susp_name")
    
    vname <- if (average) "replicate" else "analysis"
    ret <- melt(ret, measure.vars = intCols, variable.name = vname, value.name = "intensity")
    ret[, (vname) := stripADTIntSuffix(get(vname))]
    
    if (!average)
        ret[, replicate := anaInfo$replicate[match(analysis, anaInfo$analysis)]]
    
    mconcs <- melt(concs, measure.vars = concReps, variable.name = "replicate", value.name = "conc") 
    
    ret <- merge(ret, mconcs[, c("name", "replicate", "conc"), with = FALSE], by = c("name", "replicate"))
    
    ret <- removeDTColumnsIfPresent(ret, c("analysis", "replicate", "group"))
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
    charge <- rep(abs(charge), length.out = length(mobility))
    
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
    charge <- rep(abs(charge), length.out = length(ccs))
    
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

#' @export
setMethod("assignMobilities", "data.table", function(obj, from = NULL, matchFromBy = "InChIKey",
                                                     overwrite = FALSE, adducts = c("[M+H]+", "[M-H]-", "none"),
                                                     adductDef = NULL, predictAdductOnly = TRUE, CCSParams = NULL,
                                                     prepareChemProps = TRUE, prefCalcChemProps = TRUE,
                                                     neutralChemProps = FALSE, virtualenv = "patRoon-c3sdb")
{
    # NOTE: keep args in sync with other methods
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkChoice(from, c("pubchemlite", "c3sdb")),
        checkmate::checkDataFrame(from, min.rows = 1),
        checkmate::checkNull(from),
        .var.name = "from", add = ac
    )
    checkmate::assertChoice(matchFromBy, c("InChIKey", "InChIKey1", "InChI", "SMILES", "name"), add = ac)
    aapply(checkmate::assertFlag, . ~ overwrite + predictAdductOnly + prepareChemProps + prefCalcChemProps +
               neutralChemProps, fixed = list(add = ac))
    checkmate::assertCharacter(adducts, min.chars = 1, any.missing = FALSE, add = ac)
    assertCCSParams(CCSParams, null.ok = TRUE, add = ac)
    checkmate::assertString(virtualenv, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(obj[["adduct"]]))
        adducts <- union(adducts, obj[!is.na(adduct) & nzchar(adduct)]$adduct) # expand with adducts from suspect list
    if (!is.null(adductDef))
        adductDef <- checkAndToAdduct(adductDef)

    hash <- makeHash(obj, from, matchFromBy, overwrite, adducts, adductDef, predictAdductOnly, CCSParams,
                     prepareChemProps, prefCalcChemProps, neutralChemProps)
    cd <- loadCacheData("assignMobilitiesDT", hash)
    if (!is.null(cd))
        return(cd)
    
    do_c3sdb <- identical(from, "c3sdb")
    if (do_c3sdb)
    {
        checkPackage("reticulate")
        if (is.null(obj[["SMILES"]]))
            stop("No SMILES data found in the input data, this is necessary for c3sdb predictions.", call. = FALSE)
    }
    
    obj <- copy(obj)
    
    if (prepareChemProps)
    {
        if (is.null(adductDef) && (is.null(obj[["adduct"]]) || anyNA(obj$adduct) || any(!nzchar(obj$adduct))))
            stop("Please set the 'adductDef' argument so that m/z values can be calculated.", call. = FALSE)
        obj <- prepareChemTable(obj, prefCalcChemProps, neutralChemProps)
        obj <- calcSuspMZs(obj, adductDef)
    }
    
    if (do_c3sdb || !is.null(CCSParams))
    {
        msg <- paste("Mass data is necessary for c3sdb pedictions and mobility <--> CCS conversions.",
                     "This data can be automatically calculated by setting prepareChemProps=TRUE and ensuring a SMILES,",
                     "InChI or formula column is present.")
        if (is.null(obj[["mz"]]))
            stop(paste("No m/z data found in the input data (mz column).", msg), call. = FALSE)
        if (anyNA(obj$mz))
            warning(paste(sprintf("The following rows do not contain m/z data and will be ignored: %s.",
                                  paste0(obj[is.na(mz), which = TRUE], collapse = ", ")), msg), call. = FALSE)
    }
    
    if (do_c3sdb)
    {
        if (is.null(obj[["SMILES"]]))
            stop("No SMILES data found in the input data (SMILES column). This is necessary for c3sdb predictions",
                 call. = FALSE)
        if (anyNA(obj$SMILES))
            warning("The following rows do not contain SMILES data and will be excluded from c3sdb predictions: ",
                    paste0(obj[is.na(SMILES), which = TRUE], collapse = ", "), call. = FALSE)
    }
    
    if (!is.null(from))
    {
        if (matchFromBy == "InChIKey1" && is.null(obj[["InChIKey1"]]) && !is.null(obj[["InChIKey"]]))
            obj[, InChIKey1 := getInChIKey1(InChIKey)]
        
        if (is.null(obj[[matchFromBy]]))
            stop(sprintf("Column '%s' not found to match input data.", matchFromBy), call. = FALSE)
        if (anyNA(obj[[matchFromBy]]))
            warning(sprintf("The following rows in the input data are NA in the match column ('%s') and are therefore not matched: %s.",
                            matchFromBy, paste0(obj[is.na(get(matchFromBy)), which = TRUE], collapse = ", ")),
                    call. = FALSE)
        
        adductsNoNone <- setdiff(adducts, "none")
        
        if (is.data.frame(from))
        {
            from <- makeDT(from)
            from <- prepareChemTable(from, prefCalcChemProps, neutralChemProps)
        }
        else if (from == "pubchemlite")
        {
            printf("Loading PubChemLite... ")
            path <- getExtDepPath("metfragpcl")
            from <- fread(path)
            printf("Done!\n")
            setnames(from, "CompoundName", "name") # UNDONE: doc that this column is used for name matching, support Synonym col as well?
            pat <- "^pred_CCS_A2_"
            predCols <- grep(pat, names(from), value = TRUE)
            if (length(predCols) == 0)
                stop(sprintf("No CCS prediction columns found in the configured PubChemLite data (%s). Please ensure a version with CCS values.",
                             path), call. = FALSE)
            setnames(from, predCols, gsub(pat, "CCS_", predCols))
        }
        else # c3sdb
        {
            if (!is.null(virtualenv))
                reticulate::use_virtualenv(virtualenv)
            py_pickle <- reticulate::import("pickle")
            py_c3sdb <- reticulate::import("c3sdb.ml.data")
            py_bi <- reticulate::import_builtins()
            
            kmcm_svr <- with(py_bi$open(py_c3sdb$pretrained_data("c3sdb_kmcm_svr.pkl"), "rb"), as = "pf",
                             py_pickle$load(pf))
            
            # only do predictions for suspects with SMILES and mass data
            objDo <- obj[!is.na(SMILES) & nzchar(SMILES) & !is.na(mz)]
            
            doPred <- function(mzs, SMILES, add)
            {
                if (length(mzs) == 0)
                    return(numeric())
                
                if (!add %in% c3sdbAdducts())
                {
                    warning(sprintf("Skipping predictions for unsupported adduct '%s'.", add), call. = FALSE)
                    return(rep(NA_real_, length(mzs))) 
                }

                # use as.list() to handle scalars: https://github.com/rstudio/reticulate/issues/258
                dfi <- py_c3sdb$data_for_inference(as.list(mzs), as.list(rep(add, length(mzs))), as.list(SMILES),
                                                   py_c3sdb$pretrained_data("c3sdb_OHEncoder.pkl"),
                                                   py_c3sdb$pretrained_data("c3sdb_SScaler.pkl"))
                ret <- rep(NA_real_, length(mzs))
                ret[dfi[[2]]] <- kmcm_svr$predict(dfi[[1]])
                if (any(!dfi[[2]]))
                    warning(sprintf("Predictions failed for adduct '%s' for the following SMILES: %s.", add,
                                    paste0(SMILES[dfi[[2]] == FALSE], collapse = ", ")), call. = FALSE)
                return(ret)
            }
            
            from <- data.table(objDo[[matchFromBy]], mz = objDo$mz, SMILES = objDo$SMILES)
            setnames(from, 1, matchFromBy)
            
            hasAddCol <- !is.null(objDo[["adduct"]])
            
            for (add in adductsNoNone)
            {
                wh <- !hasAddCol | !predictAdductOnly | is.na(objDo$adduct) | !nzchar(objDo$adduct) | objDo$adduct == add
                wh <- rep(wh, length.out = nrow(objDo))
                if (any(wh))
                {
                    printf("Predicting %d CCS values for adduct '%s'... ", sum(wh), add)
                    from[wh, paste0("CCS_", add) := doPred(mz, SMILES, add)]
                    printf(" Done!\n")
                }
            }
        }
        
        # merge in data in input obj
        
        addMobCols <- paste0("mobility_", adductsNoNone)
        addCCSCols <- paste0("CCS_", adductsNoNone)
        if ("none" %in% adducts)
        {
            addMobCols <- c(addMobCols, "mobility")
            addCCSCols <- c(addCCSCols, "CCS")
        }
        predCols <- c(grep("^mobility", names(from), value = TRUE), grep("^CCS", names(from), value = TRUE))
        predCols <- intersect(predCols, c(addMobCols, addCCSCols))
        if (length(predCols) == 0)
            stop("No (adduct relevant) mobility/CCS columns found in the provided data.", call. = FALSE)
        
        if (matchFromBy == "InChIKey1" && is.null(from[["InChIKey1"]]) && !is.null(from[["InChIKey"]]))
            from[, InChIKey1 := getInChIKey1(InChIKey)]
        
        if (is.null(from[[matchFromBy]]))
            stop(sprintf("Column '%s' not found to match data from.", matchFromBy), call. = FALSE)

        for (col in predCols)
        {
            if (is.null(obj[[col]]))
            {
                m <- match(obj[[matchFromBy]], from[[matchFromBy]])
                set(obj, j = col, value = from[[col]][m])
            }
            else
            {
                takeVal <- function(x) overwrite | is.na(x) | (is.character(x) & !nzchar(x))
                obj[from, (col) := fifelse(takeVal(get(col)), get(paste0("i.", col)), get(col)),
                    on = matchFromBy][]
            }
        }
    }
    
    # fill in missing mobility/CCS values by conversions
    if (!is.null(CCSParams))
    {
        doConvert <- function(start, values, masses, charges)
        {
            doIt <- if (start == "mobility") convertMobilityToCCS else convertCCSToMobility

            if (!is.character(values))
                return(doIt(values, masses, CCSParams, charges))
            
            charges <- rep(charges, length.out = length(values))
            
            # handle semi-colon separated values from suspect lists
            tab <- rbindlist(lapply(seq_along(values), function(i)
            {
                list(x = as.numeric(strsplit(values[i], ";", fixed = TRUE)[[1]]), m = masses[i], ch = charges[i])
            }), idcol = TRUE)
            tab[, convx := as.character(doIt(x, m, CCSParams, ch))]
            tab[, convx := paste0(convx, collapse = ";"), by = .id]
            return(unique(tab, by = ".id")$convx)
        }
        
        printf("Calculating missing mobility/CCS values... ")
        for (addChr in adducts)
        {
            mCol <- "mobility"; cCol <- "CCS"
            if (addChr != "none")
            {
                mCol <- paste0(mCol, "_", addChr); cCol <- paste0(cCol, "_", addChr)
            }
            
            if (is.null(obj[[mCol]]) && is.null(obj[[cCol]]))
                next
            
            charges <- if (addChr == "none")
            {
                if (!is.null(obj[["adduct"]]))
                {
                    wh <- !is.na(obj$adduct) & nzchar(obj$adduct)
                    unAddChr <- unique(obj$adduct[wh])
                    unAddCharge <- sapply(unAddChr, function(uac) as.adduct(uac)@charge)
                    fifelse(wh, unAddCharge[obj$adduct], CCSParams$defaultCharge)
                }
                else
                    CCSParams$defaultCharge
            }
            else
                as.adduct(addChr)@charge
            
            if (is.null(obj[[mCol]]))
                obj[, (mCol) := doConvert("CCS", get(cCol), obj$mz, charges)]
            else if (is.null(obj[[cCol]]))
                obj[, (cCol) := doConvert("mobility", get(mCol), obj$mz, charges)]
            else # both present: only consider NAs
            {
                hasVal <- function(x) !is.na(x) & (!is.character(x) | nzchar(x))
                # temporarily add charges to handle subset assignments below
                # UNDONE: handle the (unlikely) case that the column is already present?
                obj[, .charge := charges]
                obj[!hasVal(get(mCol)) & hasVal(get(cCol)), (mCol) := doConvert("CCS", get(cCol), obj$mz, .charge)]
                obj[hasVal(get(mCol)) & !hasVal(get(cCol)), (cCol) := doConvert("mobility", get(mCol), obj$mz, .charge)]
                obj[, .charge := NULL]
            }
        }
        printf("Done!\n")
    }
    
    saveCacheData("assignMobilitiesDT", obj, hash)
    
    return(obj[])
})

#' @export
setMethod("assignMobilities", "data.frame", function(obj, ...)
{
    return(as.data.frame(assignMobilities(as.data.table(obj), ...)))
})

#' @export
installc3sdb <- function(envname = "patRoon-c3sdb", clearEnv = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(envname, min.chars = 1, null.ok = TRUE)
    checkmate::assertFlag(clearEnv, add = ac)
    checkmate::reportAssertions(ac)
    
    checkPackage("reticulate")
    
    if (clearEnv && !is.null(envname))
        reticulate::virtualenv_remove(envname)
    
    reticulate::py_install("git+https://github.com/dylanhross/c3sdb", envname = envname, ...)
}

#' @export
installTIMSCONVERT <- function(envname = "patRoon-TIMSCONVERT", clearEnv = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(envname, min.chars = 1, null.ok = TRUE)
    checkmate::assertFlag(clearEnv, add = ac)
    checkmate::reportAssertions(ac)
    
    checkPackage("reticulate")
    
    if (clearEnv && !is.null(envname))
        reticulate::virtualenv_remove(envname)
    
    reticulate::py_install("git+https://github.com/gtluu/timsconvert", envname = envname, ...)
}
