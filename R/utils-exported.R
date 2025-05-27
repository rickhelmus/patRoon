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

#' Get supported MS file types
#' @export
getMSFileTypes <- function() c("raw", "centroid", "profile", "ims")

#' Get supported MS file formats
#'
#' Returns the supported file formats of MS files in \pkg{patRoon}.
#'
#' @param fileType The type of file for which formats should be returned (see \code{\link{getMSFileTypes}}). If
#'   \code{NULL} then all file formats are returned.
#'
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
            ret <- c("mzML", "mzXML")
    }
    return(ret)
}

#' @details \code{generateAnalysisInfo} is an utility function that automatically generates analysis information. It
#'   scans given directories for analysis files, and uses this to automatically fill in the \code{analysis} and
#'   \code{path_*} columns. This function automatically groups together analyses that are stored with different file
#'   types and formats (see further details below).
#'
#' @param fromRaw,fromCentroid,fromProfile,fromIMS One or more file paths that should be used for finding analyses that
#'   are stored as raw, centroided, profile or IMS data, respectively (see details below). Set to \code{NULL} to skip
#'   file detection for a particular file type.
#' @param convCentroid,convProfile,convIMS These arguments specify the \link[=MSConversion]{MS file conversion}
#'   destination paths for centroided, profile and IMS data, respectively. These paths are used for those analyses for
#'   which no file with a particular file type could be found in the directories specified by the respective
#'   \code{from*} arguments. Set to \code{NULL} to not set any destination directory. If multiple paths are specified
#'   then these will be recycled to fill the table rows.
#' @param \dots Any other columns that should be added to the analysis information table, such as \code{replicate} and
#'   \code{blank}. The arguments specified by \dots should be named. Vectors are recycled to the number of rows of the
#'   table.
#'
#' @return \code{generateAnalysisInformation} returns a \code{data.frame} with automatically generated analysis
#'   information.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(fromRaw = NULL, fromCentroid = NULL, fromProfile = NULL, fromIMS = NULL,
                                 convCentroid = NULL, convProfile = NULL, convIMS = NULL, ...)
{
    otherArgs <- list(...)

    ac <- checkmate::makeAssertCollection()
    aapply(assertDirExists, . ~ fromRaw + fromCentroid + fromProfile + fromIMS, access = "r", null.ok = TRUE,
           fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ convCentroid + convProfile + convIMS, min.chars = 1,
           any.missing = FALSE, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    if (length(otherArgs) > 0)
        checkmate::assertNames(names(otherArgs), "unique", .var.name = "...", add = ac)
    checkmate::qassertr(otherArgs, "V", .var.name = "...")
    checkmate::reportAssertions(ac)

    addPath <- function(path, type)
    {
        if (is.null(path))
            return(NULL)
        files <- listMSFiles(path, getMSFileFormats(type))
        if (length(files) == 0)
        {
            warning(sprintf("No analyses found in %s!", path), call. = FALSE)
            return(NULL)
        }
        ret <- data.table(analysis = simplifyAnalysisNames(files))
        ret[, (paste0("path_", type)) := dirname(files)]
        return(ret)
    }
    
    # UNDONE: use getMSFileTypes() somehow?
    allPaths <- Map(list(fromRaw, fromCentroid, fromIMS, fromProfile), c("raw", "centroid", "profile", "ims"), f = addPath)
    allPaths <- pruneList(allPaths)
    anaInfo <- if (length(allPaths) > 0)
        Reduce(allPaths, f = function(x, y) merge(x, y, by = "analysis", all = TRUE))
    else
        data.table(analysis = character())
    
    for (ft in getMSFileTypes())
    {
        if (is.null(anaInfo[[paste0("path_", ft)]]))
            anaInfo[, (paste0("path_", ft)) := ""]
    }
    
    convPaths <- setNames(list(convCentroid, convProfile, convIMS), c("centroid", "profile", "ims"))
    for (ft in getMSFileTypes())
    {
        if (!is.null(convPaths[[ft]]))
        {
            col <- paste0("path_", ft)
            anaInfo[is.na(get(col)) | !nzchar(get(col)), (col) := rep(convPaths[[ft]], length.out = .N)]
        }
    }
    
    for (oa in names(otherArgs))
        anaInfo[, (oa) := rep(otherArgs[[oa]], length.out = nrow(anaInfo))]
    
    if (is.null(anaInfo[["replicate"]]))
        anaInfo[, replicate := analysis]
    else
        anaInfo[!nzchar(replicate), replicate := analysis]
    if (is.null(anaInfo[["blank"]]))
        anaInfo[, blank := ""]
    
    return(as.data.frame(anaInfo))
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

#' Obtains extracted ion chromatograms (EICs)
#'
#' This function generates one or more EIC(s) for given retention time, \emph{m/z} and optionally mobility ranges.
#'
#' @template analysisInfo-arg
#'
#' @param ranges A \code{list} with for each analysis a \code{data.frame} with \code{numeric} columns \code{"retmin"},
#'   \code{"retmax"}, \code{"mzmin"}, \code{"mzmax"} with the lower/upper ranges of the retention time and \emph{m/z}.
#'   Furthermore, columns \code{"mobmin"} and \code{"mobmax"} can be added for mobility lower/upper ranges in IMS data.
#' @param output Should be \code{"fill"}, \code{"pad"} or \code{"raw"}. Internally, EIC data is compressed by omitting
#'   any zero intensity data points. If \code{output="fill"} then the zero intensity points are re-added to obtain
#'   continuous chromatograms. If \code{output="pad"} then zero intensity points are only re-added that surround others,
#'   which is sufficient for \emph{e.g.} plotting. If \code{output="raw"} then the original compressed data is returned.
#'
#' @return A \code{list} with for each analysis a \code{list} with EIC data for each of the rows in \code{ranges}.
#' 
#'   If \code{output="raw"} then additional columns with \emph{e.g.} mean-averaged and base peak \emph{m/z} values for
#'   each data point are returned. Furthermore, the \code{allXValues} attribute is set that can be used to obtain the
#'   original retention time values to reconstruct the original complete chromatogram.
#'
#' @templateVar what \code{getEICs}
#' @template uses-msdata
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
    ret <- doGetEICs(analysisInfo, ranges, mode = "full", minIntensityIMS = minIntensityIMS, pad = output == "pad")
    if (output == "fill")
    {
        ret <- lapply(ret, function(anaEICs)
        {
            at <- attr(anaEICs, "allXValues")
            if (is.null(at))
                return(anaEICs) # no EICs
        
            return(lapply(anaEICs, function(eic) data.frame(time = at,
                                                            intensity = doFillEIXIntensities(at, eic$time, eic$intensity))))
        })
    }
    return(ret)
}

#' Background MS/MS peak detection
#'
#' Detects background MS/MS peaks by gathering frequently occurring peaks in MS/MS spectra from blanks.
#'
#' This function iterates through all MS/MS spectra of the given analyses and collects the most frequently occurring
#' peaks. It first averages all spectra within the same analyses, and retains those peaks above an abundance threshold
#' (set by \code{avgSpectraParams}). The analyses-averaged spectra are then also averaged, and again peaks with a
#' minimum abundance are kept (set by \code{avgAnalysesParams}).
#'
#' The frequent occurrence of a peak throughout MS/MS spectra, including DDA spectra different isolation m/z values, is
#' often a sign of contamination from the analytical system (\emph{e.g.} from the mobile phase, MS quadrupoles etc).
#' Hence, the analyses used for blank subtraction are typically just measurements of \emph{e.g.} ultrapure water or
#' another solvent, and not need to be treated by any extraction procedure.
#'
#' The output of this function can directly be used to the \code{\link[=filter,MSPeakLists-method]{filter method for
#' MSPeakLists}} to subsequently remove these background peaks, which possibly improves subsequent feature annotation.
#'
#' @param anaInfo The \link[=analysis-information]{analysis info} object with the analyses to be used for background
#'   peak detection.
#' @param replicates A \code{character} with names of replicates for the analyses used for background peak detection. If
#'   \code{NULL} then all the analyses defined in \code{anaInfo} will be used.
#' @param MSLevel The MS level of the spectra to be used for background peak detection. This should be \samp{1} or
#'   \samp{2}. This function is only tested with \code{MSLevel=2}. And while \code{MSLevel=1} is possible, it may be
#'   rather computationally intensive due to the relatively large number of MS peaks to iterate through.
#' @param retentionRange,mobilityRange A two-sized \code{numeric} vector with the minimum and maximum retention time and
#'   ion mobility to consider, respectively. If \code{NULL} then the full range is used.
#' @param minBPIntensity The minimum basepeak intensity of a spectrum to be considered for background peak detection.
#'   This is primarily intended to optimize the detection procedure.
#' @param avgSpectraParams,avgAnalysesParams A \code{list} with parameters used for averaging the MS spectra within and
#'   between analyses, respectively. See \code{\link{getDefAvgPListParams}} for more details.
#'
#' @return A \code{\link{data.table}} with the detected background peaks and abundance statistics. The table can
#'   directly be passed to the \code{removeMZs} argument of the \code{\link[=filter,MSPeakLists-method]{filter method for
#' MSPeakLists}}.
#'
#' @references \insertRef{Helmus2025}{patRoon}
#'
#' @export
getBGMSMSPeaks <- function(anaInfo, replicates = NULL, MSLevel = 2, retentionRange = NULL, mobilityRange = NULL,
                           minBPIntensity = 5000,
                           avgSpectraParams = getDefAvgPListParams(minAbundanceRel = 0.1, topMost = 25),
                           avgAnalysesParams = getDefAvgPListParams(minAbundanceRel = 0.8, topMost = 25))
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
    setnames(ret, c("abundance_rel", "abundance_abs", "abundance_prev_rel", "abundance_prev_abs"),
             c("abundance_rel_ana", "abundance_abs_ana", "abundance_rel_spec", "abundance_abs_spec"))
    printf("Done!\n")
    
    return(ret[])
}

#' Extracted Ion Chromatogram and Mobilogram parameters
#'
#' Parameters for creation of extracted ion chromatograms and mobilograms.
#'
#' The following parameters exist to configure the creation of extracted ion chromatograms (EICs) and extracted ion
#' mobilograms (EIMs):
#'
#' \itemize{
#'
#' \item \code{window} A value that is subtracted or added to the minimum and maximum retention time (EICs) or mobility
#' (EIMs) of the feature. Thus, setting this value to \samp{>0} will 'zoom out' on the x-axis of a chromatogram or
#' mobilogram. Defaults to \code{defaultLim("retention", "wide")} (EICs) or \code{defaultLim("mobility", "wide")} (EIMs)
#' (see \link{limits}).
#'
#' \item \code{topMost} Only create EICs/EIMs for this number of top most intense features. If \code{NULL} then
#' EICs/EIMs are created for all features.
#'
#' \item \code{topMostByReplicate} If set to \code{TRUE} and \code{topMost} is set: only create EICs/EIMs for the top
#' most features in each replicate. For instance, when \code{topMost=1} and \code{topMostByReplicate=TRUE}, then only
#' the most intense feature of each replicate is considered.
#'
#' \item \code{onlyPresent} If \code{TRUE} then EICs/EIMs are created only for analyses in which a feature was detected,
#' if \code{onlyPresent=FALSE} then data is generated for \strong{all} analyses. The latter is handy to evaluate if a
#' peak was 'missed' during peak detection or removed during \emph{e.g.} filtering.
#'
#' \item \code{mzExpIMSWindow} Specifically for IMS data: additional \emph{m/z} tolerance on top of the feature limits.
#' This is for IMS workflows where features were detected from centroided LC-MS like data, while EICs/EIMs are generated
#' from raw IMS data. In this case the feature \emph{m/z} limits were derived from centroided data, which typically has
#' smaller \emph{m/z} deviations across scans compared to IMS data. The \code{mzExpIMSWindow} parameter sets an
#' additional \emph{m/z} tolerance to specifically handle this case. Defaults to \code{defaultLim("mz", "default")}
#' (see \link{limits}).
#'
#' \item \code{minIntensityIMS} Raw intensity threshold for IMS data. This is primarily intended to speed up raw data
#' processing.
#'
#' }
#'
#' if \code{onlyPresent=FALSE} then the following parameters are also relevant: \itemize{
#'
#' \item \code{mzExpWindow},\code{mobExpWindow} To create EICs or EIMs for analyses in which no feature was found, the
#' \emph{m/z} or mobility value is derived from the min/max values of all features in the feature group. The value of
#' \code{mzExpWindow} and \code{mobExpWindow} further expands this window to allow a greater tolerance. Defaults to
#' \code{defaultLim("mz", "very_narrow")} and \code{defaultLim("mobility", "very_narrow")} (see \link{limits}).
#'
#' \item \code{setsAdductPos},\code{setsAdductNeg} \setsWF In sets workflows the adduct must be known to calculate the
#' ionized \emph{m/z}. If a feature is completely absent in a particular set then it follows no adduct annotations are
#' available and the value of \code{setsAdductPos} (positive ionization data) or \code{setsAdductNeg} (negative
#' ionization data) will be used instead.
#'
#' }
#'
#' The following additional parameters exist specifically for EIMs (\code{EIMParams}): \itemize{
#'
#'   \item \code{maxRTWindow} Maximum retention time window (seconds, +/- feature retention time) in which mobilograms
#'   are collected and averaged. Defaults to \code{defaultLim("retention", "very_narrow")} (see \link{limits}).
#'
#'   \item \code{IMSWindow},\code{clusterMethod}: IMS clustering parameters. \code{IMSWindow} defaults to
#'   \code{defaultLim("mobility", "medium")} (see \link{limits}).
#'
#' }
#'
#' These parameters are passed as a named \code{list} as the \code{EICParams} or \code{EIMParams} argument to functions
#' that work with EIC or EIM data. The \code{getDefEICParams} and \code{getDefEIMParams} functions generate such
#' parameter list with defaults.
#'
#' @param \dots optional named arguments that override defaults.
#'
#' @name EIXParams
#' @export
getDefEICParams <- function(...)
{
    def <- getDefEIXParams()
    def$window <- defaultLim("retention", "wide")
    
    mod <- list(...)
    if (!is.null(mod[["rtWindow"]]))
    {
        warning("The \"rtWindow\" parameter is deprecated and was renamed to \"window\"", call. = FALSE)
        mod$window <- mod$rtWindow
        mod$rtWindow <- NULL
    }
    
    return(modifyList(def, mod, keep.null = TRUE))
}

#' @rdname EIXParams
#' @export
getDefEIMParams <- function(...)
{
    def <- getDefEIXParams()
    def <- modifyList(def, list(
        window = defaultLim("mobility", "wide"),
        maxRTWindow = defaultLim("retention", "very_narrow"),
        IMSWindow = defaultLim("mobility", "medium"),
        clusterMethod = "distance"
    ))
    return(modifyList(def, list(...), keep.null = TRUE))
}

#' Parameters to handle TP data with structural information
#'
#' Parameters used by \code{\link{generateTPs}} with algorithms that use structural information.
#'
#' The following parameters can be configured: \itemize{
#'
#'   \item \code{calcLogP} A \code{character} specifying whether \verb{log P} values should be calculated with
#'   \code{\link[rcdk:get.xlogp]{rcdk::get.xlogp}} (\code{calcLogP="rcdk"}),
#'   \href{https://github.com/openbabel/openbabel}{OpenBabel} (\code{calcLogP="obabel"}) or not at all
#'   (\code{calcLogP="none"}). The \code{log P} are values of parents and TPs are used for \link[=retDir]{retention
#'   order calculation}.
#'
#'   \item \code{forceCalcLogP} Force calculation of \verb{Log P} values, even if already provided by the TP generation
#'   algorithm. This is primarily useful to obtain log P values that were consistently calculated with the same
#'   algorithm, as some algorithms may only partially output these values (\emph{e.g.} not for parents).
#'
#'   \item \code{forceCalcRetDir} Force calculation of \link[=retDir]{retention order directions}, even if
#'   already provided by the TP generation algorithm. This is primarily intended for re-calculation of library TP data,
#'   which may have been calculated with different log P values.
#'
#'   \item \code{minLogPDiff} The minimum difference in \verb{log P} values between a parent and its TPs to be
#'   considered eluting differently. This is used for \link[=retDir]{retention order calculation}.
#'
#'   \item \code{calcSims} If set to \code{TRUE} then structural similarities between the parent and its TPs are
#'   calculated. The \link[=filter,transformationProductsStructure-method]{filter method} can be used to threshold
#'   structural similarities. This may be useful under the assumption that parents and TPs who have a high structural
#'   similarity, also likely have a high MS/MS spectral similarity (which can be evaluated after componentization with
#'   \code{\link{generateComponentsTPs}}).
#'
#'   \item \code{fpType} The type of structural fingerprint that should be calculated. See the \code{type} argument of
#'   the \code{\link{get.fingerprint}} function of \CRANpkg{rcdk}.
#'
#'   \item \code{fpSimMethod} The method for calculating similarities (i.e. not dissimilarity!). See the \code{method}
#'   argument of the \code{\link{fp.sim.matrix}} function of the \CRANpkg{fingerprint} package.
#'
#' }
#'
#' These parameters are passed as a named \code{list} as the \code{TPStructParams} argument to functions.
#'
#' The \code{getDefTPStructParams} function generates such parameter list with defaults.
#'
#' @param \dots optional named arguments that override defaults.
#'
#' @references \addCitations{rcdk}{1} \cr\cr \insertRef{OBoyle2011}{patRoon}
#'
#' @export
getDefTPStructParams <- function(...)
{
    def <- list(
        calcLogP = "rcdk",
        forceCalcLogP = FALSE,
        forceCalcRetDir = FALSE,
        minLogPDiff = 1,
        calcSims = FALSE,
        fpType = "extended",
        fpSimMethod = "tanimoto"
    )
    
    return(modifyList(def, list(...)))
}

#' Peak detection parameters
#'
#' Algorithms and parameters for automatic detection of peaks in chromatograms and mobilograms.
#'
#' The algorithm and its parameters for peak detection should be in a named \code{list} with the format:
#'
#' \code{list(algorithm = <algorithm>, param1 = ..., param2 = ..., ...)}
#'
#' Where \code{<algorithm>} is the name of the algorithm and \code{param1}, \code{param2} etc are the parameters. The
#' \code{getDefPeakParams} function generates such parameter list with the algorithm and default parameters.
#'
#' The following algorithms are currently supported: \itemize{
#'
#'   \item \code{"openms"}: uses
#'   \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MRMTransitionGroupPicker.html}{MRMTransitionGroupPicker}
#'   tool from \href{http://www.openms.de}{OpenMS}.
#'
#'   \item \code{"xcms3"}: uses the \code{\link[xcms:peaksWithCentWave]{xcms::peaksWithCentWave}} function.
#'
#'   \item \code{"envipick"}: uses the \code{\link[enviPick:mzpick]{enviPick::mzpick}} function.
#'
#'   \item \code{"piek"}: uses an optimized peak detection algorithm derived from \insertCite{Dietrich2021}{patRoon}.
#' }
#'
#' The parameters are discussed in the next sections.
#'
#' @param type The type of parameter defaults: \code{"chrom"} for chromatograms and \code{"bruker_ims"} and
#'   \code{"agilent_ims"} for mobilograms coming from Bruker and Agilent systems, respectively.
#' @param algorithm The peak detection algorithm: \code{"openms"}, \code{"xcms3"}, \code{"envipick"} or \code{"piek"}.
#' @param \dots optional named arguments that override defaults.
#'
#' @section General parameters: These parameters are applicable to all algorithms \itemize{
#'
#'   \item \code{forcePeakRange} a two-sized \code{numeric} vector with the minimum and maximum width for a peak. Peaks
#'   that are more narrow or wide will be clamped to this range. This is especially useful for algorithms that consider
#'   an extensive part of the fronting/tailing noise as part as the peak. Set to \code{c(0, 0)} to disable.
#'
#'   \item \code{relMinIntensity} the minimum intensity threshold for a peak relative to the highest peak in the same
#'   chromatogram/mobilogram. This is \emph{e.g.} useful to exclude noise in mobilograms where normally few peaks are
#'   expected.
#'
#'   }
#'
#' @section Parameters for \code{openms}: The parameters directly map to the command line options for
#'   \code{MRMTransitionGroupPicker}, please see
#'   \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MRMTransitionGroupPicker.html}{its
#'   documentation}. \itemize{
#'
#'     \item \code{minPeakWidth} the minimum peak width, sets the \command{min_peak_width} option.
#'
#'     \item \code{backgroundSubtraction} the background subtraction method, sets the
#'     \command{-algorithm:background_subtraction} option.
#'
#'     \item \code{SGolayFrameLength} the frame length for Savitzky-Golay smoothing, sets the
#'     \command{-algorithm:PeakPickerMRM:sgolay_frame_length} option.
#'
#'     \item \code{SGolayPolyOrder} order of the polynomial, sets the
#'     \command{-algorithm:PeakPickerMRM:sgolay_polynomial_order} option.
#'
#'     \item \code{useGauss} set to \code{TRUE} to use Gaussian smoothing (instead of Savitzky-Golay, sets the
#'     \command{-algorithm:PeakPickerMRM:use_gauss} option.
#'
#'     \item \code{gauss_width} the Gaussian width, estimated peak size, sets the
#'     \command{-algorithm:PeakPickerMRM:gauss_width} option.
#'
#'     \item \code{SN} signal to noise threshold, sets the \command{-algorithm:PeakPickerMRM:signal_to_noise} option.
#'
#'     \item \code{SNWinLen} \code{SN} window length, sets the \command{-algorithm:PeakPickerMRM:sn_win_len} option.
#'
#'     \item \code{SNBinCount} \code{SN} bin count, sets the \command{-algorithm:PeakPickerMRM:sn_bin_count} option.
#'
#'     \item \code{method} peak picking method, sets the \command{-algorithm:PeakPickerMRM:method} option.
#'
#'     \item \code{integrationType} the integration technique, sets the
#'     \command{-algorithm:PeakIntegrator:integration_type} option.
#'
#'     \item \code{baselineType} the baseline type, sets the \command{-algorithm:PeakIntegrator:baseline_type} option.
#'
#'     \item \code{fitEMG} if \code{TRUE} then the EMG model is used for fitting, sets the
#'     \command{-algorithm:PeakIntegrator:fit_EMG} option.
#'
#'   }
#'
#' @section Parameters for \code{xcms3} and \code{envipick}: See the documentation for
#'   \code{\link[xcms:peaksWithCentWave]{xcms::peaksWithCentWave}} and \code{\link[enviPick:mzpick]{enviPick::mzpick}}
#'   for \code{xcms3} and \code{envipick}, respectively.
#'
#' @section Parameters for \code{piek}: \itemize{
#'
#'   \item \code{minIntensity} the minimum intensity of a peak.
#'
#'   \item \code{SN} the signal to noise ratio.
#'
#'   \item \code{peakWidth} two-sized \code{vector} with the minimum and maximum peak width (seconds)
#'
#'   \item \code{RTRange} two-sized \code{vector} with the minimum and maximum retention time range (seconds). Set the
#'   2nd element to \code{Inf} for no upper limit.
#'
#'   \item \code{maxPeaksPerSignal} upper threshold for consecutive maxima of similar size to be regarded as noise.
#'
#'   }
#'
#' @note The peak detection used by \code{algorithm="openms"} is different than that of
#'   \code{\link{findFeaturesOpenMS}}.
#'
#'   The \code{patRoon.threads} package option sets the number of threads for the \code{piek} algorithm.
#'
#' @template refs-openms
#' @references \insertAllCited{} \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
#'
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

#' Parameters for CCS calculation
#'
#' Configuration and description of parameters used for CCS calculation.
#'
#' The following parameters exist to configure the CCS calculation: \itemize{
#'
#'   \item \code{method} The CCS calculation method. Should be \code{"bruker"}, \code{"mason-schamp_k"},
#'   \code{"mason-schamp_1/k"} or \code{"agilent"}. See details below.
#'
#'   \item \code{defaultCharge} The default charge of the ions. This is used when no charge information is available.
#'
#'   \item \code{temperature},\code{massGas} The temperature (Kelvin) and exact mass of the drift gas. See calculation
#'   details below.
#'
#'   \item \code{MasonSchampConstant} The Mason-Schamp constant. See calculation details below.
#'
#'   \item \code{calibrant} If \code{method="agilent"}: the calibrant data to be used for CCS calculation. This should
#'   either be \itemize{
#'
#'     \item A path to an Agilent \file{.d} file.
#'
#'     \item A path to an \file{OverrideImsCal.xml} file (found in \file{sample.d/AcqData}).
#'
#'     \item A named \code{list} with the elements \code{massGas}, \code{TFix} and \code{beta}.
#'
#'    }
#'
#' }
#'
#' The CCS calculation depends on the \code{method} parameter: \itemize{
#'
#'   \item \code{bruker}: uses the Bruker \command{TDF-SDK} for calculations. See \link{msdata} for configuration
#'   options. Only applicable to TIMS data.
#'
#'   \item \code{mason-schamp_k}: uses the Mason-Schamp equation:
#'     \deqn{CCS = C \cdot \frac{charge}{\sqrt{u \cdot T}} \cdot \frac{1}{mobility}}
#'
#'     With \itemize{
#'
#'       \item \emph{C} the Mason-Schamp constant, can be changed by setting the \code{MasonSchampConstant} parameter.
#'       See \insertCite{George2024}{patRoon} for details.
#'
#'       \item \emph{u} the reduced mass of the drift gas and the ion:
#'       \deqn{u = \frac{m_{gas} \cdot m_{ion}}{m_{gas} + m_{ion}}}
#'       The mass of the drift gas is defined by the \code{massGas} parameter.
#'
#'       \item \emph{T} the temperature (Kelvin) as defined by the \code{temperature} parameter.
#'
#'     }
#'
#'   \item \code{mason-schamp_1/k}: as \code{mason-schamp_k} but assuming an inversed mobility (\eqn{\frac{1}{k}}).
#'   This is meant for TIMS data. Compared to \code{method="bruker"}, this doesn't rely on the \command{TDF-SDK} but may
#'   produce results with very minor differences \insertCite{George2024}{patRoon}.
#'
#'   \item \code{agilent}: uses Agilent calibration data with the following equation:
#'
#'     \deqn{CCS = (mobility - t_{fix}) \cdot \frac{charge}{\beta} \cdot \frac{1}{\sqrt{\frac{m_{ion}}{m_{ion} + m_{gas}}}}}
#'
#'     With \eqn{t_{fix}} and \eqn{\beta} the \code{TFix} and \code{beta} values from the calibration data. The
#'     \code{massGas} parameter sets the \eqn{m_{gas}} value.
#'
#' }
#'
#' The \code{getCCSParams} function generates such parameter list with defaults.
#'
#' @param method,calibrant Sets the CCS calculation method and the calibrant data (only required if
#'   \code{method="agilent"}). See details.
#' @param \dots optional named arguments that override defaults.
#'
#' @source The calculation formulas was derived from \insertCite{Haler2017}{patRoon}, \insertCite{George2024}{patRoon}
#'   and the implementation used by \href{https://systemsomicslab.github.io/compms/msdial/main.html}{MS-DIAL}
#'   \insertCite{Tsugawa2020}{patRoon}. (\code{MobilityToCrossSection} method from the \code{IonMobilityUtility} class).
#'
#' @references \insertAllCited{}
#'
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

#' Parameters to specify a IMS data range
#'
#' Parameters that define a range of mobility or \acronym{CCS} values.
#'
#' The following parameters are used to define an IMS range: \itemize{
#'
#'   \item \code{param} Should be \code{"mobility"} or \code{"CCS"} to specify a mobility or \acronym{CCS} range,
#'   respectively.
#'   
#'   \item \code{lower},\code{upper} The lower and upper range.
#'   
#'   \item \code{mzRelative} Set to \code{TRUE} to specify an IMS range that is normalized by \emph{m/z}.
#'
#' }
#'
#' These parameters are passed as a named \code{list} as the \code{IMSRangeParams} argument to functions.
#'
#' The \code{getIMSRangeParams} function generates such parameter list with defaults.
#'
#' @param param,lower,upper,mzRelative Arguments to specify the IMS range parameters, see Details.
#' 
#' @export
getIMSRangeParams <- function(param, lower, upper, mzRelative = FALSE)
{
    return(list(param = param, lower = lower, upper = upper, mzRelative = mzRelative))
}

#' Parameters for IMS matching
#'
#' Parameters that define how mobility or \acronym{CCS} values between \emph{e.g.} features and suspects should be
#' matched.
#'
#' The following parameters should be defined: \itemize{
#'
#'   \item \code{param} Should be \code{"mobility"} or \code{"CCS"} to match by mobility or \acronym{CCS}, respectively.
#'   
#'   \item \code{window},\code{relative} The \code{window} parameter sets the tolerance window size used for matching.
#'   If \code{relative=TRUE} then the tolerance is relative (\samp{0-1}). The defaults for \code{window} are
#'   (see \link{limits}): \itemize{
#'
#'     \item \code{defaultLim("mobility", "medium")} (\code{param="mobility"} and \code{relative=FALSE})
#'
#'     \item \code{defaultLim("mobility", "medium_rel")} (\code{param="mobility"} and \code{relative=TRUE})
#'
#'     \item \code{defaultLim("CCS", "medium")} (\code{param="CCS"} and \code{relative=FALSE})
#'
#'     \item \code{defaultLim("CCS", "medium_rel")} (\code{param="CCS"} and \code{relative=TRUE})
#'
#'   }
#'   
#'   \item \code{minMatches} The minimum number of mobility/\acronym{CCS} matches for a suspect hit. If the number of
#'   available reference mobility/\acronym{CCS} values in the suspect list is less than \code{minMatches}, then that
#'   number is used as threshold. Set to \code{0} to disable.
#'
#' }
#'
#' These parameters are passed as a named \code{list} as the \code{IMSMatchParams} argument to functions.
#'
#' The \code{getIMSMatchParams} function generates such parameter list with defaults.
#'
#' @param param Should be \code{"mobility"} or \code{"CCS"} to match by mobility or \acronym{CCS}, respectively.
#' @param \dots optional named arguments that override defaults.
#'
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
        if (ret$param == "mobility")
            ret$window <- if (ret$relative) defaultLim("mobility", "medium_rel") else defaultLim("mobility", "medium")
        else # CCS
            ret$window <- if (ret$relative) defaultLim("CCS", "medium_rel") else defaultLim("CCS", "medium")
    }
    
    return(ret)
}

#' @templateVar consider for generating calibrant data
#' @templateVar append Best to keep this \code{"maybe"}, as calibration typically doesn't support IMS filtered data.
#' @template IMS-arg
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
        unlist(lapply(concReps, function(r) anaInfo[replicate == r]$analysis))
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

#' Conversion between mobility and CCS
#'
#' Utility functions to convert between mobility and \acronym{CCS} data.
#'
#' @param mobility,ccs A \code{numeric} vector with mobility or \acronym{CCS} values that should be converted.
#' @param mz A \code{numeric} vector with the \emph{m/z} values that map to the input mobility or \acronym{CCS} values.
#' @param charge A \code{numeric} vector with the ion charges that map to the input mobility or \acronym{CCS} data. Will
#'   be recycled if necessary. If \code{NULL} then the charge configured in \code{CCSParams} is used.
#' 
#' @template CCSParams-arg
#'
#' @name CCS-Conversion
#' @export
convertMobilityToCCS <- function(mobility, mz, CCSParams, charge = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertMobilityConversionArgs(mobility, mz, CCSParams, charge, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(charge))
        charge <- CCSParams$defaultCharge
    charge <- rep(abs(charge), length.out = length(mobility))
    
    m <- mz * charge
    
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
        u <- ((m * CCSParams$massGas) / (m + CCSParams$massGas))
        return(CCSParams$MasonSchampConst * (charge / (sqrt(u * CCSParams$temperature))) * mobility)
    }
    else if (CCSParams$method == "agilent")
    {
        calibrant <- prepareAgilentIMSCalib(CCSParams$calibrant, CCSParams$massGas)
        return((mobility - calibrant$TFix) * charge / calibrant$beta / sqrt(m / (m + calibrant$massGas)))
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

#' @rdname CCS-Conversion
#' @export
convertCCSToMobility <- function(ccs, mz, CCSParams, charge = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertMobilityConversionArgs(ccs, mz, CCSParams, charge, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(charge))
        charge <- CCSParams$defaultCharge
    charge <- rep(abs(charge), length.out = length(ccs))
    
    m <- mz * charge
    
    if (CCSParams$method == "bruker")
    {
        if (!doInitBrukerLib())
            stop("Cannot convert mobility data with Bruker TIMS library: failed to initialize library", call. = FALSE)
        return(getBrukerMob(ccs, charge, mz))
    }
    else if (CCSParams$method %in% c("mason-schamp_k", "mason-schamp_1/k"))
    {
        u <- ((m * CCSParams$massGas) / (m + CCSParams$massGas))
        mob <- ccs / CCSParams$MasonSchampConst / (charge / (sqrt(u * CCSParams$temperature)))
        if (CCSParams$method == "mason-schamp_k")
            mob <- 1 / mob
        return(mob)
    }
    else if (CCSParams$method == "agilent")
    {
        calibrant <- prepareAgilentIMSCalib(CCSParams$calibrant, CCSParams$massGas)
        return((ccs * sqrt(m / (m + calibrant$massGas)) * calibrant$beta / charge) + calibrant$TFix)
    }
    # UNDONE: support waters later, procedure seems unclear how it affects different instruments and need test data
    # else # waters
    # {
    # }
}

#' Assign IMS data to suspects
#'
#' Adds calculated mobility and/or \acronym{CCS} data to a suspect list.
#'
#' The \code{assignMobilities} method for suspect lists is used to (1) add IMS data to suspects from predictions or
#' library data and (2) convert (previously added) mobility <--> \acronym{CCS} values. These steps are controlled by the
#' \code{from} and \code{CCSParams} arguments, respectively.
#'
#' Mobility and \acronym{CCS} values assigned in the suspect list are either adduct specific or not. Adduct specific
#' values are preferred, as the 'correct' value can be automatically selected during suspect screening based on the
#' adduct assigned to the feature (or passed as the \code{adduct} argument to \code{\link{screenSuspects}}). The
#' non-adduct specific values are typically used when the corresponding adduct for the mobility/\acronym{CCS} value is
#' unknown (or not of interest). These values get precedence over adduct specific values. The adduct specific values are
#' stored in \code{mobility_<adduct>} and \code{CCS_<adduct>} columns, where \code{<adduct>} is the adduct name
#' (\emph{e.g.} \code{[M+H]+}, \code{[M-H]-}). The \code{mobility} and \acronym{CCS} columns store any non-adduct
#' specific values. The \code{adducts} argument ultimately defines the use of adduct and non-adduct specific values.
#'
#' The mobility <--> \acronym{CCS} conversions occur both ways, \emph{i.e.} missing \acronym{CCS} values will be
#' converted from mobility values and \emph{vice versa}. If adduct specific values are converted then the charge value
#' used for these calculations is taken from the corresponding adduct. For non-adduct specific values the charge is
#' taken from the adduct specified in suspect list if present, or from the default charge specified in \code{CCSParams}
#' otherwise.
#'
#' @param obj The suspect list to which the mobility and/or CCS data should be added. Should be a \code{data.frame} or
#'   \code{data.table}.
#' @param from Specifies from where IMS data is added to the suspect list. This can be the following: \itemize{
#'
#'   \item \code{"pubchemlite"}: \acronym{CCS} data is matched from predicted values of the
#'   \href{https://pubchemlite.lcsb.uni.lu/}{PubChemLite} database. \strong{Note}: this requires a local copy of the
#'   \href{https://zenodo.org/records/15311000}{CCS amended PubChemLite database} (see the Handbook for more details),
#'   which is automatically installed by \pkg{patRoonExt}.
#'
#'   \item \code{"c3sdb"}: Uses the \href{https://github.com/dylanhross/c3sdb}{C3SDB} \command{Python} package to
#'   predict \acronym{CCS} values. This requires a local installation of \command{C3SDB}, \emph{e.g.} performed with
#'   \code{\link{installC3SDB}}.
#'
#'   \item A \code{data.table} or \code{data.frame} to which IMS data is matched. Should contain the column defined by
#'   \code{matchFromBy} and columns storing (non-)adduct specific mobility/\acronym{CCS} columns (see Details).
#'
#'   \item \code{NULL}: No IMS data is added to the suspect list.
#'
#'   }
#'
#'   Any \code{NA} values in \code{from} are ignored.
#' @param matchFromBy Which column should be used to match the IMS data from \code{from} and suspects. Valid options are
#'   \code{"InChIKey"}, \code{"InChIKey1"} (first block InChIKey), \code{"InChI"}, \code{"SMILES"}, \code{"name"}.
#'   However, this also depends on which columns are available in either of the data sources. \code{InChIKey1} values
#'   are automatically calculated from \code{InChIKey}s, if possible.
#' @param overwrite Set to \code{TRUE} to overwrite any existing suspect IMS data with data from \code{from}.
#' @param adducts A \code{character} with the adduct(s) to consider for assigning mobility data to suspects and mobility
#'   <--> \acronym{CCS} conversions. This may be limited by what is available in the data source specified by
#'   \code{from} (see \code{\link{C3SDBAdducts}} for \code{from="c3sdb"}). Inclusion of \code{NA} in \code{adducts}
#'   allows the use of non-adduct specific values (see Details).
#'
#'   The value for \code{adducts} is automatically expanded by the adducts specified in the \code{adduct} column of the
#'   suspect list. Hence, \code{adducts} can be empty (\code{character()}) if no calculations for other adducts are
#'   desired.
#' @param predictAdductOnly If \code{from="c3sdb"} and \code{predictAdductOnly=TRUE} then only predictions are performed
#'   for the adduct specified in the \code{adduct} column in the suspect list (if present).
#' @param prepareChemProps Set to \code{TRUE} to perform chemical property calculation and validation on the suspect
#'   list (described below).
#' @param virtualenv The virtual \command{Python} environment in which \command{C3SDB} is installed. This is passed to
#'   \code{\link[reticulate:use_virtualenv]{reticulate::use_virtualenv}}. Set to \code{NULL} to skip this and not setup
#'   the environment.
#'
#' @template CCSParams-arg
#'
#' @templateVar whatCP suspect data (if \code{prepareChemProps=TRUE}) and \code{from} data (if a table)
#' @template chemPropCalc
#'
#' @references \insertRef{Schymanski2021}{patRoon} \cr\cr \insertRef{Elapavalore2025}{patRoon} \cr\cr
#'   \insertRef{Ross2020}{patRoon}
#'
#' @name assignMobilities_susp
#' @aliases assignMobilities,data.table-method
#' @export
setMethod("assignMobilities", "data.table", function(obj, from = NULL, matchFromBy = "InChIKey",
                                                     overwrite = FALSE, adducts = c("[M+H]+", "[M-H]-", NA),
                                                     predictAdductOnly = TRUE, CCSParams = NULL,
                                                     prepareChemProps = TRUE, prefCalcChemProps = TRUE,
                                                     neutralChemProps = FALSE, virtualenv = "patRoon-C3SDB")
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
    checkmate::assertCharacter(adducts, min.chars = 1, any.missing = TRUE, add = ac)
    assertCCSParams(CCSParams, null.ok = TRUE, add = ac)
    checkmate::assertString(virtualenv, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(obj[["adduct"]]))
        adducts <- union(adducts, obj[!is.na(adduct) & nzchar(adduct)]$adduct) # expand with adducts from suspect list

    hash <- makeHash(obj, from, matchFromBy, overwrite, adducts, predictAdductOnly, CCSParams, prepareChemProps,
                     prefCalcChemProps, neutralChemProps)
    cd <- loadCacheData("assignMobilitiesDT", hash)
    if (!is.null(cd))
        return(cd)
    
    do_C3SDB <- identical(from, "c3sdb")
    if (do_C3SDB)
        checkPackage("reticulate")
    
    obj <- copy(obj)
    
    if (prepareChemProps)
        obj <- prepareChemTable(obj, prefCalcChemProps, neutralChemProps)
    
    adductsNoNone <- na.omit(adducts)
    needMZs <- do_C3SDB || !is.null(CCSParams)
    mzTab <- NULL
    if (needMZs)
    {
        mzTab <- copy(subsetDTColumnsIfPresent(obj, c("name", "neutralMass", "mz")))
        
        if (length(adductsNoNone) > 0)
        {
            msg <- paste("Mass data is necessary for C3SDB pedictions and mobility <--> CCS conversions.",
                         "This data can be automatically calculated by setting prepareChemProps=TRUE and ensuring a SMILES,",
                         "InChI or formula column is present.")
            collAdd <- paste0(adductsNoNone, collapse = ", ")
            if (is.null(mzTab[["neutralMass"]]))
            {
                warning(sprintf("Calculations for adducts (%s) will be skipped as there is no neutralMass data.",
                                collAdd), msg, call. = FALSE)
                adductsNoNone <- character()
                adducts <- intersect(adducts, "none")
            }
            else if (anyNA(mzTab$neutralMass))
                warning(sprintf("Calculations for adducts (%s) will be skipped for these rows with NA neutralMass values: %s. ",
                                collAdd, paste0(mzTab[is.na(neutralMass), which = TRUE], collapse = ", ")),
                        msg, call. = FALSE)
        }
        if (anyNA(adducts) && !is.null(CCSParams))
        {
            if (is.null(mzTab[["mz"]]))
                warning("Mobility <--> CCS conversions for non-adduct specific values will be skipped as there is no mz data.",
                        call. = FALSE)
            else if (anyNA(mzTab$mz))
                warning(sprintf("Mobility <--> CCS conversions for non-adduct specific values will be skipped for the following rows with NA mz values: %s.",
                        paste0(mzTab[is.na(mz), which = TRUE], collapse = ", ")), call. = FALSE)
        }
        
        if (!is.null(mzTab[["neutralMass"]]))
        {
            for (add in adductsNoNone)
            {
                col <- paste0("mz_", add)
                mzTab[, (col) := calculateMasses(neutralMass,
                                                 checkAndToAdduct(add, .var.name = sprintf("adducts ('%s')", add)), "mz")]
            }
        }
    }
    
    if (do_C3SDB)
    {
        if (is.null(obj[["neutralMass"]]) || is.null(obj[["SMILES"]]))
        {
            warning("No neutralMass or SMILES data found in the input data. C3SDB predictions will be skipped.",
                    call. = FALSE)
            do_C3SDB <- FALSE
            from <- NULL
        }
        if (!is.null(obj[["SMILES"]]) && (anyNA(obj$SMILES) || any(!nzchar(obj$SMILES))))
            warning("The following rows do not contain SMILES data and will be excluded from C3SDB predictions: ",
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
        else # C3SDB
        {
            if (!is.null(virtualenv))
                reticulate::use_virtualenv(virtualenv)
            py_pickle <- reticulate::import("pickle")
            py_C3SDB <- reticulate::import("c3sdb.ml.data")
            py_bi <- reticulate::import_builtins()
            
            kmcm_svr <- with(py_bi$open(py_C3SDB$pretrained_data("c3sdb_kmcm_svr.pkl"), "rb"), as = "pf",
                             py_pickle$load(pf))
            
            # only do predictions for suspects with SMILES and mass data
            wh <- !is.na(obj$neutralMass) & !is.na(obj$SMILES) & nzchar(obj$SMILES)
            objDo <- obj[wh]; mzTabDo <- mzTab[wh]
            
            doPred <- function(mzs, SMILES, add)
            {
                if (length(mzs) == 0)
                    return(numeric())
                
                if (!add %in% C3SDBAdducts())
                {
                    warning(sprintf("Skipping predictions for unsupported adduct '%s'.", add), call. = FALSE)
                    return(rep(NA_real_, length(mzs))) 
                }

                # use as.list() to handle scalars: https://github.com/rstudio/reticulate/issues/258
                dfi <- py_C3SDB$data_for_inference(as.list(mzs), as.list(rep(add, length(mzs))), as.list(SMILES),
                                                   py_C3SDB$pretrained_data("c3sdb_OHEncoder.pkl"),
                                                   py_C3SDB$pretrained_data("c3sdb_SScaler.pkl"))
                ret <- rep(NA_real_, length(mzs))
                ret[dfi[[2]]] <- kmcm_svr$predict(dfi[[1]])
                if (any(!dfi[[2]]))
                    warning(sprintf("Predictions failed for adduct '%s' for the following SMILES: %s.", add,
                                    paste0(SMILES[dfi[[2]] == FALSE], collapse = ", ")), call. = FALSE)
                return(ret)
            }
            
            from <- data.table(objDo[[matchFromBy]], SMILES = objDo$SMILES)
            setnames(from, 1, matchFromBy)
            
            hasAddCol <- !is.null(objDo[["adduct"]])
            
            for (add in adductsNoNone)
            {
                wh <- if (!hasAddCol || !predictAdductOnly)
                    rep(TRUE, nrow(objDo))
                else
                    is.na(objDo$adduct) | !nzchar(objDo$adduct) | objDo$adduct == add
                if (any(wh))
                {
                    printf("Predicting %d CCS values for adduct '%s'... ", sum(wh), add)
                    from[wh, paste0("CCS_", add) := doPred(mzTabDo[wh][[paste0("mz_", add)]], SMILES, add)]
                    printf(" Done!\n")
                }
            }
        }
        
        # merge in data in input obj
        
        addMobCols <- paste0("mobility_", adductsNoNone)
        addCCSCols <- paste0("CCS_", adductsNoNone)
        if (anyNA(adducts))
        {
            addMobCols <- c(addMobCols, "mobility")
            addCCSCols <- c(addCCSCols, "CCS")
        }
        predCols <- intersect(names(from), c(addMobCols, addCCSCols))
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
            mCol <- "mobility"; cCol <- "CCS"; mzCol <- "mz"
            if (!is.na(addChr))
            {
                mCol <- paste0(mCol, "_", addChr); cCol <- paste0(cCol, "_", addChr); mzCol <- paste0(mzCol, "_", addChr)
            }
            
            if (is.null(obj[[mCol]]) && is.null(obj[[cCol]]))
                next
            if (is.null(mzTab[[mzCol]]))
                next
            
            charges <- if (is.na(addChr))
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
                obj[, (mCol) := doConvert("CCS", get(cCol), mzTab[[mzCol]], charges)]
            else if (is.null(obj[[cCol]]))
                obj[, (cCol) := doConvert("mobility", get(mCol), mzTab[[mzCol]], charges)]
            else # both present: only consider NAs
            {
                hasVal <- function(x) !is.na(x) & (!is.character(x) | nzchar(x))
                # temporarily add charges to handle subset assignments below
                # UNDONE: handle the (unlikely) case that the column is already present?
                obj[, .charge := charges]
                obj[!hasVal(get(mCol)) & hasVal(get(cCol)), (mCol) := doConvert("CCS", get(cCol), mzTab[[mzCol]], .charge)]
                obj[hasVal(get(mCol)) & !hasVal(get(cCol)), (cCol) := doConvert("mobility", get(mCol), mzTab[[mzCol]], .charge)]
                obj[, .charge := NULL]
            }
        }
        printf("Done!\n")
    }
    
    saveCacheData("assignMobilitiesDT", obj, hash)
    
    return(obj[])
})

#' @rdname assignMobilities_susp
#' @aliases assignMobilities,data.frame-method
#' @export
setMethod("assignMobilities", "data.frame", function(obj, ...)
{
    return(as.data.frame(assignMobilities(as.data.table(obj), ...)))
})

#' Returns the adducts supported by C3SDB
#' 
#' Returns the adducts supported by the \href{https://github.com/dylanhross/c3sdb}{C3SDB} \command{Python} package.
#' 
#' @export
C3SDBAdducts <- function()
{
    c("[M+H]+", "[M+Na]+", "[M-H]-", "[M+NH4]+", "[M+K]+", "[M+H-H2O]+", "[M+HCOO]-", "[M+CH3COO]-", "[M+Na-2H]-")
}

#' Automatically installs C3SDB
#'
#' Automatically installs the \href{https://github.com/dylanhross/c3sdb}{C3SDB} \command{Python} package
#'
#' This function uses \CRANpkg{reticulate} to install the \href{https://github.com/dylanhross/c3sdb}{C3SDB}
#' \command{Python} package in a virtual environment.
#'
#' @param envname The name of the virtual \command{Python} environment to install \command{C3SDB} into. Passed to
#'   \link[reticulate:py_install]{reticulate::py_install}. Set to \code{NULL} to not set any virtual \command{Python}
#'   environment.
#' @param clearEnv Set to \code{TRUE} to remove the virtual environment if it already exists (using
#'   \link[reticulate:virtualenv_remove]{reticulate::virtualenv_remove}). Ignored if \code{envname} is \code{NULL}.
#' @param \dots Further arguments passed to \code{\link[reticulate:py_install]{reticulate::py_install}}.
#'
#' @references \insertRef{Ross2020}{patRoon}
#'
#' @export
installC3SDB <- function(envname = "patRoon-C3SDB", clearEnv = FALSE, ...)
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

#' Automatically installs TIMSCONVERT
#'
#' Automatically installs the \href{https://gtluu.github.io/timsconvert/}{TIMSCONVERT} \command{Python} package
#'
#' This function uses \CRANpkg{reticulate} to install the \href{https://gtluu.github.io/timsconvert/}{TIMSCONVERT}
#' Python package in a virtual environment.
#'
#' @param envname The name of the virtual Python environment to install \command{TIMSCONVERT} into. Passed to
#'   \link[reticulate:py_install]{reticulate::py_install}. Set to \code{NULL} to not set any virtual Python environment.
#' @param clearEnv Set to \code{TRUE} to remove the virtual environment if it already exists (using
#'   \link[reticulate:virtualenv_remove]{reticulate::virtualenv_remove}). Ignored if \code{envname} is \code{NULL}.
#' @param \dots Further arguments passed to \code{\link[reticulate:py_install]{reticulate::py_install}}.
#'
#' @references \insertRef{Luu2022}{patRoon}
#'
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

#' @details \code{numericIDLevel} Extracts the numeric part of a given identification level (\emph{e.g.} \code{"3a"}
#'   becomes \samp{3}).
#' @param level The identification level to be converted.
#' @rdname id-conf
#' @export
numericIDLevel <- function(level)
{
    checkmate::assertCharacter(level, any.missing = TRUE, min.chars = 1)
    ret <- integer(length(level))
    ret[is.na(level)] <- NA_integer_
    ret[!is.na(level)] <- as.integer(gsub("[[:alpha:]]*", "", level[!is.na(level)]))
    return(ret)
}

#' @details \code{genIDLevelRulesFile} Generates a template YAML file that is used to configure the rules for automatic
#'   estimation of identification levels. This file can then be used as input for \code{\link{estimateIDConfidence}}.
#' @param out The file path to the target file.
#' @param inLevels,exLevels A \link[=regex]{regular expression} for the identification levels to include or exclude,
#'   respectively. For instance, \code{exLevels="4|5"} would exclude level 4 and 5 from the output file. Set to
#'   \code{NULL} to ignore.
#' @rdname id-conf
#' @export
genIDLevelRulesFile <- function(out, inLevels = NULL, exLevels = NULL)
{
    aapply(checkmate::assertCharacter, . ~ inLevels + exLevels, null.ok = TRUE)
    checkmate::assertPathForOutput(out, overwrite = TRUE)
    
    defFile <- system.file("misc", "IDLevelRules.yml", package = "patRoon")
    
    if (is.null(inLevels) && is.null(exLevels))
        file.copy(defFile, out, overwrite = TRUE)
    else
    {
        rules <- readYAML(defFile)
        if (!is.null(inLevels))
            rules <- rules[grepl(inLevels, names(rules))]
        if (!is.null(exLevels))
            rules <- rules[!grepl(exLevels, names(rules))]
        # UNDONE: this quotes ID levels without sub-level, fix?
        writeYAML(rules, out)
    }
    invisible(NULL)
}
