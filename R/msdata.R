# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' Interface for HRMS and IMS-HRMS raw data
#'
#' Description, configuration and utilities for the raw (IMS-)HRMS data interface of \pkg{patRoon}.
#'
#' Version 3.0 of \pkg{patRoon} introduced an extensible and highly optimized interface to read raw data from HRMS and
#' IMS-HRMS instruments. This interface supports chooseable 'backends' which perform the reading of file data from
#' various formats. Subsequent steps such as the formation of extracted ion chromatograms, mobilograms and collection
#' and averaging of mass spectra are then performed in \code{patRoon}. The interface is largely coded in C++ (using
#' \CRANpkg{Rcpp}), uses \href{https://www.openmp.org/}{OpenMP} parallelization and applies several other optimization
#' strategies to make it suitable to rapidly process large amounts of raw data, \emph{e.g.} as encountered in IMS-HRMS
#' workflow.
#'
#' The following backends for reading (IMS-)HRMS data are currently available: \itemize{
#'
#'   \item \code{"opentims"}: uses the \href{https://github.com/michalsta/opentims}{OpenTIMS} library to read Bruker
#'   TIMS data. This backends supports very fast reading of raw instrument \file{.d} data files directly, and therefore
#'   does not require any file conversions. This backend only supports 64 bit Windows and Linux systems. See the
#'   \verb{Backend installation} section below installation instructions.
#'
#'   \item \code{"mzr"}: uses the \CRANpkg{mzR} package to read \file{.mzML} and \file{.mzXML} files. This package was
#'   more or less the default in \pkg{patRoon} prior to 3.0, and due to its popularity and age is a stable and well
#'   tested option. The \code{mzr} backend currently reads the complete analysis file at once, which makes it more RAM
#'   intensive compared to other backends. This backend currently does not support IMS data. Since \pkg{mzR} is a
#'   dependency of \pkg{patRoon}, no additional installation is necessary.
#'
#'   \item \code{"mstoolkit"}: Uses the \href{https://github.com/mhoopmann/mstoolkit}{MSToolKit} C++ library to read
#'   \file{.mzML} and \file{.mzXML} files, including IMS-HRMS data. The \command{MSToolKit} library has been developed
#'   for many years, and was recently updated with IMS-HRMS support. See the \verb{Backend installation} section below
#'   installation instructions.
#'
#'   \item \code{"streamcraft"}: Uses the \href{https://github.com/odea-project/streamcraft}{StreamCraft} C++ library to
#'   read \file{.mzML} and \file{.mzXML} files, including IMS-HRMS data. The \command{StreamCraft} library is still
#'   young and somewhat experimental. The library is integrated within \pkg{patRoon} and therefore does not require any
#'   further installation.
#'
#' }
#'
#' @section Configuration: The following package options influence the behavior of raw data interface: \itemize{
#'
#'   \item \code{patRoon.MS.backends}: A \code{character} vector with the names of the backends that may be choosen. The
#'   default is all backends. The first backend will be chosen that is available, is able to read at least one of the
#'   available analysis file types and formats (as configured by the \link[=analysis-information]{analysis information})
#'   and supports IMS if needed.
#'
#'   \item \code{patRoon.MS.preferIMS}: A \code{logical} value that indicates whether the IMS data should be prefferred,
#'   even if the processing step does not require IMS data and non-IMS data is also available. Setting this to
#'   \code{TRUE} probably result in some additional computational overhead, but may avoid any inconsistencies between
#'   the IMS data and non-IMS data that may have been introduced during the conversion step of the latter. This option
#'   is only relevant for the \code{mstoolkit} and \code{streamcraft} backends (and if one of these backends is actually
#'   used).
#'
#'   \item \code{patRoon.threads}: An \code{integer} value that indicates the number of threads to use for
#'   parallelization (multithreading). The default is determined from the number of physical cores of the system
#'   (obtained with the \link[parallel:detectCores]{parallel::detectCores} function).
#'
#'   \item \code{patRoon.path.BrukerTIMS}: The file path to the Bruker \command{TDF-SDK} library file. See the
#'   \verb{Backend installation} section below.
#'
#'   }
#'
#'
#' @section Backend installation: The \code{opentims} backend requires the \file{win64/timsdata.dll} (Windows) or
#'   \file{linux64/libtimsdata.so} (Linux) file from the \command{TDF-SDK} from
#'   \href{https://www.bruker.com/protected/en/services/software-downloads/mass-spectrometry/raw-data-access-libraries.html}{Bruker}
#'   (requires login). The \pkg{patRoonExt} package makes these files automatically available for \pkg{patRoon}.
#'   Otherwise the \code{patRoon.path.BrukerTIMS} option should be manually set to the file path of the
#'   \file{timsdata.dll} or \file{linux64/libtimsdata.so} file.
#'
#'   When \pkg{patRoon} is installed from source, \emph{e.g.} on Linux/macOS systems or when using
#'   \code{\link[remotes:install_github]{remotes::install_github}} for installation, then the
#'   \href{https://github.com/rickhelmus/Rmstoolkitlib}{https://github.com/rickhelmus/Rmstoolkitlib} \R package must be
#'   installed in advance.
#'
#'   The \code{availableBackends} function can be used to verify if the dependencies for these backends are met.
#'
#' @references \addCitations{Rcpp} \cr\cr \insertRef{Dagum1998}{patRoon} \cr\cr \insertRef{Lacki2021}{patRoon} \cr\cr
#'   \addCitations{mzR}
#'
#' @name msdata
NULL

Rcpp::loadModule("MSReadBackend", TRUE)

maybeGetMSFilesForOTIMS <- function(anaInfo, types, formats, needIMS)
{
    if (!"raw" %in% types || !"bruker_ims" %in% formats)
        return(NULL)
    
    ret <- getMSFilesFromAnaInfo(anaInfo, "raw", "bruker_ims", mustExist = FALSE)
    if (!is.null(ret))
    {
        # try to load the Bruker TIMS library
        if (!doInitBrukerLib())
            ret <- NULL
    }
    return(ret)
}

maybeGetMSFilesForMzR <- function(anaInfo, types, formats, needIMS)
{
    if (!"centroid" %in% types || needIMS)
        return(NULL)
    if (!requireNamespace("mzR", quietly = TRUE))
        return(NULL) # UNDONE: this will only be relevant when we actually make mzR a soft dependency
    return(getCentroidedMSFilesFromAnaInfo(anaInfo, formats = intersect(formats, c("mzML", "mzXML")),
                                           mustExist = FALSE))
}

maybeGetMSFilesForSC <- function(anaInfo, types, formats, needIMS)
{
    prefIMS <- getOption("patRoon.MS.preferIMS", FALSE)
    
    filesIMS <- if ("ims" %in% types || needIMS) getMSFilesFromAnaInfo(anaInfo, "ims", "mzML", FALSE)
    filesCentroid <- if ("centroid" %in% types && !needIMS && (is.null(filesIMS) || !prefIMS))
        getCentroidedMSFilesFromAnaInfo(anaInfo, formats = intersect(formats, c("mzML", "mzXML")), mustExist = FALSE)
    
    if (!is.null(filesIMS) && !is.null(filesCentroid))
        return(if (prefIMS) filesIMS else filesCentroid)
    if (!is.null(filesIMS))
        return(filesIMS)
    return(filesCentroid)
}

maybeGetMSFilesForMSTK <- function(anaInfo, types, formats, needIMS)
{
    # so far the same constraints as SC
    return(maybeGetMSFilesForSC(anaInfo, types, formats, needIMS))
}

maybeGetMSFiles <- function(bn, ...)
{
    return(switch(bn,
                  opentims = maybeGetMSFilesForOTIMS(...),
                  mzr = maybeGetMSFilesForMzR(...),
                  streamcraft = maybeGetMSFilesForSC(...),
                  mstoolkit = maybeGetMSFilesForMSTK(...),
                  NULL))
}
createMSBackend <- function(backendName)
{
    return(switch(backendName,
                  opentims = new(MSReadBackendOTIMS),
                  mzr = new(MSReadBackendMem),
                  streamcraft = new(MSReadBackendSC),
                  mstoolkit = new(MSReadBackendMSTK)))
}

setMSReadBackendMetadata <- function(backend, fileHash, generator, cacheDB = NULL)
{
    # NOTE: include class name in hash in case different backends generate different metadata for the same file
    hash <- makeHash(class(backend), fileHash)
    meta <- loadCacheData("MSReadBackendMetadata", hash, cacheDB)
    
    if (is.null(meta))
    {
        meta <- generator()
        saveCacheData("MSReadBackendMetadata", meta, hash, cacheDB)
    }
    
    setSpecMetadata(backend, meta$MS1, meta$MS2)
    
    if (backend$getNeedIMS())
    {
        mobilities <- loadCacheData("MSReadBackendMobilities", hash, cacheDB)
        if (is.null(mobilities))
        {
            mobilities <- backend$generateMobilities()
            saveCacheData("MSReadBackendMobilities", mobilities, hash, cacheDB)
        }
        backend$setMobilities(mobilities)
    }
    
    return(backend)
}

setMethod("initMSReadBackend", "Rcpp_MSReadBackendOTIMS", function(backend)
{
    fileHash <- getMSDataFileHash(file.path(backend$getCurrentFile(), "analysis.tdf_bin"))
    
    setMSReadBackendMetadata(backend, fileHash, function()
    {
        TIMSDB <- withr::local_db_connection(DBI::dbConnect(RSQLite::SQLite(),
                                                            file.path(backend$getCurrentFile(), "analysis.tdf")))
        
        getTIMSMetaTable <- function(name, cols)
        {
            as.data.table(DBI::dbGetQuery(TIMSDB, sprintf("SELECT %s FROM %s", paste0(cols, collapse = ","), name)))
        }
        getTIMSMetaTableMSMS <- function(type)
        {
            ret <- if (type == "MSMS")
                getTIMSMetaTable("FrameMsMsInfo", c("Frame", "TriggerMass", "IsolationWidth"))
            else # PASEF
            {
                getTIMSMetaTable("PasefFrameMsMsInfo",
                                 c("Frame", "ScanNumBegin", "ScanNumEnd", "IsolationMz", "IsolationWidth", "Precursor"))
            }
            setnames(ret, c("Frame", "IsolationMz", "TriggerMass", "ScanNumBegin", "ScanNumEnd"),
                     c("scan", "precursorMZ", "precursorMZ", "subScan", "subScanEnd"), skip_absent = TRUE)
            ret[, c("isolationRangeMin", "isolationRangeMax") := .(IsolationWidth/2, IsolationWidth/2)]
            ret[, IsolationWidth := NULL]
            
            if (type == "MSMS")
            {
                # so we have a table with consistent columns
                # NOTE: having a consistent PASEF/nonPASEF structure is handy, in case data files have both MS/MS types
                ret[, c("subScan", "subScanEnd") := 0]
            }
            else # PASEF
            {
                # Fix isolationMz and range from precursor table
                
                precTab <- getTIMSMetaTable("Precursors", c("Id", "MonoisotopicMz"))
                setnames(ret, "precursorMZ", "precursorMZ_orig")
                ret[, precursorMZ := precTab[match(ret$Precursor, Id)]$MonoisotopicMz]
                ret[is.na(precursorMZ), precursorMZ := precursorMZ_orig] # in case precursor is missing (why?? but it happens...)
                # adjust range in such a way that it covers the original range
                ret[, c("isolationRangeMin", "isolationRangeMax") := .(precursorMZ - (precursorMZ_orig - isolationRangeMin),
                                                                       (precursorMZ_orig + isolationRangeMax) - precursorMZ)]
                ret[, c("Precursor", "precursorMZ_orig") := NULL]
            }
            
            return(ret)
        }
        
        frames <- getTIMSMetaTable("Frames", c("Id", "Time", "ScanMode", "MsMsType", "SummedIntensities",
                                               "MaxIntensity", "Polarity"))
        setnames(frames, c("Id", "Time", "SummedIntensities", "MaxIntensity", "Polarity"),
                 c("scan", "time", "TIC", "BPC", "polarity"))
        frames[, polarity := fcase(polarity == "+", 1L,
                                   polarity == "-", -1L,
                                   default = 0L)]
        
        framesMS <- frames[MsMsType == 0][, MsMsType := NULL]
        framesMS2 <- frames[MsMsType != 0]
        
        if (any(!framesMS2$ScanMode %in% c(0, 4, 8)))
            warning("The OpenTIMS backend has only been tested with MS only, bbCID and PASEF data.", call. = FALSE)

        MSMSInfo <- if (2 %in% framesMS2$MsMsType) getTIMSMetaTableMSMS("MSMS")
        PASEFInfo <- if (8 %in% framesMS2$MsMsType) getTIMSMetaTableMSMS("PASEF")
        
        framesMS2 <- merge(framesMS2, rbind(MSMSInfo, PASEFInfo), by = "scan", sort = FALSE)
        # zero out isolation range for isCID/bbCID (ie DIA)
        framesMS2[ScanMode %in% c(3, 4), c("isolationRangeMin", "isolationRangeMax") := 0]
        framesMS2[, c("ScanMode", "MsMsType") := NULL]

        return(list(MS1 = framesMS, MS2 = framesMS2))
    })
})

setMethod("initMSReadBackend", "Rcpp_MSReadBackendMem", function(backend)
{
    path <- backend$getCurrentFile()
    hash <- getMSDataFileHash(path)
    db <- openCacheDBScope()
    
    hd <- NULL # might be set below, so we don't need to load it twice
    
    backend <- setMSReadBackendMetadata(backend, hash, function()
    {
        msf <- mzR::openMSfile(path)
        hd <<- as.data.table(mzR::header(msf))
        mzR::close(msf)
        
        if (!is.null(hd[["centroided"]]) && !all(hd$centroided))
            stop(sprintf("Please make sure that file '%s' is centroided!", path), call. = FALSE)
        
        setnames(hd, c("acquisitionNum", "retentionTime", "totIonCurrent", "basePeakIntensity"),
                 c("scan", "time", "TIC", "BPC"))
        hd[polarity == 0, polarity := -1] # 0 is negative mode for mzR
        setnames(hd, c("isolationWindowLowerOffset", "isolationWindowUpperOffset"),
                 c("isolationRangeMin", "isolationRangeMax"))
        
        cols <- c("scan", "time", "TIC", "BPC", "polarity")
        
        return(list(MS1 = hd[msLevel == 1, cols, with = FALSE],
                    MS2 = hd[msLevel > 1, c(cols, "isolationRangeMin", "isolationRangeMax", "precursorMZ"), with = FALSE]))
    })
    
    specs <- loadCacheData("MSReadBackendMzR", hash, db)
    
    if (is.null(specs))
    {
        msf <- mzR::openMSfile(path)
        if (is.null(hd))
            hd <- as.data.table(mzR::header(msf))
        ps <- mzR::peaks(msf)
        mzR::close(msf)
        specs <- list(MS1 = ps[hd$msLevel == 1], MS2 = ps[hd$msLevel > 1])
        saveCacheData("MSReadBackendMzR", specs, hash, db)
    }
    
    backend$setSpectra(specs$MS1, specs$MS2)
    
    return(backend)
})

setMethod("initMSReadBackend", "Rcpp_MSReadBackendSC", function(backend)
{
    setMSReadBackendMetadata(backend, getMSDataFileHash(backend$getCurrentFile()), function()
    {
        backend$generateSpecMetadata()
        return(list(MS1 = getMSMetadata(backend, 1), MS2 = getMSMetadata(backend, 2)))
    })
})

setMethod("initMSReadBackend", "Rcpp_MSReadBackendMSTK", function(backend)
{
    setMSReadBackendMetadata(backend, getMSDataFileHash(backend$getCurrentFile()), function()
    {
        backend$generateSpecMetadata()
        return(list(MS1 = getMSMetadata(backend, 1), MS2 = getMSMetadata(backend, 2)))
    })
})

openMSReadBackend <- function(backend, path)
{
    setOMPThreads()
    backend$open(path)
    backend <- initMSReadBackend(backend)
    return(backend)
}

applyMSData <- function(anaInfo, func,  ..., types = getMSFileTypes(), formats = names(MSFileExtensions()),
                        needIMS = FALSE, showProgress = TRUE)
{
    backends <- getOption("patRoon.MS.backends", character())
    
    for (bn in backends)
    {
        if (!backendAvailable(bn))
            next
        
        filePaths <- maybeGetMSFiles(bn, anaInfo, types, formats, needIMS)
        if (!is.null(filePaths))
        {
            backend <- createMSBackend(bn)
            backend$setNeedIMS(needIMS)
            printf("Using '%s' backend for reading MS data.\n", bn) # UNDONE: make printing optional (arg/option?)
            # NOTE: disable future parallelization as the backends are already OpenMP parallelized
            # NOTE: the callback can return cached data so opening the file should happen there.
            return(doMap(FALSE, prog = showProgress, stripEnv = FALSE, data = anaInfo$analysis, filePaths, ..., f = func,
                         MoreArgs = list(backend = backend)))
        }
    }

    msg <- "Could not load MS data."
    if (!setequal(types, getMSFileTypes()))
        msg <- paste(msg, sprintf("Allowed types: %s.", paste0("\"", types, "\"", collapse = ", ")))
    if (!setequal(formats, names(MSFileExtensions())))
        msg <- paste(msg, sprintf("Allowed formats: %s.", paste0("\"", formats, "\"", collapse = ", ")))
    if (needIMS)
        msg <- paste(msg, "IMS data is required.")
    stop(msg, " Please ensure all data is present and the \"patRoon.MS.backends\" option is configured properly: see ?patRoon",
         call. = FALSE)
}

#' @details The \code{availableBackends} function is used to query the available backends on the system.
#'
#' @param anaInfo Optional. If not \code{NULL} then \code{anaInfo} should be a \link[=analysis-information]{analysis
#'   information} table, and only those backends that can read each of the analyses are returned.
#' @param needIMS Only applicable if \code{anaInfo} is set: set to \code{TRUE} if IMS data is required.
#' @param verbose Set to \code{TRUE} to print the status of each backend.
#'
#' @return \code{availableBackends} returns (invisibly) a \code{character} vector with the names of the available
#'   backends.
#'
#' @rdname msdata
#' @export
availableBackends <- function(anaInfo = NULL, needIMS = FALSE, verbose = TRUE)
{
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, null.ok = TRUE)
    checkmate::assertFlag(needIMS)
    checkmate::assertFlag(verbose)
    
    allBackends <- getMSReadBackends()
    
    unselected <- setdiff(allBackends, getOption("patRoon.MS.backends", character()))
    notCompiled <- allBackends[!sapply(allBackends, backendAvailable)]
    noAnas <- if (is.null(anaInfo))
        character()
    else
        allBackends[sapply(allBackends, function(b) is.null(maybeGetMSFiles(b, anaInfo, getMSFileTypes(), names(MSFileExtensions()), needIMS)))]
    
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
