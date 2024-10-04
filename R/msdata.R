#' @include main.R
NULL

Rcpp::loadModule("MSReadBackend", TRUE)

maybeGetMSFilesForOTIMS <- function(anaInfo, types, formats)
{
    if (!"raw" %in% types || !"bruker_ims" %in% formats)
        return(NULL)
    
    ret <- getMSFilesFromAnaInfo(anaInfo, "raw", "bruker_ims", mustExist = FALSE)
    if (!is.null(ret))
    {
        # try to load the Bruker TIMS library
        libp <- getOption("patRoon.path.BrukerTIMS", "")
        if (!nzchar(libp) || !initBrukerLibrary(libp))
            ret <- NULL
    }
    return(ret)
}

maybeGetMSFilesForMzR <- function(anaInfo, types, formats)
{
    if (!"centroid" %in% types)
        return(NULL)
    if (!requireNamespace("mzR", quietly = TRUE))
        return(NULL) # UNDONE: this will only be relevant when we actually make mzR a soft dependency
    return(getCentroidedMSFilesFromAnaInfo(anaInfo, formats = intersect(formats, c("mzML", "mzXML")),
                                           mustExist = FALSE))
}

maybeGetMSFilesForSC <- function(anaInfo, types, formats)
{
    # UNDONE: mechanism to prefer IM/centroided data
    ret <- if ("ims" %in% types) getMSFilesFromAnaInfo(anaInfo, "ims", "mzML", FALSE)
    if (is.null(ret) && "centroid" %in% types)
        ret <- getCentroidedMSFilesFromAnaInfo(anaInfo, formats = intersect(formats, c("mzML", "mzXML")),
                                               mustExist = FALSE)
    return(ret)
}

maybeGetMSFilesForMSTK <- function(anaInfo, types, formats)
{
    # UNDONE: mechanism to prefer IM/centroided data
    ret <- if ("ims" %in% types) getMSFilesFromAnaInfo(anaInfo, "ims", "mzML", FALSE)
    if (is.null(ret) && "centroid" %in% types)
        ret <- getCentroidedMSFilesFromAnaInfo(anaInfo, formats = intersect(formats, c("mzML", "mzXML")),
                                               mustExist = FALSE)
    return(ret)
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
                                 c("Frame", "ScanNumBegin", "ScanNumEnd", "IsolationMz", "IsolationWidth"))
            }
            setnames(ret, c("Frame", "IsolationMz", "TriggerMass", "ScanNumBegin", "ScanNumEnd"),
                     c("scan", "precursorMZ", "precursorMZ", "subScan", "subScanEnd"), skip_absent = TRUE)
            ret[, c("isolationRangeMin", "isolationRangeMax") := .(precursorMZ - IsolationWidth/2,
                                                                   precursorMZ + IsolationWidth/2)]
            ret[, c("precursorMZ", "IsolationWidth") := NULL]
            
            if (type == "MSMS")
            {
                # so we have a table with consistent columns
                # NOTE: having a consistent PASEF/nonPASEF structure is handy, in case data files have both MS/MS types
                ret[, c("subScan", "subScanEnd") := 0]
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
    
    backend <- setMSReadBackendMetadata(backend, hash, function()
    {
        msf <- mzR::openMSfile(path)
        hd <- as.data.table(mzR::header(msf))
        mzR::close(msf)
        
        setnames(hd, c("acquisitionNum", "retentionTime", "totIonCurrent", "basePeakIntensity"),
                 c("scan", "time", "TIC", "BPC"))
        hd[polarity == 0, polarity := -1] # 0 is negative mode for mzR
        hd[, isolationRangeMin := precursorMZ - isolationWindowLowerOffset]
        hd[, isolationRangeMax := precursorMZ + isolationWindowUpperOffset]
        
        cols <- c("scan", "time", "TIC", "BPC", "polarity")
        
        return(list(MS1 = hd[msLevel == 1, cols, with = FALSE],
                    MS2 = hd[msLevel > 1, c(cols, "isolationRangeMin", "isolationRangeMax"), with = FALSE]))
    })
    
    specs <- loadCacheData("MSReadBackendMzR", hash, db)
    
    if (is.null(specs))
    {
        msf <- mzR::openMSfile(path)
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
    backend$open(path)
    backend <- initMSReadBackend(backend)
    return(backend)
}

applyMSData <- function(anaInfo, func,  ..., types = getMSFileTypes(), formats = names(MSFileExtensions()),
                        showProgress = TRUE)
{
    backends <- getOption("patRoon.MSBackends", character())
    
    for (bn in backends)
    {
        if (!backendAvailable(bn))
            next
        
        filePaths <- maybeGetMSFiles(bn, anaInfo, types, formats)
        if (!is.null(filePaths))
        {
            backend <- createMSBackend(bn)
            # NOTE: disable future parallelization as the backends are already OpenMP parallelized
            # NOTE: the callback can return cached data so opening the file should happen there.
            return(doApply("Map", doPar = FALSE, prog = showProgress, data = anaInfo$analysis, filePaths, ..., f = func,
                           MoreArgs = list(backend = backend)))
        }
    }
    
    stop("Failed to load a correct MS read backend. Please ensure patRoon.MSBackends is configured properly. See ?patRoon",
         call. = FALSE)
}
