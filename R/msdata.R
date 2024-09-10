#' @include main.R
NULL

Rcpp::loadModule("MSReadBackend", TRUE)

maybeGetMSFilesForOTIMS <- function(anaInfo)
{
    ret <- getMSFilesFromAnaInfo(anaInfo, "raw", "bruker_tims", mustExist = FALSE)
    if (!is.null(ret))
    {
        # try to load the Bruker TIMS library
        libp <- getOption("patRoon.path.BrukerTIMS", "")
        if (!nzchar(libp) || !initBrukerLibrary(libp))
            ret <- NULL
    }
    return(ret)
}

maybeGetMSFilesForMzR <- function(anaInfo)
{
    if (!requireNamespace("mzR", quietly = TRUE))
        return(NULL) # UNDONE: this will only be relevant when we actually make mzR a soft dependency
    return(getCentroidedMSFilesFromAnaInfo(anaInfo, mustExist = FALSE))
}

maybeGetMSFilesForSC <- function(anaInfo)
{
    # UNDONE: mechanism to prefer IM/centroided data
    ret <- getMSFilesFromAnaInfo(anaInfo, "ims", "mzML", FALSE)
    if (is.null(ret))
        ret <- getCentroidedMSFilesFromAnaInfo(anaInfo, mustExist = FALSE)
    return(ret)
}

maybeGetMSFilesForMSTK <- function(anaInfo)
{
    # UNDONE: mechanism to prefer IM/centroided data
    # UNDONE: check if MSTK was compiled in
    ret <- getMSFilesFromAnaInfo(anaInfo, "ims", "mzML", FALSE)
    if (is.null(ret))
        ret <- getCentroidedMSFilesFromAnaInfo(anaInfo, mustExist = FALSE)
    return(ret)
}

createMSBackend <- function(backendName)
{
    return(switch(backendName,
                  opentims = new(MSReadBackendOTIMS),
                  mzr = new(MSReadBackendMem),
                  streamcraft = new(MSReadBackendSC),
                  mstoolkit = new(MSReadBackendMSTK)))
}

setMSReadBackendMetadata <- function(backend, generator, fileHash = NULL, cacheDB = NULL)
{
    if (is.null(fileHash))
        fileHash <- getMSDataFileHash(backend$getCurrentFile())
    
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
    # UNDONE: test with non-PASEF data
    setMSReadBackendMetadata(backend, function()
    {
        TIMSDB <- withr::local_db_connection(DBI::dbConnect(RSQLite::SQLite(),
                                                            file.path(backend$getCurrentFile(), "analysis.tdf")))
        
        getTIMSMetaTable <- function(name, cols)
        {
            as.data.table(DBI::dbGetQuery(TIMSDB, sprintf("SELECT %s FROM %s", paste0(cols, collapse = ","), name)))
        }
        
        frames <- getTIMSMetaTable("Frames", c("Id", "Time", "MsMsType", "SummedIntensities", "MaxIntensity"))
        setnames(frames, c("Id", "Time", "SummedIntensities", "MaxIntensity"), c("scan", "time", "TIC", "BPC"))
        framesMS <- frames[MsMsType == 0][, MsMsType := NULL]
        framesMS2 <- frames[MsMsType != 0][, MsMsType := NULL]
        PASEFInfo <- getTIMSMetaTable("PasefFrameMsMsInfo", c("Frame", "ScanNumBegin", "ScanNumEnd",
                                                              "IsolationMz", "isolationWidth"))
        setnames(PASEFInfo, c("Frame", "ScanNumBegin", "ScanNumEnd"), c("scan", "subScan", "subScanEnd"))
        # UNDONE: do we need to divide by two?
        PASEFInfo[, c("isolationRangeMin", "isolationRangeMax") := .(IsolationMz - IsolationWidth/2,
                                                                     IsolationMz + IsolationWidth/2)]
        PASEFInfo[, c("IsolationMz", "IsolationWidth") := NULL]
        framesMS2 <- merge(framesMS2, PASEFInfo, by = "scan", sort = FALSE)
        
        return(list(MS1 = framesMS, MS2 = framesMS2))
    })
})

setMethod("initMSReadBackend", "Rcpp_MSReadBackendMem", function(backend)
{
    path <- backend$getCurrentFile()
    hash <- getMSDataFileHash(path)
    db <- openCacheDBScope()
    
    backend <- setMSReadBackendMetadata(backend, function()
    {
        msf <- mzR::openMSfile(path)
        hd <- as.data.table(mzR::header(msf))
        mzR::close(msf)
        
        setnames(hd, c("acquisitionNum", "retentionTime", "totIonCurrent", "basePeakIntensity"),
                 c("scan", "time", "TIC", "BPC"))
        hd[, isolationRangeMin := precursorMZ - isolationWindowLowerOffset]
        hd[, isolationRangeMax := precursorMZ + isolationWindowUpperOffset]
        
        cols <- c("scan", "time", "TIC", "BPC")
        
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
    setMSReadBackendMetadata(backend, function()
    {
        backend$generateSpecMetadata()
        return(list(MS1 = getMSMetadata(backend, 1), MS2 = getMSMetadata(backend, 2)))
    })
})

setMethod("initMSReadBackend", "Rcpp_MSReadBackendMSTK", function(backend)
{
    setMSReadBackendMetadata(backend, function()
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

applyMSData <- function(anaInfo, func, ...)
{
    backends <- getOption("patRoon.MSBackends", character())
    
    for (bn in backends)
    {
        filePaths <- switch(bn,
                            opentims = maybeGetMSFilesForOTIMS(anaInfo),
                            mzr = maybeGetMSFilesForMzR(anaInfo),
                            streamcraft = maybeGetMSFilesForSC(anaInfo),
                            mstoolkit = maybeGetMSFilesForMSTK(anaInfo),
                            NULL)
        if (!is.null(filePaths))
        {
            backend <- createMSBackend(bn)
            # NOTE: disable future parallelization as the backends are already OpenMP parallelized
            # NOTE: the callback can return cached data so opening the file should happen there.
            return(doApply("Map", doPar = FALSE, data = anaInfo$analysis, filePaths, f = func,
                           MoreArgs = list(backend = backend, ...)))
        }
    }
    
    stop("Failed to load a correct MS read backend. Please ensure patRoon.MSBackends is configured properly. See ?patRoon",
         call. = FALSE)
}
