# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

doSIRIUSLogin <- function(login, force, SIRIUSAPI)
{
    if (isFALSE(login))
        return(invisible(NULL)) # no need to do anything
    
    isLoggedIn <- SIRIUSAPI$login_and_account_api$IsLoggedIn()
    
    if (length(login) == 1 && login == "check" && !isLoggedIn)
        stop("There is no active SIRIUS login. Please consult the SIRIUS documentation and patRoon handbook for details.")
    
    else if (force || !isLoggedIn)
    {
        if (length(login) == 1 && login == "interactive")
        {
            if (!interactive())
                stop("Cannot perform interactive login in non-interactive R sessions!", call. = FALSE)
            
            # NOTE: if "username" is part of the prompt then the RStudio backend of getPass won't hide the text input
            login <- c(username = getPass::getPass("Please enter your SIRIUS username", noblank = TRUE),
                       password = getPass::getPass("Please enter your SIRIUS password", noblank = TRUE))
        }
        
        if (!"username" %in% names(login) || !nzchar(login["username"]))
            stop("Please provide the username of your SIRIUS account", call. = FALSE)
        if (!"password" %in% names(login) || !nzchar(login["password"]))
            stop("Please provide the password of your SIRIUS account", call. = FALSE)
        
        SIRCreds <- RSirius::AccountCredentials$new(username = login["username"], password = login["password"])
        SIRIUSAPI$login_and_account_api$Login(accept_terms = TRUE, account_credentials = SIRCreds)
    }
    invisible(NULL)
}

startSIRIUS <- function(path)
{
    checkPackage("RSirius", "sirius-ms/sirius-client-openAPI", ghSubDir = "client-api_r/generated")
    
    sdk <- RSirius::SiriusSDK$new()
    # SIRIUSAPI <- sdk$attach_or_start_sirius()
    SIRIUSAPI <- sdk$attach_to_sirius()
    shutdownSIR <- FALSE
    if (is.null(SIRIUSAPI))
    {
        SIRIUSAPI <- sdk$start_sirius(sirius_path = path)
        shutdownSIR <- TRUE
    }
    withr::defer_parent({
        tryCatch({
            if (shutdownSIR)
                sdk$shutdown_sirius()
        }, error = function(e) NULL)
    })
    
    return(SIRIUSAPI)
}

openSIRIUSProject <- function(projectPath, SIRIUSAPI, runMode)
{
    if (runMode == "read" && (is.null(projectPath) || !file.exists(projectPath)))
        stop("projectPath must be provided and exist when runMode is 'read'", call. = FALSE)
    
    projectID <- if (!is.null(projectPath) && length(names(projectPath) > 0)) names(projectPath)[1] else "patRoonProjectID"
    projectPath <- if (is.null(projectPath))
        tempfile("patRoonSIRIUS", fileext = ".sirius")
    else
        normalizePath(projectPath, mustWork = FALSE, winslash = "/")
    
    if (file.exists(projectPath) && runMode == "read")
        SIRIUSAPI$projects_api$OpenProject(projectID, projectPath)
    else
    {
        unlink(projectPath)
        SIRIUSAPI$projects_api$CreateProject(projectID, projectPath)
    }
    
    withr::defer_parent({
        tryCatch({
            SIRIUSAPI$projects_api$CloseProject(projectID)
        }, error = function(e) NULL)
    })
    
    return(projectID)
}

getSIRIUSFormulaCandidates <- function(projectID, SIRIUSAPI, SIRFeatID, adduct)
{
    formCands <- SIRIUSAPI$features_api$GetFormulaCandidates(projectID, SIRFeatID, opt_fields = "statistics")
    if (length(formCands) == 0)
        return(data.table())
    tab <- rbindlist(lapply(formCands, function(fc)
    {
        l <- fc$toList()
        l$medAbsMassDev <- l$medianMassDeviation$absolute
        l$medRelMassDev <- l$medianMassDeviation$ppm
        l$medianMassDeviation <- NULL
        if (is.null(l$medAbsMassDev))
            l$medAbsMassDev <- NA_real_
        if (is.null(l$medRelMassDev))
            l$medRelMassDev <- NA_real_
        return(l)
    }), fill = TRUE)
    # NOTE: explainedPeaks is re-caclulated later, but putting it now places the column before explainablePeaks
    setnames(tab,
             c("molecularFormula", "siriusScore", "siriusScoreNormalized", "isotopeScore", "treeScore",
               "numOfExplainedPeaks", "numOfExplainablePeaks"),
             c("neutral_formula", "score", "scoreNormalized", "isoScore", "MSMSScore", "explainedPeaks",
               "explainablePeaks"), skip_absent = TRUE)
    return(tab)
}

getSIRIUSFragInfos <- function(projectID, SIRIUSAPI, SIRFeatID, SIRFormIDs, PLMS2)
{
    # HACK: no PL IDs (yet) for feature specific peak lists
    hasPLID <- !is.null(PLMS2[["ID"]])
    
    # NOTE: there may be duplicate formIDs --> only query unique ones and then expand to all SIRFormIDs
    fragInfos <- sapply(unique(SIRFormIDs), function(fid)
    {
        as <- SIRIUSAPI$features_api$GetFormulaAnnotatedSpectrum(projectID, SIRFeatID, fid)
        ionform <- calculateIonFormula(as$spectrumAnnotation$molecularFormula,
                                       gsub(" ", "", as$spectrumAnnotation$adduct))
        rbindlist(lapply(as$peaks, function(p)
        {
            if (is.null(p$peakAnnotation) || !p$peakAnnotation$isValid()) # no annotations
            {
                return(data.table(mz = numeric(0), ion_formula_mz = character(0), formula_SIR = character(0),
                                  adduct = character(0), error_mz = numeric(0), error_ppm = numeric(0),
                                  PLID = if (hasPLID) numeric(0), ion_formula = character(0),
                                  neutral_loss = character(0)))
            }
            
            anns <- data.table(mz = p$mz, ion_formula_mz = p$peakAnnotation$exactMass,
                               formula_SIR = p$peakAnnotation$molecularFormula,
                               adduct = gsub(" ", "", p$peakAnnotation$adduct),
                               error_mz = p$peakAnnotation$massDeviationMz,
                               error_ppm = p$peakAnnotation$massDeviationPpm)
            if (hasPLID)
                anns[, PLID := PLMS2[which.min(abs(p$mz - mz))]$ID]
            
            # SIRIUS neutralizes fragments, make them ion again
            anns[, ion_formula := mapply(formula_SIR, adduct, FUN = calculateIonFormula)]
            anns[, neutral_loss := sapply(ion_formula, subtractFormula, formula1 = ionform)]
            return(anns)
        }))
    }, simplify = FALSE)
    
    return(fragInfos[SIRFormIDs])
}

getSIRIUSFingerprints <- function(projectID, SIRIUSAPI, SIRFeatID, SIRFormIDs, fingerIDData)
{
    fps <- data.table()
    for (fid in unique(SIRFormIDs))
    {
        # NOTE: GetFingerprintPrediction() throws an error if there are not FPs, so we use GetFormulaCandidate() with
        # opt_fields instead, which returns NULL if there are no FPs
        formC <- SIRIUSAPI$features_api$GetFormulaCandidate(projectID, SIRFeatID, fid, opt_fields = "predictedFingerprint")
        if (!is.null(formC$predictedFingerprint))
            set(fps, j = formC$molecularFormula, value = unlist(formC$predictedFingerprint))
    }
    
    if (nrow(fps) > 0)
        fps <- cbind(fps, fingerIDData)

    return(fps)
}

runSIRIUS <- function(runMode, fGroups, MSPeakLists, IMSSpecSims, adduct, SIRIUSAPI, SIRIUSPath, projectPath, config,
                      login, alwaysLogin, formulasOnly, calculateFeatures, cacheName, getFingerprints,
                      topMostStructures = NULL)
{
    doFGroup <- function(grp, ana = NULL)
    {
        if (!is.null(IMSSpecSims) && grp %chin% IMSSpecSims$group)
            return(FALSE)
        pl <- if (!is.null(ana)) MSPeakLists[[ana, grp]] else MSPeakLists[[grp]]
        if (is.null(pl[["MS"]]) || is.null(pl[["MSMS"]]) || !any(pl[["MS"]]$precursor))
            return(FALSE)
        return(TRUE)
    }
    
    db <- openCacheDBScope()
    baseHash <- makeHash(config, formulasOnly, calculateFeatures, getFingerprints, topMostStructures)
    
    gNamesTBD <- names(fGroups)[sapply(names(fGroups), doFGroup)]
    fgAdd <- getFGroupAdducts(gNamesTBD, annotations(fGroups)[match(gNamesTBD, group)], adduct, "sirius")
    
    hashes <- cachedData <- NULL
    if (calculateFeatures)
    {
        hashes <- pruneList(sapply(analyses(fGroups), function(ana)
        {
            anaHashes <- sapply(gNamesTBD, function(g)
            {
                if (!doFGroup(g, ana))
                    return(NA_character_)
                makeHash(baseHash, fgAdd$grpAdductsChr[[g]], MSPeakLists[[ana, g]]$MS[precursor == TRUE]$mz, MSPeakLists[[ana, g]]$MSMS)
            })
            anaHashes <- anaHashes[!is.na(anaHashes)]
        }, simplify = FALSE), checkEmptyElements = TRUE)
        
        cachedData <- pruneList(sapply(analyses(fGroups), function(ana)
        {
            pruneList(setNames(loadCacheData(cacheName, hashes[[ana]], dbArg = db, simplify = FALSE), names(hashes[[ana]])))
        }, simplify = FALSE), checkEmptyElements = TRUE)
        
        # update groups to be done: if at least one ana/fgroup pair is NULL then it needs to be done
        gNamesTBD <- sapply(gNamesTBD, function(g)
        {
            for (ana in analyses(fGroups))
            {
                if (is.null(hashes[[ana]][[g]]) || !is.null(cachedData[[ana]][[g]]))
                    return(NA_character_)
            }
            return(g)
        })
        gNamesTBD <- gNamesTBD[!is.na(gNamesTBD)]
    }
    else
    {
        hashes <- sapply(gNamesTBD, function(g)
        {
            makeHash(baseHash, fgAdd$grpAdductsChr[[g]], MSPeakLists[[g]]$MS[precursor == TRUE]$mz, MSPeakLists[[g]]$MSMS)
        })
        cachedData <- pruneList(setNames(loadCacheData(cacheName, hashes, dbArg = db, simplify = FALSE), gNamesTBD))
        gNamesTBD <- setdiff(gNamesTBD, names(cachedData))
    }
    
    results <- list()
    if (length(gNamesTBD) > 0)
    {
        if (is.null(SIRIUSAPI))
            SIRIUSAPI <- startSIRIUS(SIRIUSPath)
        
        doSIRIUSLogin(login, alwaysLogin, SIRIUSAPI)
        projectID <- openSIRIUSProject(projectPath, SIRIUSAPI, runMode)
        
        if (runMode == "execute")
        {
            makeSIRSpec <- function(pl, lev, pmz)
            {
                peaks <- Map(pl$mz, pl$intensity, f = RSirius::SimplePeak$new)
                ret <- RSirius::BasicSpectrum$new(msLevel = lev, peaks = peaks)
                if (lev > 1)
                    ret$precursorMz = pmz
                return(ret)
            }
            addSIRFeature <- function(id, pl, add, ch)
            {
                plmz <- pl$MS[precursor == TRUE]$mz
                # UNDONE: charge can be set, but only charge 1 is supported by SIRIUS docs? (https://v6.docs.sirius-ms.io/adducts/)
                if (ch < -1L)
                    ch <- -1L
                else if (ch > 2L)
                    ch <- 1L
                RSirius::FeatureImport$new(externalFeatureId = id, ionMass = plmz,
                                           detectedAdducts = list(add), charge = ch,
                                           mergedMs1 = makeSIRSpec(pl$MS, 1L, plmz),
                                           ms2Spectra = list(makeSIRSpec(pl$MSMS, 2L, plmz)))
            }
            
            SIRFeatList <- if (calculateFeatures)
            {
                sfeats <- sapply(analyses(fGroups), function(ana)
                {
                    pruneList(sapply(gNamesTBD, function(fg)
                    {
                        if (is.null(cachedData[[ana]][[fg]]) && doFGroup(fg, ana))
                            addSIRFeature(paste0(ana, "_", fg), MSPeakLists[[ana, fg]], fgAdd$grpAdductsChr[[fg]],
                                          fgAdd$grpAdducts[[fg]]@charge)
                    }, simplify = FALSE))
                }, simplify = FALSE)
                unlist(sfeats, recursive = FALSE)
            }
            else
            {
                sapply(gNamesTBD, \(fg) addSIRFeature(fg, MSPeakLists[[fg]], fgAdd$grpAdductsChr[[fg]],
                                                      fgAdd$grpAdducts[[fg]]@charge), simplify = FALSE)
            }
            
            stopifnot(length(SIRFeatList) > 0)
            
            SIRIUSAPI$features_api$AddAlignedFeatures(project_id = projectID, SIRFeatList)
            
            if (is.null(config))
                config <- getSIRIUSConfig(SIRIUSAPI = SIRIUSAPI, SIRIUSPath = SIRIUSPath, login = FALSE) # NOTE: should already be logged in
            config$spectraSearchParams$enabled <- FALSE
            config$formulaIdParams$enabled <- TRUE
            config$fingerprintPredictionParams$enabled <- getFingerprints
            config$structureDbSearchParams$enabled <- !formulasOnly
            config$canopusParams$enabled <- !formulasOnly # NOTE: must be enabled for compounds
            config$msNovelistParams$enabled <- FALSE
            
            printf("Running SIRIUS job...\n")
            job <- SIRIUSAPI$jobs_api$StartJob(projectID, config)
            
            # NOTE: maxProgress can change during the job execution, so we normalize the current progress to it at each update
            prog <- openProgBar(0, 1)
            repeat
            {
                Sys.sleep(1)
                jp <- SIRIUSAPI$jobs_api$GetJob(projectID, job$id)$progress
                if (jp$state %in% c("CANCELED", "FAILED", "DONE"))
                    break
                setTxtProgressBar(prog, jp$currentProgress / jp$maxProgress)
            }
            setTxtProgressBar(prog, 1)
            close(prog)
        }
        
        SIRFeatListImp <- SIRIUSAPI$features_api$GetAlignedFeatures(projectID)
        names(SIRFeatListImp) <- sapply(SIRFeatListImp, \(f) if (is.null(f$externalFeatureId)) NA_character_ else f$externalFeatureId)
        
        # UNDONE: this sometimes fails, why?
        fingerIDData <- tryCatch(SIRIUSAPI$projects_api$GetFingerIdData(projectID, charge = 1L),
                                 error = \(...) data.table())
        
        getResFromFeat <- function(sirFeat, PLMS2)
        {
            ret <- list()
            
            ret$formCands <- getSIRIUSFormulaCandidates(projectID, SIRIUSAPI, sirFeat$alignedFeatureId)
            if (nrow(ret$formCands) == 0)
                return(NULL)
            
            # NOTE: frag info for SIRIUS is only available from formula candidates(!)
            fragInfos <- getSIRIUSFragInfos(projectID, SIRIUSAPI, sirFeat$alignedFeatureId, ret$formCands$formulaId, PLMS2)
            set(ret$formCands, j = "fragInfo", value = fragInfos[ret$formCands$formulaId])
            set(ret$formCands, j = "explainedPeaks", value = sapply(ret$formCands$fragInfo, nrow))
            
            ret$formCands <- addMiscFormulaInfo(ret$formCands, fgAdd$grpAdducts[[sirFeat$externalFeatureId]])
            
            if (!formulasOnly)
            {
                # BUG: opt_fields doesn't seem to do anything
                # UNDONE: support opt_fields="libraryMatches"?
                ret$structCands <- SIRIUSAPI$features_api$GetStructureCandidatesPaged(projectID,
                                                                                      sirFeat$alignedFeatureId,
                                                                                      page = 0, size = topMostStructures,
                                                                                      opt_fields = c("dbLinks"))
                ret$structCands <- rbindlist(lapply(ret$structCands$content, \(sc) sc$toList()), fill = TRUE)
                if (nrow(ret$structCands) > 0)
                {
                    setnames(ret$structCands, c("inchiKey", "smiles", "structureName", "xlogP", "molecularFormula", "csiScore"),
                             c("InChIKey1", "SMILES", "compoundName", "XlogP", "neutral_formula", "score"), skip_absent = TRUE)
                    # UNDONE: add InChI and InChIKey?
                    if (is.null(ret$structCands[["compoundName"]]))
                        ret$structCands[, compoundName := InChIKey1]
                    else
                        ret$structCands[, compoundName := fifelse(nzchar(compoundName), compoundName, InChIKey1)]
                }
                
                ret$fingerprints <- if (getFingerprints && nrow(ret$structCands) > 0)
                {
                    # only get fingerprints for relavant formulae
                    getSIRIUSFingerprints(projectID, SIRIUSAPI, sirFeat$alignedFeatureId, ret$structCands$formulaId,
                                          fingerIDData)
                }
                else
                    data.table()
                
            }
            else if (getFingerprints)
            {
                # get ALL fingerprints
                ret$fingerprints <- getSIRIUSFingerprints(projectID, SIRIUSAPI, sirFeat$alignedFeatureId,
                                                          ret$formCands$formulaId, fingerIDData)
            }
            
            return(ret)
        }
        
        printf("Importing SIRIUS results...\n")
        
        if (calculateFeatures)
        {
            results <- withProg(length(analyses(fGroups)), FALSE, sapply(analyses(fGroups), function(ana)
            {
                ret <- pruneList(sapply(gNamesTBD, function(grp)
                {
                    featID <- paste0(ana, "_", grp)
                    if (is.null(SIRFeatListImp[[featID]]))
                        return(NULL)
                    return(getResFromFeat(SIRFeatListImp[[featID]], MSPeakLists[[ana, grp]]$MSMS))
                }, simplify = FALSE))
                saveCacheDataList(cacheName, ret, hashes[[ana]][names(ret)], dbArg = db)
                doProgress()
                return(ret)
            }, simplify = FALSE))
        }
        else
        {
            results <- withProg(length(gNamesTBD), FALSE, sapply(gNamesTBD, function(grp)
            {
                # NOTE: check if feature is present: already done by doFGroup(), but may be removed if we're importing data
                if (is.null(SIRFeatListImp[[grp]]))
                    return(NULL)
                ret <- getResFromFeat(SIRFeatListImp[[grp]], MSPeakLists[[grp]]$MSMS)
                doProgress()
                return(ret)
            }, simplify = FALSE))
            results <- pruneList(results)
            saveCacheDataList(cacheName, results, hashes[names(results)], dbArg = db)
        }
    }
    
    # add cached data
    if (calculateFeatures)
    {
        for (ana in names(cachedData))
        {
            results[[ana]] <- c(cachedData[[ana]], results[[ana]])
            results[[ana]] <- results[[ana]][intersect(names(fGroups), names(results[[ana]]))]
        }
    }
    else
    {
        results <- c(cachedData, results)
        results <- results[intersect(names(fGroups), names(results))]
    }
    
    return(results)
}


getMS2QTFPs <- function(featAnnSIR)
{
    # get all SIRIUS fingerprints in a table that is compatible with MS2Quant/MS2Tox
    
    FPs <- featAnnSIR@fingerprints[sapply(featAnnSIR@fingerprints, nrow) > 0]
    if (length(FPs) == 0)
        return(data.table())
    
    allFPs <- rbindlist(lapply(FPs, transpose, keep.names = "neutral_formula", make.names = "absoluteIndex"),
                        idcol = "group")
    fpColRange <- seq(3, ncol(allFPs))
    setnames(allFPs, fpColRange, paste0("Un", names(allFPs)[fpColRange]))
    allFPs[, id := as.character(seq_len(nrow(allFPs)))]
    allFPs[, ionization := mapply(group, neutral_formula, FUN = function(g, f)
    {
        # UNDONE: need to check for empty fragInfo?
        return(featAnnSIR[[g]][neutral_formula == f]$fragInfo[[1]]$ionization[1])
    })]
    
    allFPs[, predion := paste0(neutral_formula, "_", ionization)]
    
    # dummies for MS2Tox
    allFPs[, foldernumber := 0]
    allFPs[, predform := ""]
    
    return(allFPs[])
}

predictRespFactorsSIRFPs <- function(featAnnSIR, gInfo, calibrants, eluent, organicModifier, pHAq, concUnit)
{
    featAnnSIR <- featAnnSIR[gInfo$group]
    
    allFPs <- getMS2QTFPs(featAnnSIR)
    
    if (any(!allFPs$ionization %chin% c("[M+H]+", "[M]+")))
    {
        warning("One or more features are with adducts other than [M+H]+/[M]+. ",
                "These are not (yet) supported by MS2Quant and will be ignored.", call. = FALSE)
        allFPs <- allFPs[ionization %chin% c("[M+H]+", "[M]+")]
    }
    
    if (nrow(allFPs) == 0)
        return(list(RFs = data.table(), MD = list()))
    
    # NOTE: we set the area to one to easily get the response factor below
    unknowns <- data.table(identifier = allFPs$id, retention_time = gInfo[match(allFPs$group, group)]$ret,
                           SMILES = NA_character_, conc_M = NA_real_, area = 1)
    
    baseHash <- makeHash(calibrants, eluent, organicModifier, pHAq)
    splFPs <- split(allFPs[, -c("id", "group")], seq_len(nrow(allFPs)))
    hashes <- sapply(seq_len(nrow(allFPs)), function(i) makeHash(baseHash, splFPs[[i]], unknowns$retention_time[i]))
    
    cachedData <- pruneList(loadCacheData("RF_SIRFP", hashes, simplify = FALSE))
    indsTODO <- which(!hashes %in% names(cachedData))
    hashesTODO <- hashes[indsTODO]
    
    MS2QRes <- NULL
    if (length(indsTODO) > 0)
    {
        MS2QRes <- getMS2QuantRes(calibrants, unknowns[indsTODO], eluent, organicModifier, pHAq, allFPs[indsTODO])
        saveCacheData("MS2QMD", MS2QRes$MD, baseHash)
        MS2QRes$RFs <- MS2QRes$RFs[, c("identifier", "RF_pred"), with = FALSE]
        setnames(MS2QRes$RFs, "RF_pred", "RF_SIRFP")
        for (i in seq_len(nrow(MS2QRes$RFs)))
            saveCacheData("RF_SIRFP", MS2QRes$RFs$RF_SIRFP[i], hashesTODO[i])
        
        MS2QRes$RFs <- merge(allFPs[, c("group", "neutral_formula", "id"), with = FALSE], MS2QRes$RFs, by.x = "id",
                             by.y = "identifier", sort = FALSE)
    }
    
    if (length(cachedData) > 0)
    {
        cachedRFs <- rbindlist(lapply(cachedData, function(cd) data.table(RF_SIRFP = cd)), idcol = "hash")
        cachedRFs[, neutral_formula := allFPs$neutral_formula[match(hash, hashes)]]
        cachedRFs[, group := allFPs$group[match(hash, hashes)]]
        cachedRFs[, hash := NULL]
        
        if (is.null(MS2QRes))
        {
            MD <- loadCacheData("MS2QMD", baseHash)
            if (is.null(MD))
            {
                warning("Could not find cached calibration data! You may have an old cache file. ",
                        "Please clear any cached data, eg by running: clearCache(\"RF_SMILES\")", call. = FALSE)
                MD <- list()
            }
            MS2QRes <- list(RFs = cachedRFs, MD = MD)
        }
        else
        {
            MS2QRes$RFs <- rbind(MS2QRes$RFs, cachedRFs)
            # restore order
            boundHashes <- c(hashesTODO, names(cachedData))
            MS2QRes$RFs <- MS2QRes$RFs[match(hashes, boundHashes)]
        }
    }
    
    if (!is.null(MS2QRes$RFs[["id"]]))
        MS2QRes$RFs[, id := NULL]

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    # NOTE: need to take the inverse before conversion
    MS2QRes$RFs[, RF_SIRFP := 1/convertConc(1/RF_SIRFP[1], "M", concUnit, formulaMW(neutral_formula[1])),
                by = "neutral_formula"][]
    
    return(MS2QRes)
}

predictLC50SIRFPs <- function(featAnnSIR, LC50Mode, concUnit)
{
    # UNDONE: check supported adducts?
    
    allFPs <- getMS2QTFPs(featAnnSIR)
    if (nrow(allFPs) == 0)
        return(data.table())

    baseHash <- makeHash(LC50Mode)
    hashes <- sapply(split(allFPs[, -c("id", "group")], seq_len(nrow(allFPs))), function(s) makeHash(baseHash, s))
    
    cachedData <- pruneList(loadCacheData("LC50_SIRFP", hashes, simplify = FALSE))
    indsTODO <- which(!hashes %in% names(cachedData))
    hashesTODO <- hashes[indsTODO]
    
    LC50s <- NULL
    if (length(indsTODO) > 0)
    {
        allFPsTODO <- allFPs[indsTODO]
        allFPsTODO[, exactMass := sapply(neutral_formula, getFormulaMass)]
        suppressMessages(utils::capture.output(LC50s <- MS2Tox::FishLC50Prediction(allFPsTODO, LC50Mode)))
    
        LC50s <- merge(allFPsTODO[, c("group", "neutral_formula", "id"), with = FALSE],
                       LC50s[, c("id", "LC50_predicted")], by = "id", sort = FALSE)
        setnames(LC50s, "LC50_predicted", "LC50_SIRFP")
        LC50s[, id := NULL]
        
        for (i in seq_len(nrow(LC50s)))
            saveCacheData("LC50_SIRFP", LC50s$LC50_SIRFP[i], hashesTODO[i])
    }
    
    if (length(cachedData) > 0)
    {
        cachedLC50s <- rbindlist(lapply(cachedData, function(cd) data.table(LC50_SIRFP = cd)), idcol = "hash")
        cachedLC50s[, neutral_formula := allFPs$neutral_formula[match(hash, hashes)]]
        cachedLC50s[, group := allFPs$group[match(hash, hashes)]]
        cachedLC50s[, hash := NULL]
        
        if (is.null(LC50s))
            LC50s <- cachedLC50s
        else
        {
            LC50s <- rbind(LC50s, cachedLC50s)
            # restore order
            boundHashes <- c(hashesTODO, names(cachedData))
            LC50s <- LC50s[match(hashes, boundHashes)]
        }
    }

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    # NOTE: need to take the inverse before conversion
    LC50s[, LC50_SIRFP := convertConc(LC50_SIRFP[1], "log10 mM", concUnit, formulaMW(neutral_formula[1])),
        by = "neutral_formula"]

    return(LC50s[])
}

syncSIRFPs <- function(obj)
{
    # sync fingerprints
    obj@fingerprints <- obj@fingerprints[names(obj@fingerprints) %chin% groupNames(obj)]
    obj@fingerprints <- pruneList(Map(obj@fingerprints, obj@groupAnnotations[names(obj@fingerprints)], f = function(fp, ann)
    {
        fpForms <- intersect(names(fp), ann$neutral_formula)
        fp <- fp[, c(fpForms, "absoluteIndex"), with = FALSE]
        return(if (ncol(fp) == 1) NULL else fp) # nullify if no candidates left
    }))
    return(obj)
}
