# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

makeAnaInfoR <- function(anaInfo, varName = NULL)
{
    # NOTE: for constructive we convert to a data.frame to avoid need for data.table and allow read.table() to
    # construct the table
    anaInfo <- as.data.frame(anaInfo)
    
    # UNDONE: read.table() looks nicer, but need to see for larger tables. Make it optional?
    
    if (is.null(varName))
        return(constructive::construct(anaInfo, constructive::opts_data.frame("read.table"))$code)
    
    return(constructive::construct_multi(setNames(list(anaInfo), varName),
                                         constructive::opts_data.frame("read.table"))$code)
}

scriptGenerator <- setRefClass("scriptGenerator",
                               fields = list(settings = "list", txtCon = "ANY"))
scriptGenerator$methods(
    initialize = function(...)
    {
        txtCon <<- textConnection(NULL, "w")
        return(callSuper(...))
    },
    
    finalize = function()
    {
        if (!is.null(txtCon))
        {
            close(txtCon)
            txtCon <<- NULL
        }
    },
    
    getCode = function() paste0(textConnectionValue(txtCon), collapse = "\n"),
    
    doQuote = function(txt) paste0("\"", txt, "\""),
    
    addText = function(txt) writeLines(txt, txtCon),
    
    addNL = function(count = 1) addText(rep("", count)),
    
    addComment = function(txt, condition = TRUE) { if (condition) addText(paste("#", txt)) },
    
    addHeader = function(title)
    {
        hd <- strrep("-", 25)
        addNL()
        addComment(c(hd, title, hd))
        addNL()
    },
    
    addAssignment = function(var, val, quote = FALSE) addText(paste(var, "<-", if (quote) doQuote(val) else val)),
    
    addCall = function(var, func, args, condition = TRUE, indent = 0)
    {
        if (!condition)
            return(NULL)
        
        if (!is.list(args[[1]]))
            args <- list(args) # HACK: allow unnested lists for cases with one arg
        
        args <- lapply(args, function(a)
        {
            # handle args
            #   - add names if necessary
            #   - handle conditional args
            #   - handle strings
            #   - handle conditional NULLs
            
            if (!is.null(a[["condition"]]) && !a$condition)
                return(NULL)
            
            ret <- a$value
            if (is.null(a[["value"]]) || (!is.null(a[["isNULL"]]) && a$isNULL) ||
                (!is.null(a[["zeroToNULL"]]) && a$zeroToNULL && a$value == 0))
                ret <- "NULL"
            else if (!is.null(a[["quote"]]) && a$quote)
                ret <- doQuote(ret)
            
            if (length(ret) > 1)
                ret <- paste0("c(", paste0(ret, collapse = ", "), ")")
            
            if (!is.null(a[["name"]]))
                ret <- paste(a$name, "=", ret)
            
            return(ret)
        })
        args <- pruneList(args)
        
        argText <- paste0(args, collapse = ", ")
        callPrefix <- paste0(strrep(" ", indent), if (!is.null(var)) paste0(var, " <- ") else "", func, "(")
        singleLineLength <- nchar(callPrefix) + nchar(argText) + nchar(")")
        
        if (singleLineLength > 120)
        {
            argStart <- nchar(callPrefix)
            argPos <- argStart
            argText <- ""
            for (i in seq_along(args))
            {
                a <- args[[i]]
                isLast <- i == length(args)
                if (!isLast)
                    a <- paste0(a, ", ")
                na <- nchar(a)
                if (isLast)
                    na <- na + 1 # +1: for trailing ")" NOTE: the actual addition of ")" is done below
                
                if (argPos > argStart && (argPos + na) > 120)
                {
                    # add a newline and indent if adding the next argument would make the line longer than 120 chars
                    # NOTE: this will not work if the first argument is longer than 120 chars, ignoring that...
                    # NOTE: use trimws() as there might be a trailing whitespace of the previous arg from the addition
                    # of ", " above
                    argText <- paste0(trimws(argText, "right", " "), "\n", strrep(" ", argStart), a)
                    argPos <- argStart + na
                }
                else
                {
                    argText <- paste0(argText, a)
                    argPos <- argPos + na
                }
            }
        }
        
        cl <- paste0(callPrefix, argText, ")")
        addText(cl)
    },
    
    addLoadSuspCall = function(var, file)
    {
        addCall(var, "read.csv", list(
            list(value = file, quote = TRUE)
        ))
    },

    addLoadSuspList = function(ionization, suspects, varBase)
    {
        if (ionization != "both")
            addLoadSuspCall(varBase, suspects$single)
        else
        {
            varPos <- paste0(varBase, "Pos"); varNeg <- paste0(varBase, "Neg")
            addLoadSuspCall(varPos, suspects$sets$pos)
            if (nzchar(suspects$sets$neg))
                addLoadSuspCall(varNeg, suspects$sets$neg)
            else
                addAssignment(varNeg, paste0(varBase, "Pos"))
        }
    },
    
    addSuspListEx = function(ionization, exPos, exNeg, varBase)
    {
        if (ionization == "positive")
            addAssignment(varBase, exPos)
        else if (ionization == "negative")
            addAssignment(varBase, exNeg)
        else
        {
            addAssignment(paste0(varBase, "Pos"), exPos)
            addAssignment(paste0(varBase, "Neg"), exNeg)
        }
    },
    
    addScreenCall = function(susp, adductArg, oh = TRUE, am = NULL)
    {
        addCall("fGroups", "screenSuspects", list(
            list(value = "fGroups"),
            list(value = susp),
            list(name = "rtWindow", value = 12),
            list(name = "mzWindow", value = 0.005),
            adductArg,
            list(name = "onlyHits", value = oh),
            list(name = "amend", value = am, condition = !is.null(am))
        ))
    }
)

genScriptInitBlock <- function(CCSCalibrant, anaInfoData, settingsGen, settingsAna, settingsPre, generator)
{
    addAnaInfo <- function(anaInfoVarName, aid, comment, exPol)
    {
        if (settingsAna$generateAnaInfo == "table")
        {
            if (settingsAna$analysisTableFileType == "embedded")
            {
                generator$addComment("Create analysis table", condition = comment)
                generator$addText(makeAnaInfoR(aid$tab, anaInfoVarName))
            }
            else
            {
                generator$addComment("Load analysis table", condition = comment)
                if (settingsAna$analysisTableFileType == "CSV")
                    generator$addCall(anaInfoVarName, "read.csv", list(value = aid$file, quote = TRUE))
                else # R: NOTE: cannot use addCall() as there is a nested function call
                    generator$addAssignment(anaInfoVarName, sprintf("eval(parse(\"%s\"))", aid$file))
            }
        }
        else if (settingsAna$generateAnaInfo == "dynamic")
        {
            generator$addCall(anaInfoVarName, "generateAnalysisInfo", list(
                list(name = "fromRaw", value = aid$fromRaw, quote = TRUE),
                list(name = "fromCentroid", value = aid$fromCentroid, quote = TRUE),
                list(name = "fromProfile", value = aid$fromProfile, quote = TRUE),
                list(name = "fromIMS", value = aid$fromIMS, quote = TRUE),
                list(name = "convCentroid", value = settingsPre$output$centroid, quote = TRUE),
                list(name = "convProfile", value = settingsPre$output$profile, quote = TRUE),
                list(name = "convIMS", value = settingsPre$output$ims, quote = TRUE)
            ))
        }
        else if (settingsAna$generateAnaInfo == "example")
        {
            pd <- if (settingsGen$IMS$mode != "none") "patRoonDataIMS" else "patRoonData"
            generator$addComment(sprintf("Example data from %s package (triplicate solvent blank + triplicate standard)",
                                        pd), condition = comment)
            generator$addCall(anaInfoVarName, sprintf("%s::exampleAnalysisInfo", pd), list(value = exPol, quote = TRUE))
        }
        else # none
        {
            generator$addComment("NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.",
                                 condition = comment)
            generator$addCall(anaInfoVarName, "data.frame", list(
                list(name = "path_centroid", value = "character()"),
                list(name = "analysis", value = "character()"),
                list(name = "replicate", value = "character()"),
                list(name = "blank", value = "character()")
            ))
        }
    }
    addPrepBlock <- function(anaInfoVarName, DAMethodVarName)
    {
        generator$addCall(NULL, "setDAMethod", list(
            list(value = anaInfoVarName),
            list(value = settingsPre$brukerCalib[[DAMethodVarName]], quote = TRUE)
        ), condition = settingsPre$brukerCalib$enabled && nzchar(settingsPre$brukerCalib[[DAMethodVarName]]), indent = 4)
        generator$addCall(NULL, "recalibrarateDAFiles", list(value = anaInfoVarName),
                          condition = settingsPre$brukerCalib$enabled, indent = 4)
        itf <- splitConvTypeFormat(settingsPre$steps$from); otf <- splitConvTypeFormat(settingsPre$steps$to)
        for (i in seq_len(nrow(settingsPre$steps)))
        {
            generator$addCall(NULL, "convertMSFiles", list(
                list(value = anaInfoVarName),
                list(name = "typeFrom", value = itf[[i]][1], quote = TRUE),
                list(name = "formatFrom", value = itf[[i]][2], quote = TRUE),
                list(name = "typeTo", value = otf[[i]][1], quote = TRUE),
                list(name = "formatTo", value = otf[[i]][2], quote = TRUE),
                list(name = "algorithm", value = settingsPre$steps$algorithm[i], quote = TRUE),
                list(name = "overwrite", value = FALSE)
            ), indent = 4)
        }
    }
    
    generator$addHeader("initialization")
    
    generator$addAssignment("workPath", settingsGen$destination, quote = TRUE)
    generator$addCall(NULL, "setwd", list(value = "workPath"))
    generator$addNL()
    
    if (settingsGen$ionization != "both")
        addAnaInfo("anaInfo", anaInfoData, TRUE, settingsGen$ionization)
    else
    {
        addAnaInfo("anaInfoPos", anaInfoData$positive, TRUE, "positive")
        addAnaInfo("anaInfoNeg", anaInfoData$negative, FALSE, "negative")
    }
    
    if (nrow(settingsPre$steps) > 0 || settingsPre$brukerCalib$enabled)
    {
        generator$addNL()
        generator$addComment("Set to FALSE to skip data pre-treatment")
        generator$addAssignment("doDataPretreatment", TRUE)
        generator$addText("if (doDataPretreatment)\n{")
        if (settingsGen$ionization != "both")
            addPrepBlock("anaInfo", "method")
        else
        {
            addPrepBlock("anaInfoPos", "methodPos")
            generator$addNL()
            addPrepBlock("anaInfoNeg", "methodNeg")
        }
        generator$addText("}")
    }
    
    if (settingsGen$IMS$mode != "none" && settingsGen$IMS$CCSMethod != "none")
    {
        generator$addNL()
        generator$addCall("CCSParams", "getCCSParams", list(
            list(name = "method", value = settingsGen$IMS$CCSMethod, quote = TRUE),
            list(name = "calibrant", value = CCSCalibrant, quote = TRUE,
                 condition = settingsGen$IMS$CCSMethod == "agilent")
        ))
    }
}

genScriptSuspListsBlock <- function(ionization, IMSMode, settingsFeat, doSusps, doPiekSusps, doISTDs, generator)
{
    generator$addHeader("setup suspect lists")
    
    pd <- if (IMSMode != "none") "patRoonDataIMS" else "patRoonData"
    
    if (doSusps)
    {
        if (settingsFeat$exSuspList)
        {
            generator$addComment("Get example suspect list")
            generator$addSuspListEx(ionization, sprintf("%s::suspectsPos", pd), sprintf("%s::suspectsNeg", pd), "suspList")
        }
        else
        {
            generator$addComment("Load suspect list")
            generator$addLoadSuspList(ionization, settingsFeat$suspects, "suspList")
        }
        
        if (IMSMode != "none" && settingsFeat$IMSSuspCCSPred != "none")
        {
            # UNDONE: also call if IMSSuspCCSPred == "none"?
            
            addACall <- function(varName, add)
            {
                generator$addCall(varName, "assignMobilities", list(
                    list(value = varName),
                    list(name = "from", value = settingsFeat$IMSSuspCCSPred, quote = TRUE),
                    list(name = "CCSParams", value = "CCSParams"),
                    list(name = "adducts", value = c(paste0('"', add, '"'), NA)),
                    list(name = "overwrite", value = FALSE)
                ))
            }
            
            generator$addNL()
            generator$addComment("Add mobility and CCS to suspect list(s)")
            if (ionization != "both")
                addACall("suspList", if (ionization == "positive") "[M+H]+" else "[M-H]-")
            else
            {
                addACall("suspListPos", "[M+H]+")
                addACall("suspListNeg", "[M-H]-")
            }
        }
    }
    
    if (doPiekSusps)
    {
        if (doSusps)
            generator$addNL()
        generator$addComment("Init suspect list for piek")
        
        copySusps <- (ionization != "both" && !nzchar(settingsFeat$piekParams$suspects$single)) ||
            (ionization == "both" && !nzchar(settingsFeat$piekParams$suspects$sets$pos))
        if (copySusps)
        {
            if (ionization != "both")
                generator$addAssignment("suspListPiek", "suspList")
            else
            {
                generator$addAssignment("suspListPiekPos", "suspListPos")
                generator$addAssignment("suspListPiekNeg", "suspListNeg")
            }
        }
        else
            generator$addLoadSuspList(ionization, settingsFeat$piekParams$suspects, "suspListPiek")
    }
    
    if (doISTDs)
    {
        if (doSusps || doPiekSusps)
            generator$addNL()
        
        doEx <- (ionization != "both" && !nzchar(settingsFeat$fGroupsAdv$ISTDLists$single)) ||
            (ionization == "both" && !nzchar(settingsFeat$fGroupsAdv$ISTDLists$sets$pos))
        if (doEx)
        {
            generator$addComment("Get example ISTD list")
            generator$addSuspListEx(ionization, sprintf("%s::ISTDListPos", pd), sprintf("%s::ISTDListNeg", pd), "ISTDList")
        }
        else
        {
            generator$addComment("Load ISTD list")
            generator$addLoadSuspList(ionization, settingsFeat$fGroupsAdv$ISTDLists, "ISTDList")
        }
    }
}

genScriptFeaturesBlock <- function(ionization, IMS, settingsFeat, generator)
{
    addFindFeatures <- function(varName, anaInfoVarName, suspPiekL, suspPiekA)
    {
        peakParams <- if (settingsFeat$featAlgo == "piek")
            sprintf("getDefPeakParams(\"chrom\", \"%s\")", settingsFeat$piekParams$peaksAlgo)
        
        fa <- settingsFeat$featAlgo
        generator$addCall(varName, "findFeatures", list(
            list(value = anaInfoVarName),
            list(value = if (fa == "XCMS") "xcms3" else tolower(fa), quote = TRUE),
            list(name = "noiseThrInt", value = 1000, condition = fa == "OpenMS"),
            list(name = "chromSNR", value = 3, condition = fa == "OpenMS"),
            list(name = "chromFWHM", value = 5, condition = fa == "OpenMS"),
            list(name = "minFWHM", value = 1, condition = fa == "OpenMS"),
            list(name = "maxFWHM", value = 30, condition = fa == "OpenMS"),
            list(name = "kmeans", value = TRUE, condition = fa == "KPIC2"),
            list(name = "level", value = 1000, condition = fa == "KPIC2"),
            list(name = "IMS", value = IMS$mode == "direct", condition = fa == "piek" && IMS$mode != "none"),
            list(name = "genEICParams", value = "genEICParams", condition = fa == "piek"),
            list(name = "peakParams", value = peakParams, condition = fa == "piek"),
            list(name = "suspects", value = suspPiekL,
                 condition = fa == "piek" && settingsFeat$piekParams$filter == "suspects"),
            list(name = "adduct", value = suspPiekA, quote = TRUE,
                 condition = fa == "piek" && settingsFeat$piekParams$filter == "suspects"),
            list(name = "doFMF", value = TRUE, condition = fa == "Bruker")
        ))
    }
    
    generator$addHeader("features")
    
    if (settingsFeat$fGroupsAlgo != "SIRIUS") # NOTE: never the case with sets
    {
        generator$addComment("Find all features")
        generator$addComment("NOTE: see the reference manual for many more options",
                             condition = !settingsFeat$featAlgo %in% c("Bruker", "SIRIUS"))
        
        if (settingsFeat$featAlgo == "piek")
        {
            mmz <- settingsFeat$piekParams$filter; mims <- settingsFeat$piekParams$filterIMS
            def <- list(
                none = getPiekEICParams(),
                susp = getPiekEICParams(filter = "suspects", filterIMS = "suspects"),
                ms2 = getPiekEICParams(filter = "ms2", filterIMS = "ms2")
            )
            
            doDirectIMS <- IMS$mode == "direct"
            generator$addCall("genEICParams", "getPiekEICParams", list(
                list(name = "filter", value = mmz, quote = TRUE),
                list(name = "filterIMS", value = mims, quote = TRUE, condition = doDirectIMS),
                list(name = "mzRange", value = def$none$mzRange),
                list(name = "mzStep", value = def$none$mzStep),
                list(name = "mobRange", value = def$none$mobRange, condition = doDirectIMS),
                list(name = "mobStep", value = def$none$mobStep, condition = doDirectIMS),
                list(name = "rtWindow", value = def$susp$rtWindow, condition = mmz == "suspects"),
                list(name = "mzWindow", value = def$susp$mzWindow, condition = mmz == "suspects"),
                list(name = "IMSWindow", value = def$susp$IMSWindow, condition = mims == "suspects" && doDirectIMS),
                list(name = "rtWindow", value = def$ms2$rtWindow, condition = mmz == "ms2"),
                list(name = "mzWindow", value = def$ms2$mzWindow, condition = mmz == "ms2"),
                list(name = "IMSWindow", value = def$ms2$IMSWindow, condition = mims == "ms2" && doDirectIMS),
                list(name = "minTIC", value = def$ms2$minTIC, condition = mmz == "ms2")
            ))
        }
        
        if (ionization != "both")
            addFindFeatures("fList", "anaInfo", "suspListPiek", if (ionization == "positive") "[M+H]+" else "[M-H]-")
        else
        {
            addFindFeatures("fListPos", "anaInfoPos", "suspListPiekPos", "[M+H]+")
            addFindFeatures("fListNeg", "anaInfoNeg", "suspListPiekNeg", "[M-H]-")
            generator$addCall("fList", "makeSet", list(
                list(value = "fListPos"),
                list(value = "fListNeg"),
                list(name = "adducts", value = c("[M+H]+", "[M-H]-"), quote = TRUE)
            ))
        }
        generator$addNL()
    }
    
    generator$addComment("Group and align features between analyses")
    doRTAlign <- !settingsFeat$fGroupsAlgo %in% c("SIRIUS", "Greedy") &&
        (ionization != "both" || settingsFeat$fGroupsAlgo == "OpenMS")
    generator$addCall("fGroups", "groupFeatures", list(
        list(value = "fList", condition = settingsFeat$fGroupsAlgo != "SIRIUS"),
        list(value = "anaInfo", condition = settingsFeat$fGroupsAlgo == "SIRIUS"),
        list(value = if (settingsFeat$fGroupsAlgo == "XCMS") "xcms3" else tolower(settingsFeat$fGroupsAlgo), quote = TRUE),
        list(name = "rtalign", value = TRUE, condition = doRTAlign),
        list(name = "groupParam", value = "xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$replicate)",
             condition = settingsFeat$fGroupsAlgo == "XCMS"),
        list(name = "retAlignParam", value = "xcms::ObiwarpParam()", condition = doRTAlign && settingsFeat$fGroupsAlgo == "XCMS"),
        list(name = "scoreWeights", value = "c(retention = 1, mz = 1, mobility = 1)", condition = settingsFeat$fGroupsAlgo == "Greedy")
    ))
    
    retRange <- settingsFeat$fGroupsAdv$retention
    if (all(retRange == 0))
        retRange <- NULL
    else if (retRange[2] == 0)
        retRange[2] <- Inf
    mzRange <- settingsFeat$fGroupsAdv$mz
    if (all(mzRange == 0))
        mzRange <- NULL
    else if (mzRange[2] == 0)
        mzRange[2] <- Inf
    generator$addNL()
    generator$addComment("Basic rule based filtering")
    generator$addCall("fGroups", "filter", list(
        list(value = "fGroups"),
        list(name = "preAbsMinIntensity", value = settingsFeat$fGroupsAdv$preIntThr, zeroToNULL = TRUE),
        list(name = "absMinIntensity", value = settingsFeat$fGroupsAdv$intThr, zeroToNULL = TRUE),
        list(name = "relMinReplicateAbundance", value = settingsFeat$fGroupsAdv$repAbundance, zeroToNULL = TRUE),
        list(name = "maxReplicateIntRSD", value = settingsFeat$fGroupsAdv$maxRepRSD, zeroToNULL = TRUE),
        list(name = "blankThreshold", value = settingsFeat$fGroupsAdv$blankThr, zeroToNULL = TRUE),
        list(name = "removeBlanks", value = settingsFeat$fGroupsAdv$removeBlanks),
        list(name = "retentionRange", value = retRange),
        list(name = "mzRange", value = mzRange)
    ))
    
    generator$addNL()
    generator$addComment("Update group properties")
    generator$addCall("fGroups", "updateGroups", list(
        list(value = "fGroups"),
        list(name = "what", value = c("ret", "mz", "mobility"), quote = TRUE),
        list(name = "intWeight", value = FALSE)
    ))
    
    # NOTE: for direct we only want to add CCS calculations, for post there is also the mobility assignment, which can
    # be done better later
    if (IMS$mode == "direct" && IMS$CCSMethod != "none")
    {
        generator$addNL()
        generator$addComment("Assign CCS values")
        generator$addCall("fGroups", "assignMobilities", list(
            list(value = "fGroups"),
            list(name = "CCSParams", value = "CCSParams")
        ))
    }
}

genScriptComponBlock <- function(ionization, settingsAnnon, generator)
{
    generator$addHeader("componentization")
    
    generator$addComment("Perform automatic generation of components")
    generator$addCall("components", "generateComponents", list(
        list(value = "fGroups"),
        list(value = tolower(settingsAnnon$componAlgo), quote = TRUE),
        list(name = "ionization", value = ionization, quote = TRUE, condition = ionization != "both"),
        list(name = "rtRange", value = c(-120, 120), condition = settingsAnnon$componAlgo == "nontarget"),
        list(name = "mzRange", value = c(5, 120), condition = settingsAnnon$componAlgo == "nontarget"),
        list(name = "elements", value = c("C", "H", "O"), quote = TRUE, condition = settingsAnnon$componAlgo == "nontarget"),
        list(name = "rtDev", value = defaultLim("retention", "wide"), condition = settingsAnnon$componAlgo == "nontarget"),
        list(name = "absMzDev", value = 0.002, condition = settingsAnnon$componAlgo == "nontarget")
    ))
    
    if (settingsAnnon$selectIons && settingsAnnon$componAlgo != "nontarget")
    {
        pa <- switch(ionization,
                     positive = "[M+H]+",
                     negative = "[M-H]-",
                     both = c("[M+H]+", "[M-H]-"))
        generator$addCall("fGroups", "selectIons", list(
            list(value = "fGroups"),
            list(value = "components"),
            list(name = "prefAdduct", value = pa, quote = TRUE),
            list(name = "onlyMonoIso", value = TRUE)
        ))
    }
}

genScriptFGNormBlock <- function(ionization, settingsFeat, adductArg, generator)
{
    adductArg$condition <- adductArg$condition && settingsFeat$fGroupsAdv$featNorm == "istd"
    
    generator$addHeader("normalization")
    
    if (settingsFeat$fGroupsAdv$featNorm == "istd")
    {
        generator$addNL()
        if (ionization != "both")
            generator$addComment("Set adduct to NULL if ISTD list contains an adduct column")
        
        standards <- if (ionization != "both") "ISTDList" else c("ISTDListPos", "ISTDListNeg")
    }
    else
        standards <- NULL
    
    generator$addCall("fGroups", "normInts", list(
        list(value = "fGroups"),
        list(name = "featNorm", value = settingsFeat$fGroupsAdv$featNorm, quote = TRUE),
        list(name = "groupNorm", value = settingsFeat$fGroupsAdv$groupNorm),
        list(name = "normFunc", value = "max"),
        list(name = "standards", value = standards, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        list(name = "ISTDRTWindow", value = 120, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        list(name = "ISTDMZWindow", value = 300, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        list(name = "minISTDs", value = 3, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        list(name = "rtWindow", value = 12, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        list(name = "mzWindow", value = 0.005, condition = settingsFeat$fGroupsAdv$featNorm == "istd"),
        adductArg
    ))
}

genScriptScreenBlock <- function(ionization, settingsFeat, adductArg, amendScrForTPs, generator)
{
    generator$addHeader("suspect screening")
    
    if (amendScrForTPs)
        generator$addComment("NOTE: onlyHits is set to FALSE to ensure TPs can be found below")
    else
        generator$addComment("Set onlyHits to FALSE to retain features without suspects (eg for full NTA)")
    if (ionization != "both")
        generator$addComment("Set adduct to NULL if suspect list contains an adduct column")
    
    if (ionization != "both")
        generator$addScreenCall("suspList", adductArg, !amendScrForTPs)
    else
        generator$addScreenCall("list(suspListPos, suspListNeg)", adductArg, !amendScrForTPs)
}

genScriptTPInitBlock <- function(ionization, settingsTP, adductArg, amendScrForTPs, generator)
{
    adductArg$condition <- adductArg$condition && settingsTP$TPsAlgo == "logic"
    
    generator$addHeader("transformation products")
    
    if (settingsTP$TPsAlgo != "logic" && settingsTP$TPGenInput == "suspects")
    {
        generator$addComment("Load parent suspect list")
        generator$addLoadSuspCall("suspListParents", settingsTP$TPSuspectList)
        generator$addNL()
    }
    
    generator$addComment("Obtain TPs")
    generator$addCall("TPs", "generateTPs", list(
        list(value = settingsTP$TPsAlgo, quote = TRUE),
        list(name = "parents", value = switch(settingsTP$TPGenInput,
                                              suspects = "suspListParents",
                                              screening = "fGroups",
                                              all = "NULL"),
             condition = settingsTP$TPsAlgo != "Logic"),
        list(value = "fGroups", condition = settingsTP$TPsAlgo == "logic"),
        adductArg,
        list(name = "type", value = "env", quote = TRUE, condition = settingsTP$TPsAlgo == "biotransformer"),
        list(name = "transLibrary", value = "photolysis_ranked", quote = TRUE, condition = settingsTP$TPsAlgo == "cts"),
        list(name = "generations", value = if (settingsTP$TPsAlgo == "biotransformer") 2 else 1,
             condition = settingsTP$TPsAlgo %in% c("biotransformer", "cts", "lobrary")),
        list(name = "formulas", value = "formulas", condition = settingsTP$TPsAlgo == "ann_form"),
        list(name = "compounds", value = "compounds", condition = settingsTP$TPsAlgo == "ann_comp"),
        list(name = "TPsRef", value = "NULL", condition = settingsTP$TPsAlgo == "ann_comp"),
        list(name = "fGroupsComps", value = "fGroups", condition = settingsTP$TPsAlgo == "ann_comp"),
        list(name = "TPStructParams", value = "getDefTPStructParams()",
             condition = settingsTP$TPsAlgo %in% c("biotransformer", "cts", "library", "ann_comp"))
    ))
    
    if (!settingsTP$TPsAlgo %in% c("ann_form", "ann_comp"))
    {
        generator$addNL()
        generator$addComment("Screen TPs")
        generator$addCall("suspListTPs", "convertToSuspects", list(
            list(value = "TPs"),
            list(name = "includeParents", value = FALSE)
        ))
        if (amendScrForTPs)
            generator$addComment("Amend with TP screening")
        generator$addScreenCall("suspListTPs", adductArg, am = if (amendScrForTPs) TRUE else NULL)
    }
}

genScriptFeatMobBlock <- function(IMS, settingsFeat, doSusps, generator)
{
    defIMM <- getIMSMatchParams("mobility")
    
    generator$addHeader("mobility assignment")
    generator$addCall("fGroups", "assignMobilities", list(
        list(value = "fGroups"),
        list(name = "mobPeakParams", value = sprintf("getDefPeakParams(type = \"%s_ims\", algorithm = \"%s\")",
                                                    IMS$limits, settingsFeat$IMSPeaksMob)),
        list(name = "chromPeakParams", value = sprintf("getDefPeakParams(type = \"chrom\", algorithm = \"%s\")",
                                                     settingsFeat$IMSPeaksChrom)),
        list(name = "fallbackEIC", value = TRUE),
        list(name = "calcArea", value = "integrate", quote = TRUE),
        list(name = "CCSParams", value = "CCSParams", condition = IMS$CCSMethod != "none"),
        list(name = "fromSuspects", value = FALSE, condition = doSusps),
        list(name = "IMSMatchParams",
             value = sprintf("getIMSMatchParams(\"mobility\", relative = %s, window = %.2f)", defIMM$relative,
                             defIMM$window), condition = doSusps)
    ))
}

genScriptAnnBlock <- function(ionization, IMS, settingsAnn, adductArg, doSusps, addMFDB, generator)
{
    doForms <- nzchar(settingsAnn$formulasAlgo); doComps <- nzchar(settingsAnn$compoundsAlgo)
    
    generator$addHeader("annotation")
    
    generator$addComment("Retrieve MS peak lists")
    generator$addCall("avgMSListParams", "getDefAvgPListParams", list(name = "clusterMzWindow", value = defaultLim("mz", "narrow")))
    generator$addCall("mslists", "generateMSPeakLists", list(
        list(value = "fGroups"),
        list(name = "maxMSRtWindow", value = defaultLim("retention", "narrow")),
        list(name = "avgFeatParams", value = "avgMSListParams"),
        list(name = "avgFGroupParams", value = "avgMSListParams")
    ))
    generator$addComment("Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.")
    generator$addCall("mslists", "filter", list(
        list(value = "mslists"),
        list(name = "MSLevel", value = 2),
        list(name = "absMinIntensity", value = "NULL"),
        list(name = "relMinIntensity", value = 0.05),
        list(name = "topMostPeaks", value = 25),
        list(name = "maxMZOverPrec", value = 4)
    ))
    
    if (doForms)
    {
        generator$addNL()
        generator$addComment("Calculate formula candidates")
        generator$addCall("formulas", "generateFormulas", list(
            list(value = "fGroups"),
            list(value = "mslists"),
            list(value = tolower(settingsAnn$formulasAlgo), quote = TRUE),
            list(name = "relMzDev", value = 5, condition = settingsAnn$formulasAlgo != "Bruker"),
            list(name = "precursorMzSearchWindow", value = defaultLim("mz", "narrow"),
                 condition = settingsAnn$formulasAlgo == "Bruker"),
            adductArg,
            list(name = "elements", value = "CHNOP", quote = TRUE, condition = settingsAnn$formulasAlgo != "Bruker"),
            list(name = "oc", value = FALSE, condition = settingsAnn$formulasAlgo == "GenForm"),
            list(name = "profile", value = "qtof", quote = TRUE, condition = settingsAnn$formulasAlgo == "SIRIUS"),
            list(name = "calculateFeatures", value = "FALSE", condition = settingsAnn$formulasAlgo != "Bruker"),
            list(name = "featThresholdAnn", value = 0.75),
            list(name = "setThresholdAnn", value = 0, condition = ionization == "both")
        ))
        if ("formulas" %in% settingsAnn$estIDConf)
        {
            generator$addCall("formulas", "estimateIDConfidence", list(
                list(value = "formulas"),
                list(name = "IDFile", value = "idlevelrules.yml", quote = TRUE)
            ))
        }
    }
    
    if (doComps)
    {
        generator$addNL()
        
        if (settingsAnn$compoundsAlgo == "Library")
        {
            generator$addComment("Load MS library. You may want to filter it, please see the manuals for more details.")
            generator$addCall("mslibrary", "loadMSLibrary", list(
                list(value = settingsAnn$MSLibraryPath, quote = TRUE),
                list(value = settingsAnn$MSLibraryFormat, quote = TRUE),
                list(name = "absMzDev", value = 0.002)
            ))
        }
        
        addMFDB <- settingsAnn$compoundsAlgo == "MetFrag" && addMFDB
        if (addMFDB)
        {
            generator$addCall(NULL, "convertToMFDB", list(
                list(value = "TPs"),
                list(value = "TP-database.csv", quote = TRUE),
                list(name = "includeParents", value = TRUE)))
        }
        
        generator$addComment("Calculate compound structure candidates")
        
        if (settingsAnn$compoundsAlgo == "SIRIUS")
            generator$addComment("Please see the handbook for SIRIUS login options. If you want to disable automatic login set login=FALSE")
        
        generator$addCall("compounds", "generateCompounds", list(
            list(value = "fGroups"),
            list(value = "mslists"),
            list(value = tolower(settingsAnn$compoundsAlgo), quote = TRUE),
            list(name = "dbRelMzDev", value = 5, condition = settingsAnn$compoundsAlgo == "MetFrag"),
            list(name = "fragRelMzDev", value = 5, condition = settingsAnn$compoundsAlgo == "MetFrag"),
            list(name = "fragAbsMzDev", value = 0.002, condition = settingsAnn$compoundsAlgo == "MetFrag"),
            list(name = "relMzDev", value = 5, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            adductArg,
            list(name = "database", value = "pubchemlite", quote = TRUE, condition = settingsAnn$compoundsAlgo == "MetFrag" && !addMFDB),
            list(name = "database", value = "csv", quote = TRUE, condition = addMFDB),
            list(name = "extraOpts", value = "list(LocalDatabasePath = \"TP-database.csv\")", condition = addMFDB),
            list(name = "maxCandidatesToStop", value = 2500, condition = settingsAnn$compoundsAlgo == "MetFrag"),
            list(name = "fingerIDDatabase", value = "pubchem", quote = TRUE, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            list(name = "elements", value = "CHNOP", quote = TRUE, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            list(name = "profile", value = "qtof", quote = TRUE, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            list(name = "login", value = "interactive", quote = TRUE, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            list(name = "alwaysLogin", value = FALSE, condition = settingsAnn$compoundsAlgo == "SIRIUS"),
            list(name = "MSLibrary", value = "mslibrary", condition = settingsAnn$compoundsAlgo == "Library"),
            list(name = "minSim", value = 0.75, condition = settingsAnn$compoundsAlgo == "Library"),
            list(name = "absMzDev", value = 0.002, condition = settingsAnn$compoundsAlgo == "Library"),
            list(name = "specSimParams", value = "getDefSpecSimParams()", condition = settingsAnn$compoundsAlgo == "Library"),
            list(name = "setThresholdAnn", value = 0, condition = ionization == "both")
        ))
        if (settingsAnn$compoundsAlgo == "MetFrag" && doForms)
        {
            generator$addCall("compounds", "addFormulaScoring", list(
                list(value = "compounds"),
                list(value = "formulas"),
                list(name = "updateScore", value = TRUE)
            ))
        }
        if (IMS$mode != "none")
        {
            doConv <- IMS$CCSMethod != "none"; doPred <- settingsAnn$compCCSPred != "none"
            if (doConv || doPred)
            {
                generator$addNL()
                if (doConv && doPred)
                    generator$addComment("Predict and convert mobility and CCS data")
                else if (doConv)
                    generator$addComment("Convert mobility and CCS data (if possible)")
                else
                    generator$addComment("Predict CCS data")
                generator$addCall("compounds", "assignMobilities", list(
                    list(value = "compounds"),
                    list(name = "from", value = settingsAnn$compCCSPred, quote = TRUE, condition = doPred),
                    list(name = "overwrite", value = FALSE),
                    adductArg,
                    list(name = "CCSParams", value = "CCSParams", condition = doConv)
                ))
            }
        }
        if ("compounds" %in% settingsAnn$estIDConf)
        {
            generator$addNL()
            generator$addCall("compounds", "estimateIDConfidence", list(
                list(value = "compounds"),
                list(name = "MSPeakLists", value = "mslists"),
                list(name = "formulas", value = "formulas", isNULL = !doForms),
                list(name = "IDFile", value = "idlevelrules.yml", quote = TRUE)
            ))
        }
    }
    
    if (doSusps && "suspects" %in% settingsAnn$estIDConf && (doForms || doComps))
    {
        generator$addNL()
        generator$addComment("Annotate suspects")
        generator$addCall("fGroups", "estimateIDConfidence", list(
            list(value = "fGroups"),
            list(name = "formulas", value = "formulas", isNULL = !doForms),
            list(name = "compounds", value = "compounds", isNULL = !doComps),
            list(name = "MSPeakLists", value = "mslists"),
            list(name = "IDFile", value = "idlevelrules.yml", quote = TRUE)
        ))
    }
}

genScriptTPCompBlock <- function(doFormAnn, doCompAnn, TPsAlgo, generator)
{
    generator$addHeader("Parent and TP linkage")
    
    generator$addComment("You probably want to prioritize the data before componentization. Please see the handbook for more info.")
    generator$addCall("componentsTP", "generateComponents", list(
        list(value = "fGroups"),
        list(value = "tp", quote = TRUE),
        list(name = "fGroupsTPs", value = "fGroups"),
        list(name = "TPs", value = "TPs"),
        list(name = "MSPeakLists", value = "mslists", condition = doFormAnn || doCompAnn),
        list(name = "formulas", value = "formulas", condition = doFormAnn),
        list(name = "compounds", value = "compounds", condition = doCompAnn)
    ))
    
    generator$addNL()
    generator$addComment("You may want to configure the filtering step below. See the manuals for more details.")
    generator$addCall("componentsTP", "filter", list(
        list(value = "componentsTP"),
        list(name = "retDirMatch", value = FALSE),
        list(name = "minSpecSim", value = "NULL"),
        list(name = "minSpecSimPrec", value = "NULL"),
        list(name = "minSpecSimBoth", value = "NULL"),
        list(name = "minFragMatches", value = "NULL"),
        list(name = "minNLMatches", value = "NULL"),
        list(name = "formulas", value = "formulas", isNULL = !doFormAnn, TPsAlgo == "logic")
    ))
    
    generator$addNL()
    generator$addComment("Only keep linked feature groups")
    generator$addAssignment("fGroups", "fGroups[results = componentsTP]")
}

genScriptReportBlock <- function(settingsAnn, settingsReport, doTPs, generator)
{
    # UNDONE: for now TP components are always reported instead of others
    componVal <- if (doTPs)
        "componentsTP"
    else if (nzchar(settingsAnn$componAlgo) && (!settingsAnn$selectIons || settingsAnn$componAlgo == "nontarget"))
        "components"
    else
        "NULL"
    doForms <- nzchar(settingsAnn$formulasAlgo); doComps <- nzchar(settingsAnn$compoundsAlgo)
    
    generator$addHeader("reporting")
    
    if ("HTML" %in% settingsReport$reportGen)
    {
        generator$addComment("Advanced report settings can be edited in the report.yml file.")
        generator$addCall(NULL, "report", list(
            list(value = "fGroups"),
            list(name = "MSPeakLists", value = "mslists", isNULL = !doForms && !doComps),
            list(name = "formulas", value = "formulas", isNULL = !doForms),
            list(name = "compounds", value = "compounds", isNULL = !doComps),
            list(name = "components", value = componVal),
            list(name = "TPs", value = "TPs", condition = doTPs),
            list(name = "settingsFile", value = "report.yml", quote = TRUE),
            list(name = "openReport", value = TRUE)
        ))
    }
    
    if ("legacy" %in% settingsReport$reportGen)
    {
        if ("HTML" %in% settingsReport$reportGen)
            generator$addNL()
        
        generator$addComment("Generate reports with legacy interface.")
        generator$addCall(NULL, "reportCSV", condition = "CSV" %in% settingsReport$reportLegacy, list(
            list(value = "fGroups"),
            list(name = "path", value = "report", quote = TRUE),
            list(name = "formulas", value = "formulas", isNULL = !doForms),
            list(name = "compounds", value = "compounds", isNULL = !doComps),
            list(name = "components", value = componVal)
        ))
        generator$addCall(NULL, "reportPDF", condition = "PDF" %in% settingsReport$reportLegacy, list(
            list(value = "fGroups"),
            list(name = "path", value = "report", quote = TRUE),
            list(name = "formulas", value = "formulas", isNULL = !doForms),
            list(name = "compounds", value = "compounds", isNULL = !doComps),
            list(name = "MSPeakLists", value = "mslists", condition = doForms || doComps),
            list(name = "components", value = componVal)
        ))
    }
}

getScriptCode <- function(CCSCalibrant, anaInfoData, settings, noDate)
{
    ionization <- settings$general$ionization
    IMSMode <- settings$general$IMS$mode
    doSusps <- settings$features$exSuspList || (ionization != "both" && nzchar(settings$features$suspects$single)) ||
        (ionization == "both" && nzchar(settings$features$suspects$sets$pos))
    doPiekSusps <- settings$features$featAlgo == "piek" && settings$features$piekParams$filter == "suspects"
    doFGNorm <- settings$features$fGroupsAdv$featNorm != "none" || settings$features$fGroupsAdv$groupNorm
    doISTDs <- doFGNorm && settings$features$fGroupsAdv$featNorm == "istd"
    doTPs <- nzchar(settings$TP$TPsAlgo)
    doTPsStructScr <- doTPs && settings$TP$TPsAlgo %in% c("biotransformer", "cts", "library")
    doTPsAnn <- doTPs && settings$TP$TPsAlgo %in% c("ann_form", "ann_comp")
    amendScrForTPs <- doSusps && doTPsStructScr && settings$TP$TPGenInput == "screening"
    
    # This will be passed to script generators for suspect screening, annotation etc. We do the conditional here, so we
    # don't need to pass all the conditional variables to the generator functions.
    adductArg <- list(
        name = "adduct", value = if (ionization == "positive") "[M+H]+" else "[M-H]-",
        quote = TRUE,
        condition = ionization != "both" &&
            (!nzchar(settings$annotation$componAlgo) || settings$annotation$componAlgo == "nontarget" || !settings$annotation$selectIons)
    )
    
    generator <- scriptGenerator$new()
    
    if (!noDate)
        generator$addComment(paste("Script automatically generated on", date()))
    generator$addNL()
    generator$addCall(NULL, "library", list(value = "patRoon"))
    
    genScriptInitBlock(CCSCalibrant, anaInfoData, settings$general, settings$analyses, settings$preTreatment, generator)

    if (doSusps || doPiekSusps || doISTDs)
        genScriptSuspListsBlock(ionization, IMSMode, settings$features, doSusps, doPiekSusps, doISTDs, generator)
        
    genScriptFeaturesBlock(ionization, settings$general$IMS, settings$features, generator)

    if (nzchar(settings$annotation$componAlgo))
        genScriptComponBlock(ionization, settings$annotation, generator)
    
    if (doFGNorm)
        genScriptFGNormBlock(ionization, settings$features, adductArg, generator)
    
    if (doSusps)
        genScriptScreenBlock(ionization, settings$features, adductArg, amendScrForTPs, generator)
    
    # NOTE: we do TPs from ann_form/ann_comp later as it needs feature annotations (and doesn't need screening of TP
    # suspects)
    if (doTPs && !doTPsAnn)
        genScriptTPInitBlock(ionization, settings$TP, adductArg, amendScrForTPs, generator)
    
    if (IMSMode == "post")
        genScriptFeatMobBlock(settings$general$IMS, settings$features, doSusps, generator)
    
    if (nzchar(settings$annotation$formulasAlgo) || nzchar(settings$annotation$compoundsAlgo))
    {
        addMFDB <- doTPsStructScr && settings$TP$TPDoMFDB
        genScriptAnnBlock(ionization, settings$general$IMS, settings$annotation, adductArg, doSusps, addMFDB, generator)
    }
    
    if (doTPs)
    {
        if (doTPsAnn)
            genScriptTPInitBlock(ionization, settings$TP, adductArg, amendScrForTPs, generator)
        genScriptTPCompBlock(nzchar(settings$annotation$formulasAlgo), nzchar(settings$annotation$compoundsAlgo),
                             settings$TP$TPsAlgo, generator)
    }
    
    if (length(settings$report$reportGen) > 0)
        genScriptReportBlock(settings$annotation, settings$report, doTPs, generator)
    
    generator$addNL()
    
    return(generator$getCode())
}

doCreateProject <- function(CCSCalibrant, anaInfoTabs, settings, noDate = FALSE)
{
    mkdirp(settings$general$destination)
    
    prepareAnaInfoData <- function(pol)
    {
        ret <- list()
        
        # Make analysis table
        if (settings$analyses$generateAnaInfo == "table")
        {
            aTab <- copy(anaInfoTabs[[pol]])
            aTab[is.na(replicate) | !nzchar(replicate), replicate := analysis]
            aTab <- aTab[, -"type"]
            
            aTab[is.na(path_centroid) | !nzchar(path_centroid), path_centroid := settings$preTreatment$output$centroid]
            aTab[is.na(path_profile) | !nzchar(path_profile), path_profile := settings$preTreatment$output$profile]
            aTab[is.na(path_ims) | !nzchar(path_ims), path_ims := settings$preTreatment$output$ims]
            
            if (settings$analyses$analysisTableFileType != "embedded")
            {
                n <- paste0("analysisTableFile", settings$analyses$analysisTableFileType)
                if (settings$general$ionization == "both")
                    n <- paste0(n, if (pol == "positive") "Pos" else "Neg")
                fp <- file.path(settings$general$destination, settings$analyses[[n]])
                
                if (settings$analyses$analysisTableFileType == "CSV")
                    write.csv(aTab, fp, row.names = FALSE)
                else # "R"
                    cat(makeAnaInfoR(aTab), file = fp, sep = "\n")
                ret$file <- fp
            }
            
            ret$tab <- aTab
        }
        else if (settings$analyses$generateAnaInfo == "dynamic")
        {
            pf <- if (settings$general$ionization != "both") "" else if (pol == "positive") "Pos" else "Neg"
            ret <- list(
                fromRaw = settings$analyses[[paste0("genAnaInfoDynRaw", pf)]],
                fromCentroid = settings$analyses[[paste0("genAnaInfoDynCentroid", pf)]],
                fromIMS = settings$analyses[[paste0("genAnaInfoDynIMS", pf)]],
                fromProfile = settings$analyses[[paste0("genAnaInfoDynProfile", pf)]]
            )
            ret[!sapply(ret, nzchar)] <- list(NULL)
        }
        return(ret)
    }
    
    if (settings$general$ionization != "both")
        aid <- prepareAnaInfoData(settings$general$ionization)
    else
        aid <- list(positive = prepareAnaInfoData("positive"), negative = prepareAnaInfoData("negative"))
    
    for (t in names(settings$preTreatment$output))
    {
        if (!nzchar(settings$preTreatment$output[[t]]))
            settings$preTreatment$output[[t]] <- "."
    }
    
    genLimitsFile(file.path(settings$general$destination, "limits.yml"), settings$general$IMS$limits)
    
    doSusps <- settings$features$exSuspList || (settings$general$ionization != "both" && nzchar(settings$features$suspects$single)) ||
        (settings$general$ionization == "both" && nzchar(settings$features$suspects$sets$pos))
    doFormAnn <- nzchar(settings$annotation$formulasAlgo); doCompAnn <- nzchar(settings$annotation$compoundsAlgo)
    if ((doSusps && "suspects" %in% settings$annotation$estIDConf && (doFormAnn || doCompAnn)) ||
        (doFormAnn && "formulas" %in% settings$annotation$estIDConf) ||
        (doCompAnn && "compounds" %in% settings$annotation$estIDConf))
    {
        genIDLevelRulesFile(file.path(settings$general$destination, "idlevelrules.yml"))
    }
    
    if ("HTML" %in% settings$report$reportGen)
        genReportSettingsFile(file.path(settings$general$destination, "report.yml"))
    
    code <- getScriptCode(CCSCalibrant, aid, settings, noDate)
    if (!nzchar(settings$general$scriptFile))
    {
        # insert at end of current document
        rstudioapi::insertText(Inf, code, rstudioapi::getSourceEditorContext()$id)
    }
    else
    {
        sp <- file.path(settings$general$destination, settings$general$scriptFile)
        cat(code, file = sp, sep = "")
        
        if (settings$general$createRStudioProj)
        {
            rstudioapi::initializeProject(settings$general$destination)
            rstudioapi::openProject(settings$general$destination)
        }
        else
            rstudioapi::navigateToFile(sp)
    }
}
