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

getScriptCode <- function(anaInfoData, settings)
{
    txtCon <- withr::local_connection(textConnection(NULL, "w"))
    
    # helper functions
    doQuote <- function(txt) paste0("\"", txt, "\"")
    addText <- function(txt) writeLines(txt, txtCon)
    addNL <- function(count = 1) addText(rep("", count))
    addComment <- function(txt, condition = TRUE) { if (condition) addText(paste("#", txt)) }
    addHeader <- function(title)
    {
        hd <- strrep("-", 25)
        addNL()
        addComment(c(hd, title, hd))
        addNL()
    }
    addAssignment <- function(var, val, quote = FALSE) addText(paste(var, "<-", if (quote) doQuote(val) else val))
    addCall <- function(var, func, args, condition = TRUE, indent = 0)
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
            # HACK: temporary replace " = " with "|=|" in order to keep var assignments on the same line
            args <- lapply(args, function(a) 
            {
                a <- gsub(", ", ",|", a, fixed = TRUE) # items in a list
                a <- gsub(" = ", "|=|", a, fixed = TRUE) # named arg assignments
                return(a)
            })
            argText <- paste0(args, collapse = ", ")
            
            argText <- paste0(strwrap(argText, 120 - nchar(callPrefix) - 1, # -1: trailing ")"
                                      exdent = nchar(callPrefix)),
                              collapse = "\n")
            
            argText <- gsub("|", " ", argText, fixed = TRUE)
        }
        
        cl <- paste0(callPrefix, argText, ")")
        addText(cl)
    }
    addAnaInfo <- function(anaInfoVarName, aid, comment)
    {
        if (settings$analyses$generateAnaInfo == "table")
        {
            if (settings$analyses$analysisTableFileType == "embedded")
            {
                addComment("Create analysis table", condition = comment)
                addText(makeAnaInfoR(aid$tab, anaInfoVarName))
            }
            else
            {
                addComment("Load analysis table", condition = comment)
                if (settings$analyses$analysisTableFileType == "CSV")
                    addCall(anaInfoVarName, "read.csv", list(value = aid$file, quote = TRUE))
                else # R
                    addCall(anaInfoVarName, "source", list(value = aid$file, quote = TRUE))
            }
        }
        else if (settings$analyses$generateAnaInfo == "dynamic")
        {
            addCall(anaInfoVarName, "generateAnalysisInfo", list(
                list(name = "fromRaw", value = aid$fromRaw, quote = TRUE),
                list(name = "fromCentroid", value = aid$fromCentroid, quote = TRUE),
                list(name = "fromProfile", value = aid$fromProfile, quote = TRUE),
                list(name = "fromIMS", value = aid$fromIMS, quote = TRUE),
                list(name = "convCentroid", value = settings$preTreatment$output$centroid, quote = TRUE),
                list(name = "convProfile", value = settings$preTreatment$output$profile, quote = TRUE),
                list(name = "convIMS", value = settings$preTreatment$output$ims, quote = TRUE)
            ))
        }
        else if (settings$analyses$generateAnaInfo == "example")
        {
            addComment("Example data from patRoonData package (triplicate solvent blank + triplicate standard)",
                       condition = comment)
            addCall(anaInfoVarName, "patRoonData::exampleAnalysisInfo", list(value = ionization, quote = TRUE))
        }
        else # none
        {
            addComment("NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.",
                       condition = comment)
            addCall(anaInfoVarName, "data.frame", list(
                list(name = "path_centroid", value = "character()"),
                list(name = "analysis", value = "character()"),
                list(name = "replicate", value = "character()"),
                list(name = "blank", value = "character()")
            ))
        }
    }
    addPrepBlock <- function(anaInfoVarName, DAMethodVarName)
    {
        addCall(NULL, "setDAMethod", list(
            list(value = anaInfoVarName),
            list(value = settings$preTreatment$brukerCalib[[DAMethodVarName]], quote = TRUE)
        ), condition = settings$preTreatment$brukerCalib$enabled && nzchar(settings$preTreatment$brukerCalib[[DAMethodVarName]]), indent = 4)
        addCall(NULL, "recalibrarateDAFiles", list(value = anaInfoVarName),
                condition = settings$preTreatment$brukerCalib$enabled, indent = 4)
        itf <- splitConvTypeFormat(settings$preTreatment$steps$from); otf <- splitConvTypeFormat(settings$preTreatment$steps$to)
        for (i in seq_len(nrow(settings$preTreatment$steps)))
        {
            addCall(NULL, "convertMSFilesAnaInfo", list(
                list(value = anaInfoVarName),
                list(name = "typeFrom", value = itf[[i]][1], quote = TRUE),
                list(name = "formatFrom", value = itf[[i]][2], quote = TRUE),
                list(name = "typeTo", value = otf[[i]][1], quote = TRUE),
                list(name = "formatTo", value = otf[[i]][2], quote = TRUE),
                list(name = "algorithm", value = settings$preTreatment$steps$algorithm[i], quote = TRUE),
                list(name = "overWrite", value = FALSE)
            ), indent = 4)
        }
    }
    addFindFeatures <- function(varName, anaInfoVarName)
    {
        addCall(varName, "findFeatures", list(
            list(value = anaInfoVarName),
            list(value = if (settings$features$featAlgo == "XCMS") "xcms3" else tolower(settings$features$featAlgo), quote = TRUE),
            list(name = "noiseThrInt", value = 1000, condition = settings$features$featAlgo == "OpenMS"),
            list(name = "chromSNR", value = 3, condition = settings$features$featAlgo == "OpenMS"),
            list(name = "chromFWHM", value = 5, condition = settings$features$featAlgo == "OpenMS"),
            list(name = "minFWHM", value = 1, condition = settings$features$featAlgo == "OpenMS"),
            list(name = "maxFWHM", value = 30, condition = settings$features$featAlgo == "OpenMS"),
            list(name = "kmeans", value = TRUE, condition = settings$features$featAlgo == "KPIC2"),
            list(name = "level", value = 1000, condition = settings$features$featAlgo == "KPIC2"),
            list(name = "doFMF", value = TRUE, condition = settings$features$featAlgo == "Bruker")
        ))
    }
    getAdductArg <- function(cond = TRUE) list(name = "adduct", value = if (settings$general$ionization == "positive") "[M+H]+" else "[M-H]-",
                                               quote = TRUE,
                                               condition = settings$general$ionization != "both" &&
                                                   (!nzchar(settings$annotation$componAlgo) || settings$annotation$componAlgo == "nontarget" || !settings$annotation$selectIons) &&
                                                   cond)
    addLoadSuspCall <- function(var, file)
    {
        addCall(var, "read.csv", list(
            list(value = file, quote = TRUE),
            list(name = "stringsAsFactors", value = FALSE)
        ))
    }
    addLoadSuspList <- function(name, ionization, doEx, featSettings, varBase, exVarBase, inpVarBase)
    {
        # Generalized for suspect screening and ISTD normalization
        
        varBasePos <- paste0(varBase, "Pos"); varBaseNeg <- paste0(varBase, "Neg")
        
        if (doEx)
        {
            addComment(paste("Get example", name))
            
            exVarBasePos <- paste0("patRoonData::", exVarBase, "Pos")
            exVarBaseNeg <- paste0("patRoonData::", exVarBase, "Neg")
            
            if (ionization == "positive")
                addAssignment(varBase, exVarBasePos)
            else if (ionization == "negative")
                addAssignment(varBase, exVarBaseNeg)
            else
            {
                addAssignment(varBasePos, exVarBasePos)
                addAssignment(varBaseNeg, exVarBaseNeg)
            }
        }
        else
        {
            addComment(paste("Load", name))
            
            if (ionization != "both")
                addLoadSuspCall(varBase, featSettings[[inpVarBase]])
            else
            {
                inpVarBasePos <- paste0(inpVarBase, "Pos"); inpVarBaseNeg <- paste0(inpVarBase, "Neg")
                if (nzchar(featSettings[[inpVarBaseNeg]]))
                {
                    addLoadSuspCall(varBasePos, featSettings[[inpVarBasePos]])
                    addLoadSuspCall(varBaseNeg, featSettings[[inpVarBaseNeg]])
                }
                else
                    addLoadSuspCall(varBase, featSettings[[inpVarBasePos]])
            }
        }
    }
    addScreenCall <- function(susp, oh = TRUE, am = NULL)
    {
        addCall("fGroups", "screenSuspects", list(
            list(value = "fGroups"),
            list(value = susp),
            list(name = "rtWindow", value = 12),
            list(name = "mzWindow", value = 0.005),
            getAdductArg(),
            list(name = "onlyHits", value = oh),
            list(name = "amend", value = am, condition = !is.null(am))
        ))
    }
    
    
    addComment(paste("Script automatically generated on", date()))
    addNL()
    addCall(NULL, "library", list(value = "patRoon"))
    
    addHeader("initialization")
    
    addAssignment("workPath", settings$general$destination, quote = TRUE)
    addCall(NULL, "setwd", list(value = "workPath"))
    addNL()
    
    if (settings$general$ionization != "both")
        addAnaInfo("anaInfo", anaInfoData, TRUE)
    else
    {
        addAnaInfo("anaInfoPos", anaInfoData$positive, TRUE)
        addAnaInfo("anaInfoNeg", anaInfoData$negative, FALSE)
    }
    
    if (nrow(settings$preTreatment$steps) > 0 || settings$preTreatment$brukerCalib$enabled)
    {
        addNL()
        addComment("Set to FALSE to skip data pre-treatment")
        addAssignment("doDataPretreatment", TRUE)
        addText("if (doDataPretreatment)\n{")
        if (settings$general$ionization != "both")
            addPrepBlock("anaInfo", "method")
        else
        {
            addPrepBlock("anaInfoPos", "methodPos")
            addNL()
            addPrepBlock("anaInfoNeg", "methodNeg")
        }
        addText("}")
    }
    
    addHeader("features")
    
    if (settings$features$fGroupsAlgo != "SIRIUS") # NOTE: never the case with sets
    {
        addComment("Find all features")
        addComment(sprintf("NOTE: see the %s manual for many more options",
                           if (settings$features$featAlgo == "OpenMS") "reference" else settings$features$featAlgo),
                   condition = !settings$features$featAlgo %in% c("Bruker", "SIRIUS"))
        if (settings$general$ionization != "both")
            addFindFeatures("fList", "anaInfo")
        else
        {
            addFindFeatures("fListPos", "anaInfoPos")
            addFindFeatures("fListNeg", "anaInfoNeg")
            addCall("fList", "makeSet", list(
                list(value = "fListPos"),
                list(value = "fListNeg"),
                list(name = "adducts", value = c("[M+H]+", "[M-H]-"), quote = TRUE)
            ))
        }
        
        addNL()
    }
    
    addComment("Group and align features between analyses")
    addCall("fGroups", "groupFeatures", list(
        list(value = "fList", condition = settings$features$fGroupsAlgo != "SIRIUS"),
        list(value = "anaInfo", condition = settings$features$fGroupsAlgo == "SIRIUS"),
        list(value = if (settings$features$fGroupsAlgo == "XCMS") "xcms3" else tolower(settings$features$fGroupsAlgo), quote = TRUE),
        list(name = "rtalign", value = TRUE, condition = settings$features$fGroupsAlgo != "SIRIUS"),
        list(name = "groupParam", value = "xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$replicate)",
             condition = settings$features$fGroupsAlgo == "XCMS"),
        list(name = "retAlignParam", value = "xcms::ObiwarpParam()", condition = settings$features$fGroupsAlgo == "XCMS")
    ))
    
    retRange <- settings$features$advanced$retention
    if (all(retRange == 0))
        retRange <- NULL
    else if (retRange[2] == 0)
        retRange[2] <- Inf
    mzRange <- settings$features$advanced$mz
    if (all(mzRange == 0))
        mzRange <- NULL
    else if (mzRange[2] == 0)
        mzRange[2] <- Inf
    addNL()
    addComment("Basic rule based filtering")
    addCall("fGroups", "filter", list(
        list(value = "fGroups"),
        list(name = "preAbsMinIntensity", value = settings$features$advanced$preIntThr, zeroToNULL = TRUE),
        list(name = "absMinIntensity", value = settings$features$advanced$intThr, zeroToNULL = TRUE),
        list(name = "relMinReplicateAbundance", value = settings$features$advanced$repAbundance, zeroToNULL = TRUE),
        list(name = "maxReplicateIntRSD", value = settings$features$advanced$maxRepRSD, zeroToNULL = TRUE),
        list(name = "blankThreshold", value = settings$features$advanced$blankThr, zeroToNULL = TRUE),
        list(name = "removeBlanks", value = settings$features$advanced$removeBlanks),
        list(name = "retentionRange", value = retRange),
        list(name = "mzRange", value = mzRange)
    ))
    
    if (nzchar(settings$annotation$componAlgo))
    {
        addHeader("componentization")
        
        addComment("Perform automatic generation of components")
        addCall("components", "generateComponents", list(
            list(value = "fGroups"),
            list(value = tolower(settings$annotation$componAlgo), quote = TRUE),
            list(name = "ionization", value = settings$general$ionization, quote = TRUE, condition = settings$general$ionization != "both"),
            list(name = "rtRange", value = c(-120, 120), condition = settings$annotation$componAlgo == "nontarget"),
            list(name = "mzRange", value = c(5, 120), condition = settings$annotation$componAlgo == "nontarget"),
            list(name = "elements", value = c("C", "H", "O"), quote = TRUE, condition = settings$annotation$componAlgo == "nontarget"),
            list(name = "rtDev", value = defaultLim("retention", "wide"), condition = settings$annotation$componAlgo == "nontarget"),
            list(name = "absMzDev", value = 0.002, condition = settings$annotation$componAlgo == "nontarget")
        ))
        
        if (settings$annotation$selectIons && settings$annotation$componAlgo != "nontarget")
        {
            pa <- switch(settings$general$ionization,
                         positive = "[M+H]+",
                         negative = "[M-H]-",
                         both = c("[M+H]+", "[M-H]-"))
            addCall("fGroups", "selectIons", list(
                list(value = "fGroups"),
                list(value = "components"),
                list(name = "prefAdduct", value = pa, quote = TRUE),
                list(name = "onlyMonoIso", value = TRUE)
            ))
        }
    }
    
    if (settings$features$advanced$featNorm != "none" || settings$features$advanced$groupNorm)
    {
        addHeader("normalization")
        
        if (settings$features$advanced$featNorm == "istd")
        {
            addLoadSuspList("ISTD list", settings$general$ionization,
                            (settings$general$ionization != "both" && !nzchar(settings$features$advanced$ISTDList)) ||
                                (settings$general$ionization == "both" && !nzchar(settings$features$advanced$ISTDListPos)),
                            settings$features$advanced, "ISTDList", "ISTDList", "ISTDList")
            addNL()
            if (settings$general$ionization != "both")
                addComment("Set adduct to NULL if ISTD list contains an adduct column")
            
            standards <- if (settings$general$ionization != "both" ||
                             (!settings$features$exSuspList && !nzchar(settings$features$advanced$ISTDListNeg)))
                "ISTDList"
            else
                c("ISTDListPos", "ISTDListNeg")
        }
        else
            standards <- NULL
        
        addCall("fGroups", "normInts", list(
            list(value = "fGroups"),
            list(name = "featNorm", value = settings$features$advanced$featNorm, quote = TRUE),
            list(name = "groupNorm", value = settings$features$advanced$groupNorm),
            list(name = "normFunc", value = "max"),
            list(name = "standards", value = standards, condition = settings$features$advanced$featNorm == "istd"),
            list(name = "ISTDRTWindow", value = 120, condition = settings$features$advanced$featNorm == "istd"),
            list(name = "ISTDMZWindow", value = 300, condition = settings$features$advanced$featNorm == "istd"),
            list(name = "minISTDs", value = 3, condition = settings$features$advanced$featNorm == "istd"),
            list(name = "rtWindow", value = 12, condition = settings$features$advanced$featNorm == "istd"),
            list(name = "mzWindow", value = 0.005, condition = settings$features$advanced$featNorm == "istd"),
            getAdductArg(settings$features$advanced$featNorm == "istd")
        ))
    }
    
    doSusps <- settings$features$exSuspList || (settings$general$ionization != "both" && nzchar(settings$features$suspectList)) ||
        (settings$general$ionization == "both" && nzchar(settings$features$suspectListPos))
    
    if (doSusps)
    {
        addHeader("suspect screening")
        
        addLoadSuspList("suspect list", settings$ionization, settings$features$exSuspList, settings$features,
                        "suspList", "suspects", "suspectList")
        
        useScrForTPScreening <- settings$TP$doTPs && settings$TP$TPGen != "Logic" && settings$TP$TPGenInput == "screening"
        
        addNL()
        
        if (useScrForTPScreening)
            addComment("NOTE: onlyHits is set to FALSE to ensure TPs can be found below")
        else
            addComment("Set onlyHits to FALSE to retain features without suspects (eg for full NTA)")
        if (settings$general$ionization != "both")
            addComment("Set adduct to NULL if suspect list contains an adduct column")
        
        if (settings$general$ionization != "both" || (!settings$features$exSuspList && !nzchar(settings$features$suspectListNeg)))
            addScreenCall("suspList", !useScrForTPScreening)
        else
            addScreenCall("list(suspListPos, suspListNeg)", !useScrForTPScreening)
    }
    else
        useScrForTPScreening <- FALSE
    
    if (settings$TP$doTPs)
    {
        addHeader("transformation products")
        
        if (settings$TP$TPGen != "Logic" && settings$TP$TPGenInput == "suspects")
        {
            addComment("Load parent suspect list")
            addLoadSuspCall("suspListParents", settings$TP$TPSuspectList)
            addNL()
        }
        
        addComment("Obtain TPs")
        addCall("TPs", "generateTPs", list(
            list(value = tolower(settings$TP$TPGen), quote = TRUE),
            list(name = "parents", value = switch(settings$TP$TPGenInput,
                                                  suspects = "suspListParents",
                                                  screening = "fGroups",
                                                  all = "NULL"),
                 condition = settings$TP$TPGen != "Logic"),
            list(value = "fGroups", condition = settings$TP$TPGen == "Logic"),
            getAdductArg(settings$TP$TPGen == "Logic"),
            list(name = "type", value = "env", quote = TRUE, condition = settings$TP$TPGen == "BioTransformer"),
            list(name = "transLibrary", value = "photolysis_ranked", quote = TRUE, condition = settings$TP$TPGen == "CTS"),
            list(name = "generations", value = if (settings$TP$TPGen == "BioTransformer") 2 else 1,
                 condition = settings$TP$TPGen != "Logic"),
            list(name = "calcSims", value = FALSE, condition = settings$TP$TPGen != "Logic")
        ))
        
        addNL()
        addComment("Screen TPs")
        addCall("suspListTPs", "convertToSuspects", list(
            list(value = "TPs"),
            list(name = "includeParents", value = FALSE)
        ))
        if (useScrForTPScreening)
            addComment("Amend with TP screening")
        addScreenCall("suspListTPs", am = if (useScrForTPScreening) TRUE else NULL)
    }
    
    doMSPL <- nzchar(settings$annotation$formulasAlgo) || nzchar(settings$annotation$compoundsAlgo)
    
    if (nzchar(settings$annotation$formulasAlgo) || nzchar(settings$annotation$compoundsAlgo))
    {
        addHeader("annotation")
        
        useFMF <- settings$features$featAlgo == "Bruker" && settings$annotation$peakListGen == "Bruker"
        addComment("Retrieve MS peak lists")
        addCall("avgMSListParams", "getDefAvgPListParams", list(name = "clusterMzWindow", value = 0.005))
        addCall("mslists", "generateMSPeakLists", list(
            list(value = "fGroups"),
            list(value = "mzr", quote = TRUE, condition = settings$annotation$peakListGen == "mzR"),
            list(value = if (useFMF) "brukerfmf" else "bruker", quote = TRUE, condition = settings$annotation$peakListGen == "Bruker"),
            list(name = "maxMSRtWindow", value = 5, condition = !useFMF),
            list(name = "precursorMzWindow", value = if (settings$annotation$DIA) "NULL" else settings$annotation$precursorMzWindow,
                 settings$annotation$peakListGen == "mzR"),
            list(name = "bgsubtr", value = TRUE, condition = !useFMF && settings$annotation$peakListGen == "Bruker"),
            list(name = "MSMSType", value = if (settings$annotation$DIA) "BBCID" else "MSMS", quote = TRUE,
                 condition = !useFMF && settings$annotation$peakListGen == "Bruker"),
            list(name = "avgFeatParams", value = "avgMSListParams", condition = settings$annotation$peakListGen == "mzR"),
            list(name = "avgFGroupParams", value = "avgMSListParams")
        ))
        addComment("Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.")
        addCall("mslists", "filter", list(
            list(value = "mslists"),
            list(name = "absMSIntThr", value = "NULL"),
            list(name = "absMSMSIntThr", value = "NULL"),
            list(name = "relMSIntThr", value = "NULL"),
            list(name = "relMSMSIntThr", value = 0.05),
            list(name = "topMSPeaks", value = "NULL"),
            list(name = "topMSMSPeaks", value = 25)
        ))
        
        if (nzchar(settings$annotation$formulasAlgo))
        {
            addNL()
            addComment("Calculate formula candidates")
            addCall("formulas", "generateFormulas", list(
                list(value = "fGroups"),
                list(value = "mslists"),
                list(value = tolower(settings$annotation$formulasAlgo), quote = TRUE),
                list(name = "relMzDev", value = 5, condition = settings$annotation$formulasAlgo != "Bruker"),
                list(name = "precursorMzSearchWindow", value = defaultLim("mz", "narrow"),
                     condition = settings$annotation$formulasAlgo == "Bruker"),
                getAdductArg(),
                list(name = "elements", value = "CHNOP", quote = TRUE, condition = settings$annotation$formulasAlgo != "Bruker"),
                list(name = "oc", value = FALSE, condition = settings$annotation$formulasAlgo == "GenForm"),
                list(name = "profile", value = "qtof", quote = TRUE, condition = settings$annotation$formulasAlgo == "SIRIUS"),
                list(name = "calculateFeatures", value = "TRUE", condition = settings$annotation$formulasAlgo != "Bruker"),
                list(name = "featThresholdAnn", value = 0.75),
                list(name = "setThresholdAnn", value = 0, condition = settings$general$ionization == "both")
            ))
        }
        
        if (nzchar(settings$annotation$compoundsAlgo))
        {
            addNL()
            
            if (settings$annotation$compoundsAlgo == "Library")
            {
                addComment("Load MS library. You may want to filter it, please see the manuals for more details.")
                addCall("mslibrary", "loadMSLibrary", list(
                    list(value = settings$annotation$MSLibraryPath, quote = TRUE),
                    list(value = settings$annotation$MSLibraryFormat, quote = TRUE),
                    list(name = "absMzDev", value = 0.002)
                ))
            }
            
            doTPDB <- settings$annotation$compoundsAlgo == "MetFrag" && settings$TP$doTPs && settings$TP$TPGen != "Logic" && settings$TP$TPDoMFDB
            if (doTPDB)
                addCall(NULL, "convertToMFDB", list(
                    list(value = "TPs"),
                    list(value = "TP-database.csv", quote = TRUE),
                    list(name = "includeParents", value = TRUE)))
            
            
            addComment("Calculate compound structure candidates")
            
            if (settings$annotation$compoundsAlgo == "SIRIUS")
                addComment("Please see the handbook for SIRIUS login options. If you want to disable automatic login set login=FALSE")
            
            addCall("compounds", "generateCompounds", list(
                list(value = "fGroups"),
                list(value = "mslists"),
                list(value = tolower(settings$annotation$compoundsAlgo), quote = TRUE),
                list(name = "dbRelMzDev", value = 5, condition = settings$annotation$compoundsAlgo == "MetFrag"),
                list(name = "fragRelMzDev", value = 5, condition = settings$annotation$compoundsAlgo == "MetFrag"),
                list(name = "fragAbsMzDev", value = 0.002, condition = settings$annotation$compoundsAlgo == "MetFrag"),
                list(name = "relMzDev", value = 5, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                getAdductArg(),
                list(name = "database", value = "pubchem", quote = TRUE, condition = settings$annotation$compoundsAlgo == "MetFrag" && !doTPDB),
                list(name = "database", value = "csv", quote = TRUE, condition = doTPDB),
                list(name = "extraOpts", value = "list(LocalDatabasePath = \"TP-database.csv\")", condition = doTPDB),
                list(name = "maxCandidatesToStop", value = 2500, condition = settings$annotation$compoundsAlgo == "MetFrag"),
                list(name = "fingerIDDatabase", value = "pubchem", quote = TRUE, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                list(name = "elements", value = "CHNOP", quote = TRUE, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                list(name = "profile", value = "qtof", quote = TRUE, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                list(name = "login", value = "interactive", quote = TRUE, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                list(name = "alwaysLogin", value = FALSE, condition = settings$annotation$compoundsAlgo == "SIRIUS"),
                list(name = "MSLibrary", value = "mslibrary", condition = settings$annotation$compoundsAlgo == "Library"),
                list(name = "minSim", value = 0.75, condition = settings$annotation$compoundsAlgo == "Library"),
                list(name = "absMzDev", value = 0.002, condition = settings$annotation$compoundsAlgo == "Library"),
                list(name = "specSimParams", value = "getDefSpecSimParams()", condition = settings$annotation$compoundsAlgo == "Library"),
                list(name = "setThresholdAnn", value = 0, condition = settings$general$ionization == "both")
            ))
            if (settings$annotation$compoundsAlgo == "MetFrag")
            {
                addCall("compounds", "addFormulaScoring", list(
                    list(value = "compounds"),
                    list(value = "formulas"),
                    list(name = "updateScore", value = TRUE)
                ))
            }
        }
        
        if (doSusps && settings$annotation$annotateSus && (nzchar(settings$annotation$formulasAlgo) || nzchar(settings$annotation$compoundsAlgo)))
        {
            addNL()
            addComment("Annotate suspects")
            addCall("fGroups", "annotateSuspects", list(
                list(value = "fGroups"),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(settings$annotation$formulasAlgo)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(settings$annotation$compoundsAlgo)),
                list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
                list(name = "IDFile", value = "idlevelrules.yml", quote = TRUE, condition = settings$annotation$genIDLevelFile)
            ))
        }
    }
    
    if (settings$TP$doTPs)
    {
        addHeader("Parent and TP linkage")
        
        addComment("You probably want to prioritize the data before componentization. Please see the handbook for more info.")
        addCall("componentsTP", "generateComponents", list(
            list(value = "fGroups"),
            list(value = "tp", quote = TRUE),
            list(name = "fGroupsTPs", value = "fGroups"),
            list(name = "TPs", value = "TPs"),
            list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
            list(name = "formulas", value = "formulas", condition = nzchar(settings$annotation$formulasAlgo)),
            list(name = "compounds", value = "compounds", condition = nzchar(settings$annotation$compoundsAlgo))
        ))
        
        addNL()
        addComment("You may want to configure the filtering step below. See the manuals for more details.")
        addCall("componentsTP", "filter", list(
            list(value = "componentsTP"),
            list(name = "retDirMatch", value = FALSE),
            list(name = "minSpecSim", value = "NULL"),
            list(name = "minSpecSimPrec", value = "NULL"),
            list(name = "minSpecSimBoth", value = "NULL"),
            list(name = "minFragMatches", value = "NULL"),
            list(name = "minNLMatches", value = "NULL"),
            list(name = "formulas", value = "formulas", isNULL = !nzchar(settings$annotation$formulasAlgo),
                 condition = settings$TP$TPGen == "Logic")
        ))
        
        addNL()
        addComment("Only keep linked feature groups")
        addAssignment("fGroups", "fGroups[results = componentsTP]")
    }
    
    if (length(settings$report$reportGen) > 0)
    {
        # UNDONE: for now TP components are always reported instead of others
        componVal <- if (settings$TP$doTPs)
            "componentsTP"
        else if (nzchar(settings$annotation$componAlgo) && (!settings$annotation$selectIons || settings$annotation$componAlgo == "nontarget"))
            "components"
        else
            "NULL"
        
        addHeader("reporting")
        
        if ("HTML" %in% settings$report$reportGen)
        {
            addComment("Advanced report settings can be edited in the report.yml file.")
            addCall(NULL, "report", list(
                list(value = "fGroups"),
                list(name = "MSPeakLists", value = "mslists", isNULL = !doMSPL),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(settings$annotation$formulasAlgo)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(settings$annotation$compoundsAlgo)),
                list(name = "components", value = componVal),
                list(name = "TPs", value = "TPs", condition = settings$TP$doTPs),
                list(name = "settingsFile", value = "report.yml", quote = TRUE),
                list(name = "openReport", value = TRUE)
            ))
        }
        
        if ("legacy" %in% settings$report$reportGen)
        {
            if ("HTML" %in% settings$report$reportGen)
                addNL()
            
            addComment("Generate reports with legacy interface.")
            addCall(NULL, "reportCSV", condition = "CSV" %in% settings$report$reportLegacy, list(
                list(value = "fGroups"),
                list(name = "path", value = "report", quote = TRUE),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(settings$annotation$formulasAlgo)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(settings$annotation$compoundsAlgo)),
                list(name = "components", value = componVal)
            ))
            addCall(NULL, "reportPDF", condition = "PDF" %in% settings$report$reportLegacy, list(
                list(value = "fGroups"),
                list(name = "path", value = "report", quote = TRUE),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(settings$annotation$formulasAlgo)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(settings$annotation$compoundsAlgo)),
                list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
                list(name = "components", value = componVal)
            ))
        }
    }
    
    addNL()
    
    return(paste0(textConnectionValue(txtCon), collapse = "\n"))
}

doCreateProject <- function(anaInfoTabs, settings)
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
                    n <- paste0(n, if (settings$general$ionization == "positive") "Pos" else "Neg")
                fp <- file.path(settings$general$destination, settings$analyses[[n]])
                
                if (settings$analyses$analysisTableFileType == "CSV")
                    write.csv(aTab, fp, row.names = FALSE)
                else # "R"
                    cat(makeAnaInfoR(aTab), file = fp, sep = "\n")
                ret$file <- fp
            }
            
            ret$tab <- aTab
        }
        else (settings$analyses$generateAnaInfo == "dynamic")
        {
            pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
            ret <- list(
                fromRaw = settings$analyses[[paste0("genAnaInfoDynRaw", pol)]],
                fromCentroid = settings$analyses[[paste0("genAnaInfoDynCentroid", pol)]],
                fromIMS = settings$analyses[[paste0("genAnaInfoDynIMS", pol)]],
                fromProfile = settings$analyses[[paste0("genAnaInfoDynProfile", pol)]]
            )
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
    
    doSusps <- settings$features$exSuspList || (settings$general$ionization != "both" && nzchar(settings$features$suspectList)) ||
        (settings$general$ionization == "both" && nzchar(settings$features$suspectListPos))
    if (doSusps && settings$annotation$annotateSus && settings$annotation$genIDLevelFile)
        genIDLevelRulesFile(file.path(settings$general$destination, "idlevelrules.yml"))
    
    if ("HTML" %in% settings$report$reportGen)
        genReportSettingsFile(file.path(settings$general$destination, "report.yml"))
    
    code <- getScriptCode(aid, settings)
    if (settings$general$outputScriptTo == "curFile")
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
