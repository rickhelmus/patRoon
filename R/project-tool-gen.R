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

getScriptCode <- function(input, anaInfoData, MSConvSettings)
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
        if (input$generateAnaInfo == "table")
        {
            if (input$analysisTableFileType == "embedded")
            {
                addComment("Create analysis table", condition = comment)
                addText(makeAnaInfoR(aid$tab, anaInfoVarName))
            }
            else
            {
                addComment("Load analysis table", condition = comment)
                if (input$analysisTableFileType == "CSV")
                    addCall(anaInfoVarName, "read.csv", list(value = aid$file, quote = TRUE))
                else # R
                    addCall(anaInfoVarName, "source", list(value = aid$file, quote = TRUE))
            }
        }
        else if (input$generateAnaInfo == "dynamic")
        {
            addCall(anaInfoVarName, "generateAnalysisInfo", list(
                list(name = "fromRaw", value = aid$fromRaw, quote = TRUE),
                list(name = "fromCentroid", value = aid$fromCentroid, quote = TRUE),
                list(name = "fromProfile", value = aid$fromProfile, quote = TRUE),
                list(name = "fromIMS", value = aid$fromIMS, quote = TRUE),
                list(name = "convCentroid", value = MSConvSettings$output$centroid, quote = TRUE),
                list(name = "convProfile", value = MSConvSettings$output$profile, quote = TRUE),
                list(name = "convIMS", value = MSConvSettings$output$ims, quote = TRUE)
            ))
        }
        else if (input$generateAnaInfo == "example")
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
            list(value = MSConvSettings$brukerCalib[[DAMethodVarName]], quote = TRUE)
        ), condition = MSConvSettings$brukerCalib$enabled && nzchar(MSConvSettings$brukerCalib[[DAMethodVarName]]), indent = 4)
        addCall(NULL, "recalibrarateDAFiles", list(value = anaInfoVarName),
                condition = MSConvSettings$brukerCalib$enabled, indent = 4)
        itf <- splitConvTypeFormat(MSConvSettings$steps$from); otf <- splitConvTypeFormat(MSConvSettings$steps$to)
        for (i in seq_len(nrow(MSConvSettings$steps)))
        {
            addCall(NULL, "convertMSFilesAnaInfo", list(
                list(value = anaInfoVarName),
                list(name = "typeFrom", value = itf[[i]][1], quote = TRUE),
                list(name = "formatFrom", value = itf[[i]][2], quote = TRUE),
                list(name = "typeTo", value = otf[[i]][1], quote = TRUE),
                list(name = "formatTo", value = otf[[i]][2], quote = TRUE),
                list(name = "algorithm", value = MSConvSettings$steps$algorithm[i], quote = TRUE),
                list(name = "overWrite", value = FALSE)
            ), indent = 4)
        }
    }
    addFindFeatures <- function(varName, anaInfoVarName)
    {
        addCall(varName, "findFeatures", list(
            list(value = anaInfoVarName),
            list(value = if (input$featFinder == "XCMS") "xcms3" else tolower(input$featFinder), quote = TRUE),
            list(name = "noiseThrInt", value = 1000, condition = input$featFinder == "OpenMS"),
            list(name = "chromSNR", value = 3, condition = input$featFinder == "OpenMS"),
            list(name = "chromFWHM", value = 5, condition = input$featFinder == "OpenMS"),
            list(name = "minFWHM", value = 1, condition = input$featFinder == "OpenMS"),
            list(name = "maxFWHM", value = 30, condition = input$featFinder == "OpenMS"),
            list(name = "kmeans", value = TRUE, condition = input$featFinder == "KPIC2"),
            list(name = "level", value = 1000, condition = input$featFinder == "KPIC2"),
            list(name = "doFMF", value = TRUE, condition = input$featFinder == "Bruker")
        ))
    }
    getAdductArg <- function(cond = TRUE) list(name = "adduct", value = if (input$ionization == "positive") "[M+H]+" else "[M-H]-",
                                               quote = TRUE,
                                               condition = input$ionization != "both" &&
                                                   (!nzchar(input$components) || input$components == "nontarget" || !input$selectIons) &&
                                                   cond)
    addLoadSuspCall <- function(var, file)
    {
        addCall(var, "read.csv", list(
            list(value = file, quote = TRUE),
            list(name = "stringsAsFactors", value = FALSE)
        ))
    }
    addLoadSuspList <- function(name, doEx, varBase, exVarBase, inpVarBase)
    {
        # Generalized for suspect screening and ISTD normalization
        
        varBasePos <- paste0(varBase, "Pos"); varBaseNeg <- paste0(varBase, "Neg")
        
        if (doEx)
        {
            addComment(paste("Get example", name))
            
            exVarBasePos <- paste0("patRoonData::", exVarBase, "Pos")
            exVarBaseNeg <- paste0("patRoonData::", exVarBase, "Neg")
            
            if (input$ionization == "positive")
                addAssignment(varBase, exVarBasePos)
            else if (input$ionization == "negative")
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
            
            if (input$ionization != "both")
                addLoadSuspCall(varBase, input[[inpVarBase]])
            else
            {
                if (nzchar(input$suspectListNeg))
                {
                    addLoadSuspCall(varBasePos, input[[paste0(inpVarBase, "Pos")]])
                    addLoadSuspCall(varBaseNeg, input[[paste0(inpVarBase, "Neg")]])
                }
                else
                    addLoadSuspCall(varBase, input[[paste0(inpVarBase, "Pos")]])
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
    
    addAssignment("workPath", input$destinationPath, quote = TRUE)
    addCall(NULL, "setwd", list(value = "workPath"))
    addNL()
    
    if (input$ionization != "both")
        addAnaInfo("anaInfo", anaInfoData, TRUE)
    else
    {
        addAnaInfo("anaInfoPos", anaInfoData$positive, TRUE)
        addAnaInfo("anaInfoNeg", anaInfoData$negative, FALSE)
    }
    
    if (nrow(MSConvSettings$steps) > 0 || MSConvSettings$brukerCalib$enabled)
    {
        addNL()
        addComment("Set to FALSE to skip data pre-treatment")
        addAssignment("doDataPretreatment", TRUE)
        addText("if (doDataPretreatment)\n{")
        if (input$ionization != "both")
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
    
    if (input$featGrouper != "SIRIUS") # NOTE: never the case with sets
    {
        addComment("Find all features")
        addComment(sprintf("NOTE: see the %s manual for many more options",
                           if (input$featFinder == "OpenMS") "reference" else input$featFinder),
                   condition = !input$featFinder %in% c("Bruker", "SIRIUS"))
        if (input$ionization != "both")
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
        list(value = "fList", condition = input$featGrouper != "SIRIUS"),
        list(value = "anaInfo", condition = input$featGrouper == "SIRIUS"),
        list(value = if (input$featGrouper == "XCMS") "xcms3" else tolower(input$featGrouper), quote = TRUE),
        list(name = "rtalign", value = TRUE, condition = input$featGrouper != "SIRIUS"),
        list(name = "groupParam", value = "xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$replicate)",
             condition = input$featGrouper == "XCMS"),
        list(name = "retAlignParam", value = "xcms::ObiwarpParam()", condition = input$featGrouper == "XCMS")
    ))
    
    retRange <- c(input[["retention-min"]], input[["retention-max"]])
    if (all(retRange == 0))
        retRange <- NULL
    else if (retRange[2] == 0)
        retRange[2] <- Inf
    mzRange <- c(input[["mz-min"]], input[["mz-max"]])
    if (all(mzRange == 0))
        mzRange <- NULL
    else if (mzRange[2] == 0)
        mzRange[2] <- Inf
    addNL()
    addComment("Basic rule based filtering")
    addCall("fGroups", "filter", list(
        list(value = "fGroups"),
        list(name = "preAbsMinIntensity", value = input$preIntThr, zeroToNULL = TRUE),
        list(name = "absMinIntensity", value = input$intThr, zeroToNULL = TRUE),
        list(name = "relMinReplicateAbundance", value = input$repAbundance, zeroToNULL = TRUE),
        list(name = "maxReplicateIntRSD", value = input$maxRepRSD, zeroToNULL = TRUE),
        list(name = "blankThreshold", value = input$blankThr, zeroToNULL = TRUE),
        list(name = "removeBlanks", value = input$removeBlanks),
        list(name = "retentionRange", value = retRange),
        list(name = "mzRange", value = mzRange)
    ))
    
    if (nzchar(input$components))
    {
        addHeader("componentization")
        
        addComment("Perform automatic generation of components")
        addCall("components", "generateComponents", list(
            list(value = "fGroups"),
            list(value = tolower(input$components), quote = TRUE),
            list(name = "ionization", value = input$ionization, quote = TRUE, condition = input$ionization != "both"),
            list(name = "rtRange", value = c(-120, 120), condition = input$components == "nontarget"),
            list(name = "mzRange", value = c(5, 120), condition = input$components == "nontarget"),
            list(name = "elements", value = c("C", "H", "O"), quote = TRUE, condition = input$components == "nontarget"),
            list(name = "rtDev", value = defaultLim("retention", "wide"), condition = input$components == "nontarget"),
            list(name = "absMzDev", value = 0.002, condition = input$components == "nontarget")
        ))
        
        if (input$selectIons && input$components != "nontarget")
        {
            pa <- switch(input$ionization,
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
    
    if (input$featNorm != "none" || input$groupNorm)
    {
        addHeader("normalization")
        
        if (input$featNorm == "istd")
        {
            addLoadSuspList("ISTD list",
                            (input$ionization != "both" && !nzchar(input$ISTDList)) ||
                                (input$ionization == "both" && !nzchar(input$ISTDListPos)),
                            "ISTDList", "ISTDList", "ISTDList")
            addNL()
            if (input$ionization != "both")
                addComment("Set adduct to NULL if ISTD list contains an adduct column")
            
            standards <- if (input$ionization != "both" || (!input$exSuspList && !nzchar(input$suspectListNeg)))
                "ISTDList"
            else
                c("ISTDListPos", "ISTDListNeg")
        }
        else
            standards <- NULL
        
        addCall("fGroups", "normInts", list(
            list(value = "fGroups"),
            list(name = "featNorm", value = input$featNorm, quote = TRUE),
            list(name = "groupNorm", value = input$groupNorm),
            list(name = "normFunc", value = "max"),
            list(name = "standards", value = standards, condition = input$featNorm == "istd"),
            list(name = "ISTDRTWindow", value = 120, condition = input$featNorm == "istd"),
            list(name = "ISTDMZWindow", value = 300, condition = input$featNorm == "istd"),
            list(name = "minISTDs", value = 3, condition = input$featNorm == "istd"),
            list(name = "rtWindow", value = 12, condition = input$featNorm == "istd"),
            list(name = "mzWindow", value = 0.005, condition = input$featNorm == "istd"),
            getAdductArg(input$featNorm == "istd")
        ))
    }
    
    doSusps <- input$exSuspList || (input$ionization != "both" && nzchar(input$suspectList)) ||
        (input$ionization == "both" && nzchar(input$suspectListPos))
    
    if (doSusps)
    {
        addHeader("suspect screening")
        
        addLoadSuspList("suspect list", input$exSuspList, "suspList", "suspects", "suspectList")
        
        useScrForTPScreening <- input$doTPs && input$TPGen != "Logic" && input$TPGenInput == "screening"
        
        addNL()
        
        if (useScrForTPScreening)
            addComment("NOTE: onlyHits is set to FALSE to ensure TPs can be found below")
        else
            addComment("Set onlyHits to FALSE to retain features without suspects (eg for full NTA)")
        if (input$ionization != "both")
            addComment("Set adduct to NULL if suspect list contains an adduct column")
        
        if (input$ionization != "both" || (!input$exSuspList && !nzchar(input$suspectListNeg)))
            addScreenCall("suspList", !useScrForTPScreening)
        else
            addScreenCall("list(suspListPos, suspListNeg)", !useScrForTPScreening)
    }
    else
        useScrForTPScreening <- FALSE
    
    if (input$doTPs)
    {
        addHeader("transformation products")
        
        if (input$TPGen != "Logic" && input$TPGenInput == "suspects")
        {
            addComment("Load parent suspect list")
            addLoadSuspCall("suspListParents", input$TPSuspectList)
            addNL()
        }
        
        addComment("Obtain TPs")
        addCall("TPs", "generateTPs", list(
            list(value = tolower(input$TPGen), quote = TRUE),
            list(name = "parents", value = switch(input$TPGenInput,
                                                  suspects = "suspListParents",
                                                  screening = "fGroups",
                                                  all = "NULL"),
                 condition = input$TPGen != "Logic"),
            list(value = "fGroups", condition = input$TPGen == "Logic"),
            getAdductArg(input$TPGen == "Logic"),
            list(name = "type", value = "env", quote = TRUE, condition = input$TPGen == "BioTransformer"),
            list(name = "transLibrary", value = "photolysis_ranked", quote = TRUE, condition = input$TPGen == "CTS"),
            list(name = "generations", value = if (input$TPGen == "BioTransformer") 2 else 1,
                 condition = input$TPGen != "Logic"),
            list(name = "calcSims", value = FALSE, condition = input$TPGen != "Logic")
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
    
    doMSPL <- nzchar(input$formulaGen) || nzchar(input$compIdent)
    
    if (nzchar(input$formulaGen) || nzchar(input$compIdent))
    {
        addHeader("annotation")
        
        useFMF <- input$featFinder == "Bruker" && input$peakListGen == "Bruker"
        addComment("Retrieve MS peak lists")
        addCall("avgMSListParams", "getDefAvgPListParams", list(name = "clusterMzWindow", value = 0.005))
        addCall("mslists", "generateMSPeakLists", list(
            list(value = "fGroups"),
            list(value = "mzr", quote = TRUE, condition = input$peakListGen == "mzR"),
            list(value = if (useFMF) "brukerfmf" else "bruker", quote = TRUE, condition = input$peakListGen == "Bruker"),
            list(name = "maxMSRtWindow", value = 5, condition = !useFMF),
            list(name = "precursorMzWindow", value = if (input$DIA) "NULL" else input$precursorMzWindow,
                 input$peakListGen == "mzR"),
            list(name = "bgsubtr", value = TRUE, condition = !useFMF && input$peakListGen == "Bruker"),
            list(name = "MSMSType", value = if (input$DIA) "BBCID" else "MSMS", quote = TRUE,
                 condition = !useFMF && input$peakListGen == "Bruker"),
            list(name = "avgFeatParams", value = "avgMSListParams", condition = input$peakListGen == "mzR"),
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
        
        if (nzchar(input$formulaGen))
        {
            addNL()
            addComment("Calculate formula candidates")
            addCall("formulas", "generateFormulas", list(
                list(value = "fGroups"),
                list(value = "mslists"),
                list(value = tolower(input$formulaGen), quote = TRUE),
                list(name = "relMzDev", value = 5, condition = input$formulaGen != "Bruker"),
                list(name = "precursorMzSearchWindow", value = defaultLim("mz", "narrow"),
                     condition = input$formulaGen == "Bruker"),
                getAdductArg(),
                list(name = "elements", value = "CHNOP", quote = TRUE, condition = input$formulaGen != "Bruker"),
                list(name = "oc", value = FALSE, condition = input$formulaGen == "GenForm"),
                list(name = "profile", value = "qtof", quote = TRUE, condition = input$formulaGen == "SIRIUS"),
                list(name = "calculateFeatures", value = "TRUE", condition = input$formulaGen != "Bruker"),
                list(name = "featThresholdAnn", value = 0.75),
                list(name = "setThresholdAnn", value = 0, condition = input$ionization == "both")
            ))
        }
        
        if (nzchar(input$compIdent))
        {
            addNL()
            
            if (input$compIdent == "Library")
            {
                addComment("Load MS library. You may want to filter it, please see the manuals for more details.")
                addCall("mslibrary", "loadMSLibrary", list(
                    list(value = input$MSLibraryPath, quote = TRUE),
                    list(value = input$MSLibraryFormat, quote = TRUE),
                    list(name = "absMzDev", value = 0.002)
                ))
            }
            
            doTPDB <- input$compIdent == "MetFrag" && input$doTPs && input$TPGen != "Logic" && input$TPDoMFDB
            if (doTPDB)
                addCall(NULL, "convertToMFDB", list(
                    list(value = "TPs"),
                    list(value = "TP-database.csv", quote = TRUE),
                    list(name = "includeParents", value = TRUE)))
            
            
            addComment("Calculate compound structure candidates")
            
            if (input$compIdent == "SIRIUS")
                addComment("Please see the handbook for SIRIUS login options. If you want to disable automatic login set login=FALSE")
            
            addCall("compounds", "generateCompounds", list(
                list(value = "fGroups"),
                list(value = "mslists"),
                list(value = tolower(input$compIdent), quote = TRUE),
                list(name = "dbRelMzDev", value = 5, condition = input$compIdent == "MetFrag"),
                list(name = "fragRelMzDev", value = 5, condition = input$compIdent == "MetFrag"),
                list(name = "fragAbsMzDev", value = 0.002, condition = input$compIdent == "MetFrag"),
                list(name = "relMzDev", value = 5, condition = input$compIdent == "SIRIUS"),
                getAdductArg(),
                list(name = "database", value = "pubchem", quote = TRUE, condition = input$compIdent == "MetFrag" && !doTPDB),
                list(name = "database", value = "csv", quote = TRUE, condition = doTPDB),
                list(name = "extraOpts", value = "list(LocalDatabasePath = \"TP-database.csv\")", condition = doTPDB),
                list(name = "maxCandidatesToStop", value = 2500, condition = input$compIdent == "MetFrag"),
                list(name = "fingerIDDatabase", value = "pubchem", quote = TRUE, condition = input$compIdent == "SIRIUS"),
                list(name = "elements", value = "CHNOP", quote = TRUE, condition = input$compIdent == "SIRIUS"),
                list(name = "profile", value = "qtof", quote = TRUE, condition = input$compIdent == "SIRIUS"),
                list(name = "login", value = "interactive", quote = TRUE, condition = input$compIdent == "SIRIUS"),
                list(name = "alwaysLogin", value = FALSE, condition = input$compIdent == "SIRIUS"),
                list(name = "MSLibrary", value = "mslibrary", condition = input$compIdent == "Library"),
                list(name = "minSim", value = 0.75, condition = input$compIdent == "Library"),
                list(name = "absMzDev", value = 0.002, condition = input$compIdent == "Library"),
                list(name = "specSimParams", value = "getDefSpecSimParams()", condition = input$compIdent == "Library"),
                list(name = "setThresholdAnn", value = 0, condition = input$ionization == "both")
            ))
            if (input$compIdent == "MetFrag")
            {
                addCall("compounds", "addFormulaScoring", list(
                    list(value = "compounds"),
                    list(value = "formulas"),
                    list(name = "updateScore", value = TRUE)
                ))
            }
        }
        
        if (doSusps && input$annotateSus && (nzchar(input$formulaGen) || nzchar(input$compIdent)))
        {
            addNL()
            addComment("Annotate suspects")
            addCall("fGroups", "annotateSuspects", list(
                list(value = "fGroups"),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(input$compIdent)),
                list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
                list(name = "IDFile", value = "idlevelrules.yml", quote = TRUE, condition = input$genIDLevelFile)
            ))
        }
    }
    
    if (input$doTPs)
    {
        addHeader("Parent and TP linkage")
        
        addComment("You probably want to prioritize the data before componentization. Please see the handbook for more info.")
        addCall("componentsTP", "generateComponents", list(
            list(value = "fGroups"),
            list(value = "tp", quote = TRUE),
            list(name = "fGroupsTPs", value = "fGroups"),
            list(name = "TPs", value = "TPs"),
            list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
            list(name = "formulas", value = "formulas", condition = nzchar(input$formulaGen)),
            list(name = "compounds", value = "compounds", condition = nzchar(input$compIdent))
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
            list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen),
                 condition = input$TPGen == "Logic")
        ))
        
        addNL()
        addComment("Only keep linked feature groups")
        addAssignment("fGroups", "fGroups[results = componentsTP]")
    }
    
    if (length(input$reportGen) > 0)
    {
        # UNDONE: for now TP components are always reported instead of others
        componVal <- if (input$doTPs)
            "componentsTP"
        else if (nzchar(input$components) && (!input$selectIons || input$components == "nontarget"))
            "components"
        else
            "NULL"
        
        addHeader("reporting")
        
        if ("HTML" %in% input$reportGen)
        {
            addComment("Advanced report settings can be edited in the report.yml file.")
            addCall(NULL, "report", list(
                list(value = "fGroups"),
                list(name = "MSPeakLists", value = "mslists", isNULL = !doMSPL),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(input$compIdent)),
                list(name = "components", value = componVal),
                list(name = "TPs", value = "TPs", condition = input$doTPs),
                list(name = "settingsFile", value = "report.yml", quote = TRUE),
                list(name = "openReport", value = TRUE)
            ))
        }
        
        if ("legacy" %in% input$reportGen)
        {
            if ("HTML" %in% input$reportGen)
                addNL()
            
            addComment("Generate reports with legacy interface.")
            addCall(NULL, "reportCSV", condition = "CSV" %in% input$reportLegacy, list(
                list(value = "fGroups"),
                list(name = "path", value = "report", quote = TRUE),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(input$compIdent)),
                list(name = "components", value = componVal)
            ))
            addCall(NULL, "reportPDF", condition = "PDF" %in% input$reportLegacy, list(
                list(value = "fGroups"),
                list(name = "path", value = "report", quote = TRUE),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(input$compIdent)),
                list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
                list(name = "components", value = componVal)
            ))
        }
    }
    
    addNL()
    
    return(paste0(textConnectionValue(txtCon), collapse = "\n"))
}

doCreateProject <- function(input, anaInfoTabs, MSConvSettings)
{
    mkdirp(input$destinationPath)
    
    prepareAnaInfoData <- function(pol)
    {
        ret <- list()
        
        # Make analysis table
        if (input$generateAnaInfo == "table")
        {
            aTab <- copy(anaInfoTabs[[pol]])
            aTab[is.na(replicate) | !nzchar(replicate), replicate := analysis]
            aTab <- aTab[, -"type"]
            
            aTab[is.na(path_centroid) | !nzchar(path_centroid), path_centroid := MSConvSettings$output$centroid]
            aTab[is.na(path_profile) | !nzchar(path_profile), path_profile := MSConvSettings$output$profile]
            aTab[is.na(path_ims) | !nzchar(path_ims), path_ims := MSConvSettings$output$ims]
            
            if (input$analysisTableFileType != "embedded")
            {
                n <- paste0("analysisTableFile", input$analysisTableFileType)
                if (input$ionization == "both")
                    n <- paste0(n, if (input$ionization == "positive") "Pos" else "Neg")
                fp <- file.path(input$destinationPath, input[[n]])
                
                if (input$analysisTableFileType == "CSV")
                    write.csv(aTab, fp, row.names = FALSE)
                else # "R"
                    cat(makeAnaInfoR(aTab), file = fp, sep = "\n")
                ret$file <- fp
            }
            
            ret$tab <- aTab
        }
        else (input$generateAnaInfo == "dynamic")
        {
            pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
            ret <- list(
                fromRaw = input[[paste0("genAnaInfoDynRaw", pol)]],
                fromCentroid = input[[paste0("genAnaInfoDynCentroid", pol)]],
                fromIMS = input[[paste0("genAnaInfoDynIMS", pol)]],
                fromProfile = input[[paste0("genAnaInfoDynProfile", pol)]]
            )
        }
        return(ret)
    }
    
    if (input$ionization != "both")
        aid <- prepareAnaInfoData(input$ionization)
    else
        aid <- list(positive = prepareAnaInfoData("positive"), negative = prepareAnaInfoData("negative"))
    
    for (t in getMSFileTypes())
    {
        if (!nzchar(MSConvSettings$output[[t]]))
            MSConvSettings$output[[t]] <- "."
    }
    
    doSusps <- input$exSuspList || (input$ionization != "both" && nzchar(input$suspectList)) ||
        (input$ionization == "both" && nzchar(input$suspectListPos))
    if (doSusps && input$annotateSus && input$genIDLevelFile)
        genIDLevelRulesFile(file.path(input$destinationPath, "idlevelrules.yml"))
    
    if ("HTML" %in% input$reportGen)
        genReportSettingsFile(file.path(input$destinationPath, "report.yml"))
    
    code <- getScriptCode(input, aid, MSConvSettings$steps)
    if (input$outputScriptTo == "curFile")
    {
        # insert at end of current document
        rstudioapi::insertText(Inf, code, rstudioapi::getSourceEditorContext()$id)
    }
    else
    {
        sp <- file.path(input$destinationPath, input$scriptFile)
        cat(code, file = sp, sep = "")
        
        if (input$createRStudioProj)
        {
            rstudioapi::initializeProject(input$destinationPath)
            rstudioapi::openProject(input$destinationPath)
        }
        else
            rstudioapi::navigateToFile(sp)
    }
}
