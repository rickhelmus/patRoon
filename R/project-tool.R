# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

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

getTPGenInputs <- function(isLib)
{
    ret <- c("Suspect list" = "suspects",
             "Screening results" = "screening")
    if (isLib)
        ret <- c(ret, "All library parents" = "all")
    return(ret)
}

fileSelect <- function(idText, idButton, label, value = "", ...)
{
    fillRow(
        flex = c(1, NA),
        textInput(idText, label, value, width = "100%", ...),
        actionButton(idButton, "", icon("folder-open"), style = "margin: 25px 0 0 15px")
    )
}

getAnaTableFileUI <- function(pol)
{
    # NOTE: pol may be NULL
    basef <- "analyses"
    if (!is.null(pol))
        basef <- paste0(basef, if (pol == "positive") "-pos" else "-neg")
    pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
    
    fillRow(
        flex = NA,
        conditionalPanel(
            condition = "input.analysisTableFileType == \"CSV\"",
            textInput(paste0("analysisTableFileCSV", pf), "File name", paste0(basef, ".csv"), width = "100%")
        ),
        conditionalPanel(
            condition = "input.analysisTableFileType == \"R\"",
            textInput(paste0("analysisTableFileR", pf), "File name", paste0(basef, ".R"), width = "100%")
        )
    )
}

genAnaDynUI <- function(pol)
{
    pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
    
    htmltools::tagList(
        fileSelect(paste0("genAnaInfoDynRaw", pf), "genAnaInfoDynRawButton",
                   "raw analyses", "", placeholder = "Leave empty for none"),
        fileSelect(paste0("genAnaInfoDynCentroid", pf), "genAnaInfoDynCentroidButton",
                   "centroided analyses", "", placeholder = "Leave empty for none"),
        fileSelect(paste0("genAnaInfoDynIMS", pf), "genAnaInfoDynIMSButton",
                   "IMS analyses", "", placeholder = "Leave empty for none"),
        fileSelect(paste0("genAnaInfoDynProfile", pf), "genAnaInfoDynProfileButton",
                   "profile analyses", "", placeholder = "Leave empty for none")
    )
}

textNote <- function(txt) div(style = "margin: 8px 0 12px; font-size: small", txt)

getNewProjectUI <- function(destPath)
{
    rangeNumeric <- function(id, label, minVal = 0, maxVal = 0, ...)
    {
        fillRow(
            numericInput(paste0(id, "-min"), paste("Min.", label), value = minVal, ..., width = "95%"),
            numericInput(paste0(id, "-max"), paste("Max.", label), value = maxVal, ..., width = "100%")
        )
    }
    
    anaTableCondition <- function(pol, wh)
    {
        ret <- "input.generateAnaInfo == \"table\""
        if (!missing(pol))
        {
            if (is.null(pol))
                ret <- paste(ret, "&& input.ionization != \"both\"")
            else
            {
                ret <- paste(ret, "&& input.ionization == \"both\"")
                if (pol != "both")
                    ret <- paste0(ret, " && input.currentSetStatic == \"", pol, "\"")
            }
        }
        return(ret)
    }
    
    anaDynCondition <- function(pol, wh)
    {
        ret <- "input.generateAnaInfo == \"dynamic\""
        if (is.null(pol))
            ret <- paste(ret, "&& input.ionization != \"both\"")
        else
        {
            ret <- paste(ret, "&& input.ionization == \"both\"")
            if (pol != "both")
                ret <- paste0(ret, " && input.currentSetDynamic == \"", pol, "\"")
        }
        return(ret)
    }
    
    miniUI::miniPage(
        shinyjs::useShinyjs(),
        
        # based on https://stackoverflow.com/a/63882648
        tags$style(
            type = 'text/css',
            '.modal-dialog { width: fit-content !important; height: fit-content !important; }'
        ),
        
        tags$script(htmlwidgets::JS('
                Shiny.addCustomMessageHandler("selectHotRows", function(data)
                {
                    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
                    var ht = HTMLWidgets.getInstance(window[data.tab]).hot;
                    ht.selectRows(data.range[0] - 1, data.range[1] - 1);
                });
            ')
        ),
        
        miniUI::gadgetTitleBar("Create project tool", right = miniUI::miniTitleBarButton("create", "Create", TRUE)),

        miniUI::miniTabstripPanel(
            miniUI::miniTabPanel(
                "General", icon = icon("save"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        fileSelect("destinationPath", "projectDestButton", "Project destination",
                                   if (is.null(destPath)) "~/" else destPath),
                        br(),
                        fillRow(
                            height = 100,
                            radioButtons("outputScriptTo", "Insert code into", c("New file" = "newFile",
                                                                                 "Current file" = "curFile")),
                            conditionalPanel(
                                condition = "input.outputScriptTo == \"newFile\"",
                                textInput("scriptFile", "Script file", "process.R", width = "80%"),
                                checkboxInput("createRStudioProj", "Create (and open) RStudio project", value = TRUE)
                            )
                        ),
                        br(),
                        fillRow(
                            height = 100,
                            radioButtons("ionization", "Ionization", c("positive", "negative", "both (sets)" = "both"))
                        )
                    )
                ),
                miniUI::miniButtonBlock(
                    actionButton("loadParams", "Load parameters"),
                    actionButton("saveParams", "Save parameters")
                )
            ),
            
            miniUI::miniTabPanel(
                "Analyses", icon = icon("folder-open"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        fillRow(
                            height = 180,
                            fillCol(
                                flex = NA,
                                radioButtons("generateAnaInfo", "Generate analysis information",
                                             c("None" = "none", "Analyses table" = "table",
                                               "Dynamically in script" = "dynamic", "Example data" = "example")),
                                conditionalPanel(
                                    condition = anaTableCondition(),
                                    radioButtons("analysisTableFileType", "Save analyses table",
                                                 c("as CSV file" = "CSV", "as R file" = "R", "embedded in script" = "embedded"),
                                                 inline = TRUE)
                                ),
                                conditionalPanel(
                                    condition = "input.generateAnaInfo == \"example\"",
                                    fillRow(
                                        height = 30,
                                        textNote("Make sure that the patRoonData package is installed.")
                                    )
                                )
                            ),
                            fillCol(
                                flex = NA,
                                conditionalPanel(
                                    condition = anaTableCondition("both"),
                                    radioButtons("currentSetStatic", "Define table for set", c("positive", "negative"),
                                                 inline = TRUE)
                                ),
                                conditionalPanel(
                                    condition = anaTableCondition(NULL),
                                    getAnaTableFileUI(NULL)
                                ),
                                conditionalPanel(
                                    condition = anaTableCondition("positive"),
                                    getAnaTableFileUI("positive")
                                ),
                                conditionalPanel(
                                    condition = anaTableCondition("negative"),
                                    getAnaTableFileUI("negative")
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"dynamic\"",
                            fillCol(
                                flex = NA,
                                conditionalPanel(
                                    condition = anaDynCondition("both"),
                                    radioButtons("currentSetDynamic", "Define table for set", c("positive", "negative"),
                                                 inline = TRUE)
                                ),
                                strong(em("Generate analysis information from the following directories")),
                                br(),
                                conditionalPanel(
                                    condition = anaDynCondition(NULL),
                                    genAnaDynUI(NULL)
                                ),
                                conditionalPanel(
                                    condition = anaDynCondition("positive"),
                                    genAnaDynUI("positive")
                                ),
                                conditionalPanel(
                                    condition = anaDynCondition("negative"),
                                    genAnaDynUI("negative")
                                )
                            )
                        ),
                        conditionalPanel(
                            height = 30,
                            condition = paste(anaTableCondition(), "|| input.generateAnaInfo == \"dynamic\""),
                            textNote(paste("Raw analyses are instrument files (.raw, .d etc);",
                                           "others are mzXML/mzML files.",
                                           "Exporting data is configured in Data pre-treatment."))
                        ),
                        conditionalPanel(
                            condition = anaTableCondition(),
                            rhandsontable::rHandsontableOutput("analysesHot")
                        )
                    )
                ),
                
                conditionalPanel(
                    condition = anaTableCondition(),
                    miniUI::miniButtonBlock(
                        actionButton("setAnaInfoFiles", label = "Files", icon = icon("file-import"),
                                     title = "Set analysis files"),
                        shinyjs::disabled(actionButton("removeAnaInfoRows", label = "Remove", icon = icon("trash"),
                                                       title = "Remove selected row(s)")),
                        shinyjs::disabled(actionButton("anaInfoRowsUp", label = "Move up", icon = icon("arrow-up"),
                                                       title = "Move selected row(s) up")),
                        shinyjs::disabled(actionButton("anaInfoRowsDown", label = "Move down", icon = icon("arrow-down"),
                                                       title = "Move selected row(s) down")),
                        actionButton("importAnaInfoCSV", label = "Import CSV", icon = icon("file-csv"),
                                     title = "Import previously generated analyses information from a csv file")
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Data pre-treatment", icon = icon("upload"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = c(NA, 1),
                        fillRow(
                            height = 30,
                            strong("MS data conversion steps")
                        ),
                        fillRow(
                            rhandsontable::rHandsontableOutput("MSConversionHot")
                        )
                    )
                ),
                miniUI::miniButtonBlock(
                    actionButton("addMSConversion", label = "Add", icon = icon("plus"),
                                 title = "Add MS conversion step"),
                    shinyjs::disabled(actionButton("removeMSConversionRows", label = "Remove", icon = icon("trash"),
                                                   title = "Remove selected row(s)")),
                    shinyjs::disabled(actionButton("MSConversionRowsUp", label = "Move up", icon = icon("arrow-up"),
                                                   title = "Move selected row(s) up")),
                    shinyjs::disabled(actionButton("MSConversionRowsDown", label = "Move down", icon = icon("arrow-down"),
                                                   title = "Move selected row(s) down")),
                    actionButton("configConvPaths", label = "Output paths", icon = icon("folder-open")),
                    actionButton("configBrukerCalib", label = "Bruker m/z calibration", icon = icon("cog"),
                                 title = "Configure Bruker DataAnalysis m/z calibration")
                )
            ),
            miniUI::miniTabPanel(
                "Features", icon = icon("chart-area"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        fillCol(
                            height = 90,
                            flex = NA,
                            fillRow(
                                selectInput("featFinder", "Feature finder", c("OpenMS", "XCMS", "enviPick", "SIRIUS", "KPIC2",
                                                                              "Bruker DataAnalysis" = "Bruker"),
                                            "OpenMS", FALSE, width = "95%"),
                                fillCol(
                                    flex = c(1, NA),
                                    height = 90,
                                    selectInput("featGrouper", "Feature grouper", c("OpenMS", "XCMS", "KPIC2", "SIRIUS"),
                                                "OpenMS", FALSE, width = "100%"),
                                    conditionalPanel(
                                        condition = "input.featGrouper == \"SIRIUS\"",
                                        textNote(HTML("This will always find <b>and</b> group features with SIRIUS."))
                                    )
                                )
                            )
                        ),
                        fillCol(
                            height = 125,
                            flex = NA,
                            conditionalPanel(
                                condition = "input.ionization != \"both\"",
                                fileSelect("suspectList", "suspectListButton", "Suspect list", placeholder = "Leave empty for no suspect screening")
                            ),
                            conditionalPanel(
                                condition = "input.ionization == \"both\"",
                                fillRow(
                                    height = 60,
                                    fillCol(
                                        width = "95%",
                                        fileSelect("suspectListPos", "suspectListButtonPos", "Suspect list (positive)",
                                                   placeholder = "Leave empty for no suspect screening")
                                    ),
                                    fillCol(
                                        width = "95%",
                                        fileSelect("suspectListNeg", "suspectListButtonNeg", "Suspect list (negative)",
                                                   placeholder = "Leave empty if same as positive")
                                    )
                                )
                            ),
                            fillRow(
                                height = 50,
                                checkboxInput("exSuspList", "Example suspect list(s)")
                            )
                        ),
                        hr(),
                        fillCol(
                            flex = NA,
                            height = 70,
                            strong("Post-Filtering of feature groups"),
                            textNote("Set below values to zero to disable a particular filter.")
                        ),
                        fillCol(
                            height = 325,
                            fillRow(
                                numericInput("preIntThr", "Pre-Intensity threshold", 1E2, 0, step = 100, width = "95%"),
                                numericInput("intThr", "Intensity threshold", 1E4, 0, step = 1000, width = "100%")
                            ),
                            fillRow(
                                numericInput("repAbundance", "Min. replicate abundance (relative)", 1, 0, 1.0, 0.1, width = "95%"),
                                numericInput("maxRepRSD", "Max. replicate intensity RSD", 0.75, 0, step = 0.1, width = "100%")
                            ),
                            fillRow(
                                numericInput("blankThr", "Min. blank threshold", 5, 0, step = 1, width = "95%"),
                                checkboxInput("removeBlanks", "Discard blanks after filtering", TRUE)
                            ),
                            rangeNumeric("retention", "retention time (s)", step = 10),
                            rangeNumeric("mz", "m/z", step = 10)
                        ),
                        hr(),
                        fillCol(
                            height = 125,
                            flex = NA,
                            
                            selectInput("featNorm", "Feature normalization", c("None" = "none",
                                                                               "Internal standard" = "istd",
                                                                               "Internal standard concentration" = "conc",
                                                                               "TIC" = "tic")),
                            checkboxInput("groupNorm", "Group normalization"),

                            conditionalPanel(
                                condition = "input.featNorm == \"istd\" && input.ionization != \"both\"",
                                fileSelect("ISTDList", "ISTDListButton", "Internal standard list",
                                           placeholder = "Leave empty for example list")
                            ),
                            conditionalPanel(
                                condition = "input.featNorm == \"istd\" && input.ionization == \"both\"",
                                fillRow(
                                    height = 60,
                                    fillCol(
                                        width = "95%",
                                        fileSelect("ISTDListPos", "ISTDListButtonPos", "Internal standard list (positive)",
                                                   placeholder = "Leave empty for example list")
                                    ),
                                    fillCol(
                                        width = "95%",
                                        fileSelect("ISTDListNeg", "ISTDListButtonNeg", "Internal standard list (negative)",
                                                   placeholder = "Leave empty if same as positive")
                                    )
                                )
                            ),
                            conditionalPanel(
                                condition = "input.featNorm == \"istd\" || input.featNorm == \"conc\"",
                                textNote(HTML("Please make sure that the <i>norm_conc</i> column of the analysis information is set."))
                            )
                        )
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Annotation", icon = icon("chart-bar"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        
                        fillCol(
                            height = 100,
                            selectInput("components", "Component generation",
                                        c("None" = "", "RAMClustR", "CAMERA", "OpenMS", "CliqueMS", "nontarget"),
                                        multiple = FALSE, width = "100%"),
                            conditionalPanel(
                                condition = "input.components != \"nontarget\" && input.components != \"\"",
                                checkboxInput("selectIons", "Select feature adduct ions", value = TRUE)
                            )
                        ),
                        fillRow(
                            height = 90,
                            fillCol(
                                flex = c(1, NA),
                                selectInput("formulaGen", "Formula generation",
                                            c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                                            multiple = FALSE, width = "95%"),
                                textNote("DataAnalysis only works with features from DataAnalysis")
                            ),
                            selectInput("compIdent", "Compound identification",
                                        c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag", "Library"),
                                        multiple = FALSE, width = "100%")
                        ),
                        conditionalPanel(
                            condition = "input.formulaGen != \"\" || input.compIdent != \"\"",
                            fillRow(
                                height = 110,
                                fillCol(
                                    selectInput("peakListGen", "Peak list generator",
                                                c("mzR", "Bruker DataAnalysis" = "Bruker"),
                                                "mzR", multiple = FALSE, width = "95%"),
                                    checkboxInput("DIA", "Data Independent Acquisition (DIA) MS/MS")
                                ),
                                conditionalPanel(
                                    condition = "input.peakListGen == \"mzR\" && !input.DIA",
                                    numericInput("precursorMzWindow", "MS/MS precursor m/z search window", 4, width = "100%"),
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.compIdent == \"Library\"",
                            fillRow(
                                height = 90,
                                fillCol(
                                    width = "95%",
                                    selectInput("MSLibraryFormat", "Library format",
                                                c("MSP" = "msp", "MoNA JSON" = "json"), multiple = FALSE)
                                ),
                                fillCol(
                                    width = "95%",
                                    fileSelect("MSLibraryPath", "MSLibraryPathButton", "MS library path")
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "(input.exSuspList || (input.ionization != \"both\" && input.suspectList != \"\") || (input.ionization == \"both\" && input.suspectListPos != \"\")) && (input.formulaGen != \"\" || input.compIdent != \"\")",
                            fillRow(
                                height = 90,
                                fillCol(
                                    strong("Suspect annotation"),
                                    checkboxInput("annotateSus", "Annotate suspects", TRUE, width = "100%"),
                                    conditionalPanel(
                                        condition = "input.annotateSus",
                                        checkboxInput("genIDLevelFile", "Generate template file with configurable identification levels",
                                                      TRUE, width = "100%")
                                    ),
                                    textNote("Suspect annotation is currently only optimized for GenForm/MetFrag")
                                )
                            )
                        )
                    )
                )
            ),
            miniUI::miniTabPanel(
                "TP screening", icon = icon("react"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        checkboxInput("doTPs", "Perform transformation product screening"),
                        conditionalPanel(
                            condition = "input.doTPs",
                            selectInput("TPGen", "TP algorithm", c("BioTransformer", "CTS", "Library", "Logic")),
                            conditionalPanel(
                                condition = "input.TPGen != \"Logic\"",
                                selectInput("TPGenInput", "Parent input", getTPGenInputs(FALSE)),
                                conditionalPanel(
                                    condition = "input.TPGenInput == \"suspects\"",
                                    fileSelect("TPSuspectList", "TPSuspButton", "Parent suspect list",
                                               placeholder = "Please specify parent suspect list")
                                ),
                                checkboxInput("TPDoMFDB", "Generate TP MetFrag database", TRUE)
                            )
                        )
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Reporting", icon = icon("file-medical-alt"),
                miniUI::miniContentPanel(
                    fillRow(
                        # flex = NA,
                        checkboxGroupInput("reportGen", "Report generation", c("HTML reports" = "HTML",
                                                                               "Legacy interface" = "legacy"),
                                           "HTML", width = "100%"),
                        conditionalPanel(
                            condition = "input.reportGen.includes('legacy')",
                            checkboxGroupInput("reportLegacy", "Legacy report formats", c("CSV", "PDF"),
                                               "CSV", width = "100%")
                        )
                    )
                )
            )
        )
    )
}

getNewProjectWidgetTypes <- function()
{
    list(
        outputScriptTo = "radio",
        scriptFile = "text",
        createRStudioProj = "check",
        ionization = "radio",
        analysisTableFile = "text",
        analysisTableFilePos = "text",
        analysisTableFileNeg = "text",
        convAlgo = "select",
        convFrom = "select",
        convTo = "select",
        peakPicking = "check",
        peakPickingVendor = "check",
        DAMethod = "text", # UNDONE: remove
        DAMethodPos = "text", # UNDONE: remove
        DAMethodNeg = "text", # UNDONE: remove
        doBrukerCalib = "check", # UNDONE: remove
        featFinder = "select",
        featGrouper = "select",
        suspectList = "text",
        suspectListPos = "text",
        suspectListNeg = "text",
        exSuspList = "check",
        preIntThr = "numeric",
        intThr = "numeric",
        repAbundance = "numeric",
        maxRepRSD = "numeric",
        blankThr = "numeric",
        "retention-min" = "numeric",
        "retention-max" = "numeric",
        "mz-min" = "numeric",
        "mz-max" = "numeric",
        removeBlanks = "check",
        featNorm = "select",
        groupNorm = "check",
        ISTDList = "text",
        ISTDListPos = "text",
        ISTDListNeg = "text",
        components = "select",
        selectIons = "check",
        formulaGen = "select",
        compIdent = "select",
        peakListGen = "select",
        DIA = "check",
        precursorMzWindow = "numeric",
        MSLibraryFormat = "select",
        MSLibraryPath = "text",
        annotateSus = "check",
        genIDLevelFile = "check",
        doTPs = "check",
        TPGen = "select",
        TPGenInput = "select",
        TPSuspectList = "text",
        TPDoMFDB = "check",
        reportGen = "checkGroup",
        reportLegacy = "checkGroup"
    )
}

loadNewProjectParams <- function(file, input, session)
{
    wtypes <- getNewProjectWidgetTypes()
    values <- readYAML(file)
    for (param in names(values))
    {
        if (wtypes[[param]] == "radio")
            updateRadioButtons(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "text")
            updateTextInput(session, param, value = values[[param]])
        else if (wtypes[[param]] == "check")
            updateCheckboxInput(session, param, value = values[[param]])
        else if (wtypes[[param]] == "checkGroup")
            updateCheckboxGroupInput(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "select")
            updateSelectInput(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "numeric")
            updateNumericInput(session, param, value = values[[param]])
    }
}

saveNewProjectParams <- function(file, input)
{
    values <- getNewProjectWidgetTypes()
    values <- isolate(reactiveValuesToList(input))[names(values)]
    writeYAML(values, file)
}

anaFilesToAnaInfo <- function(anaFiles)
{
    ret <- dcast(anaFiles, analysis ~ type, value.var = "path")
    setnames(ret, getMSFileTypes(), paste0("path_", getMSFileTypes()), skip_absent = TRUE)
    for (ft in getMSFileTypes())
    {
        col <- paste0("path_", ft)
        if (is.null(ret[[col]]))
            set(ret, j = col, value = character(0))
    }
    return(ret)
}

anaInfoToAnaFiles <- function(anaInfo)
{
    ret <- melt(anaInfo, id.vars = "analysis", measure.vars = getAnaInfoPathCols(anaInfo), variable.name = "type",
                value.name = "path", na.rm = TRUE)
    ret[, type := gsub("^path_", "", type)]
    return(ret[])
}


groupConvTypesFormats <- function(algo, inOut)
{
    fTypes <- intersect(getMSFileTypes(), if (inOut == "input") validConvFromTypes(algo) else validConvToTypes(algo))
    return(sapply(fTypes, function(ft)
    {
        formats <- if (inOut == "input") getMSInConversionFormats(algo, ft) else getMSOutConversionFormats(algo, ft)
        return(setNames(paste0(ft, "_", formats), paste0(formats, " (", ft, ")")))
    }, simplify = FALSE))
}

splitConvTypeFormat <- function(typeFormat) strsplit(typeFormat, "_")

convTypeFormatToLabel <- function(typeFormat)
{
    tf <- splitConvTypeFormat(typeFormat)
    return(paste0(tf[[1]][1], " (", tf[[1]][2], ")"))
}

#' Easily create new \pkg{patRoon} projects
#' 
#' The \code{newProject} function is used to quickly generate a processing R script. This tool allows the user to
#' quickly select the targeted analyses, workflow steps and configuring some of their common parameters. This function
#' requires to be run within a \href{https://www.rstudio.com/}{RStudio} session. The resulting script is either added to
#' the current open file or to a new file. The \link[=analysis-information]{analysis information} will be written to a
#' \file{.csv} file so that it can easily be modified afterwards.
#'
#' @param destPath Set destination path value to this value (useful for debugging). Set to \code{NULL} for a default
#'   value.
#'
#' @export
newProject <- function(destPath = NULL)
{
    rstudioapi::verifyAvailable()

    # UNDONE: warning/message about empty groups

    # NOTE: disable column sorting so we don't have to worry about getting correct row index
    # (https://github.com/jrowen/rhandsontable/issues/257)
    # NOTE: preventOverflow should not be set when columnSorting = FALSE
    # https://github.com/handsontable/handsontable/issues/4303
    # NOTE: set selectionMode to range as only row series can currently be queried
    # (https://github.com/jrowen/rhandsontable/issues/313)
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE,
                    sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    selectionMode = "range", outsideClickDeselects = FALSE,
                    contextMenu = FALSE, manualColumnResize = TRUE)

    emptyAnaTable <- function()
    {
        ai <- data.table(analysis = character(0), type = character(0), replicate = character(0), blank = character(0),
                         conc = numeric(0), norm_conc = numeric(0))
        ai[, paste0("path_", getMSFileTypes()) := character(0)]
        return(ai)
    }
    
    anaInfoTabs <- list(positive = emptyAnaTable(), negative = emptyAnaTable())
    MSConvSettings <- list(
        steps = data.table(algorithm = character(0), from = character(0), to = character(0)),
        output = list(centroid = file.path("converted", "centroid"), profile = file.path("converted", "profile"),
                      ims = file.path("converted", "ims")),
        brukerCalib <- list(enabled = FALSE, method = "", methodPos = "", methodNeg = "")
    )
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(triggerAnaInfoHotUpdate = 0,
                                  analysisFiles = data.table(analysis = character(), type = character(), path = character()),
                                  triggerMSConvHotUpdate = 0)

        doObserveSelDir <- function(textID, buttonID)
        {
            observeEvent(input[[buttonID]], {
                d <- rstudioapi::selectDirectory("Select directory", path = input[[textID]])
                if (!is.null(d))
                    updateTextInput(session, textID, value = d)
            })
        }
        makeAnalysesHot <- function()
        {
            rValues$triggerAnaInfoHotUpdate
            
            dt <- copy(anaInfoTabs[[getCurPolarityAnaTab()]])
            pcols <- getAnaInfoPathCols(dt)
            dt[, type := {
                p <- unlist(mget(pcols))
                paste0(sub("^path_", "", names(p)[!is.na(p)]), collapse = ", ")
            }, by = .I]
            dt[, (pcols) := NULL]
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(dt, height = 250, maxRows = nrow(dt), columnSorting = FALSE),
                             hotOpts)) %>%
                rhandsontable::hot_col(c("replicate", "blank"), readOnly = FALSE, type = "text") %>%
                rhandsontable::hot_col(c("conc", "norm_conc"), readOnly = FALSE, type = "numeric")
            
            return(hot)
        }
        triggerAnaInfoHotUpdate <- function() rValues$triggerAnaInfoHotUpdate <- rValues$triggerAnaInfoHotUpdate + 1
        makeAnalysisFilesHot <- function()
        {
            tab <- copy(rValues$analysisFiles)
            tab[, format := mapply(analysis, path, type, FUN = function(a, p, t) {
                exts <- unique(MSFileExtensions()[getMSFileFormats(t)])
                paths <- file.path(p, paste0(a, ".", exts))
                exts <- exts[file.exists(paths)]
                if (length(exts) == 0)
                    "not found!"
                else
                    paste0(exts, collapse = ", ")
            })]
            setcolorder(tab, c("analysis", "type", "format", "path"))
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(tab, height = 300, maxRows = nrow(rValues$analysisFiles), rowHeaders = NULL),
                             hotOpts))
            return(hot)
        }
        getCurPolarityAnaTab <- function()
        {
            if (input$ionization == "positive" || (input$ionization == "both" && input$currentSetStatic == "positive"))
                return("positive")
            return("negative")
        }
        moveSelectedRows <- function(dir, sel, tab)
        {
            # NOTE: assume selection is a block range
            mv <- if (dir == "up")
                seq(min(sel) - 1, max(sel))
            else
                seq(min(sel), max(sel) + 1)
            
            tab <- copy(tab)
            tab[, index := .I]
            tab[mv, index := shift(index, n = if (dir == "up") - 1L else 1L, type = "cyclic")]
            tab <- tab[(index), -"index"]
            return(tab)
        }
        moveTabSel <- function(dir, sel, id)
        {
            newRange <- c(min(sel), max(sel)) + if (dir == "up") -1 else 1
            session$sendCustomMessage("selectHotRows", list(tab = id, range = newRange))
        }
        moveSelectedAnalyses <- function(dir)
        {
            # NOTE: assume selection is a block range
            sel <- input$analysesHot_select$select$rAll
            pol <- getCurPolarityAnaTab()
            anaInfoTabs[[pol]] <<- moveSelectedRows(dir, sel, anaInfoTabs[[pol]])
            triggerAnaInfoHotUpdate()
            moveTabSel(dir, sel, "analysesHot")
        }
        doObserveGenAnaInfoDynSelDir <- function(textID, buttonID)
        {
            observeEvent(input[[buttonID]], {
                pf <- if (input$ionization != "both") "" else if (input$currentSetDynamic == "positive") "Pos" else "Neg"
                textID <- paste0(textID, pf)
                d <- rstudioapi::selectDirectory("Select directory", path = input[[textID]])
                if (!is.null(d))
                    updateTextInput(session, textID, value = d)
            })
        }
        verifyAnalysesOK <- function()
        {
            verifyAny <- function(pol)
            {
                if (nrow(anaInfoTabs[[pol]]) == 0 && input$ionization %in% c(pol, "both"))
                {
                    rstudioapi::showDialog("No analyses", paste0("Please select some analyses for ", pol, " mode!"), "")
                    return(FALSE)
                }
            }
            verifyAny("positive"); verifyAny("negative")

            if (input$generateAnaInfo == "table")
            {
                n <- paste0("analysisTableFile", input$analysisTableFileType)
                checkAnas <- if (input$ionization != "both")
                    input[[n]]
                else
                    input[paste0(n, c("Pos", "Neg"))]
                for (f in checkAnas)
                {
                    p <- file.path(input$destinationPath, f)
                    if (file.exists(p))
                    {
                        ov <- rstudioapi::showQuestion("Analysis table file already exists",
                                                       sprintf("Analysis table file already exists: '%s'.\nOverwrite?", p),
                                                       "Yes", "No")
                        if (!ov)
                            return(FALSE)
                    }
                }
            }
            
            return(TRUE)
        }
        makeMSConversionHot <- function()
        {
            rValues$triggerMSConvHotUpdate
            tab <- copy(MSConvSettings$steps)
            tab[, from := sapply(from, convTypeFormatToLabel)]
            tab[, to := sapply(to, convTypeFormatToLabel)]
            do.call(rhandsontable::rhandsontable, c(list(tab, height = 250, maxRows = nrow(tab), columnSorting = FALSE),
                                                    hotOpts))
        }
        moveSelectedConversions <- function(dir)
        {
            # NOTE: assume selection is a block range
            sel <- input$MSConversionHot_select$select$rAll
            MSConvSettings$steps <<- moveSelectedRows(dir, sel, MSConvSettings$steps)
            triggerMSConvHotUpdate()
            moveTabSel(dir, sel, "MSConversionHot")
        }
        selectDAMethod <- function(inputName)
        {
            dm <- rstudioapi::selectDirectory("Select DataAnalysis method")
            if (!is.null(dm))
            {
                if (!file.exists(file.path(dm, "DataAnalysis.method")))
                {
                    rstudioapi::showDialog("Invalid DataAnalysis method", "Please select a valid DataAnalysis method!", "")
                    dm <- NULL
                }
            }
            
            if (!is.null(dm))
                updateTextInput(session, inputName, value = dm)
        }
        triggerMSConvHotUpdate <- function() rValues$triggerMSConvHotUpdate <- rValues$triggerMSConvHotUpdate + 1
        selectSuspList <- function(inputName)
        {
            sl <- rstudioapi::selectFile("Select suspect list", filter = "csv files (*.csv)")
            if (!is.null(sl))
            {
                csvTab <- tryCatch(fread(sl), error = function(e) FALSE, warning = function(w) FALSE)
                cols <- names(csvTab)
                massCols <- c("mz", "neutralMass", "formula", "SMILES", "InChI")
                
                err <- NULL
                if (is.logical(csvTab))
                    err <- "Failed to open/parse selected csv file!"
                else if (nrow(csvTab) == 0)
                    err <- "The selected files seems to be empty."
                else if (!"name" %in% cols)
                    err <- "The selected file does not have a name column"
                else if (!any(massCols %in% cols))
                    err <- paste("The selected file should have at least one of the columns:",
                                 paste(massCols, collapse = ", "))
                
                if (!is.null(err))                
                    rstudioapi::showDialog("Error", err, "")
                else
                    updateTextInput(session, inputName, value = sl)
            }
        }
        
        observeEvent(input$create, {
            if (!nzchar(input$destinationPath))
                rstudioapi::showDialog("Invalid destination", "Please select a destination path!", "")
            else if (input$outputScriptTo != "curFile" && !nzchar(input$scriptFile))
                rstudioapi::showDialog("No script file", "Please select a destination script file!", "")
            else if (input$generateAnaInfo == "table" && !verifyAnalysesOK())
            {}
            else if (input$outputScriptTo != "curFile" && file.exists(file.path(input$destinationPath, input$scriptFile)) &&
                     !rstudioapi::showQuestion("Script file already exists",
                                               sprintf("Script file already exists: '%s'.\nOverwrite?",
                                                       file.path(input$destinationPath, input$scriptFile)),
                                               "Yes", "No"))
            {}
            else if (input$compIdent == "Library" && !nzchar(input$MSLibraryPath))
                rstudioapi::showDialog("No MS library", "Please select an MS library!", "")
            else if (input$doTPs && input$TPGen != "Logic" && input$TPGenInput == "suspects" &&
                     !nzchar(input$TPSuspectList))
                rstudioapi::showDialog("No parent suspect list", "Please select a parent suspect list!", "")
            else if ("legacy" %in% input$reportGen && length(input$reportLegacy) == 0)
                rstudioapi::showDialog("No legacy format", "Please select at least one legacy reporting format!", "")
            else
            {
                doCreateProject(input, anaInfoTabs, MSConvSettings$steps)
                stopApp(TRUE)
            }
        })

        doObserveSelDir("destinationPath", "projectDestButton")
        
        observeEvent(input$loadParams, {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)")
            if (!is.null(sl))
                loadNewProjectParams(sl, input, session)
        })
        
        observeEvent(input$saveParams, {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)", existing = FALSE)
            if (!is.null(sl))
                saveNewProjectParams(sl, input)
        })
        
        observeEvent(input$analysesHot, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (!is.null(input$analysesHot) && input$analysesHot$params$maxRows > 0)
            {
                # sync from table edits
                dt <- rhandsontable::hot_to_r(input$analysesHot)
                pol <- getCurPolarityAnaTab()
                ai <- anaInfoTabs[[pol]]
                mcols <- c("replicate", "blank", "conc", "norm_conc")
                ai[, (mcols) := dt[, mcols, with = FALSE]]
                anaInfoTabs[[pol]] <<- ai
            }
            else
            {
                shinyjs::disable("removeAnaInfoRows")
                shinyjs::disable("anaInfoRowsUp")
                shinyjs::disable("anaInfoRowsDown")
            }
        })
        
        observeEvent(input$analysesHot_select$select$r, {
            pol <- getCurPolarityAnaTab()
            sel <- input$analysesHot_select$select$rAll
            e <- nrow(anaInfoTabs[[pol]]) > 0 && length(sel) > 0
            shinyjs::toggleState("removeAnaInfoRows", condition = e)
            shinyjs::toggleState("anaInfoRowsUp", condition = e && min(sel) > 1)
            shinyjs::toggleState("anaInfoRowsDown", condition = e && max(sel) < nrow(anaInfoTabs[[pol]]))
        })
        
        observeEvent(input$setAnaInfoFiles, {
            pol <- getCurPolarityAnaTab()
            rValues$analysisFiles <- if (nrow(anaInfoTabs[[pol]]) == 0)
                rValues$analysisFiles[0, ]
            else
                anaInfoToAnaFiles(anaInfoTabs[[pol]])
            
            showModal(modalDialog(
                title = "Select analyses",
                rhandsontable::rHandsontableOutput("analysisFilesHot", height = 300, width = 700),
                easyClose = TRUE,
                fade = FALSE, # TRUE messes up HOT
                footer = tagList(
                    # use column() to (1) make sure selectInput only occupies a single row and (2) a button can
                    # be placed right next to it and (3) left alignment
                    column(
                        width = 2,
                        style = "padding-right: 0px;",
                        htmltools::tagAppendAttributes(
                            selectInput("analysisFilesAddType", NULL, c("raw", "centroid", "profile", "ims"),
                                        selectize = FALSE, width = "100%"),
                            style = "margin-bottom: 0px;")
                    ),
                    column(
                        width = 1,
                        style = "padding-left: 0px;",
                        title = "Add all analysis files with the selected type within a directory",
                        actionButton("analysisFilesAdd", label = NULL, icon = icon("plus"))
                    ),
                    column(
                        width = 2,
                        style = "padding-right: 0px;",
                        htmltools::tagAppendAttributes(
                            selectInput("analysisFilesChangeType", NULL, c("centroid", "profile", "ims"),
                                        selectize = FALSE, width = "100%"),
                            style = "margin-bottom: 0px;")
                    ),
                    column(
                        width = 1,
                        style = "padding-left: 0px;",
                        title = "Change type of selected analysis files",
                        shinyjs::disabled(actionButton("analysisFilesChange", label = NULL, icon = icon("pen-to-square")))
                    ),
                    column(
                        width = 1,
                        shinyjs::disabled(actionButton("analysisFilesRemove", label = NULL, icon = icon("trash")))
                    ),
                    modalButton("Cancel"),
                    actionButton("analysisFilesOK", "OK")
                )
            ))
        })
        
        observeEvent(input$removeAnaInfoRows, {
            sel <- input$analysesHot_select$select$rAll
            pol <- getCurPolarityAnaTab()
            anaInfoTabs[[pol]] <<- anaInfoTabs[[pol]][-sel]
            triggerAnaInfoHotUpdate()
        })

        observeEvent(input$anaInfoRowsUp, { moveSelectedAnalyses("up") })
        observeEvent(input$anaInfoRowsDown, { moveSelectedAnalyses("down") })
        
        observeEvent(input$analysisFilesHot, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (is.null(input$analysesHot) || input$analysisFilesHot$params$maxRows == 0)
            {
                shinyjs::disable("analysisFilesChangeType")
                shinyjs::disable("analysisFilesChange")
                shinyjs::disable("analysisFilesRemove")
            }
        })
        observeEvent(input$analysisFilesHot_select$select$r, {
            if (nrow(rValues$analysisFiles) > 0 && length(input$analysisFilesHot_select$select$rAll) > 0)
            {
                sel <- input$analysisFilesHot_select$select$rAll
                rawSelected <- "raw" %in% rValues$analysisFiles[sel, "type"]
                shinyjs::toggleState("analysisFilesChangeType", condition = !rawSelected)
                shinyjs::toggleState("analysisFilesChange", condition = !rawSelected)
                
                # Disallow IMS selection for mzXML files (not supported by format)
                mzXMLSelected <- "mzxml" %in% tolower(tools::file_ext(rValues$analysisFiles[sel ,"path"]))
                ch <- c("centroid", "profile", if (!mzXMLSelected) "ims")
                s <- rValues$analysisFiles[sel ,"type"]
                if (uniqueN(s) > 1)
                    s <- ch[1] # just default to first choice if selection isn't homogenous
                updateSelectInput(session, "analysisFilesChangeType", selected = s, choices = ch)
                
                shinyjs::enable("analysisFilesRemove")
            }
            else
            {
                shinyjs::disable("analysisFilesChange")
                shinyjs::disable("analysisFilesRemove")
            }
        })
        
        observeEvent(input$analysisFilesAdd, {
            anaDir <- rstudioapi::selectDirectory(path = "~/")
            if (!is.null(anaDir))
            {
                files <- listMSFiles(anaDir, getMSFileFormats(input$analysisFilesAddType))
                af <- data.table(analysis = tools::file_path_sans_ext(basename(files)),
                                 type = if (length(files) > 0) input$analysisFilesAddType else character(),
                                 path = dirname(files))
                af <- unique(af) # in case of multiple files with the same analysis name and path (eg mzXML+mzML)
                rValues$analysisFiles <- rbind(rValues$analysisFiles, af)
            }
        })
        
        observeEvent(input$analysisFilesChange, {
            sel <- input$analysisFilesHot_select$select$rAll
            rValues$analysisFiles$type[sel] <- input$analysisFilesChangeType
        })
        
        observeEvent(input$analysisFilesRemove, {
            sel <- input$analysisFilesHot_select$select$rAll
            rValues$analysisFiles <- rValues$analysisFiles[-sel]
        })
        
        observeEvent(input$analysisFilesOK, {
            removeModal()
            pol <- getCurPolarityAnaTab()
            
            if (nrow(rValues$analysisFiles) > 0)
            {
                ai <- anaFilesToAnaInfo(rValues$analysisFiles)
                
                # remove the analyses for which all file paths were removed
                anaInfoTabs[[pol]] <<- anaInfoTabs[[pol]][analysis %in% rValues$analysisFiles$analysis]
    
                # overlap: update paths
                aiOv <- ai[analysis %in% anaInfoTabs[[pol]]$analysis]
                if (nrow(aiOv) > 0)
                {
                    pcols <- getAnaInfoPathCols(anaInfoTabs[[pol]])
                    anaInfoTabs[[pol]][aiOv, (pcols) := mget(paste0("i.", pcols)), on = "analysis"]
                }

                # add analyses for new files                
                aiNew <- ai[!analysis %in% anaInfoTabs[[pol]]$analysis]
                if (nrow(aiNew) > 0)
                    anaInfoTabs[[pol]] <<- rbind(anaInfoTabs[[pol]], aiNew, fill = TRUE)
            }
            else
                anaInfoTabs[[pol]] <<- anaInfoTabs[[pol]][0]
            
            triggerAnaInfoHotUpdate()
        })
        
        observeEvent(input$importAnaInfoCSV, {
            csvFile <- rstudioapi::selectFile(path = "~/", filter = "csv files (*.csv)")
            if (!is.null(csvFile))
            {
                ai <- tryCatch(assertAndPrepareAnaInfo(fread(csvFile)), error = function(e)
                {
                    rstudioapi::showDialog("Invalid CSV file",
                                           sprintf("The selected CSV file is not a valid analysis information file: '%s'", e))
                    return(FALSE)
                })
                
                if (!isFALSE(ai) && nrow(ai) > 0)
                {
                    pol <- getCurPolarityAnaTab()
                    
                    # check for overlap
                    if (any(ai$analysis %in% anaInfoTabs[[pol]]$analysis))
                    {
                        overwrite <- rstudioapi::showQuestion("Overwrite analysis information",
                                                              "One or more analyses are already defined. Do you want to overwrite or keep the current information for these analyses?",
                                                              "Overwrite", "Keep", timeout = NULL)
                        if (overwrite) # UNDONE: try to keep original order?
                            anaInfoTabs[[pol]] <<- anaInfoTabs[[pol]][!analysis %in% ai$analysis]
                        else
                            ai <- ai[!analysis %in% anaInfoTabs[[pol]]$analysis]
                    }
                    
                    anaInfoTabs[[pol]] <<- rbind(anaInfoTabs[[pol]], ai, fill = TRUE)
                    triggerAnaInfoHotUpdate()
                }
            }
        })
        
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynRaw", "genAnaInfoDynRawButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynCentroid", "genAnaInfoDynCentroidButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynIMS", "genAnaInfoDynIMSButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynProfile", "genAnaInfoDynProfileButton")
        
        observeEvent(input$addMSConversion, {
            showModal(modalDialog(
                title = "Add MS data conversion step",
                fillCol(
                    flex = 1,
                    height = 150,
                    width = 600,
                    fillCol(
                        selectInput("convAlgo", "Algorithm",
                                    c("ProteoWizard" = "pwiz", "OpenMS" = "openms", "Bruker" = "bruker",
                                      "IM collapse" = "im_collapse", "TIMSCONVERT" = "timsconvert"),
                                    multiple = FALSE, width = "100%")
                    ),
                    fillRow(
                        flex = c(1, NA, 1),
                        selectInput("convFrom", "From", groupConvTypesFormats("pwiz", "input"), multiple = FALSE,
                                    width = "100%"),
                        fillCol(width = "20px"),
                        selectInput("convTo", "To", groupConvTypesFormats("pwiz", "output"), multiple = FALSE,
                                    width = "100%")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton("addMSConversionOK", "OK")
                )
            ))
        })
        doObserveSelDir("MSConversionPath", "MSConversionPathButton")
        observeEvent(input$configConvPaths, {
            showModal(modalDialog(
                title = "Set MS conversion output paths",
                fillCol(
                    flex = c(1, 1, 1, NA),
                    height = 275,
                    width = 400,
                    fileSelect("MSConvOutputCentroid", "MSConvOutputCentroidButton", "Centroid", MSConvSettings$output$centroid),
                    fileSelect("MSConvOutputProfile", "MSConvOutputProfileButton", "Profile", MSConvSettings$output$profile),
                    fileSelect("MSConvOutputIMS", "MSConvOutputIMSButton", "IMS", MSConvSettings$output$ims),
                    fillCol(
                        height = 25,
                        textNote("Relative paths are relative to the project directory.")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton("MSConvOutputOK", "OK")
                )
            ))
        })
        doObserveSelDir("MSConvOutputRaw", "MSConvOutputRawButton")
        doObserveSelDir("MSConvOutputCentroid", "MSConvOutputCentroidButton")
        doObserveSelDir("MSConvOutputProfile", "MSConvOutputProfileButton")
        doObserveSelDir("MSConvOutputIMS", "MSConvOutputIMSButton")
        
        observeEvent(input$configBrukerCalib, {
            showModal(modalDialog(
                title = "Configure Bruker DataAnalysis m/z re-calibration",
                fillCol(
                    flex = NA,
                    width = 600,
                    fillCol(
                        height = 50,
                        checkboxInput("doBrukerCalib", "Perform m/z re-calibration", value = MSConvSettings$brukerCalib$enabled)
                    ),
                    conditionalPanel(
                        condition = "input.doBrukerCalib",
                        fillCol(
                            flex = NA,
                            height = 90,
                            fillCol(
                                flex = NA,
                                height = 50,
                                conditionalPanel(
                                    condition = "input.ionization != \"both\"",
                                    fileSelect("DAMethod", "DAMethodButton", "DataAnalysis method",
                                               MSConvSettings$brukerCalib$method),
                                ),
                                conditionalPanel(
                                    condition = "input.ionization == \"both\"",
                                    fillRow(
                                        fillCol(
                                            width = "95%",
                                            fileSelect("DAMethodPos", "DAMethodButtonPos",
                                                       "DataAnalysis method (positive)",
                                                       MSConvSettings$brukerCalib$methodPos)
                                        ),
                                        fillCol(
                                            width = "95%",
                                            fileSelect("DAMethodNeg", "DAMethodButtonNeg",
                                                       "DataAnalysis method (negative)",
                                                       MSConvSettings$brukerCalib$methodNeg)
                                        )
                                    )
                                )
                            ),
                            textNote("Leaving this blank will not set any method.")
                        )
                    ),
                    fillCol(
                        height = 40,
                        textNote(HTML("Only supported with bruker data and if DataAnalysis is installed.<br>This step is always executed before other conversion steps."))
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton("brukerCalibOK", "OK")
                )
            ))
        })
        
        observeEvent(input$convAlgo, {
            f <- groupConvTypesFormats(input$convAlgo, "input")
            updateSelectInput(session, "convFrom", choices = f, selected = f[[1]][1])
            f <- groupConvTypesFormats(input$convAlgo, "output")
            updateSelectInput(session, "convTo", choices = f, selected = f[[1]][1])
        })
        observeEvent(input$DAMethodButton, selectDAMethod("DAMethod"))
        observeEvent(input$DAMethodButtonPos, selectDAMethod("DAMethodPos"))
        observeEvent(input$DAMethodButtonNeg, selectDAMethod("DAMethodNeg"))
        observeEvent(input$addMSConversionOK, {
            removeModal()
            MSConvSettings$steps <<- rbind(MSConvSettings$steps, data.table(algorithm = input$convAlgo, from = input$convFrom,
                                                                            to = input$convTo))
            triggerMSConvHotUpdate()
        })
        observeEvent(input$MSConversionHot, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (is.null(input$MSConversionHot) || input$MSConversionHot$params$maxRows == 0)
            {
                shinyjs::disable("removeMSConversionRows")
                shinyjs::disable("MSConversionRowsUp")
                shinyjs::disable("MSConversionRowsDown")
            }
        })
        observeEvent(input$MSConversionHot_select$select$r, {
            sel <- input$MSConversionHot_select$select$rAll
            e <- nrow(MSConvSettings$steps) > 0 && length(sel) > 0
            shinyjs::toggleState("removeMSConversionRows", condition = e)
            shinyjs::toggleState("MSConversionRowsUp", condition = e && min(sel) > 1)
            shinyjs::toggleState("MSConversionRowsDown", condition = e && max(sel) < nrow(MSConvSettings$steps))
        })
        observeEvent(input$removeMSConversionRows, {
            sel <- input$MSConversionHot_select$select$rAll
            MSConvSettings$steps <<- MSConvSettings$steps[-sel]
            triggerMSConvHotUpdate()
        })
        observeEvent(input$MSConversionRowsUp, { moveSelectedConversions("up") })
        observeEvent(input$MSConversionRowsDown, { moveSelectedConversions("down") })
        observeEvent(input$MSConvOutputOK, {
            removeModal()
            MSConvSettings$output$centroid <<- input$MSConvOutputCentroid
            MSConvSettings$output$profile <<- input$MSConvOutputProfile
            MSConvSettings$output$ims <<- input$MSConvOutputIMS
        })
        observeEvent(input$brukerCalibOK, {
            removeModal()
            MSConvSettings$brukerCalib <<- list(enabled = input$doBrukerCalib, method = input$DAMethod,
                                                methodPos = input$DAMethodPos, methodNeg = input$DAMethodNeg)
        })

        observeEvent(input$suspectListButton, selectSuspList("suspectList"))
        observeEvent(input$suspectListButtonPos, selectSuspList("suspectListPos"))
        observeEvent(input$suspectListButtonNeg, selectSuspList("suspectListNeg"))
        observeEvent(input$exSuspList, {
            for (id in c("suspectList", "suspectListButton", "suspectListPos", "suspectListButtonPos",
                         "suspectListNeg", "suspectListButtonNeg"))
                shinyjs::toggleState(id, !input$exSuspList)
        })
        
        observeEvent(input$ISTDListButton, selectSuspList("ISTDList"))
        observeEvent(input$ISTDListButtonPos, selectSuspList("ISTDListPos"))
        observeEvent(input$ISTDListButtonNeg, selectSuspList("ISTDListNeg"))
        
        output$analysesHot <- rhandsontable::renderRHandsontable(makeAnalysesHot())
        output$analysisFilesHot <- rhandsontable::renderRHandsontable(makeAnalysisFilesHot())
        
        output$MSConversionHot <- rhandsontable::renderRHandsontable(makeMSConversionHot())
        
        observeEvent(input$featGrouper, {
            if (input$ionization == "both" && input$featGrouper == "SIRIUS")
            {
                rstudioapi::showDialog("Not supported", "Grouping features with SIRIUS is currently not supported with sets", "")
                updateSelectInput(inputId = "featGrouper", selected = "OpenMS")
            }
            else
                shinyjs::toggleState("featFinder", input$featGrouper != "SIRIUS")
        })
        
        observeEvent(input$TPGen, {
            updateSelectInput(inputId = "TPGenInput", choices = getTPGenInputs(input$TPGen == "Library"))
        })
        
        observeEvent(input$TPGenInput, {
            if (input$TPGenInput == "screening" &&
                (!input$exSuspList && (input$ionization == "both" || !nzchar(input$suspectList)) &&
                 (input$ionization != "both" || !nzchar(input$suspectListPos))))
            {
                rstudioapi::showDialog("Enable suspect screening",
                                       "This requires a workflow with suspect screening. Please configure in the Features tab.", "")
                updateSelectInput(inputId = "TPGenInput", selected = "suspects")
            }
        })
        
        observeEvent(input$TPSuspButton, selectSuspList("TPSuspectList"))
        
        observeEvent(input$MSLibraryPathButton, {
            MSFile <- rstudioapi::selectFile(path = "~/")
            if (!is.null(MSFile))
                updateTextInput(session, "MSLibraryPath", value = MSFile)
        })
    }

    runGadget(getNewProjectUI(destPath), server, viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
