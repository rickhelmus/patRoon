# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getScriptCode <- function(input, analyses)
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
    addAnaInfo <- function(anaInfoVarName, anaTable, anaTableFile, comment, ionization)
    {
        if (input$generateAnaInfo == "table")
        {
            addComment("Load analysis table", condition = comment)
            addCall(anaInfoVarName, "read.csv", list(value = anaTableFile, quote = TRUE))
        }
        else if (input$generateAnaInfo == "script")
        {
            addCall(anaInfoVarName, "generateAnalysisInfo", list(
                list(name = "paths", value = unique(anaTable$path), quote = TRUE),
                list(name = "groups", value = anaTable$group, quote = TRUE),
                list(name = "blanks", value = anaTable$blank, quote = TRUE),
                list(name = "concs", value = anaTable$conc),
                list(name = "norm_concs", value = anaTable$norm_conc)
            ))
            if (any(anaTable$exclude))
            {
                toExclude <- paste0('c(', paste0(which(anaTable$exclude), collapse = ', '), ')')
                addAssignment(anaInfoVarName, paste0(anaInfoVarName, '[-', toExclude, ', ]'))
                addNL()
            }
        }
        else if (input$generateAnaInfo == "example")
        {
            addComment("Example data from patRoonData package (triplicate solvent blank + triplicate standard)",
                       condition = comment)
            addCall(anaInfoVarName, "patRoonData::exampleAnalysisInfo", list(value = ionization, quote = TRUE))
        }
        else # none
        {
            addComment("NOTE: please set anaInfo to a valid data.frame with analysis information. See ?`analysis-information` for more details.",
                       condition = comment)
            addCall(anaInfoVarName, "data.frame", list(
                list(name = "path", value = "character()"),
                list(name = "analysis", value = "character()"),
                list(name = "groups", value = "character()"),
                list(name = "blanks", value = "character()")
            ))
        }
    }
    addPrepBlock <- function(anaInfoVarName, DAMethodVarName)
    {
        addCall(NULL, "setDAMethod", list(
            list(value = anaInfoVarName),
            list(value = input[[DAMethodVarName]], quote = TRUE)
        ), condition = nzchar(input[[DAMethodVarName]]), indent = 4)
        addCall(NULL, "recalibrarateDAFiles", list(value = anaInfoVarName), condition = input$doDACalib, indent = 4)
        if (nzchar(input$convAlgo))
        {
            if (input$peakPicking)
                centroid <- if (input$peakPickingVendor && input$convAlgo == "pwiz") "\"vendor\"" else TRUE
            else
                centroid <- FALSE
            
            for (of in input$convTo)
            {
                addCall(NULL, "convertMSFiles", list(
                    list(name = "anaInfo", value = anaInfoVarName),
                    list(name = "from", value = input$convFrom, quote = TRUE),
                    list(name = "to", value = of, quote = TRUE),
                    list(name = "algorithm", value = input$convAlgo, quote = TRUE),
                    list(name = "centroid", value = centroid)),
                    indent = 4)
            }
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
        addAnaInfo("anaInfo", analyses, input$analysisTableFile, TRUE, input$ionization)
    else
    {
        addAnaInfo("anaInfoPos", analyses$pos, input$analysisTableFilePos, TRUE, "positive")
        addAnaInfo("anaInfoNeg", analyses$neg, input$analysisTableFileNeg, FALSE, "negative")
    }
    
    if (nzchar(input$convAlgo) || nzchar(input$DAMethod) || input$doDACalib)
    {
        addNL()
        addComment("Set to FALSE to skip data pre-treatment")
        addAssignment("doDataPretreatment", TRUE)
        addText("if (doDataPretreatment)\n{")
        if (input$ionization != "both")
            addPrepBlock("anaInfo", "DAMethod")
        else
        {
            addPrepBlock("anaInfoPos", "DAMethodPos")
            addNL()
            addPrepBlock("anaInfoNeg", "DAMethodNeg")
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
        list(name = "groupParam", value = "xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$group)",
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
            list(name = "rtDev", value = 30, condition = input$components == "nontarget"),
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
                list(name = "precursorMzSearchWindow", value = 0.002, condition = input$formulaGen == "Bruker"),
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
            addCall(NULL, "reportHTML", condition = "HTML" %in% input$reportLegacy, list(
                list(value = "fGroups"),
                list(name = "path", value = "report", quote = TRUE),
                list(name = "formulas", value = "formulas", isNULL = !nzchar(input$formulaGen)),
                list(name = "compounds", value = "compounds", isNULL = !nzchar(input$compIdent)),
                list(name = "MSPeakLists", value = "mslists", condition = doMSPL),
                list(name = "components", value = componVal),
                list(name = "TPs", value = "TPs", condition = input$doTPs),
                list(name = "reportPlots", value = c("chord", "venn", "upset", "eics", "formulas"), quote = TRUE),
                list(name = "selfContained", value = FALSE),
                list(name = "openReport", value = TRUE)
            ))
        }
    }
    
    addNL()

    return(paste0(textConnectionValue(txtCon), collapse = "\n"))
}

doCreateProject <- function(input, analyses)
{
    mkdirp(input$destinationPath)

    prepareAnas <- function(anas, tableFile)
    {
        anas <- copy(anas)
        anas[, group := ifelse(!nzchar(group), analysis, group)]
        
        # Make analysis table
        if (input$generateAnaInfo == "table")
            write.csv(anas[!anas$exclude, c("path", "analysis", "group", "blank", "conc", "norm_conc")],
                      file.path(input$destinationPath, tableFile), row.names = FALSE)
        
        return(anas)
    }

    if (input$ionization != "both")
        analyses <- prepareAnas(analyses, input$analysisTableFile)
    else
        analyses <- Map(prepareAnas, analyses, list(input$analysisTableFilePos, input$analysisTableFileNeg))
    
    doSusps <- input$exSuspList || (input$ionization != "both" && nzchar(input$suspectList)) ||
        (input$ionization == "both" && nzchar(input$suspectListPos))
    if (doSusps && input$annotateSus && input$genIDLevelFile)
        genIDLevelRulesFile(file.path(input$destinationPath, "idlevelrules.yml"))
    
    if ("HTML" %in% input$reportGen)
        genReportSettingsFile(file.path(input$destinationPath, "report.yml"))
    
    code <- getScriptCode(input, analyses)
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

getNewProjectUI <- function(destPath)
{
    textNote <- function(txt) div(style = "margin: 8px 0 12px; font-size: small", txt)
    fileSelect <- function(idText, idButton, label, value = "", ...)
    {
        fillRow(
            flex = c(1, NA),
            textInput(idText, label, value, width = "100%", ...),
            actionButton(idButton, "", icon("folder-open"), style = "margin: 25px 0 0 15px")
        )
    }
    rangeNumeric <- function(id, label, minVal = 0, maxVal = 0, ...)
    {
        fillRow(
            numericInput(paste0(id, "-min"), paste("Min.", label), value = minVal, ..., width = "95%"),
            numericInput(paste0(id, "-max"), paste("Max.", label), value = maxVal, ..., width = "100%")
        )
    }

    miniUI::miniPage(
        shinyjs::useShinyjs(),
        
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
                        flex = c(NA, NA, 1, NA),
                        fillRow(
                            height = 120,
                            radioButtons("generateAnaInfo", "Generate analysis information",
                                         c("None" = "none", "From new csv file" = "table",
                                           "Load in script" = "script", "Example data" = "example")),
                            conditionalPanel(
                                condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                                fillCol(
                                    flex = NA,
                                    conditionalPanel(
                                        condition = "input.ionization == \"both\"",
                                        radioButtons("currentSet", "Selected set", c("positive", "negative"),
                                                     inline = TRUE)
                                    ),
                                    conditionalPanel(
                                        condition = "input.ionization != \"both\" && input.generateAnaInfo == \"table\"",
                                        textInput("analysisTableFile", "Analysis table output file", "analyses.csv")
                                    ),
                                    conditionalPanel(
                                        condition = "input.ionization == \"both\" && input.generateAnaInfo == \"table\" && input.currentSet == \"positive\"",
                                        textInput("analysisTableFilePos",
                                                  "Analysis table output file", "analyses-pos.csv")
                                    ),
                                    conditionalPanel(
                                        condition = "input.ionization == \"both\" && input.generateAnaInfo == \"table\" && input.currentSet == \"negative\"",
                                        textInput("analysisTableFileNeg",
                                                  "Analysis table output file", "analyses-neg.csv")
                                    )
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                            fillRow(
                                height = 30,
                                textNote("Make sure to consider data conversion if data files are not yet in mzXML/mzML format.")
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"example\"",
                            fillRow(
                                height = 30,
                                textNote("Make sure that the patRoonData package is installed.")
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                            flex = NA,
                            conditionalPanel(
                                condition = "input.ionization != \"both\"",
                                rhandsontable::rHandsontableOutput("analysesHot")
                            ),
                            conditionalPanel(
                                condition = "input.ionization == \"both\" && input.currentSet == \"positive\"",
                                rhandsontable::rHandsontableOutput("analysesHotPos")
                            ),
                            conditionalPanel(
                                condition = "input.ionization == \"both\" && input.currentSet == \"negative\"",
                                rhandsontable::rHandsontableOutput("analysesHotNeg")
                            )
                        )
                    )
                ),
                
                conditionalPanel(
                    condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                    miniUI::miniButtonBlock(
                        actionButton("addAnalysesDir", "Add analyses from directory ..."),
                        actionButton("addAnalysesCSV", "Add analyses from csv file ...")
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Data pre-treatment", icon = icon("upload"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        height = 220,
                        fillRow(
                            height = 75,
                            selectInput("convAlgo", "Data Conversion Algorithm", c("None" = "", "ProteoWizard" = "pwiz",
                                                                                   "Bruker DataAnalysis" = "bruker",
                                                                                   "OpenMS" = "openms"),
                                        width = "100%")
                        ),
                        conditionalPanel(
                            "input.convAlgo != \"\"",
                            fillRow(
                                height = 90,
                                selectInput("convFrom", "Input format", MSFileFormats(), multiple = FALSE,
                                            width = "95%"),
                                fillCol(
                                    flex = c(1, NA),
                                    selectInput("convTo", "Output format", "mzML", multiple = FALSE,
                                                selected = "mzML", width = "100%"),
                                    textNote("enviPick/XCMS support mzXML, XCMS/OpenMS support mzML")
                                )
                            )
                        ),
                        conditionalPanel(
                            "input.convAlgo == \"pwiz\" || input.convAlgo == \"bruker\"",
                            fillRow(
                                height = 60,
                                checkboxInput("peakPicking", "Perform peak picking (line spectra)", value = TRUE),
                                conditionalPanel(
                                    condition = "input.convAlgo == \"pwiz\"",
                                    checkboxInput("peakPickingVendor", "Use vendor algorithm for peak picking", value = TRUE)
                                )
                            )
                        )
                    ),
                    hr(),
                    fillCol(
                        flex = NA,
                        height = 70,
                        strong("Bruker DataAnalysis options"),
                        textNote("Only supported with bruker data and if DataAnalysis is installed.")
                    ),
                    fillCol(
                        flex = NA,
                        height = 120,
                        fillCol(
                            flex = NA,
                            height = 50,
                            conditionalPanel(
                                condition = "input.ionization != \"both\"",
                                fileSelect("DAMethod", "DAMethodButton", "DataAnalysis method"),
                            ),
                            conditionalPanel(
                                condition = "input.ionization == \"both\"",
                                fillRow(
                                    fillCol(
                                        width = "95%",
                                        fileSelect("DAMethodPos", "DAMethodButtonPos", "DataAnalysis method (positive)")
                                    ),
                                    fillCol(
                                        width = "95%",
                                        fileSelect("DAMethodNeg", "DAMethodButtonNeg", "DataAnalysis method (negative)")
                                    )
                                )
                            )
                        ),
                        textNote("Leaving this blank will not set any method"),
                        checkboxInput("doDACalib", "Perform m/z re-calibration")
                    )
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
                            checkboxGroupInput("reportLegacy", "Legacy report formats", c("CSV", "PDF", "HTML"),
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
        DAMethod = "text",
        DAMethodPos = "text",
        DAMethodNeg = "text",
        doDACalib = "check",
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
                    columnSorting = FALSE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    selectionMode = "range", outsideClickDeselects = FALSE,
                    contextMenu = FALSE, manualColumnResize = TRUE)

    emptyAnaTable <- function() data.table(exclude = logical(0), analysis = character(0), format = character(0),
                                           group = character(0), blank = character(0), conc = numeric(0), 
                                           norm_conc = numeric(0), path = character(0))
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(analyses = emptyAnaTable(),
                                  analysesPos = emptyAnaTable(),
                                  analysesNeg = emptyAnaTable())

        makeAnalysesHot <- function(rvName)
        {
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(rValues[[rvName]], height = 250, maxRows = nrow(rValues[[rvName]])),
                             hotOpts)) %>%
                rhandsontable::hot_col(c("group", "blank"), readOnly = FALSE, type = "text") %>%
                rhandsontable::hot_col(c("conc", "norm_conc"), readOnly = FALSE, type = "numeric") %>%
                rhandsontable::hot_col("exclude", readOnly = FALSE, type = "checkbox")
            
            return(hot)
        }
        getCurAnaHotName <- function()
        {
            if (input$ionization != "both")
                return("analysesHot")
            else if (input$currentSet == "positive")
                return("analysesHotPos")
            return("analysesHotNeg")
        }
        getCurAnaRVName <- function()
        {
            if (input$ionization != "both")
                return("analyses")
            else if (input$currentSet == "positive")
                return("analysesPos")
            return("analysesNeg")
        }
        verifyAnalysesOK <- function()
        {
            noAnas <- input$ionization != "both" && nrow(rValues$analyses) == 0
            noAnasPos <- input$ionization == "both" && nrow(rValues$analysesPos) == 0
            noAnasNeg <- input$ionization == "both" && nrow(rValues$analysesNeg) == 0
            
            if (noAnas || noAnasPos || noAnasNeg)
            {
                msg <- "Please select some analyses"
                if (noAnasPos && !noAnasNeg)
                    msg <- paste(msg, "for positive mode")
                else if (!noAnasPos && noAnasNeg)
                    msg <- paste(msg, "for negative mode")
                rstudioapi::showDialog("No analyses selected", paste0(msg, "!"), "")
                return(FALSE)
            }
            
            if (input$generateAnaInfo == "table")
            {
                checkAnas <- if (input$ionization != "both")
                    input$analysisTableFile
                else
                    c(input$analysisTableFilePos, input$analysisTableFileNeg)
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
            else if (input$generateAnaInfo %in% c("table", "script") && !verifyAnalysesOK())
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
                anas <- if (input$ionization != "both") rValues$analyses else list(pos = rValues$analysesPos,
                                                                                   neg = rValues$analysesNeg)
                doCreateProject(input, anas)
                stopApp(TRUE)
            }
        })

        observeEvent(input$projectDestButton, {
            dest <- rstudioapi::selectDirectory("Select destination directory", path = input$destinationPath)
            if (!is.null(dest))
                updateTextInput(session, "destinationPath", value = dest)
        })
        
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
        
        doObserveAnaHot <- function(name, rvName)
        {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (!is.null(input[[name]]) && input[[name]]$params$maxRows > 0)
            {
                df <- rhandsontable::hot_to_r(input[[name]])
                rValues[[rvName]][, c("exclude", "group", "blank", "conc", "norm_conc") := .(df$exclude, df$group, df$blank, df$conc, df$norm_conc)]
            }
        }
        observeEvent(input$analysesHot, doObserveAnaHot("analysesHot", "analyses"))
        observeEvent(input$analysesHotPos, doObserveAnaHot("analysesHotPos", "analysesPos"))
        observeEvent(input$analysesHotNeg, doObserveAnaHot("analysesHotNeg", "analysesNeg"))
        
        observeEvent(input$addAnalysesDir, {
            anaDir <- rstudioapi::selectDirectory(path = "~/")
            if (!is.null(anaDir))
            {
                files <- listMSFiles(anaDir, MSFileFormats())

                if (length(files) > 0)
                {
                    dt <- data.table(exclude = FALSE, path = dirname(files), analysis = simplifyAnalysisNames(files),
                                     group = "", blank = "", conc = NA_real_, norm_conc = NA_real_)

                    msExts <- MSFileExtensions()
                    dt[, format := sapply(tolower(tools::file_ext(files)), function(ext)
                    {
                        paste0(names(msExts)[sapply(msExts, function(e) ext %in% tolower(e))], collapse = "/")
                    })]

                    dt[, format := paste0(.SD$format, collapse = ", "), by = .(path, analysis)]
                    dt <- unique(dt, by = c("analysis", "path"))
         
                    setcolorder(dt, c("exclude", "analysis", "format", "group", "blank", "conc", "norm_conc", "path"))

                    rvName <- getCurAnaRVName()
                    rValues[[rvName]] <- rbind(rValues[[rvName]], dt)
                }
            }
        })

        observeEvent(input$addAnalysesCSV, {
            csvFile <- rstudioapi::selectFile(path = "~/", filter = "csv files (*.csv)")
            if (!is.null(csvFile))
            {
                csvTab <- tryCatch(fread(csvFile, select = c("path", "analysis", "group", "blank"),
                                         colClasses = "character"),
                                   error = function(e) FALSE, warning = function(w) FALSE)
                if (is.logical(csvTab))
                    rstudioapi::showDialog("Error", "Failed to open/parse selected csv file!", "")
                else if (nrow(csvTab) > 0)
                {
                    msExts <- MSFileExtensions()
                    msFiles <- normalizePath(listMSFiles(csvTab$path, MSFileFormats()), winslash = "/")
                    msFilesNoExt <- tools::file_path_sans_ext(msFiles)
                    formats <- mapply(csvTab$analysis, csvTab$path, FUN = function(ana, path)
                    {
                        fps <- msFiles[msFilesNoExt == file.path(path, ana)]
                        if (length(fps) == 0)
                            return("")
                        ret <- sapply(tolower(tools::file_ext(fps)), function(ext)
                        {
                            paste0(names(msExts)[sapply(msExts, function(e) ext %in% tolower(e))], collapse = "/")
                        })
                        return(paste0(ret, collapse = ", "))
                    })
                    
                    csvTab[, format := formats]
                    csvTab <- csvTab[nzchar(format)] # prune unknown files (might have been removed?)
                    
                    rvName <- getCurAnaRVName()
                    rValues[[rvName]] <- rbind(rValues[[rvName]], csvTab, fill = TRUE)
                }
            }
        })
        
        observeEvent(input$convAlgo, {
            from <- switch(input$convAlgo,
                           pwiz = MSFileFormats("pwiz"),
                           bruker = "bruker",
                           openms = MSFileFormats("openms"),
                           ""
                   )
            sel <- ""
            if (nzchar(input$convAlgo))
                sel <- MSFileFormats(input$convAlgo, input$convAlgo != "openms")[1]

            updateSelectInput(session, "convFrom", choices = from, selected = sel)
        })

        observeEvent(input$DAMethodButton, selectDAMethod("DAMethod"))
        observeEvent(input$DAMethodButtonPos, selectDAMethod("DAMethodPos"))
        observeEvent(input$DAMethodButtonNeg, selectDAMethod("DAMethodNeg"))

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
        
        output$analysesHot <- rhandsontable::renderRHandsontable(makeAnalysesHot("analyses"))
        output$analysesHotPos <- rhandsontable::renderRHandsontable(makeAnalysesHot("analysesPos"))
        output$analysesHotNeg <- rhandsontable::renderRHandsontable(makeAnalysesHot("analysesNeg"))
        
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
