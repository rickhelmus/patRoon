# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

local_edition(3) # for snapshots

testBaseDir <- getWorkPath("test-np")
defaultTestDir <- file.path(testBaseDir, "default")
suspectsPaths <- list(
    pos = getWorkPath("suspects_pos.csv"),
    neg = getWorkPath("suspects_neg.csv"),
    istd_pos = getWorkPath("istds_pos.csv"),
    istd_neg = getWorkPath("istds_neg.csv")
)

defaultSettings <- list(
    general = defaultGeneralSettings(defaultTestDir),
    analyses = defaultAnalysesSettings(),
    preTreatment = defaultPreTreatSettings(),
    features = defaultFeaturesSettings(),
    annotations = defaultAnnotationSettings(),
    TPs = defaultTPSettings(),
    report = defaultReportSettings()
)

# HACK: make sure these following Shiny update functions correctly set their data with testServer()
doShinyUpdate <- function(session, inputId, value)
{
    # printf("Updating %s to %s\n", inputId, value)
    if (!is.null(value)) # NULL if something else than the value was updated
        do.call(session$setInputs, setNames(list(value), inputId))
}
local_mocked_bindings(
    updateCheckboxGroupInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...)
    {
        doShinyUpdate(session, inputId, selected)
    },
    updateCheckboxInput = function(session, inputId, label = NULL, value = NULL)
    {
        doShinyUpdate(session, inputId, value)
    },
    updateSelectInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL)
    {
        doShinyUpdate(session, inputId, selected)
    },
    updateRadioButtons = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...)
    {
        doShinyUpdate(session, inputId, selected)
    },
    updateTextInput = function(session, inputId, label = NULL, value = NULL, ...)
    {
        doShinyUpdate(session, inputId, value)
    }
)

settingsFile <- getWorkPath("project_settings.yml")
local_mocked_bindings(
    selectFile = function(...) settingsFile,
    initializeProject = function(...) NULL,
    openProject = function(...) NULL,
    .package = "rstudioapi"
)

makeNewProj <- function(settings, CCSCalib = "", aid = list())
{
    doCreateProject(CCSCalib, aid, settings, noDate = TRUE)
    return(readAllFile(file.path(settings$general$destination, settings$general$scriptFile)))
}

defaultCode <- makeNewProj(defaultSettings)

modifyDefSettings <- function(general = NULL, analyses = NULL, preTreatment = NULL,
                              features = NULL, annotations = NULL, TPs = NULL, report = NULL)
{
    settings <- defaultSettings
    if (!is.null(general))
        settings$general <- modifyList(settings$general, general)
    if (!is.null(analyses))
        settings$analyses <- modifyList(settings$analyses, analyses)
    if (!is.null(preTreatment))
    {
        # HACK: don't try to change steps with modifyList(), gives errors
        ptNoSteps <- preTreatment[names(preTreatment) != "steps"]
        settings$preTreatment <- modifyList(settings$preTreatment, ptNoSteps)
        if (!is.null(preTreatment[["steps"]]))
            settings$preTreatment$steps <- copy(preTreatment$steps)
    }
    if (!is.null(features))
        settings$features <- modifyList(settings$features, features)
    if (!is.null(annotations))
        settings$annotations <- modifyList(settings$annotations, annotations)
    if (!is.null(TPs))
        settings$TPs <- modifyList(settings$TPs, TPs)
    if (!is.null(report))
        settings$report <- modifyList(settings$report, report)
    
    return(settings)
}

testNewProj <- function(..., name, CCSCalib = "", aid = list())
{
    settings <- modifyDefSettings(...)
    path <- file.path(testBaseDir, name)
    diffp <- file.path(path, "process.diff")
    
    announce_snapshot_file(path, name = name)
    
    # verify that the settings are correctly stored/loaded
    testServer(newProjectServer(path), {
        writeProjectSettings(settings, settingsFile)
        session$setInputs("general-loadParams" = 1)
        session$flushReact()
        
        unlink(settingsFile)
        session$setInputs("general-saveParams" = 1)
        session$flushReact()
        expect_equal(readProjectSettings(settingsFile, path), settings)
    })
    
    unlink(path, recursive = TRUE)
    settings$general$destination <- path
    code <- makeNewProj(settings, CCSCalib, aid)
    
    # diffobj package is used to create diffs so we don't need snapshot whole scripts for each test
    cat(as.character(diffobj::diffFile(file.path(defaultTestDir, "process.R"), file.path(path, settings$general$scriptFile),
                                       pager = "off", format = "raw", mode = "unified", rds = FALSE, disp.width = 200)),
                     file = diffp, sep = "\n")
    expect_snapshot_file(diffp, name = name, cran = TRUE)
    
    if (settings$analyses$generateAnaInfo == "table")
    {
        checkCols <- c("analysis", "replicate", "blank", "path_centroid")
        if (settings$analyses$analysisTableFileType == "CSV")
        {
            checkAnaF <- function(f, t)
            {
                expect_equal(fread(file.path(path, f), select = checkCols), t[, checkCols, with = FALSE])
            }
            if (settings$general$ionization == "positive")
                checkAnaF(settings$analyses$analysisTableFileCSV, aid$positive)
            else if (settings$general$ionization == "negative")
                checkAnaF(settings$analyses$analysisTableFileCSV, aid$negative)
            else if (settings$general$ionization == "both")
            {
                checkAnaF(settings$analyses$analysisTableFileCSVPos, aid$positive)
                checkAnaF(settings$analyses$analysisTableFileCSVNeg, aid$negative)
            }
        }
        else if (settings$analyses$analysisTableFileType == "R")
        {
            checkAnaF <- function(f, t)
            {
                expect_equal(eval(parse(file.path(path, f)))[, checkCols], as.data.frame(t)[, checkCols])
            }
            if (settings$general$ionization == "positive")
                checkAnaF(settings$analyses$analysisTableFileR, aid$positive)
            else if (settings$general$ionization == "negative")
                checkAnaF(settings$analyses$analysisTableFileR, aid$negative)
            else if (settings$general$ionization == "both")
            {
                checkAnaF(settings$analyses$analysisTableFileRPos, aid$positive)
                checkAnaF(settings$analyses$analysisTableFileRNeg, aid$negative)
            }
        }
    }
}

test_that("Default settings", {
    expect_snapshot_file(file.path(defaultTestDir, "process.R"), name = "default_process.R", cran = TRUE)
    expect_snapshot_file(file.path(defaultTestDir, "report.yml"), name = "default_report.yml", cran = TRUE)
    expect_snapshot_file(file.path(defaultTestDir, "limits.yml"), name = "default_limits.yml", cran = TRUE)
})

test_that("General settings", {
    # NOTE: ionization pos/neg will not change the default script --> test it in later cases where it does
    testNewProj(general = list(ionization = "both"), name = "general-sets")
    testNewProj(general = list(IMS = list(mode = "direct")), features = list(featAlgo = "EIC"),
                name = "general-ims_direct")
    testNewProj(general = list(IMS = list(mode = "post")), name = "general-ims_post")
    
    testNewProj(general = list(IMS = list(mode = "post", CCSMethod = "agilent", limits = "agilent")),
                name = "general-ims_post_agilent")
    expect_equal(readYAML(file.path(testBaseDir, "general-ims_post_agilent", "limits.yml"))$general$IMS, "agilent")
})

test_that("Analysis settings", {
    aid <- list(positive = makeDT(patRoonData::exampleAnalysisInfo("positive")),
                negative = makeDT(patRoonData::exampleAnalysisInfo("negative")))
    aid$positive$type <- aid$negative$type <- "centroid"
    
    testNewProj(analyses = list(generateAnaInfo = "example"), name = "analysis-example")
    
    testNewProj(analyses = list(generateAnaInfo = "table", analysisTableFileType = "CSV"), aid = aid,
                name = "analysis-tab_csv_pos")
    testNewProj(general = list(ionization = "negative"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "CSV"), aid = aid,
                name = "analysis-tab_csv_neg")
    testNewProj(general = list(ionization = "both"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "CSV"), aid = aid,
                name = "analysis-tab_csv_both")
    
    testNewProj(analyses = list(generateAnaInfo = "table", analysisTableFileType = "R"), aid = aid,
                name = "analysis-tab_R_pos")
    testNewProj(general = list(ionization = "negative"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "R"), aid = aid,
                name = "analysis-tab_R_neg")
    testNewProj(general = list(ionization = "both"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "R"), aid = aid,
                name = "analysis-tab_R_both")
    
    testNewProj(analyses = list(generateAnaInfo = "table", analysisTableFileType = "embedded"), aid = aid,
                name = "analysis-tab_emb_pos")
    testNewProj(general = list(ionization = "negative"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "embedded"), aid = aid,
                name = "analysis-tab_emb_neg")
    testNewProj(general = list(ionization = "both"),
                analyses = list(generateAnaInfo = "table", analysisTableFileType = "embedded"), aid = aid,
                name = "analysis-tab_emb_both")
    
    testNewProj(analyses = list(generateAnaInfo = "dynamic", genAnaInfoDynRaw = "raw"), name = "analysis-dyn_raw_only")
    testNewProj(analyses = list(generateAnaInfo = "dynamic", genAnaInfoDynRaw = "raw",
                                genAnaInfoDynCentroid = "centroid", genAnaInfoDynIMS = "ims",
                                genAnaInfoDynProfile = "profile"),
                name = "analysis-dyn")
    testNewProj(general = list(ionization = "both"),
                analyses = list(generateAnaInfo = "dynamic",
                                genAnaInfoDynRawPos = "rawpos",
                                genAnaInfoDynRawNeg = "rawneg",
                                genAnaInfoDynCentroidPos = "centroidpos",
                                genAnaInfoDynCentroidNeg = "centroidneg",
                                genAnaInfoDynIMSPos = "imspos",
                                genAnaInfoDynIMSNeg = "imsneg",
                                genAnaInfoDynProfilePos = "profilepos",
                                genAnaInfoDynProfileNeg = "profileneg"),
                name = "analysis-dyn_sets")
})

test_that("Pre-treatment settings", {
    testNewProj(preTreatment = list(
        steps = rbindlist(list(
            # just try to cover most algos with some formats
            data.table(algorithm = "pwiz", from = joinConvTypeFormat("raw", "thermo"),
                       to = joinConvTypeFormat("centroid", "mzML")),
            data.table(algorithm = "openms", from = joinConvTypeFormat("centroid", "mzML"),
                       to = joinConvTypeFormat("centroid", "mzXML")),
            data.table(algorithm = "bruker", from = joinConvTypeFormat("raw", "bruker"),
                       to = joinConvTypeFormat("profile", "mzML")),
            data.table(algorithm = "im_collapse", from = joinConvTypeFormat("ims", "mzML"),
                       to = joinConvTypeFormat("centroid", "mzML")),
            data.table(algorithm = "timsconvert", from = joinConvTypeFormat("raw", "bruker_ims"),
                       to = joinConvTypeFormat("profile", "mzML"))
        ))),
        name = "pretreatment-steps"
    )
    testNewProj(preTreatment = list(brukerCalib = list(enabled = TRUE, method = "method")),
                name = "pretreatment-bruker_calib")
    testNewProj(general = list(ionization = "both"),
                preTreatment = list(brukerCalib = list(enabled = TRUE, methodPos = "method-pos",
                                                       methodNeg = "method-neg")),
                name = "pretreatment-bruker_calib_sets")
})

test_that("Feature settings", {
    testNewProj(features = list(featAlgo = "XCMS"), name = "features-xcms")
    testNewProj(features = list(featAlgo = "enviPick"), name = "features-envipick")
    testNewProj(features = list(featAlgo = "KPIC2"), name = "features-kpic2")
    testNewProj(features = list(featAlgo = "SIRIUS"), name = "features-sirius")
    testNewProj(features = list(featAlgo = "Bruker"), name = "features-bruker")
    
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "bins")), name = "features-eic_bins")
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "bins", peaksAlgo = "openms")),
                name = "features-eic_bins_openms")
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "suspects",
                                                                       suspects = list(single = "suspects"))),
                name = "features-eic_suspects")
    testNewProj(general = list(ionization = "both"),
                features = list(featAlgo = "EIC",
                                featEICParams = list(methodMZ = "suspects",
                                                     suspects = list(sets = list(pos = "suspects-pos",
                                                                                 neg = "suspects-neg")))),
                name = "features-eic_suspects_sets")
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "ms2")), name = "features-eic_ms2")
    
    testNewProj(features = list(fGroupsAlgo = "XCMS"), name = "features-fgroups-xcms")
    testNewProj(features = list(fGroupsAlgo = "KPIC2"), name = "features-fgroups-kpic2")
    testNewProj(features = list(fGroupsAlgo = "SIRIUS"), name = "features-fgroups-sirius")
    
    testNewProj(general = list(IMS = list(mode = "post")),
                features = list(IMSPeaksMob = "xcms", IMSPeaksChrom = "envipick"), name = "features-ims_peak")
    
    testNewProj(features = list(suspects = list(single = "suspects")), name = "features-suspects")
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "suspects",
                                                                       suspects = list(single = "")),
                                suspects = list(single = "suspects")), name = "features-suspects_eic")
    testNewProj(features = list(featAlgo = "EIC", featEICParams = list(methodMZ = "suspects",
                                                                       suspects = list(single = "suspectsEIC")),
                                suspects = list(single = "suspects")), name = "features-suspects_eic_suspects")
    testNewProj(general = list(ionization = "both"),
                features = list(suspects = list(sets = list(pos = "suspectsPos", neg = "suspectsNeg"))),
                name = "features-suspects_sets")
    testNewProj(general = list(ionization = "both"),
                features = list(featAlgo = "EIC",
                                featEICParams = list(methodMZ = "suspects",
                                                     suspects = list(sets = list(pos = "suspectsEICPos", 
                                                                                 neg = "suspectsEICNeg"))),
                                suspects = list(sets = list(pos = "suspectsPos", neg = "suspectsNeg"))),
                name = "features-suspects_sets_eic_suspects")
    testNewProj(features = list(exSuspList = TRUE), name = "features-suspects_ex")
    testNewProj(general = list(ionization = "both"), features = list(exSuspList = TRUE),
                name = "features-suspects_sets_ex")
    
    testNewProj(general = list(IMS = list(mode = "single")),
                features = list(exSuspList = TRUE, IMSSuspCCSPred = "c3sdb"),
                name = "features-ims_susp_pred")
    testNewProj(general = list(ionization = "both", IMS = list(mode = "single")),
                features = list(exSuspList = TRUE, IMSSuspCCSPred = "c3sdb"),
                name = "features-ims_susp_pred_sets")
    
    testNewProj(features = list(fGroupsAdv = list(preIntThr = 3E5, intThr = 3E6, repAbundance = 0.5,
                                                  maxRepRSD = 0.25, blankThr = 10, removeBlanks = FALSE,
                                                  retention = c(60, 10000), mz = c(10, 1000))),
                name = "features-filters")
    
    testNewProj(features = list(fGroupsAdv = list(featNorm = "istd", ISTDLists = list(single = "istd"))),
                name = "features-norm_istd")
    testNewProj(general = list(ionization = "both"),
                features = list(fGroupsAdv = list(featNorm = "istd",
                                                  ISTDLists = list(sets = list(pos = "istd-pos", neg = "istd-neg")))),
                name = "features-norm_istd_sets")
    testNewProj(features = list(exSuspList = TRUE, fGroupsAdv = list(featNorm = "istd")),
                name = "features-norm_istd_ex")
    testNewProj(features = list(fGroupsAdv = list(featNorm = "tic", groupNorm = TRUE)),
                name = "features-norm_tic_group")
})

test_that("annotation settings", {
    testNewProj(annotations = list(componAlgo = "RAMClustR"), name = "annotations-compon_rc")
    testNewProj(annotations = list(componAlgo = "CAMERA", selectIons = FALSE),
                name = "annotations-compon_camera_nosel")
    testNewProj(annotations = list(componAlgo = "nontarget"), name = "annotations-compon_nt")
    
    testNewProj(annotations = list(formulasAlgo = "GenForm"), name = "annotations-formulas_genform")
    testNewProj(annotations = list(formulasAlgo = "SIRIUS", estIDConf = character()),
                name = "annotations-formulas_sirius_noann")
    
    testNewProj(annotations = list(compoundsAlgo = "MetFrag"), name = "annotations-compounds_metfrag")
    testNewProj(annotations = list(compoundsAlgo = "SIRIUS", estIDConf = character()),
                name = "annotations-compounds_sirius_noann")
    testNewProj(annotations = list(compoundsAlgo = "Library", MSLibraryPath = "lib", MSLibraryFormat = "json"),
                name = "annotations-compounds_library")
    
    testNewProj(general = list(IMS = list(mode = "post", CCSMethod = "bruker")),
                annotations = list(compoundsAlgo = "MetFrag", compCCSPred = "c3sdb"),
                name = "annotations-compounds_ims_conv_pred")
    
    testNewProj(features = list(exSuspList = TRUE),
                annotations = list(formulasAlgo = "GenForm", compoundsAlgo = "MetFrag"),
                name = "annotations-susp_ann")
})

test_that("TP settings", {
    testNewProj(TPs = list(TPsAlgo = "logic"), name = "tp-logic")
    testNewProj(TPs = list(TPsAlgo = "biotransformer", TPGenInput = "suspects", TPSuspectList = "suspects"),
                name = "tp-biotransformer_suspects")
    testNewProj(annotations = list(compoundsAlgo = "MetFrag"),
                TPs = list(TPsAlgo = "biotransformer", TPGenInput = "suspects", TPSuspectList = "suspects",
                           TPDoMFDB = TRUE), name = "tp-biotransformer_mfdb")
    testNewProj(annotations = list(compoundsAlgo = "MetFrag"),
                TPs = list(TPsAlgo = "biotransformer", TPGenInput = "suspects", TPSuspectList = "suspects",
                           TPDoMFDB = FALSE), name = "tp-biotransformer_no_mfdb")
    testNewProj(features = list(exSuspList = TRUE),
                TPs = list(TPsAlgo = "cts", TPGenInput = "screening"), name = "tp-cts_scr")
    testNewProj(TPs = list(TPsAlgo = "library", TPGenInput = "all"), name = "tp-library")
    testNewProj(annotations = list(formulasAlgo = "GenForm"),
                TPs = list(TPsAlgo = "ann_form", TPGenInput = "suspects", TPSuspectList = "suspects"),
                name = "tp-ann_form")
    testNewProj(annotations = list(compoundsAlgo = "MetFrag"),
                TPs = list(TPsAlgo = "ann_comp", TPGenInput = "suspects", TPSuspectList = "suspects"),
                name = "tp-ann_comp")
})

test_that("Report settings", {
    testNewProj(report = list(reportGen = c("HTML", "legacy"), reportLegacy = c("CSV", "PDF")), 
                name = "report-html_legacy_csv_pdf")
})

test_that("Old setting conversion", {
    expect_snapshot(readProjectSettings(file.path(getTestDataPath(), "newProject-23.yml"), defaultTestDir),
                    cran = TRUE)
})
