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
        settings$preTreatment <- modifyList(settings$preTreatment, preTreatment)
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
    # UNDONE: verify limits.yml instrument
    testNewProj(general = list(IMS = list(mode = "post", CCSMethod = "agilent")), name = "general-ims_post_agilent")
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
