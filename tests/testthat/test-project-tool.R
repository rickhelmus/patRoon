# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

testDir <- getWorkPath("test-np")
suspectsPaths <- list(
    pos = getWorkPath("suspects_pos.csv"),
    neg = getWorkPath("suspects_neg.csv"),
    istd_pos = getWorkPath("istds_pos.csv"),
    istd_neg = getWorkPath("istds_neg.csv")
)

defaultSettings <- list(
    general = defaultGeneralSettings(testDir),
    analyses = defaultAnalysesSettings(),
    preTreatment = defaultPreTreatSettings(),
    features = defaultFeaturesSettings(),
    annotations = defaultAnnotationSettings(),
    TPs = defaultTPSettings(),
    report = defaultReportSettings()
)

complexSettings <- list(
    general = modifyList(
        defaultGeneralSettings(testDir),
        list(
            ionization = "both",
            IMS = list(mode = "post", CCSMethod = "agilent")
        )
    ),
    analyses = modifyList(
        defaultAnalysesSettings(),
        list(
            generateAnaInfo = "example"
        )
    ),
    preTreatment = defaultPreTreatSettings(),
    features = modifyList(
        defaultFeaturesSettings(),
        list(
            featAlgo = "EIC",
            featEICParams = list(
                methodMZ = "suspects",
                suspects = list(sets = suspectsPaths[c("pos", "neg")])
            ),
            suspects = list(sets = suspectsPaths[c("pos", "neg")]),
            IMSSuspCCSPred = "pubchemlite",
            fGroupsAdv = list(
                ISTDLists = list(sets = suspectsPaths[c("istd_pos", "istd_neg")]),
                featNorm = "istd"
            )
        )
    ),
    annotations = modifyList(
        defaultAnnotationSettings(),
        list(
            formulasAlgo = "GenForm",
            compoundsAlgo = "MetFrag",
            componAlgo = "nontarget",
            compCCSPred = "pubchemlite"
        )
    ),
    TPs = modifyList(
        defaultTPSettings(),
        list(
            TPsAlgo = "biotransformer",
            TPGenInput = "suspects",
            TPSuspectList = suspectsPaths$pos
        )
    ), 
    report = modifyList(defaultReportSettings(), list(reportLegacy = c("CSV", "PDF")))
)

# HACK: make sure these following Shiny uupdate functions correctly set their data with testServer()
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

test_that("Setting serialization", {
    settingsFile <- getWorkPath("project_settings.yml")
    local_mocked_bindings(
        selectFile = function(...) settingsFile,
        .package = "rstudioapi"
    )
    
    testServer(newProjectServer(testDir), {
        session$setInputs("general-saveParams" = 1)
        session$flushReact()
        checkmate::expect_file_exists(settingsFile)
        expect_equal(readProjectSettings(settingsFile, testDir), defaultSettings)
    })

    # change a setting, and see it is correctly loaded. There is no easy way to access the loaded settings, so we
    # simply write/load it from a YML settings file
    testServer(newProjectServer(testDir), {
        writeProjectSettings(complexSettings, settingsFile)
        
        session$setInputs("general-loadParams" = 1)
        session$flushReact()

        unlink(settingsFile)
        session$setInputs("general-saveParams" = 1)
        session$flushReact()
        
        expect_equal(readProjectSettings(settingsFile, testDir), complexSettings)
    })
    
    local_edition(3)
    expect_snapshot(cat(getScriptCode("calibrant.d", list(), complexSettings, noDate = TRUE)))
})
