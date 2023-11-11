context("compounds")

fGroups <- getCompFGroups()

plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- plists[FALSE, reAverage = TRUE]
plistsEmptyMS <- removeMSPlists(plists, "MS")
fGroupsEmpty <- getEmptyTestFGroups()

doMetFrag <- TRUE # !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))
doSIRIUS <- TRUE # !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

if (doMetFrag)
{
    compsMF <- callMF(fGroups, plists)
    ct <- annotations(compsMF)
    compsMFEmpty <- callMF(fGroupsEmpty, plists = plistsEmpty)
    compsMFEmptyPL <- callMF(fGroups, plists = plistsEmpty)
}

if (doSIRIUS)
{
    if (FALSE)
        updateSIRIUSCompProj(fGroups, plists)
    
    compsSIR <- doGenCompsSIR(fGroups, plists)
    compsSIREmpty <- doGenComps(fGroupsEmpty, plistsEmpty, "sirius")
    compsSIREmptyPL <- doGenComps(fGroups, plistsEmpty, "sirius")
}

mslibrary <- loadMSLibrary(getMSLibJSONPath(), "json")
compsLib <- doGenComps(fGroups, plists, "library", MSLibrary = mslibrary, minSim = 0) # set minSim to 0 to get more 'hits'
compsLibEmpty <- doGenComps(fGroupsEmpty, plistsEmpty, "library", MSLibrary = mslibrary)
compsLibEmptyPL <- doGenComps(fGroups, plistsEmpty, "library", MSLibrary = mslibrary)

test_that("verify MetFragCL compound generation", {
    skip_if_not(doMetFrag)
    expect_known_value(compsMF, testFile("compounds-mf"))
    expect_known_show(compsMF, testFile("compounds-mf", text = TRUE))
    expect_length(compsMF, length(groupNames(compsMF))) # should be one compound per feature group
    # make sure that all suspect names correspond to identified compounds
    expect_true(all(sapply(names(ct), function(grp) nrow(ct[[grp]]) == 1 && ct[[grp]]$identifier == screenInfo(fGroups)[group == grp]$name)))
    expect_length(compsMFEmpty, 0)
    expect_length(compsMFEmptyPL, 0)
})

test_that("verify SIRIUS compound generation", {
    skip_if_not(doSIRIUS)
    expect_known_value(compsSIR, testFile("compounds-sir"))
    expect_known_show(compsSIR, testFile("compounds-sir", text = TRUE))
    expect_length(compsSIREmpty, 0)
    expect_length(compsSIREmptyPL, 0)
    expect_length(doGenComps(fGroups, plistsEmptyMS, "sirius"), 0)
})

test_that("verify Library compound generation", {
    expect_known_value(compsLib, testFile("compounds-lib"))
    expect_known_show(compsLib, testFile("compounds-lib", text = TRUE))
    expect_length(compsLibEmpty, 0)
    expect_length(compsLibEmptyPL, 0)
    expect_length(doGenComps(fGroups, plistsEmptyMS, "library", MSLibrary = mslibrary), 0)
    expect_length(doGenComps(fGroups, plists, "library", MSLibrary = mslibrary["none"]), 0)
    lmCol <- if (testWithSets()) "libMatch-positive" else "libMatch"
    expect_gte(min(as.data.table(doGenComps(fGroups, plists, "library", MSLibrary = mslibrary, minSim = 0.25))[[lmCol]],
                   na.rm = TRUE), 0.25)
    expect_true(all(!is.na(as.data.table(compsLib, fragments = TRUE)$frag_ion_formula))) # check presence annotations
})

test_that("verify fingerprints", {
    skip_if(!doSIRIUS || testWithSets())
    expect_gt(length(compsSIR), 0)
    testSIRFPSubset(compsSIR)
})

hasCompounds <- doMetFrag || doSIRIUS

if (doMetFrag)
{
    # include some isomers to test filtering... (sirius should already have multiple compounds for feature groups)
    compsMFIso <- callMF(fGroups, plists, db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
}

# continue with one or another...
comps <- if (doMetFrag) compsMFIso else if (doSIRIUS) compsSIR
compsEmpty <- if (doMetFrag) compsMFEmptyPL else if (doSIRIUS) compsSIREmptyPL
compsExplained <- filter(comps, minExplainedPeaks = 1)

if (!testWithSets())
{
    compsT1 <- filter(comps, topMost = 1)
    compsT1N <- filter(comps, topMost = 1, negate = TRUE)
}

test_that("delete and filter", {
    skip_if_not(hasCompounds)

    checkmate::expect_names(groupNames(delete(comps, i = 1)), disjunct.from = groupNames(comps)[1])
    checkmate::expect_names(groupNames(delete(comps, i = groupNames(comps)[1])), disjunct.from = groupNames(comps)[1])
    expect_length(delete(comps, i = groupNames(comps)), 0)
    expect_false(delete(comps, j = 1)[[1]]$UID[1] == comps[[1]]$UID[1])
    expect_length(delete(comps, j = 1), length(comps) - length(groupNames(comps)))
    expect_false(delete(comps, j = function(...) 1)[[1]]$UID[1] == comps[[1]]$UID[1])
    expect_length(delete(comps, j = function(...) 1), length(comps) - length(groupNames(comps)))
    expect_length(delete(comps, j = function(...) TRUE), 0)
    expect_equal(delete(comps, i = character()), comps)
    expect_equal(delete(comps, j = integer()), comps)
    expect_length(delete(comps), 0)
    expect_length(delete(compsEmpty), 0)
    
    expect_lte(length(filter(comps, topMost = 1)), length(fGroups))
    expect_lte(length(filter(comps, topMost = 5)), 5 * length(fGroups))
    expect_lte(length(filter(comps, topMost = 1, negate = TRUE)), length(fGroups))
    expect_lte(length(filter(comps, topMost = 5, negate = TRUE)), 5 * length(fGroups))

    expect_lte(length(filter(comps, minExplainedPeaks = 2)), length(comps))
    expect_length(filter(comps, minExplainedPeaks = 1E6), 0)
    expect_length(filter(compsEmpty, minExplainedPeaks = 2, topMost = 1), 0)
    expect_equivalent(filter(comps, scoreLimits = list(fragScore = c(-Inf, Inf))), comps)
    expect_lte(length(filter(comps, minExplainedPeaks = 2, negate = TRUE)), length(comps))
    expect_equivalent(filter(comps, minExplainedPeaks = 1E6, negate = TRUE), comps)
    expect_length(filter(compsEmpty, minExplainedPeaks = 2, topMost = 1, negate = TRUE), 0)
    expect_length(filter(comps, scoreLimits = list(fragScore = c(-Inf, Inf)), negate = TRUE), 0)

    expect_length(filter(comps, elements = "C1-100"), length(comps)) # all should contain carbon
    expect_length(filter(comps, elements = c("Na1-100", "C1-100")), length(comps)) # no sodium, but carbon should be there
    expect_length(filter(comps, elements = c("H1-100", "C1-100")), length(comps)) # presence of both shouldn't affect results
    expect_length(filter(comps, elements = "Na1-100"), 0) # no sodium
    expect_length(filter(comps, elements = "Na0-100"), length(comps)) # no sodium, but optional
    expect_length(filter(comps, elements = "C1-100", negate = TRUE), 0)
    expect_length(filter(comps, elements = c("Na1-100", "C1-100"), negate = TRUE), length(comps))
    expect_length(filter(comps, elements = c("H1-100", "C1-100"), negate = TRUE), 0)
    expect_length(filter(comps, elements = "Na1-100", negate = TRUE), length(comps))
    expect_length(filter(comps, elements = "Na0-100", negate = TRUE), 0)

    expect_lte(length(filter(compsExplained, fragElements = "C1-100")), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "C"), 0) # fragments may contain only single carbon
    expect_length(filter(compsExplained, fragElements = "Na1-100"), 0)
    expect_length(filter(compsExplained, fragElements = "Na0-100"), length(compsExplained))
    expect_lte(length(filter(compsExplained, fragElements = "C1-100", negate = TRUE)), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "C", negate = TRUE), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "Na1-100", negate = TRUE), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "Na0-100", negate = TRUE), 0)

    expect_length(filter(compsExplained, lossElements = "C0-100"), length(compsExplained))
    expect_length(filter(compsExplained, lossElements = "Na0-100"), length(compsExplained))
    expect_gt(length(filter(compsExplained, lossElements = "C1-100")), 0) # NL might be empty, at least some should contain carbon though!
    expect_length(filter(compsExplained, lossElements = "Na1-100"), 0) # no sodium
    expect_length(filter(compsExplained, lossElements = "C0-100", negate = TRUE), 0)
    expect_length(filter(compsExplained, lossElements = "Na0-100", negate = TRUE), 0)
    expect_gt(length(filter(compsExplained, lossElements = "C1-100", negate = TRUE)), 0)
    expect_length(filter(compsExplained, lossElements = "Na1-100", negate = TRUE), length(compsExplained))

    expect_lt(length(filter(compsLib, minScore = 0.25)), length(compsLib))
    expect_lt(length(filter(compsLib, minScore = 0.25, negate = TRUE)), length(compsLib))
    expect_equal(filter(compsLib, minScore = 0.25), filter(compsLib, scoreLimits = list(libMatch = c(0.25, Inf))))
    
    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIso, minScore = 0.75)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minScore = 0.75, negate = TRUE)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200, negate = TRUE)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, scoreLimits = list(fragScore = c(200, Inf)), negate = TRUE)),
              length(compsMFIso))

    skip_if_not(doSIRIUS)
    expect_lt(length(filter(compsSIR, minScore = -200)), length(compsSIR))
    expect_lt(length(filter(compsSIR, minScore = -200, negate = TRUE)), length(compsSIR))
    expect_lt(length(filter(compsSIR, scoreLimits = list(score = c(-200, Inf)))), length(compsSIR))
    expect_lt(length(filter(compsSIR, scoreLimits = list(score = c(-200, Inf)), negate = TRUE)),
              length(compsSIR))
    
    skip_if(testWithSets())
    
    # in case of ties between pos/neg the score is sometimes not the highest --> skip test with sets for now
    expect_true(all(sapply(annotations(compsT1[groupNames(compsT1N)]), function(a) max(a$score)) >=
                        sapply(annotations(compsT1N), function(a) max(a$score))))
})

test_that("basic subsetting", {
    skip_if_not(hasCompounds)

    expect_length(comps["nope"], 0)
    expect_equivalent(groupNames(comps[1:2]), groupNames(comps)[1:2])
    expect_equivalent(groupNames(comps[groupNames(comps)[2:3]]), groupNames(comps)[2:3])
    expect_equivalent(groupNames(comps[c(FALSE, TRUE)]), groupNames(comps)[c(FALSE, TRUE)])
    expect_equal(length(comps[FALSE]), 0)
    expect_length(compsEmpty[1:5], 0)

    expect_equivalent(comps[[1]], annotations(comps)[[1]])
    expect_equivalent(comps[[groupNames(comps)[1]]], annotations(comps)[[1]])
    expect_equivalent(callDollar(comps, groupNames(comps)[1]), comps[[1]])
})

test_that("as.data.table() works", {
    testFeatAnnADT(comps)
    normScName <- if (testWithSets()) "fragScore-positive" else "fragScore"
    expect_range(as.data.table(comps, normalizeScores = "max")[[normScName]], c(0, 1))
})

if (doMetFrag)
{
    forms <- doGenForms(fGroups, plists, "genform")
    compsMFIsoF <- addFormulaScoring(compsMFIso, forms)
}

test_that("formula scoring works", {
    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIsoF, minFormulaScore = 0.5)), length(compsMFIsoF))
    expect_error(addFormulaScoring(compsMFEmptyPL, forms), NA)
})

verifyMSPLAnFilter <- function(mspl, obj1, obj2 = NULL, negate = FALSE)
{
    allObj <- list(obj1)
    if (!is.null(obj2))
        allObj <- c(allObj, list(obj2))
    
    msplF <- filter(mspl, annotatedBy = allObj, retainPrecursorMSMS = FALSE, negate = negate)
    if (length(msplF) == 0)
        return(succeed("Empty MSPL"))
    
    doVerify <- function(pl, ao)
    {
        for (fg in groupNames(msplF))
        {
            if (is.null(pl[[fg]][["MSMS"]]))
                next
            
            objIDs <- unlist(lapply(ao, function(o)
            {
                if (is.null(o[[fg]]))
                    return(integer())
                return(lapply(o[[fg]]$fragInfo, "[[", "PLID"))
            }))
            
            if (negate)
                expect_true(!any(pl[[fg]]$MSMS$ID %in% objIDs))
            else
                expect_true(all(pl[[fg]]$MSMS$ID %in% objIDs))
        }
    }
    
    if (testWithSets())
    {
        for (s in sets(msplF))
            doVerify(msplF[, sets = s], lapply(allObj, function(o) o[, sets = s]))
    }
    else
        doVerify(msplF, allObj)
}

test_that("annotatedBy and results filters", {
    # NOTE: do these tests here, since we have all the objects ready: peaklists, compounds & formulas
    
    verifyMSPLAnFilter(plists, comps)
    verifyMSPLAnFilter(plists, forms)
    verifyMSPLAnFilter(plists, comps, forms)
    verifyMSPLAnFilter(plistsEmpty, compsEmpty)
    
    verifyMSPLAnFilter(plists, comps, negate = TRUE)
    verifyMSPLAnFilter(plists, forms, negate = TRUE)
    verifyMSPLAnFilter(plists, comps, forms, negate = TRUE)
    
    expect_setequal(groupNames(comps), names(filter(fGroups, results = comps)))
    expect_setequal(groupNames(comps), names(fGroups[, results = comps]))
    expect_setequal(union(groupNames(forms), groupNames(comps)), names(filter(fGroups, results = list(forms, comps))))
    expect_length(intersect(groupNames(comps), names(filter(fGroups, results = comps, negate = TRUE))), 0)
    
})

# on a clean system, i.e. where ~/.jnati/repo/jnniinchi is not yet initialized, starting multiple
# MetFrag processes in parallel (i.e. when maxProcAmount>1) may result in errors. This should now
# be fixed by setting a small delay between starting up processes (delayBetweenProc arg of executeMultiProcess())
jnatiTestDir <- file.path(tempdir(), "jnati-test")
test_that("MetFrag uninitialized jniinchi workaround", {
    # this test has been working for a very long time and eats up quite some time, disabled for now
    skip("not needed anymore")
    
    skip_if_not(doMetFrag)

    # temporarily change jnati workdir so it can be safely wiped
    withr::with_envvar(c(JAVA_TOOL_OPTIONS = sprintf("-Djnati.dir=%s", jnatiTestDir)), {
        withr::with_options(c(patRoon.cache.mode = "none"), {
            for (n in seq_len(5))
            {
                info = sprintf("iter: %d", n)
                unlink(jnatiTestDir, recursive = TRUE)
                expect_warning(compsJNI <- callMF(fGroups, plists), NA, info = info)
                expect_equal(compsJNI, compsMF, info = info)
            }
        })
    })
})

test_that("MetFrag can timeout", {
    # this test is probably obsolete, skip for now to save time
    skip("already tested with multiProc tests")
    
    skip_if_not(doMetFrag)
    withr::with_options(c(patRoon.cache.mode = "none"), {
        # call with unreasonably short timeout...
        expect_warning(compsTO <- callMF(fGroups, plists, to = 1))
        expect_lt(length(compsTO), length(compsMF))
    })
})

if (doMetFrag && doSIRIUS)
    compsCons <- doCompCons(compsMF, compsSIR, MSPeakLists = plists)

test_that("consensus works", {
    skip_if_not(hasCompounds)
    expect_length(doCompCons(comps, compsEmpty, MSPeakLists = plists), length(comps))

    skip_if_not(doMetFrag && doSIRIUS)
    expect_known_value(compsCons, testFile("compounds-cons"))
    expect_known_show(compsCons, testFile("compounds-cons", text = TRUE))
    expect_lt(length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, relMinAbundance = 1)), length(compsCons))
    expect_length(doCompCons(compsMFEmptyPL, compsSIREmptyPL, MSPeakLists = plists), 0)

    expect_equal(length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, uniqueFrom = 1)) +
                 length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, uniqueFrom = 2)) +
                 length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, relMinAbundance = 1)), length(compsCons))
    expect_equal(length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, uniqueFrom = 1:2, uniqueOuter = TRUE)) +
                 length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, relMinAbundance = 1)), length(compsCons))
    expect_length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, uniqueFrom = 1:2), length(compsCons))
    expect_lt(length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, uniqueFrom = 1:2, uniqueOuter = TRUE)),
              length(compsCons))
    expect_length(doCompCons(compsMFEmptyPL, compsSIREmptyPL, MSPeakLists = plists, uniqueFrom = 1), 0)
    expect_length(doCompCons(compsMFEmptyPL, compsSIREmptyPL, MSPeakLists = plists, uniqueFrom = 1, uniqueOuter = TRUE), 0)
})

anPLGroup <- screenInfo(fGroups)[name == "1H-benzotriazole"]$group
if (doMetFrag && doSIRIUS)
{
    anPL <- annotatedPeakList(compsCons, index = 1, groupName = anPLGroup, MSPeakLists = plists, formulas = forms)
    anPLOnly <- annotatedPeakList(compsCons, index = 1, groupName = anPLGroup, MSPeakLists = plists, formulas = forms,
                                  onlyAnnotated = TRUE)
}

test_that("annotation works", {
    skip_if_not(doMetFrag && doSIRIUS)

    expect_lt(nrow(anPLOnly), nrow(anPL))
    expect_true(any(is.na(anPL$ion_formula)))
    expect_false(any(is.na(anPLOnly$ion_formula)))
    expect_true(all(compsCons[[1]]$fragInfo[[1]]$ion_formula %in% anPLOnly$ion_formula))
    expect_true(any(grepl("metfrag", anPLOnly$mergedBy)))
    expect_true(any(grepl("sirius", anPLOnly$mergedBy)))
    expect_true(any(grepl("genform", anPLOnly$mergedBy)))
})

test_that("reporting works", {
    skip_if_not(hasCompounds)

    expect_error(reportCSV(fGroups, getWorkPath(), compounds = comps), NA)
    for (grp in names(annotations(comps)))
        checkmate::expect_file_exists(getWorkPath("compounds", sprintf("%s-%s.csv", class(fGroups), grp)))

    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = comps,
                           MSPeakLists = plists), NA)
    for (grp in names(annotations(comps)))
        checkmate::expect_file_exists(getWorkPath("compounds", sprintf("%s-%s.pdf", class(fGroups), grp)))

    expect_reportHTML(makeReportHTML(fGroups, compounds = comps, MSPeakLists = plists))
})

test_that("reporting empty objects works", {
    skip_if_not(hasCompounds)
    expect_error(reportCSV(fGroups, getWorkPath(), compounds = compsEmpty), NA)
    expect_error(reportCSV(fGroupsEmpty, getWorkPath(), compounds = comps), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = compsEmpty,
                           MSPeakLists = plistsEmpty), NA)
    expect_error(reportPDF(fGroupsEmpty, getWorkPath(), reportFGroups = FALSE, compounds = comps,
                           MSPeakLists = plists), NA)
    expect_reportHTML(makeReportHTML(fGroups, compounds = compsEmpty, MSPeakLists = plistsEmpty))
    expect_error(makeReportHTML(fGroupsEmpty, compounds = comps, MSPeakLists = plists), NA)
})

anPLGroup2 <- screenInfo(fGroups)[name == "DEET"]$group

test_that("plotting works", {
    skip_if_not(doMetFrag)

    # plotting structure seems to be difficult to do reproducible between systems, so disable for vdiffr now...
    expect_doppel("compound-spec", function() plotSpectrum(compsMFIso, 1, anPLGroup, plists, plotStruct = FALSE))
    expect_doppel("compound-spec_sim", function() plotSpectrum(compsMFIso, index = c(1, 1), c(anPLGroup, anPLGroup2),
                                                               MSPeakLists = plists, plotStruct = FALSE))
    expect_plot(plotSpectrum(compsMFIso, 1, anPLGroup, plists, plotStruct = TRUE))

    # plotStructure gives an empty plot??
    # expect_doppel("struct", function() plotStructure(compsMFIso, 1, names(annotations(compsMFIso))[1]))
    expect_plot(plotStructure(compsMFIso, 1, anPLGroup))
    expect_plot(plotStructure(compsMFIso, 1:2, anPLGroup))
    expect_doppel("scores", function() plotScores(compsMFIso, 1, anPLGroup))

    skip_if_not(doSIRIUS)
    expect_doppel("venn", function() plotVenn(compsMF, compsSIR))
    expect_error(plotVenn(compsMFEmpty, compsSIREmpty))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIR))$areas[2], length(compsSIR))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIREmpty))$areas[1], length(compsMF))
    expect_equal(expect_plot(plotVenn(compsMFEmpty, compsSIR))$areas[2], length(compsSIR))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIR))$intersectionCounts,
                 length(doCompCons(compsMF, compsSIR, MSPeakLists = plists, relMinAbundance = 1)))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIREmpty))$intersectionCounts, 0)

    expect_ggplot(plotUpSet(compsMF, compsSIR))
    expect_error(plotUpSet(compsMFEmpty, compsSIREmpty))
    expect_error(plotUpSet(compsMF, compsSIREmpty))
})

if (testWithSets())
{
    fgOneEmptySet <- makeOneEmptySetFGroups(fGroups)
    compsOneEmptySet <- callMF(fgOneEmptySet, plists)
    compsAvgSpecCols <- callMF(fGroups, plists, setAvgSpecificScores = TRUE)
}

test_that("sets functionality", {
    skip_if_not(testWithSets())
    skip_if_not(doMetFrag)
    
    expect_equal(comps, comps[, sets = sets(comps)])
    expect_length(comps[, sets = character()], 0)
    expect_equal(sets(filter(comps, sets = "positive", negate = TRUE)), "negative")
    expect_setequal(groupNames(comps), unique(unlist(lapply(setObjects(comps), groupNames))))
    expect_setequal(groupNames(unset(comps, "positive")), groupNames(setObjects(comps)[[1]]))
    expect_setequal(groupNames(unset(compsOneEmptySet, "positive")), groupNames(setObjects(compsOneEmptySet)[[1]]))
    expect_length(unset(compsOneEmptySet, "negative"), 0)
    
    expect_lt(length(callMF(fgOneEmptySet, plists, setThreshold = 1)), length(compsOneEmptySet))
    expect_length(callMF(fgOneEmptySet, plists, setThresholdAnn = 0), length(compsOneEmptySet))
    
    # setAvgSpecificScores=FALSE (default)
    checkmate::expect_names(names(as.data.table(comps)),
                            must.include = paste0("fragScore-", sets(compsAvgSpecCols)[1]))
    checkmate::expect_names(names(as.data.table(compsAvgSpecCols)), must.include = "fragScore")
    
    expect_doppel("compound-spec-set", function() plotSpectrum(compsMFIso, 1, anPLGroup, plists, plotStruct = FALSE,
                                                               perSet = FALSE))
    expect_doppel("compound-spec-set-perset", function() plotSpectrum(compsMFIso, 1, anPLGroup, plists,
                                                                      plotStruct = FALSE, perSet = TRUE, mirror = FALSE))
    expect_doppel("compound-spec-set-mirror", function() plotSpectrum(compsMFIso, 1, anPLGroup, plists,
                                                                      plotStruct = FALSE, perSet = TRUE, mirror = TRUE))
    
    skip_if_not(doSIRIUS)
    
    expect_gt(length(setObjects(compsSIR)[[1]]@fingerprints), 0)
    expect_gt(length(setObjects(compsSIR)[[2]]@fingerprints), 0)
    testSIRFPSubset(setObjects(compsSIR)[[1]])
})
