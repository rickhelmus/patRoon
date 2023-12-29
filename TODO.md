# Priority

- better explain behavior of collapseSuspects=NULL and merging of predicted concs/tox
    - check if ref docs are OK
    - add to handbook
    - add NEWS

## TC

- annotations slot with sets workflows: update eg if features of one set are all removed
- annSims
    - let lib compounds use the same specsimParams default? Then annSim == libMatch by default.
    - remove support for annSimBoth in annotateSuspects? Or optionally calc for feat annotations and copy that?
        - either of these is preferred, as annSimBoth now needs separate specSimParams to annotateSuspects
- anaInfo
    - fGroups/feat subset
        - deprecate rGroups subset/filter param?
    - rename group col to replicate and warn for deprecation (like blank/ref column)
    - see if MSPL sets anaInfo can be replaced by non-exported slot with named vector with sets for each ana
    - properly handle NA values for custom cols?
    - plotVenn sets method: deprecate sets arg?
    - as.data.table()
        - avoid duplicate column names for intensities (eg if values from average_group clash)
            - always prefix? will break quite a bit...
            - auto prefix? bit inconsistent
            - or check upfront and give error if there are clashes
            - or let the user choose
        - include ion mz column for sets workflow?
            - split columns per set if features==F?
            - also if averaging with features==T? or only if average!="group"? or add ion_mz column in features?
            --> add ion_mz in features and set specific ion_mz columns in fGroups@annotations
        - adduct column
            - rename to group_adduct with features==T
            - split for sets instead of collapsing (to be consistent with ion_mz)?
            - add feature specific adducts if present
        - add specific columns from anaInfo?
            - prefix to avoid potential duplicates?
            - only if average==F? Or check if numeric?
    - FC: select multiple groups and/or anaInfo cols?

- tests
    - IDL filter
    - annSim: jaccard (as was done for suspects)
    - as.data.table()
        - changed/new average, regression, regressionBy args
    
- docs
    - specSimParams --> specSimParamsMatch
    - ES contributions for IDLs
    - explicitly mention annSim can be filtered with scoreLimits?
    - analysisInfo slot/accessor is now data.table()
    - reorderAnalyses(): doc that XCMS/XCMS3/KPIC2 fGroups internal slot is not updated, maybe also improve general docs for what is updated and for XCMS what it means for exporting data
    - analysisInfo<-()
        - doc methods (eg add limitations)
        - doc that anaInfo shouldn't be changed by reference?
    - new subset args (ni, reorder, reorder sets)
    - plotVenn(): average arg
    - plotUpSet()
        - average arg
        - nsets can be NULL
    - as.data.table()
        - average arg
            - different for features==T&&average==T
        - regressionBy arg
            - mention that "set" can be used for sets workflows
            - regressionBy column with features=T
            - mention that regressionBy values for average groups should be equal
        - fix: mention how pred results are merged with collapseSuspects=NULL (ie all non-suspect type values are removed if fGroup has suspect result)
        - updates for changed regression arg and conc_reg --> x_reg

- NEWS
    - specSimParamsMatch --> specSimParams
    - annSuspects should be faster now (no need to calc annSims, and estIDLevel is faster)
    - analysisInfo slot/accessor is now data.table()
    - no analysisInfo slot in fGroups anymore
    - analysisInfo()<-
    - feature subsetting: ni and reorder args (also reorder sets)
    - plotVenn(): average arg
    - plotUpSet():
        - average arg
        - nsets can be/defaults to NULL
    - as.data.table()
        - average arg
            - .all replaces average==T if features==T and also supports features==F
        - regressionBy arg
        - regression for features==F&&average==T: use average conc instead of first of each replicate
        - FC now possible with features==T
        - adduct column now split for features
        - updates for changed regression arg and conc_reg --> x_reg
    - sets workflows/screening: fixed getAllSuspCols() --> may affect formRank/compRank, filter() and reporting

## Ext

- R bundle
    - steps
        - Create a shortcut/batch file?
            - use R.utils::createWindowsShortcut() to link Rgui.exe
- update GHPreRel ref in handbook
- cleanup
    - remove appveyor files/scripts from patRoon/patRoonDeps

## Predict

- general
    - also score formulas?
        - general score mechanism currently not there
    - re-rank after addCompoundScore()
        - would make sense, but might break eg suspect annotation ranks
- quant
    - target RFs: those that are specified by user in a suspect list, to be used directly by calculateConcs()
    - filter suspects/annotations from conc results, ie to remove candidates with very low concs

## General

- add showProgress option for future MP
- ppm spectral averaging (Ricardo)
- SIRIUS5: doc that tool should be present even if host is not a worker. Or another approach?
- No need for R.utils dep if normalizePath() is used everywhere instead of getAbsolutePath()?
- get rid of magrittr dependency?

## Reporting

- get rid of pngquant
    - fully remove remainder when old reportHTML is fully removed (installation script, Dockerfile, options, ...)


## Features

- Screening
    - sets: warn that re-grouping (adduct()<- / selectIons) will discard screening results
        - or somehow handle this? re-screen? although original suspects will be lost...
    - function that copies suspect adduct to annotations slot?
        - can only work with exactly one match
            - could throw a warning when this happens
        - tricky with sets
            - doc that screenSuspects need to be re-run?


## TPs

- newProject()
    - selector for CTS transLibrary?


## MS library

- prepareChemTable()
    - skip/warn if obabel is not available?


# Low priority


## General

- add 'keep.rownames = FALSE' to all as.data.table methods (or see if there is a work-around)
- remove mz column from patRoonData suspects?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
- delete() for other classes
- generalize makeLegend() to new plot util


## Features

- misc
    - OpenMS: alignment may yield negative RTs...
    - OpenMS MapAligner exception
        - seems to be related when little overlap between sets --> add note in doc?
- adduct annotations
    - selectIons(): prefer adducts based on MS/MS? eg handy for Na/K adducts
    - what to do with unsupported adducts for annotation?
	    - default selectIons() to only consider 'common' adducts? or change default adducts for componentization algos?
	    - check better for what is supported by SIRIUS?
- import XCMS features: verify anaInfo (or remove necessity, eg with importAnaInfo func)
- getEICsForFeatures method for kpic2?
- load OpenMS intensities in parallel
    - either with futures or with MP and cache intensities afterwards
- XCMS: multiple features grouped in same analysis?
    - can be, but now handled by default method="medret" param. Make this configurable?
- updatePICSet(): also sync peaks list? otherwise doc
- Somehow integrate XCMS::fillChromPeaks
- Normalization
    - maybe: allow usage of normalized intensities with filter()?
    - find another way to assign close/far ISTDs: if there are multiple close ones available, it makes more sense to not consider those that are a bit far away.


## Annotation

- SusDat MF support
- parallel MSPeakLists generation?
- somehow handle different fragment formula annotations when making a consensus between formula/compounds objects
- DA formulas: also rank formula results like GF/SIRIUS?
- plotSpectrum/spectrumSimilarity: allow separate MSLevel for comparisons
- Support multiple MS/MS formula annotation candidates (ie same MS/MS peak annotated with different formulas)
    - mainly relevant for GenForm
    

## Components

- feature components
    - cliqueMS
        - change checkPackage GH link once PRs are merged
        - maxCharge --> chargeMax (same as OpenMS)? update docs
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?
- int and others: also use calculateComponentIntensities() for intensities?
- intclust
    - optionally take areas instead of intensities
    - cache results
- import check sessions?
    - needs way to match component names
    

## tests

- checkComponents() / checkFeatures()
    - server tests?
    - import
- MSPeakLists and others?: also test object that is fully empty (now still has analyses)
- ensure peaklists are sorted
- features
    - new multiple blank filtering
    - syncing of XCMS/KPIC2 objects
    - check if featindex and groups slots are in sync with features
    - subsetting and groupScores
    - plotInt sets arg?
- components: somehow verify adductConflictsUsePref

## MS library

- library
    - fixup annotation formulas too?
        - for now a little bit in C++ fixAnnFormula(), maybe implement verifyFormulas() in C++ someday?
    - future
        - plotSpectrum, plotVenn, same data format for MSPeakLists
- compounds
    - show mirror spectrum in report?
        - Would need library data somehow --> perhaps compoundsLibrary can include averaged lib spectra


# Future


## General

- test negative subset indices
- convertMSFiles(): Support OpenMS vendor conversion? (eg thermo)
- newProject()
    - also allow suspect annotation with only peak lists? currently only selectable if formulas/compounds selected
- Reduce non-exported class only methods
- future MP
    - delayBetweenProc?
    - batch mode
- msPurity integration
- algorithmObject() generic: for xset, xsa, rc, ...
- more withr wrapping? (par)
- newProject()
    - concentration column for anaInfo
    - generate more detailed script with e.g. commented examples of subsetting, extraction etc
	- import Bruker seq file?
    - fix multi line delete (when possible)


## Features

- makeSet(): also support fGroups method via comparison?
- feature optim:
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- OpenMS: Support KD grouper?
- Integration of mzMine features (package pending...), MS-DIAL and peakonly?
- suspect screening
    - automatic suspect list name assignment if that's lacking? might be handy for some NORMAN lists
- topMost filter that accepts rGroups, either as AND or OR

## Annotation

- MSPeakLists
    - isotope tagging is lost after averaging
    - test avg params
    - metadata() generic?
    - DA
        - generateMSPeakListsDA: find precursor masses with larger window
        - tests
            - utils? EICs with export/vdiffr?
            - test MS peak lists deisotoping?
- metadata for Bruker peaklists?
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- do something about negative H explained fragments by MF?
- MetFrag: auto-include suspect results if suspectListScore is selected?
- do something with sirius fingerprints? --> comparison?
- fix compoundViewer
- add new MF HD scorings and make sure default normalization equals that of MF web
- CFM-ID and MS-FINDER integration
- utility functions to make custom DBs for MetFrag and SIRIUS and support to use them with the latter
- DBE calculation for SIRIUS?
- OM reporting
- as.data.table: option to average per replicate group?
- ID levels for non-suspects
    - function to calculate ID levels from suspect list (to take RTs/MSMS if available), formulas, compounds
    - store in compounds?
    - does it make sense for formula candidates?
    - add into reporting
        - also mark if in suspect list

## Suspects

- ID level rules: add scorings for SIRIUS/DA
- interface
    - also convert TASQ?
    - annotateSuspects()
        - check why it's is sometimes slow
            - seems to be logging, disable by default? --> only slow with testthat?
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
- misc
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
            - could also fix column names, replace "-" with NAs etc
        - if yes, mention in ref docs for screenSuspects()


## components
- mass defect components
- split peak correlation and adduct etc annotation? would allow better non-target integration
- fillPeaks for CAMERA (and RAMClustR?)
- feature components
    - cliqueMS
        - current adduct conversion to this format doesn't mimic Cat and 2H/2Na etc
            - Perhaps just document limitation?
    - minimal annotation abundance across analyses (eg adduct must be annotated in >=X analyses)?
    - prefAdducts: also include eg Na by default?


## Sets
- compound/formula set consensus
    - weights for ranking (like compound consensus)?


## TPs

- filter on stability/persistence/toxicity of TP?

