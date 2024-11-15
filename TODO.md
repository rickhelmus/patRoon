# Priority

## general

- newProject: don't break lines of text strings, e.g. long paths with spaces in them passed to generateAnalysisInfo()
- lossElements: handbook suggests it takes a whole formula --> maybe make fragFormula/lossFormula filters?
- BUG: annotatedBy filter incorrectly filters feature level peak lists
    - record group IDs in feature tabs and use them to remove
    - or work with tolerances, so that reAveraged data can work (but may be less accurate)
- better explain behavior of collapseSuspects=NULL and merging of predicted concs/tox
    - check if ref docs are OK
    - add to handbook
    - add NEWS
- update sysdata.rda
- start using saveCacheDataList() on more places?
- basic and default error handling for executeCommand()?
- maxMSRtWindow --> maxMSRTWindow
- BUG: unescaped set names are now used regular expressions (problem if eg dot is present when make.unique was used)


## TC

- annSims
    - let lib compounds use the same specsimParams default? Then annSim == libMatch by default.
    - remove support for annSimBoth in annotateSuspects? Or optionally calc for feat annotations and copy that?
        - either of these is preferred, as annSimBoth now needs separate specSimParams to annotateSuspects
    - add annSims, IDLs etc to report
- anaInfo
    - fGroups/feat subset
        - deprecate rGroups subset/filter param?
    - rename group col to replicate and warn for deprecation (like blank/ref column)
    - see if MSPL sets anaInfo can be replaced by non-exported slot with named vector with sets for each ana
    - properly handle NA values for custom cols?
    - plotVenn/unique sets methods: deprecate/remove sets arg?
    - plotVenn()/overlap/plotUpSet(): deprecate/remove list arg for which?
    - FC: select multiple groups and/or anaInfo cols?
    - plotChord/report(): report sets as outer groups? make it optional?

- Misc
    - bg MS/MS subtraction
        - default abundance filtering while averaging? or apply filter in newProject()
    - don't try to do MetFrag/PubChem with MP (by default) or warn, also add docs why
    - MSPL filter to remove mass peaks with X Da higher m/z than precursor
        - do by default in newProject()?

- TP components
    - getTPParents(): remove unused calcLogP arg?
    - update patRoonInst for new deps?
    - hashing for all algos seem to take unnecessary variables (eg parent names, things for prepareChemTable() etc)
    - further limit candidate columns? Can get quite excessive
    - update filter for candidate specific frag/NL matches
    - add set specific frag/NL matches for candidates
        - report them?
    - sync topMost filter for ann TPs
    - collapseComponents(): not needed anymore? otherwise update
    - report()
        - BUG: JS errors when fGroups are subset/do not contain TPs anymore?
        - refactor/redesign
            - prefix cmp, tp etc columns
            - move main_columns to sysdata
            - move all JS callbacks etc to js file(s)?
            - improve code distribution over files
            - move flexdashboard to Suggests (and remove completely if legacy reportHTML is removed)
            - finish up
                - make title/href configurable?
        - TP structures with fmcsR/depict (maybe later?)
        - delta sign sometimes shown as N ?

- tests
    - IDL filter
    - annSim: jaccard (as was done for suspects)
    - as.data.table()
        - changed/new average, regression, regressionBy args
    - delete() for screening
        - check if SOs are properly synced
        - test if k=NA works
    - TP components
        - filter for candidate specific frag/NL matches
    
- docs
    - specSimParams --> specSimParamsMatch
    - ES contributions for IDLs
    - explicitly mention annSim can be filtered with scoreLimits?
    - annSim for sets is calculated as max
    - analysisInfo slot/accessor is now data.table()
    - reorderAnalyses(): doc that XCMS/XCMS3/KPIC2 fGroups internal slot is not updated, maybe also improve general docs for what is updated and for XCMS what it means for exporting data
    - analysisInfo<-()
        - doc methods (eg add limitations)
        - doc that anaInfo shouldn't be changed by reference?
    - new subset args (ni, reorder, reorder sets)
    - plotVenn()/unique()/overlap(): aggregate arg
    - overlap(): which can be NULL
    - plotUpSet()
        - aggregate arg
        - nsets can be NULL
    - plotChord()
        - aggregate arg (replaces average)
        - outerGroups --> groupBy and now expects anaInfo column name
            - also update (incorrect) example in handbook
    - plot()/plotChroms()
        - colourBy --> groupBy
            - "none" --> NULL
            - anaInfo col possible
    - as.data.table()
        - average arg
            - different for features==T&&average==T --> fGroups
            - rGroups also supported
        - regressionBy arg
            - mention that "set" can be used for sets workflows
            - regressionBy column with features=T
            - mention that regressionBy values for average groups should be equal
        - fix: mention how pred results are merged with collapseSuspects=NULL (ie all non-suspect type values are removed if fGroup has suspect result)
        - updates for changed regression arg and conc_reg --> x_reg
        - anaInfoCols arg
            - only works with features==T and cols must be numeric if average==T
        - make new page with all as.data.table() methods
    - minimum max intensity feature filter
    - maxMZOverPrecMS/maxMZOverPrecMSMS MSPL filters
    - plotInt()
        - xnames --> xNames
        - new args: areas, xBy, groupBy, averageFunc, regression
    - filter() for screening: k arg, including NA to clearout results for fGroups specified by j (or all if j==NULL)
    - MS/MS bg subtraction
        - abundance filter args
        - abundance columns in mspl
        - abundance parameter in average params
        - getBGMSMSPeaks()
        - new/shortened filter() args
        - removeMZs and mzWindow filter() args
            - removeMZs list for sets
    - generateComponentsTPs()
        - new args
        - as.data.table candidates arg
        - candidate specific frag/NL matches
        - new slots
    - generateTPs()
        - ann_comp/ann_form algos
        - forceCalcRetDir: only relevant for library atm, clearly mention its use in genTPsLib docs
        - TPStructParams: all parameters, explain when logPs are calculated (ie if both parent/TP is absent)
        - TPsComp/TPsForm: filter()
        - parallel=T is only useful with many candidates


- NEWS
    - specSimParamsMatch --> specSimParams
    - annSuspects should be faster now (no need to calc annSims, and estIDLevel is faster)
    - analysisInfo slot/accessor is now data.table()
    - no analysisInfo slot in fGroups anymore
    - analysisInfo()<-
    - feature subsetting: ni and reorder args (also reorder sets)
    - plot()/plotChroms()
        - colourBy --> groupBy
            - "none" --> NULL
            - anaInfo col possible
    - plotVenn()/unique()/overlap(): aggregate arg
    - overlap(): which can be NULL
    - plotUpSet():
        - aggregate arg
        - nsets can be/defaults to NULL
    - plotChord()
        - aggregate arg (replaces average)
        - outerGroups --> groupBy and now expects anaInfo column name
    - as.data.table()
        - average arg
            - "fGroups" replaces average==T if features==T and also supports features==F
        - regressionBy arg
        - regression for features==F&&average==T: use average conc instead of first of each replicate
        - FC now possible with features==T
        - adduct column now split for features
        - updates for changed regression arg and conc_reg --> x_reg
        - anaInfoCols arg
        - intensity cols are suffixed
            - clarify that same columns are used for areas?
        - adduct column: now split per set and renamed to group_adduct
    - plotInt()
        - xnames --> xNames
        - new args: areas, xBy, groupBy, averageFunc, regression
        - removed sets method, replaced by groupBy
    - sets workflows/screening: fixed getAllSuspCols() --> may affect formRank/compRank, filter() and reporting
    - sets features: ion_mz column
    - sets fGroups: ion_mz column in @annotations
    - minimum max intensity feature filter
    - maxMZOverPrecMS/maxMZOverPrecMSMS MSPL filters
    - The default metabolic logic transformations are now accessible through the `TPLogicTransformations()` function.
    - filter() for screening: k arg
    - MS/MS bg subtraction
        - abundance filter args
        - abundance columns in mspl
        - abundance parameter in average params
        - getBGMSMSPeaks()
        - new/shortened filter() args
        - removeMZs and mzWindow filter() args
    - generateComponentsTPs()
        - new args
        - as.data.table candidates arg
        - candidate specific frag/NL matches
        - new slots
    - generateTPs()
        - consistent "logP" naming, also for compounds
        - consistent calculation and configuration of logP/retDir calculation for all TP structure algos
        - new ann_form/ann_comp args
        - TPStructParams
        - log P tolerance for retDir calculation
    - report()
        - fixed: TP graphs were generated for components with absent (parent) fGroups

## msdata

- see if it makes sense to return lists instead of dataframes for other functions besides getEICList()
- minIntensityIMS
    - use it for all other functions that get EICs
    - default=25 is OK?
    - store it in an option? or otherwise in EICParams for fGroups?
    - remove default in doGetEICs() when finished
    - also use it for MSPL and EIMs?
- getEICList(): do mob assignment from base peak instead of average
- MT
    - test ThreadExceptionHandler a bit more, eg if subsequent run calls cancels well if an exception is caught
- cleanup
    - remove verifyDataCentroided()
    - remove now unused avgParams (eg averaging func)
    - get fully rid of precursorMzWindow (eg newProject()) and update NEWS/docs
    - move old MSPL generators to deprecate.R and add notifications
    - clean old tims files
- OpenTIMS doesn't seem to want to initialize the Bruker libs sometimes when repeatedly calling devtools::load_all()
    - forcing setup_bruker just before obtaining a handle seems to be a workaround
    - for now just throw an error if it's not loaded, could try above workaround if it seems to occur more frequently
- replace makeFileHash() with getMSDataFileHash() on more places?
    - for now leave, but can be done when things appear slow
- perform centroid checks in backends?
- Further test multithreading stability with OTIMS
- check assumptions about data being sorted
    - meatadata: scans, RTs, mobilities
    - spectra: m/z, mobilities (see if this can be used to speedup eg MSTK)
- MSPL
    - add mobility column?
    - store metadata?
    - more extensively test filters/summing/averaging
    - check defaults for new averaging params
    - withPrecursor: now only applied prior to other filtering systems, change?
        - if yes: adjust R and C++ code
        - doc in any case
    - further test bbCID data
        - also mixed bbCID/PASEF file (ie with different segments)
- the default backends now put mzR as last to ensure mobility data can be used, change this somehow?
- getEICs()
    - switch to anaInfo param (or optional?), otherwise add file type arg
    - allow multiple files?
    - see when we update/add similar raw data functions
- embed TIMS-SDK? --> in patRoonExt
- optionally link to MSTK depending on if it's available
    - somehow clearly state if unavailable during installation
- only enable OTIMS on Win/Lin x64
    - somehow clearly state if unavailable during installation?
    - should work after tims.cpp is removed --> test when done
- collapseIMSFiles(): see if this makes sense as an alternative to convertMSFiles()
    - pros: doesn't need pwiz, slightly faster per file conversion (eg 23 vs 30 secs) and produces smaller files (eg ~180 vs 350 mb)
    - cons: can't parallelize over multiple files, so probably not really faster, files are potentially with incomplete metadata
    - could also use MSTK and/or SC for writing
        - both need to write a bit more metadata (instrumentConfiguration and dataProcessing), as OpenMS fails otherwise
        - MSTK: load polarities
        - best would be to do write spectra while reading, so we can do MP
            - or skip MP and do complete batch in C++
        - no support for mzXML
    - initial tests of pwiz, mzR and MSTK show quite different numbers of OpenMS features...
    - if we stick with convertMSFiles(), somehow recommend skipping MS2+gzip?
- update generateAnalysisInfo() and newProject()
    - see if getAllMSFilesFromAnaInfo() can be used, otherwise remove it
- convertMSFilesXXX()
    - change function names?
    - add option for zlib compression, and enable by default? Esp useful for IMS/profile data
- more clearly state which backend is used, e.g. for general logging and citing

- test
    - new verifyFileForFormat() usage in convertMSFiles()
    - new spec averaging params

- NEWS
    - MSFileFormats() --> getMSFileConversionFormats() / getMSFileFormats()
    - new behavior of getMSFilesFromAnaInfo(): file types are checked one by one to avoid mixes and always checked to be present (mustExist was set a bit randomly...)
    - better checking of analysis file directory checking (verifyFileForFormat())
    - new/changed PListParams
    - (optimized getEICFGroupInfo() substantially)
    - loading OpenMS peak intensities is much faster, removed now unneeded intSearchRTWindow arg
    - getPICSet() isn't limited by centroided data anymore
    - updated convertMSFilesXXX() functions, including changed args
    - generateMSPeakLists()
        - now uses backends, old methods still available but deprecated

- docs
    - getMSFileFormats()
    - patRoon.threads, patRoon.MSBackends and patRoon.path.BrukerTIMS options
    - update PListParams
    - updated convertMSFilesXXX() functions, including changed args
    - update generateMSPeakLists()
    - availableBackends()

## IMS

- general ratio of mass/CCS to filter results --> eg derived from CCS DB
- be consistent in mobility vs IMS and mobilograms and EIMs
    - mobWindow and IMSWindow now randomly used
- make sure mobility ranges are set when getting EIMs
- findPeaks()
    - export? If yes, add checkmate's, documentation etc
        - if not, still need to doc args somehow
        - clashes with xcms::findPeaks()...
    - add smoothing for XCMS3/enviPick? e.g. with signal::sgolayfilt()
        - might also be nice for plotting chroms?
    - OpenMS
        - remove old function
        - do we need scaleTimeFactor? then add to params
        - check if we can get compute_peak_shape_metrics to work (maybe check OMS version?)
    - Dietrich
        - Disabled noise removal for reported intensities/areas
            - keep doing this?
            - report both?
            - document
    - Params
        - see if current chrom defaults are fine
        - add sensible IMS defaults
- findMobilities()
    - assert that instrument data has actually IMS data
    - better names for ims_parent_ID/ims_parent_group?
    - handle cases when there are already IMS assignments (or just throw an error?)
    - limit mobmin/mobmax --> both min and max, ie to prevent excessively wrong peak range assignments
    - SC seems to hang?
    - suppress XCMS warnings (or at least if no peaks are found)
    - remove mobilities slot
    - ID re-assignment works with e.g. XCMS IDs that may be numeric?
    - fGroups method
        - support updating groupQualities/groupScores/ISTDs/ISTDAssignments/annotations/concentrations/toxicities/screenInfo?
            - currently cleared with a warning
            - duplicate values where possible
        - sets methods
    - minMobilityMatches: move to filter()?
        - then needs to store nr of mobilities in suspect list (eg like MS2 peaks)
    - somehow handle eg XCMS objects which cannot do mobilities (eg clearout?)
    - test for fGroups from screenInfo(), eg for fGroups with >1 suspect assigned
    - clearly doc what IMSWindow is used for
    - see if mzWindow for assignFeatureMobilities() is now not needed anymore since getMobilograms() is fixed
    - store mobility assignment type in feature table? (eg mobilogram, suspect, PASEF)
    - method for fGroupsScreeningSet
    - clusterFGroupMobilities(): don't clear out annotation table (sets!)
        - copy data from IMS parents
        - maybe also for other slots? (with a warning/message?)
    - doc the use for fromSuspects, eg
        - doesn't rely on mobility peak detection, so might be less prone to false negatives with eg low intensities
        - scenario 1: we know the mobility very well, eg from a database --> use a narrow IMSWindow
        - scenario 2: we only have the mobility from eg a prediction and don't care so much about identification by CCS match --> use wide IMSWindow
        - somehow only do this as fallback when mobility peak detection fails?
            - doesn't seem very useful for approach 1&2 though (we already know the mobility or are not interested in an accurate one)
            - adds extra complexity when grouping mobility features, eg if different methods were used across the features
- splitMobilities()
    - keep? then update
    - remove unassigned features?
- plotMobilogram()
    - Finalize function: consistent arg names, checkmates, more functionality from plotChroms(), areas? --> see UNDONEs
    - by default only allow 1 fg to be plotted?
    - finish plotMobilogramHash()
    - tests, docs
- Suspect features
    - XCMS/XCMS3/KPIC2: doc and/or default min fractions to zero as these probably don't make a lot of sense otherwise
    - Handle mobilities
        - Does the mobmin/mobmax range make sense how it is computed now?
    - remove mobility assignment?
        - if not, support >1 mobilities in suspect list and do mobility assignment directly from suspect list like findMobilities()
- ISTDs
    - update slots in findMobilities()? Doc otherwise that normalization should be done after mob assignment
    - how do we assign ISTDs after mobilities are assigned?
        - need to make sure that minimum number of ISTDs is satisfied by number of unique IMS parents
        - for ISTDs with no suspect mobility and >1 mobility features, simply assign all+parent?
- reporting
    - convert reactable cell img functions to JS versions, so img paths are set dynamically which can save quite some space for selfContained reports
        --> see WIP changes made for chromsLarge
    - mobilograms
        - make setting (and use it instead of large chroms)
        - have small/big mobilograms? or put them in a separate tab instead of the table?

- Tests
    - more verification that fGroupsScreeningSets still works fine after removal of setObjects
    
- Docs
    - Dietrich features
    - getDefPeakParams()
    - getDefEICParams() / getDefEIMParams()
        - update for retWindow --> window
        - add docs for getDefEIMParams()
        - update handbook
    - findMobilities()
        - mention that feature properties (except intensity, rt, area) are simply copied from parent
- NEWS
    - Dietrich features
    - groupInfo is now a DT
    - FIXED: suspect list assertion only checked part of the columns
    - getDefEICParams(): retWindow --> window
    - fGroupsScreeningSet doesn't have setObjects anymore
        - makes common operations faster
    - FIXED: fGroupsScreeningSet now correctly stores RF_SMILES per set


## Features


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
- getFeatureEIXs method for kpic2?
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

