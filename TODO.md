# Priority

## general

- lossElements: handbook suggests it takes a whole formula --> maybe make fragFormula/lossFormula filters?
- BUG: annotatedBy filter incorrectly filters feature level peak lists
    - record group IDs in feature tabs and use them to remove
    - or work with tolerances, so that reAveraged data can work (but may be less accurate)
- better explain behavior of collapseSuspects=NULL and merging of predicted concs/tox
    - check if ref docs are OK
    - add to handbook
    - add NEWS
- update patRoonInst/patRoonExt/patRoonData for new deps
    - install timsconvert and c3sdb in Docker and bundle
    - embed TIMS-SDK in patRoonExt
    - add MRMTransitionGroupPicker to patRoonExt
    - switch to PCL w/ CCSbase in patRoonExt
- BUG: annSim.1 column (in formulas?)?
- fix default IDLs
- update REUSE
- Harmonize/consistent names
    - maxMSRtWindow --> maxMSRTWindow
    - be consistent in mobility vs IMS and mobilograms and EIMs
        - mobWindow and IMSWindow now randomly used
        - IMS arg but affects mobility features
        - withMobility vs withIMS
        - methodMZ and methodIMS --> methodMS+methodIMS or methodMZ+methodMob?
    - mob_area/mob_intensity vs mobmin/mobmax
    - naming: im_collapse and IMSCollapse
- FIX: IMSRangeParams: in sets workflows, the mz column is taken so relative filtering is done to neutralMasses
    - maybe add _mz rel column like compounds?
    - actually, the same is the case for the mz filter (also IMS and mz for features) --> maybe just document difference? the actual results will not change much in general


## Maybe

- see if MSPL sets anaInfo can be replaced by non-exported slot with named vector with sets for each ana
    - could be taken from analyses of sets objects
- FC: select multiple groups and/or anaInfo cols?
- groupBy: fGroups --> fGroup for consistency with anaInfo cols? (set/analysis/replicate)
- report
    - TP structures with fmcsR/depict
    - refactor/redesign
        - prefix cmp, tp etc columns
        - move main_columns to sysdata
        - move all JS callbacks etc to js file(s)?
        - improve code distribution over files
        - finish up
            - make title/href configurable?
    - overview tab of all EICs?
        - needs JS to set imgs from plot table and way to attach titles to imgs
    - comps-clust: don't have imgs double in reportPlots
- TPs
    - add set specific frag/NL matches for candidates
        - report them?
- Possible optimizations
    - Rcpp: see if it makes sense to return lists instead of dataframes, eg to further optimize MS library loading
    - replace makeFileHash() with getMSDataFileHash() on more places?
- msdata
    - store MSPL metadata?
        - better have a util that gets metadata from MS data
        --> remove metadata slot once we have that
    - collapseIMSFiles(): could use MSTK and/or SC for writing
        - both need to write a bit more metadata (instrumentConfiguration and dataProcessing), as OpenMS fails otherwise
        - best would be to do write spectra while reading, so we can do MP
            - or skip MP and do complete batch in C++
        - no support for mzXML
- add smoothing for XCMS3/enviPick peaks? e.g. with signal::sgolayfilt()
    - might also be nice for plotting chroms?
- add IMSMatch filter for suspects?
- convertMobilityToCCS() / convertCCSToMobility(): handle Waters data?
- start using saveCacheDataList() on more places?
- basic and default error handling for executeCommand()?
- getEIMs() function?
- add TP-formula filter() method, so parent-duplicates can be removed (ie removeParentIsomers filter)? If so, link to it in ann_form docs like ann_comp
- checkFeatures()
    - optionally disable mobilograms?
    - mobilogram previews?


## newProject()

- anaInfo.R
    - conc/norm conc are set to "NA_real_"
- remove conc column?
- remove default limits that are in limits.yml
- IMS changes
    - Codegen
        - Select right example anaInfo
        - Add mobility filters for suspects & compounds?
        - is assignMobilities() before TP componentization OK?
            - IMS="maybe", so probably fine?
        - predict CCS for TPs?

## Param defaults

- don't use defaults for eg feature finding (and grouping)?
- tweak agilent defaults --> after application?

## TC

- minAbundanceFeatAbs/Rel: does it make sense to also filter fGroup peak lists? If not, remove and update docs.

## msdata

- Agilent
    - SC doesn't recognize IM --> fixed now?

## IMS

- findPeaks()
    - Dietrich
        - Disabled noise removal for reported intensities/areas
            - keep doing this?
            - report both?
            - document
    - xcms3 for mobilities params seem not working
- assignMobilities()
    - better names for ims_parent_ID/ims_parent_group?
    - susps/compounds: overwrite doesn't overwrite converted mobilities, change?
- feat EICs
    - tweak default EIC filtering params, more strict for binning?
        - current defaults way to strict for patRoonData neg samples (for filtering of EICs _and_ peaks)
    - mob range: doc that defaults are for TIMS? Or just not give any defaults? Give instrument? Put it in limits?
    - split IMSWindow arg for setting mobmin/mobmax? this currently sets very narrow ranges. Same for mz and ms2
    - somehow better handle overlapping EICs for suspects? Ie for isobaric compounds, with m/z differences > very_narrow
- update featAnn consensus?
- greedy
    - set default weights
- tweak default EIMParams: clustering width should be much smaller (~0.01) --> sortof smooths data, mzExpIMSWindow can be reduced too?
    - check agilent data, if mobilities are equalish between frames too then just always use distance clustering with a small number (still configurable, remove cluster opt)


## Tests

- Features
    - SAFD?
- Feat Ann
    - MSPL
        - new spec averaging params: add IMS
        - more extensively test filters/summing/averaging
        - further test bbCID data
            - also mixed bbCID/PASEF file (ie with different segments)
            - are m/z ranges always present?
        - see if current default abundance threshold is fine
- TPs
    - verify TP retDirs?
- msdata
    - new verifyFileForFormat() usage in convertMSFiles()
    - file conversion?
        - with anaInfo function to go from eg IMS to centroid
        - testing with pwiz tricky for CI
        - would need raw data for other conversions
        --> check after we have IMS data
    - backends
        - availableBackends()
        - use eg mzR, MSTK, SC to get eg EICs or PLs
    - fixedIsolationWidth?
- IMS
    - IMSRangeParams/withIMSParent/applyIMS filters
    - more verification that normInts() works before/after assignMobilities()
    - IMS arg for [, filter(), ADT, and plotting functions
    - applyMS for filter()
    - IMSRangeFilter
    - expandForIMS()
        - verify things are copied or error is thrown for unsupported algos
    - convertMobilityToCCS() / convertCCSToMobility()
        - eg convert back and forth
    - suspects
        - test order of data selection for mobility and CCS columns, missing data etc
    - assignMobilities()
        - DT method
            - robustness with missing data in input/from
            - check warnings
        - fGroups methods
            - verify relevant slots are copied
            - from screenInfo(), eg for fGroups with >1 suspect assigned
            - mobPeakParams, fromSuspects, chromPeakParams, calcArea, fallbackEIC, CCSParams, IMSMatchParams
        - compounds method: IMS, from overwrite, CCSParams args
    - selectIMS filter: ensure that objects is fully reverted with IMS=FALSE
    - minMobSpecSim
        - verify that everything is copied, including fingerprints and scoreRanges
    - IMSRangeParams (incl mobility_mz/CCS_mz) and IMSMatchParams filters for suspects (screenSuspects() and assignMobilities()) and compounds
    - piek
        - tests with ims data
        - verify removal duplicate suspects/features?
    - groupFeaturesGreedy()
    - plotChroms(): mob annotation
    - plotMobilograms()
    - XCMS features/fGroups, eg with subsetting & exporting
    - normInts(): verify results are copied
    - CCS/Mob conversion utilities? Or rely on assignMobilities() tests?
    - spectrumSimilarityMobility()
    - susp lists with >1 mobility/CCS values. Test if mob/CCS counts differ.
    - feat comparison() and consensus()
- misc
    - genLimitsFile() and verify that it overrides
    - generateAnalysisInfo()?
- newProject()
    - see where testing is slow, possibly disable testServer() on CI?
    - add snapshots
- manual checks
    - checkFeatures()/checkComponents() verify if things still work
    - DA features and formulas?
    - manually check all HTML reporting functionality at the end

## Docs

- move all plotBPC(), getEIC() etc methods into one doc file?
- update relativeTo arg (not just replicates)
- mention that mzR backend caches as loads all at once
- convert
    - add patRoon 3.0 citation for im_collapse
- handbook
    - installation
        - add additional external deps, eg Python, TIMS-SDK, c3sdb, PCL/CCS etc
    - workflow concepts: maybe mention mobility as possible feature property
    - workflow gen
        - intro
            - update convertMSFiles() in table
        - annotation
            - update MSPL section
                - add back precursorMzWindow if we re-add it
            - estIDLevels
                - update thresholds when we change these
    - IMS chapter
        - intro
            - make sure to update step numbers of post assignment throughout the chapter if this is changed
            - table
                - add more rows, eg for componentization?
        - workflows
            - see which knitr chunks to enable for execution
            - details for IMS grouping
            - note that fallbackEIC only works for EIM mob assignments and setting chromPeakParams is mandatory when fromSuspects=TRUE
            - don't say that fGroups makeSet() is only supported for direct workflows
        - processing
            - evaluate the examples once they can be run
                - maybe add or improve plotMobilogram examples
    - Advanced
        - say something about consensus/comparison?
            - first finish (ref) docs for IMS grouping
    - mention somewhere that IM collapse doesn't export MS/MS spectra by default and why
- update tutorial
- IMS
    - piek
        - refs to paper
        - mention that bins are overlapping and twice as much
        - summarize code changes and fixes
        - mention that only mobility column for suspects is used so assignMobilities must be used for CCS conversion
    - assignMobilities()
        - refs to paper
        - fGroups method
            - add link to piek and import for direct IMS assignment
        - give some notes about sets workflows?
        - mention that fromSuspects only takes suspects assigned to a single fGroup (also in handbook?)
    - peakParams: refs to paper
    - Greedy
        - mention that links are removed for mobility features
        - mention that and why rtalign=FALSE by default
    - feat consensus
        - doc that all feats must have IMS data. For post workflows, IMS subset could be used
    - note that fallbackEIC only works for EIM mob assignments and setting chromPeakParams is mandatory when fromSuspects=TRUE

## NEWS

- formulas: calculateFeatures by default now FALSE
- TC
    - specSimParamsMatch --> specSimParams
    - annotateSuspects() --> estimateIDConfidence()
    - annotateSuspects()/estimateIDConfidence() now copies annSims from feat annotations instead of calculating
        - change in func args
        - annSimBoth is copied, so needs to be present in compounds (ie from IDLs)
    - annSuspects should be faster now (no need to calc annSims, and estIDLevel is faster)
    - sets annSims are now max value (as is done now for forms/comps)
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
        - nsets can be/defaults to NULL (all methods)
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
        - added susp_bestEstIDLevel column
        - moved docs to new page with all as.data.table() methods
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
    - delete() for screening: k arg
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
        - candidate specific frag/NL matches, also filter
        - new slots
    - generateTPs()
        - consistent "logP" naming, also for compounds
        - consistent calculation and configuration of logP/retDir calculation for all TP structure algos
        - new ann_form/ann_comp args
        - TPStructParams
        - log P tolerance for retDir calculation --> default enabled and has a large effect vs prev results --> doc how to revert to old behavior?
    - report()
        - fixed: TP graphs were generated for components with absent (parent) fGroups
    - compoundsLibrary: specSimParamsLib now defaults to specSimParams, and the latter now defaults to removing the precursor
    - replicate groups --> replicates
        - topMostByRGroup, replicateGroups(), replicateGroupSubtract()
        - NT rGroup column
        - rGroups arg for filters
    - plotVenn()/overlap()/unique(): removal of sets arg
    - plotVenn()/overlap/plotUpSet(): removed list arg possibility for which
    - FIXED: plotVenn() now warns when there are too many groups, instead of erroring
    - MetFrag: don't do MP for non-local databases to avoid connection errors
    - componentsNT: renamed "rt" to "ret" for consistency
- msdata
    - MSFileFormats() --> getMSConversionFormats() / getMSFileFormats()
    - new behavior of getMSFilesFromAnaInfo(): file types are checked one by one to avoid mixes and always checked to be present (mustExist was set a bit randomly...)
    - better checking of analysis file directory checking (verifyFileForFormat())
    - new/changed PListParams
    - (optimized getEICFGroupInfo() substantially)
    - loading OpenMS peak intensities is much faster, removed now unneeded intSearchRTWindow arg
    - getPICSet() isn't limited by centroided data anymore
    - updated convertMSFilesXXX() functions, including changed args
        - algo specific functions are now exported
    - generateMSPeakLists()
        - now uses backends, old methods still available but deprecated
        - no more precursorMzWindow, and avg params were changed
        - averaging is now intensity weighted
        - spec avg params now default to distance --> faster and suitable for large datasets (eg IMS)
    - EICParams for getPICSet() and calculatePeakQualities() (needed for m/z IMS expansion)
    - mzExpIMSWindow EIXParam
    - (subsetDTColumnsIfPresent: use order of requested cols instead of original) --> may change column order in some places
    - minIntensity arg for pwiz conversion, set by default
    - getBPCs(): doen't return m/z anymore
    - SAFD: fileType replaces profPath argument
    - patRoon.checkCentroided option removed (not needed anymore, checks are faster now and part of regular reading operations)
- IMS
    - hasMobilities slot for features
    - Dietrich features
    - groupInfo is now a DT
    - FIXED: suspect list assertion only checked part of the columns
    - getDefEICParams(): retWindow --> window
    - fGroupsScreeningSet doesn't have setObjects anymore
        - makes common operations faster
    - FIXED: fGroupsScreeningSet now correctly stores RF_SMILES per set
    - plotChroms(): annotate now has mob option
    - IMS arg for [, filter(), plotting functions, ADT, getQuantCalibFromScreening(), export(), getXCMS...()
    - applyIMS arg for all fGroups filter methods
    - filter: removeISTDs and onlyHits now use doFGroupFilter() --> caches and prints message
    - withIMSParent arg for filter()
    - ADT
        - character columns are now collapsed instead of removed with features=T and averaging
    - feature ID column is now always of character type
    - significantly optimized doScreenSuspects (emptyResults)
    - reporting
        - size optimizations, mainly for self contained (lzstring, no duplicate images)
        - candidate column in CSV of pred tables now doesn't contain images
- limits
    - lowered clusterMzWindow
- sets names are now checked to not contain any special characters (besides underscores). Automatic labels are now separated by underscores instead of dots.
- FIXED: SIRIUS with calculateFeatures=T may sometimes fail due to file name truncation
- generateAnalysisInfo()
- MS2QuantMeta slots for fGroupsScreeningSet, formulasSet and compoundsSet (latter already present but wasn;t filled in)
- Changed default db for MetFrag to pubchemlite
- Fixed: codegen for functions: don't break lines of args with long text strings
- Fixed: feature filter didn't hash qualityRange
- Fixed: newProject(): XCMS3/KPIC2 groupFeatures tried to do RT alignment in sets worfklow
- importFeatures() / importFeatureGroups()
    - made interface if both function more consistent: same (order) of arguments
    - added XCMS, XCMS3 and KPIC2 to importFeatureGroups()
    - renamed `path` to more generic `input`
    - removed `analysisInfo` argument from importFeatures()
    - added `table` import (importFeaturesTable() and importFeatureGroupsTable())
- overWrite --> overwrite for consistency, applies to predictCheckFeaturesSession() and file conversion (clearly doc in NEWS)

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

