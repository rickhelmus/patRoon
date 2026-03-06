# Priority

## general


- update patRoonInst/patRoonExt/patRoonData for new deps
    - install timsconvert and c3sdb in Docker and bundle
        - Docker: finish WIP
        - Windows: for now don't bundle, but install in GHA for testing: https://docs.github.com/en/actions/tutorials/build-and-test-code/python
    - update PCLiteCCS and update relevant README
    - update patRoonDeps dependencies for Rmstoolkitlib, reticulate, fmcsR etc
    - remove current hack to get in new deps for GHA and CircleCI
- BUG: annSim.1 column (in formulas?)?
- Harmonize/consistent names
    - better names for ims_parent_ID/ims_parent_group?
    - settle for terminology for IMS parents/IMS features
- See why parallel reporting is so slow (on Windows?)
- reticulate: use py_require()? doesn't seem to be able handle multiple py versions, so no?
- MSTK: PR for centroid status buf for mzXML?


## Maybe

- see if MSPL sets anaInfo can be replaced by non-exported slot with named vector with sets for each ana
    - could be taken from analyses of sets objects
- support filtering feature level peak lists by annotatedBy filter?
    - record group IDs in feature tabs and use them to remove
    - or work with tolerances, so that reAveraged data can work (but may be less accurate)
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
- convertMobilityToCCS() / convertCCSToMobility(): handle Waters data?
- start using saveCacheDataList() on more places?
- basic and default error handling for executeCommand()?
- getEIMs() function?
- add TP-formula filter() method, so parent-duplicates can be removed (ie removeParentIsomers filter)? If so, link to it in ann_form docs like ann_comp
- checkFeatures()
    - optionally disable mobilograms?
    - mobilogram previews?
- EIC optims
    - move filling and calcStats to C++? --> return as attributes?
- plotChroms() etc: mention that only precursors are plot? Now seems confusing as number of fGroups seems less when reporting
- greedy
    - somehow handle cases where...
        1. there is >=1 adjacent feature to the most intense
        2. the most intense groups poorer (more deviation, smaller in size) compared to the adjacent, because the adjacent is closer to the others
    - maybe a comprehensive mode?
        1. make groups for all features: do as now and take best scoring
        2. score groups from high to low
        3. iterate through groups, remove all features from current group that are also in others
        4. go to step 1, but ignore any features still in groups 
- mzRange, mzDefectRange and IMSRangeParams use neutral masses instead of m/z values in sets workflows --> change? (now documented)


## newProject()

- remove conc column?
- IMS changes
    - Codegen
        - Add mobility filters for suspects & compounds?
        - is assignMobilities() before TP componentization OK?
            - IMS="maybe", so probably fine?
        - predict CCS for TPs?
        - don't do CCS prediction with PCL+MetFrag?


## IMS

- rename and split assignMobilities()?

## Tests

- Features
    - SAFD?
- Feat Ann
    - MSPL
        - more extensively test filters/summing/averaging
        - further test bbCID data
            - also mixed bbCID/PASEF file (ie with different segments)
            - are m/z ranges always present?
- msdata
    - new verifyFileForFormat() usage in convertMSFiles()
    - file conversion?
        - with anaInfo function to go from eg IMS to centroid
        - testing with pwiz tricky for CI
        - would need raw data for other conversions
        --> check after we have IMS data
- misc
    - generateAnalysisInfo()?
- newProject()
    - see where testing is slow, possibly disable testServer() on CI?
- manual checks
    - checkFeatures()/checkComponents() verify if things still work
    - manually check all HTML reporting functionality at the end
    - check updates for mobilities to checkFeatures(), checkComponents() and importCheckFeaturesSession()
    - test updated verifyDependencies()


## Docs

- move all plotBPC(), getEIC() etc methods into one doc file?
- handbook
    - installation
        - verify R installation table when new deps in patRoonDeps are added
        - mention that for dev PKG_BUILD_EXTRA_FLAGS/pkg.build_extra_flags should be disabled
    - check: are TP docs updated for ADT(candidates=T)? And getDefTPStructParams()?
- add refs to paper: piek, assignMobilities(), peakParams, greedy, IMSCollapse
    - also in in NEWS (IMS, msdata etc)?


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

