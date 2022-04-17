# Priority

## General

- add showProgress option for future MP


## Features

- Normalization
    - maybe: allow usage of normalized intensities with filter()?
    - find another way to assign close/far ISTDs: if there are multiple close ones available, it makes more sense to not consider those that are a bit far away.
    - reportHTML: report normalized intensities in featInfo?
    - Example ISTD lists in patRoonData
        - merge into master
    - normInts(): default normFunc OK? and others?
        - update newProject() for any changes


## TPs

- CTS
    - update SMILES to InChI/InChIKey/formula/neutralMass conversion after merging with msident
- Reporting
    - Split formulaDiff like plotGraph()?
- newProject()
    - selector for CTS transLibrary?
- plotGraph()
    - keep duplicate selector?


## MS library

- library
    - fixup annotation formulas too?
    - Always sanitize SMILES/InChIs/Formulas? Or make it optional and default enabled? Now only formulas are done...
    - also verify formulas of adducts?
    - include lib adducts optionally for adduct guessing
    - don't merge synonyms? or split them somehow on export? Otherwise doc
    - naming OK ('records')?
        - metaData?
    - Column name harmonization
        - More column names?
        - Make option in as.data.table() to enable/disable name harmonization?
    - sanitize input
        - trim spaces
        - Verify SMILES, InChI, InChIKey, Formula, neutralMass (from any of the previous, if available)
        - MoNA: take calculated SMILES? See e.g. "PR308903" and "PR308904"
- compounds
    - show mirror spectrum in report? Would need library data somehow
    - separate specSimParams for lib? E.g. to assume that lib spectra are cleaner and don't need intensity cleaning
    - collapse with IK1?
        - yes for now, as e.g. sets and consensus expects this (UID column)
    - support hits without matched peaks (i.e. like MF)? Would interfere with min sim score though
- suspects
    - anMSMSSimilarity doesn't make too much sense... either take (max if consensus?) libMatch or clearly doc
    - amend suspect lists with fragments from MSLibrary
    - update ID levels
- prepareChemTable(): for suspects and status messages/progress bars?


## NEWS

- prepareChemTable() etc: faster and more thorough calculations, fixes for labeled compounds, ...
- MP: more workarounds to handle NA exit codes on Linux



## Tests


## Docs

- TPs
    - handbook
        - update for transformationProductsStructure
- Features
    - Update handbook for normFunc --> normalized
    - update sets workflows section(s)
        - standards argument for normInts
        - set argument for plotGraph
        - normInts works per set


## NEWS

- TPs
    - sim calc doesn't happen by default anymore for BT
    - TP structure base class
        - TPLibrary gains filter() method and parent structure sim calculation
    - filter() for TP base class
    - Fixed: convertMFDB() now always collapses duplicates, not just for BT
    - BT
        - Simplified/Harmonized several columns
        - Converted ID and parent IDs to integer values
        - Removed several unnecessary `parent_` columns (parent_SMILES, etc)
        - equal TPs formed from different routes (but from the same initial parent) are now named the same
        - Fixed: retDir is now derived from _original_ parent, i.e. not its direct parent
        - TP parent expansion
        - steps --> generations for consistency with other algos
    - Lib
        - IDs
        - generations
        - name --> name_lib, new name column for consistency
        - caching
    - convertToSuspects()/generateComponentsTPs(): only include TPs that are unique from each parent
    - convertToMFDB()/generateComponentsTPs(): don't include columns that are specific to a parent/TP pair
    - additional args for generateTPs()
- Removed onlyLinked argument from plotGraph generic (not components methods)
- reportHTML
    - features menu with submenus
    - improvements TPs
        - simplified parent table with separated plots
        - plotGraph() (needs TP arg)
    - Fixed: properly subscript negative element counts in formulas
- normalization changes
    - Added slots to fGroups --> cache should be cleared
    - normFunc --> normalized for plotInt(), generateComponentsIntClust(), as.data.table()
        - with auto normalization for compat
    - as.data.table() can now report normalized values for avaraged feature data (features && average && normalized) == TRUE
    - removeISTDs arg for fGroups filter()
    - plotGraph() for fGroups
    - normInts() method for fGroups
        - new normalization methods: istd, tic, conc amongst features. Normalization across groups as before.
    - normalized argument for groupTable() and plotVolcano()
    - internalStandards() / internalStandardAssignments() methods (and slots)
    - norm_conc column in analysisInfo
    - newProject() changes
- Fixed: annotatedPeakList() for compounds: avoid _unset suffixes in mergedBy column
- Fixed: newProject(): loading analysis info from CSV now works on Windows
- norm_concs argument to generateAnalysisInfo()


# Lower priority


## General

- add 'keep.rownames = FALSE' to all as.data.table methods (or see if there is a work-around)
- remove mz column from patRoonData suspects?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
- allow specifying average function in other places where as.data.table() is used (eg clustering, plotting etc)
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
- optimize hashing? Or further avoid hashing big objects/objects with lists?
- load OpenMS intensities in parallel
    - either with futures or with MP and cache intensities afterwards
- XCMS: multiple features grouped in same analysis?
    - can be, but now handled by default method="medret" param. Make this configurable?
- updatePICSet(): also sync peaks list? otherwise doc
- Somehow integrate XCMS::fillChromPeaks
- Export generic EIC generation, i.e. without the need of feature data


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
        - apply sort fix? https://github.com/osenan/cliqueMS/issues/8
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?
- int and others: also use calculateComponentIntensities() for intensities?
- plot doesnt work for componentsReduced that originates from cluster components
    - maybe drop reduced mechanism?
- intclust
    - optionally take areas instead of intensities
    - cache results
- import check sessions?
    - needs way to match component names
    

## TPs

- predictTPsBioTransformer()
    - do we still need to check for non-calculated formulae?
- spectrumSimilarity
    - defaults OK for sim params?
        - precursor FALSE?
        - thresholds not really handy for formulas/compounds
            - at least doc that annotation results may disappear


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
    - as.data.table: normalization?
    - plotInt sets arg?
- components: somehow verify adductConflictsUsePref


## docs

- ref docs
    - delete()
        - mention j=DT for fGroups method?
    - order of sets sections
- handbook
    - TPs
        - add MSPL filtering of annotated peaks (in examples)?
    - update
        - introduction
            - mention changes of patRoon 2.0?
- tutorial
    - check if all is OK
- TPs
    - doc that merging TPs (same fGroup/TP) could be done with suspect screening
    - logic: mention results can be filtered with TP components?
    - mention how sets filter work for componentsTPs?


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
- expand reporting
    - eg include suspect name in EICs
        - already now in featInfo
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects



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


## Reporting

- add more options to reportPlots argument of reportHTML()?
- onlyAnnotated argument

