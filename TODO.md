# Priority

## General

- ggplot: remove?
    - Otherwise fix plotSpectrum() (sets) for title and ticks and possibly others
- newProject
    - sets: add specific options for annotation and others?
    - TP workflow
- add 'keep.rownames = FALSE' to all as.data.table methods (or see if there is a work-around)


## Features

- neutralizing / ionization
    - selectIons(): prefer adducts based on MS/MS? eg handy for Na/K adducts
    - what to do with unsupported adducts for annotation?
	    - skip calculation with a warning?
		    - default to M+H/M-H for now with warning...
	    - default selectIons() to only consider 'common' adducts? or change default adducts for componentization algos?
	    - check better for what is supported by SIRIUS?
- importFeaturesXCMS/importFeaturesXCMS3/importFeatureGroupsXCMS: get rid of anaInfo arg requirement? (or make import func?)
- Check: units of plotChord() rt/mz graphs seems off
- misc
    - topMostByRGroup: make default? or only for reporting?
    - quality/score filters
    - updatePICSet(): also sync peaks list? otherwise doc
    - getEICsForFeatures method for kpic2?
    - optimize hashing? Or further avoid hashing big objects/objects with lists?
    - delete for features and fGroups
        - XCMS: also update groups data? --> use new function once it hits BC
    - load OpenMS intensities in parallel
        - either with futures or with MP and cache intensities afterwards
    - OpenMS: alignment may yield negative RTs...

## Annotation

- GenForm oc is FALSE by default. OK? Add to newProject?
- ion_formula_mz --> ionMass?
- min score filters: take maximum value of all merged columns?
    - changed for now, update docs if keep
- formulas filter: MSMSScore doesn't remove non-MS/MS peaks anymore, OK? If yes, doc. Otherwise re-add tests
- default setThresholdAnn=0?
- update GenForm


## Components

- feature components
    - cliqueMS
        - change checkPackage GH link once PRs are merged
        - current adduct conversion to this format doesn't mimic Cat and 2H/2Na etc
            - Perhaps just document limitation?
    - minimal annotation abundance across analyses (eg adduct must be annotated in >=X analyses)?
    - OpenMS: handle potentialAdducts per set
    - prefAdducts: also include eg Na by default?


## sets

- misc
    - formulas/compounds: update set column when subsetting on sets?
    - fix MapAligner exception with test-components
    - compound/formula set consensus
        - weights for ranking (like compound consensus)?
    - disallow commas, minus and algo consensus labels in set names
    - screening: filter on set specific rank cols? Or all?
        - for now all if no sets columns are present (ie not ID level and ranks)


## TPs

- componentsTP
    - precursor_rt etc are from suspect list, not really clear --> also for reporting
    - show status after finishing
- predictTPsBioTransformer()
    - Include BT in installation script and verifyDependencies()
    - do we still need to check for non-calculated formulae?
    - disable parallellization (by default)?
- metabolic logic
    - cite 10.1021/acs.analchem.5b02905 and possibly others if more is included
- predictTPsComponents
    - fix if empty cTab for MSMS components
    - doc that precursor should not occur in multiple components (is this relevant for users?)
    - caching
    - remove? doesn't seem useful anymore
- log2fc
    - P values are calculated properly?
- spectrumSimilarity
    - plotting
        - formulas/compounds
            - doesn't work with structures at the moment, either fix or doc and set default to FALSE
            - put similarity in title (like for MSPL)?
    - defaults OK for sim params?
        - precursor FALSE?
        - thresholds not really handy for formulas/compounds
            - at least doc that annotation results may disappear
- Consistency
    - precursor diff: 1-2 or 2-1? --> verify all
- misc
    - show method for new components classes
    - truncate MP logfiles like with suspects, eg for long suspect names with BT
    - retDir: convert to character?
    - doConvertToMFDB/getCompInfoList: works fully for both BT and Lib columns?
    - only keep relevant suspect list columns when input is screening object for lib/BT TPs
- reporting
    - padding between two tables?
    - default TP columns OK?
    - fragMatches/NLMatches: doc that it's _not_ candidate specific
        --> also add for candidate specific if possible? although this could be taken from suspect annotations
        - otherwise maybe rename to eg allFragmentMatches
    - plotInt: not normalized now? y axis isn't...

## Reporting

- featInfo: finished?


## tests

- checkComponents() / checkFeatures()
    - server tests?
    - import
    - MC integration (features)
- test DA algorithms
    - formulas: check if fragInfo etc is correct
- MSPeakLists and others?: also test object that is fully empty (now still has analyses)
- sets
    - thoroughly test consensus for compounds/formulas and ranking with both set and algo consensus
- new ranking for formulas/compounds consensus (and sets)
    - filterSets?
- ensure peaklists are sorted
- features
    - new multiple blank filtering
    - syncing of XCMS/KPIC2 objects
    - check if featindex and groups slots are in sync with features
    - subsetting and groupScores
    - peak qualities
        - filter (once implemented)
    - as.data.table: normalization?



## docs

- ref docs
    - delete()
        - mention j=DT for fGroups method?
    - parallelize args
        - feat opt: update that it doesn't support parallelization atm --> now it does, but via futures
    - finish featureAnnotations class
- handbook
    - KPIC2, cliqueMS installation
    - update MP section on what uses it (e.g. SIRIUS features, SAFD, ...)
- ffOpenMS etc: analyses also needs to be available for hashing
- update/add aliases
- remove some aliases for generics now not unique to one class (eg unique/overlap)
- check UIs and import functions
- features
    - selectIons: chargeMismatch --> note that OpenMS findFeatures removes isotopes, hence, adducts more reliable
    - topMostByRGroup: handbook?
    - check features
        - session filter: argument and its order
        - MC import/export
        - session filters: TRUE will use default yml file name
    - reportHTML: EICs made if annotations, even if not specified in reportPlots
    - progressr
- suspect screening
    - explain three mass matching methods (see comments doScreenSuspects())
    - mention mz column can now be NA
    - form/compRanks: not updated when subsetting on sets
    - logPath for annotateSuspects()
- annotation
    - minMSMSPeaks and annotatedBy MSPeakLists filters
    - clearly mention (refs, handbook) that MSPeakLists should not be filtered/subset after annotation
    - absAlignMzDev
    - compounds vs formulas formats?
        - formulas: each line is the best candidate from an analysis/set
    - set ranking: same as compounds consensus (but no weights (yet))
    - as.data.table(formulas, average=T): now removes most cols
    - formulasDA: update that MSPeakLists now is necessary (obtained via FFM PLs)
    - new/changed as.data.table args
        - maxFormulas/maxFragFormulas removed (handbook)
    - generateFormulas: MSPeakLists argument
        - also update handbook, tutorial...
    - as.data.table(): countFragElements implies fragments=TRUE
    - MSPL
        - reAverage for subset/filter
        - ID column, used by fragInfo
- sets
    - setObjects() can be used for specific slots such as algo objects and MF settings
    - filter() for features/fGroups: apply to neutral masses
    - CAMERA/RAMClustR/nontarget components: clearly mention it is simply a merge between sets
    - intclust/specclust/TPs is not a componentsSet
    - XCMS(3) grouping: exportedData/rtalign/retcorArgs not supported
    - find nice way to re-use docs
    - mention that setObjects are _not_ filtered by setThreshold for formulas/compounds
    - mention that new consensus for formulas/compounds is made after filter() and addFormulaScoring()
        - this could mean really different results if subsetting on sets is done prior to filtering
        - put message() in sync function?
    - mention that set coverage/formula feature coverages do not consider sets/analyses without any results
        - put message() in sync function?
        - although subsetting doesn't sync anymore --> clearly document all this!
    - document for every object how consensus/merge is done
    - update/check version nr mentioned in filter() for MSPeakLists
    - explain xlim/ylim behavior for annotations/mols for plotSpec()
    - SIRIUS: batch calculations done per adduct
    - clearly document that selectIons() for sets re-group --> new group names! (eg workflow objects now incompatible)
    - update/check version nr mentioned in filter() for MSPeakLists
    - explain xlim/ylim behavior for annotations/mols for plotSpec()
    - update for featThreshold and featThresholdAnn
        - put featThreshold in common arguments handbook table?
    - overlap/unique/plotVenn for sets: sets arg overrides `which`
    - adducts<-
        - for sets: make clear that order/names are taken from annTab[set == s]
        - for sets: example with reGroup
    - as.data.table() for formulas: formula column removed if average=T
    - unset(fGroups) will get adduct annotated fGroups
    - annotation
        - subsetting on sets with updateConsensus=FALSE: set coverage columns are not updated
        - @export multi methods?
        - algo consensus: done on setObjects
    - suspects
        - different rank columns
        - estIDLevel: best case, sublevels stripped if not the same
    - features
        - groupFeatures: align doesn't work for sets with KPIC/XCMS
        - xcms-conv
            - set/... arg
- components
    - OpenMS: qTry == "feature" currently not supported
    - OpenMS: adduct specification: molMult must be one, multiple additions (eg Na2) is controlled by chargeMin/max
    - OpenMS/cliqueMS adducts?
    - checkComponents/delete: for NT/TP components all fGroups are removed
    - update docs for reduced removal: only clust affected, remove note about reduced
    - extraOpts --> `...` for nontarget
- TPs
    - mention Bas as author for log2fc, spec similarity/shift etc
        - cosine based on OrgMassSpecR
    - predictTPsBioTransformer
        - use identifier as fallback for naming when no compoundName is present
        - citations, also EnviPath
        - compound similarities
    - FCParams/plotVolcano
    - update spectrumSimilarity()
    - document that relative intensity and min peaks filter for spec sim is applied after removing precursors
        - min peaks always applied lastly
    - mention that retDirMatch filter ignores any zero values for TP_retDir/retDir (ie to be safe)
    - doc that merging TPs (same fGroup/TP) could be done with suspect screening
    - doc somewhere plotInt order with sets
    - refs for PC transformations
    - generateComponentsTPs()
        - if pred is NULL using filter() is crucial
        - shift ignored from specSimParams
    

## NEWS

- progressr
- newProject
    - DA mslists: MSMSType now correctly bbCID if DIA is requested
        - added DIA checkbox because of this
    - added some default findFeatureOpenMS args
- Features
    - as.data.table(fGroups): normalization, FC, averageFunc
    - topMostByRGroup/EICTopMostByRGroup
    - as.data.table: qualities argument (and potentially faster now with features=T?)
    - optimized feature group filters
    - reportHTML: EICs shared amongst EIC and annotation tab, annotation EIC now scaled
    - ... for findFeaturesXCMS3
    - XCMS3 grouping/import with exportedData and comparison() supports xcms3
    - don't subtract blanks from each other
    - syncing XCMS objects
    - print feature counts in show(fGroups) and filter()
    - noDataPlot() for empty plots, eg by plot(), plotChroms()...
    - mzWindow --> mzExpWindow
    - groupFeatures: feat arg --> obj
    - XMCS(3): `exportedData` --> `loadRawData`
    - OpenMS: minFWHM/maxFWHM defaults lowered for findFeatures and feat opt
    - clarify reportCSV() now only reports remaining features?
- Annotation
    - MSPL
        - reAverage = FALSE for subset/filter and implications of filtering when setting it to TRUE
        - PLIndex -- > PLID and different meaning
    - Fixed Hill ordering: H wasn't alphabetical if no C is present
    - minMSMSPeaks and annotatedBy MSPeakLists filters
    - fixed: withMSMS now applied after all other filters
    - fixed: topX peaks for MSPeakLists filter would re-order peaklists
    - MetFrag
        - formula_MF in MF FragInfo
        - useSmiles=true
    - SIRIUS
        - Fixed: as.data.table() and possibly others didn't handle empty fragInfos from SIRIUS
        - More fixes to correctly handle 'adduct fragments'
    - Fixed: as.data.table(compounds, fragments=TRUE) not outputting anything with missing fragInfo
    - changed formula feature consensus
        - featThreshold and featThresholdAnn
        - former follows original definition but is now actually working properly, as pruning before could lead too high abundances
        - latter only takes annotated features into account
        - defaults changed: featThreshold=0, replaced by featThresholdAnn=0.75
        - renamed analysis column from feature consensus formulae results to analysis_from
        - added analyses column in feature consensus formulae results
        - formulas as.data.table(average=TRUE): remove analysis_from column
    - annotatedPeakList(): also add annotation columns for missing results (for consistency)
    - Fixed: `plotSpectrum()` if `xlim` is set and this yields no data then an empty plot is shown
    - Fixed: `plotSpectrum()` automatic `ylim` determination was incorrect if only one peak is shown
    - Fixed: consensus from feature formulas possibly could have fragment m/zs not in group MS/MS peaklists
    - Fixed: consensus from feature formulas possibly could have fragment m/zs that deviated from those in in group MS/MS peaklists
    - Fixed: formula algorithm consensus wrongly ranked candidates not ubiquitously present in all algorithms
    - formula/compound annotation consensus: ranking is properly scaled
    - Fixed: the `scoreLimits` filter for formulas could ignore results not obtained with MS/MS data
    - Bruker formulas: MSPeakLists argument now required
    - frag mSigma/score now also averaged
    - formula --> ion_formula/neutral_formula (formulas/compounds)
    - intensity removed from fragInfo
    - generateFormulas now includes mandatory MSPeakLists argument
    - as.data.table(): maxFormulas/maxFragFormulas removed
    - elements filter: neutral_formula used for formulas (but not fragments)
    - Fixed: MF was using wrong/inconsistent cache name
    - compounds consensus: removed minMaxNormalization param (wasn't used anyway, old left-over?)
    - formula_mz --> ion_formula_mz
    - MF now uses peaklist precursor mz instead of fGroup mz
- adducts
    - GenForm/MetFragAdducts()
        - now report generic format
        - load from cached R data --> faster as.adduct()
    - adduct as.character: err argument
    - fix: conversion of adducts with multiple additions/subtractions to GenForm/MetFrag format failed
    - Standardized GenForm/MetFrag addition/subtraction columns from their adduct tables to fix some conversions (eg NH4 --> H4N)
    - OpenMS and cliqueMS formats
- suspects
    - susp_ prefix in as.data.table
    - suspFormRank/suspCompRank --> formRank/compRank
    - Updated TASQ support
    - numericIDLevel: handle NAs
    - logPath for annotateSuspects()
- components
    - RC components: ensure that columns are the right type if all values are NA
    - changed "rt" to "ret" for component columns for consistency
    - show(): show unique fGroup counts
    - CAMERA components: handle cases when `minSize` filter results in zero components
    - CAMERA components: ensure consistent component names
    - filter(): allow negative rtIncrement values
    - componentsReduced removed, intclust checks instead
    - extraOpts --> `...` for nontarget
- plotGraph: better error if object is empty
- Fixed: cache parallelization issues (thanks to https://blog.r-hub.io/2021/03/13/rsqlite-parallel/)
- xnames/showLegend args for plotInt
- newProject: switch to new system and tweaks
- as.data.table(formulas, average=T): now removes most cols


# Lower priority

## General

- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
- allow specifying average function in other places where as.data.table() is used (eg clustering, plotting etc)
- delete() for other classes
- generalize makeLegend() to new plot util

## Features

- import XCMS features: verify anaInfo (or remove necessity)


## Annotation

- SusDat MF support
- parallel MSPeakLists generation?
- somehow handle different fragment formula annotations when making a consensus between formula/compounds objects
- DA formulas: also rank formula results like GF/SIRIUS?


## Components

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


## TPs

- filter on stability/persistence/toxicity of TP?


## Reporting

- add more options to reportPlots argument of reportHTML()?
- onlyAnnotated argument

