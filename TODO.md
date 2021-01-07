# Release

## general
- test negative subset indices
- refs to OpenBabel?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
    - Support OpenMS vendor conversion? (eg thermo)
- somehow handle different fragment formula annotations when making a consensus between formula/compounds objects


## AutoID

- ID level rules
    - add scorings for SIRIUS/DA
- interface
    - also convert TASQ?
    - newProject()
        - also allow suspect annotation with only peak lists? currently only selectable if formulas/compounds selected
    - annotateSuspects()
        - check why it's is sometimes slow
            - seems to be logging, disable by default? --> only slow with testthat?
    - filter():
        - cache?
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
    - spec similarity:
        - port from TPs someday
- misc
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
            - could also fix column names, replace "-" with NAs etc
        - if yes, mention in ref docs for screenSuspects()
- expand reporting
    - eg include suspect name in EICs
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects
- update docs and handbook
    - mention that components should be done prior to onlyHits=T?


## docs
- improve instructions for MF and SIRIUS installation?
- ref docs and exports for getXCMSnSet
- ffOpenMS etc: analyses also needs to be available for hashing
- improve docs for areas (only affects when features=FALSE) and average (different behavior when features=TRUE/FALSE) for as.data.table() of featureGroups
- update/check version nr mentioned in filter() for MSPeakLists
- explain xlim/ylim behavior for annotations/mols for plotSpec()


## sets
- methods to implement
    - consensus / comparison()
        - could first make consensus of setObjects and then make new set object from that
        - for compounds (and formulas?) need proper way to average scorings
        - or for later?
    - more sub class specific methods
        - compoundsSetMF sub-class (for settings slot)?
    - provide methods for non-implemented functionality
        - consensus()? (see above)
        - compoundViewer?
        - groupFeaturesXCMS3?
- misc
    - as.data.table() for formulas: average=T will now produce strange averaged ionized formula, for now simply remove this column...
        - also give a note in docs?
        - or maybe only remove is not all adducts are equal?
    - handle errors when object has <=1 set
        - groupFeaturesScreening()
        - mergeScreeningSetInfos()
    - fix empty MS(MS) peaklists if unavailable during merging sets
        - already fixed?
- merging setObjects
    - check if more has to be cached and may need status messages
    - filter features, annotation results on minimum abundance amongst different setObjects
        - (unlike setThreshold also take objects without results in to account)
        - add new filter()?
    - compound set consensus: scoreRanges should be re-determined from annotation results?
- setThreshold
    - filter() argument
    - keep/change default setThreshold?
        - or remove argument from generators?
- suspect screening
    - implement TASQ?
    - consensus?
        - or just new set filter described above?
    - annotation columns not in report, fine? (there are many) If yes document
    - as.data.table(fGroupsScrSets, collapseSuspects=NULL): omits sets column
        - still true? check
- neutralizing / ionization
    - mergeIons()
        - other name?
        - makes sense to not choose monoisotopic mass? not for annotation at least
		- also method for fGroupsSets?
			- do it per component set
			- used to update adducts
			- needs to re-group afterwards
			    - make general reGroup() method, that uses stored settings from groupFeatures()?
	    - prefer adducts based on MS/MS? eg handy for Na/K adducts
	- makeSet()
	    - fGroups: also support method via comparison?
	- formula/compounds: get adduct from gInfo if present
		- also for screening? could (optionally) look for matches with neutralized fGroup/susp masses
		- fragment annotation with MF OK?
		    - Also check GF
		    - report 'neutral' fragments? Check what SIRIUS does.
		- what to do with unsupported adducts?
		    - skip calculation with a warning?
		        - default to M+H/M-H for now with warning...
		    - enforce adduct argument for these cases?
		    - default mergeIons() to only consider 'common' adducts?
		    - check better for what is supported by SIRIUS?
	- suspect screening
	    - update for old adducts() calls
	    - update interface
	        - split screening for suspects with ionized mass (mz) column and the rest?
	            - mz: match suspects as usual
	            - rest: match by neutral mass
	                - get adduct from arg --> susp list --> annotations?
                - make optional to always match by ionized mass?
	- add annotations to as.data.table()
	- cliqueMS components
	- component selection tool/function
	    - otherwise perhaps make a fGroup remover function to help subsetting
			- similarly as for feature remover in fGroups...
			- del()/delete()/rmResult()/delResult() generic? could also be for other classes
			    - or use regular R set like methods (<-())
		- similarly: set() like method to change data, such as adduct annotations
- NEWS
    - [..., reAverage = FALSE] and implications of filtering when setting it to TRUE
    - Fixed Hill ordering: H wasn't alphabetical if no C is present
    - MetFrag
        - formula_MF in MF FragInfo
        - useSmiles=true
    - SIRIUS
        - Fixed: as.data.table() and possibly others didn't handle empty fragInfos from SIRIUS
        - More fixes to correctly handle 'adduct fragments'
    - Fixed: as.data.table(compounds, fragments=TRUE) not output anything with missing fragInfo
    - adducts
        - GenForm/MetFragAdducts()
            - now report generic format
            - load from cached R data --> faster as.adduct()
        - adduct as.character: err argument
        - fix: conversion of adducts with multiple additions/subtractions to GenForm/MetFrag format failed
        - Standardized GenForm/MetFrag addition/subtraction columns from their adduct tables to fix some conversions (eg NH4 --> H4N)
- docs
    - filter() for features/fGroups: apply to neutral masses
    - CAMERA/RAMClustR/nontarget components: clearly mention it is simply a merge between sets
    - intclust is not a componentsSet
    - find nice way to re-use docs
    - mention that setObjects are _not_ filtered by setThreshold for formulas/compounds
    - mention that new consensus for formulas/compounds is made after filter() and addFormulaScoring()
        - this could mean really different results if subsetting on sets is done prior to filtering
        - put message() in sync function?
    - mention that set coverage/formula feature coverages do not consider sets/analyses without any results
        - put message() in sync function?
    - document for every object how consensus/merge is done
    - improve docs for areas (only affects when features=FALSE) and average (different behavior when features=TRUE/FALSE) for as.data.table() of featureGroups
    - update/check version nr mentioned in filter() for MSPeakLists
    - explain xlim/ylim behavior for annotations/mols for plotSpec()
    - update/add aliases
    - SIRIUS: batch calculations done per adduct
    - suspect screening
        - explain three mass matching methods (see comments doScreenSuspects())
        - mention mz column can now be NA
- tests
    - handle/test empty objects
    - test DA algorithms
    


## features
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- suspect screening
    - rename patRoonData::targets?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- importFeaturesXCMS/importFeaturesXCMS3/importFeatureGroupsXCMS: get rid of anaInfo arg requirement? (or make import func?)
- Check: units of plotChord() rt/mz graphs seems off
- checkFeatures()
    - filter negate: don't remove kept fGroups if features were changed?
    - update filter for delete()
- misc
    - topMostByRGroup: make default? or only for reporting?
    - quality/score filters
    - updatePICSet(): also sync peaks list? otherwise doc
    - update docs and checkChromatograms() for mzWindow --> mzExpWindow
    - newProject(): new algorithms
    - getEICsForFeatures method for kpic2?
    - optimize hashing? Or further avoid hashing big objects/objects with lists?
    - rename exportedData, update docs (also handbook)
    - remove featuresOpenMS method for getXCMSSet(), update docs
    - groupFeatures(algorithm="sirius")
    - don't key quality/score tables?
    - print feature counts in show(fGroups) and filter()
    - delete for features and fGroups
        - XCMS: also update groups data?
- NEWS
    - topMostByRGroup/EICTopMostByRGroup
    - as.data.table: qualities argument (and potentially faster now with features=T?)
    - optimized feature group filters
    - reportHTML: EICs shared amongst EIC and annotation tab, annotation EIC now scaled
    - ... for findFeaturesXCMS3
    - XCMS3 grouping/import with exportedData and comparison() supports xcms3
    - progressr
    - don't subtract blanks from each other
    - syncing XCMS objects
- tests
    - topMostByRGroup
    - xcms3 comparison
    - new multiple blank filtering
    - syncing of XCMS/KPIC2 objects
    - delete()
- docs
    - clarify reportCSV() now only reports remaining features?
    - topMostByRGroup: handbook?
    - session filter: argument and its order
    - importCheckFeaturesSession()
    - groupQualities/Scores slots
    - GaussianSimilarity: NAs are made zero
    - as.data.table: qualities argument
    - calculatePeakQualities(), also generic
    - reportHTML: EICs made if annotations, even if not specified in reportPlots
    - handbook
        - checkFeatures()
    - update IPO docs for kpic2 (and mention min/max_width split and others)
    - ... for findFeaturesXCMS3
    - progressr

## MSPeakLists
- isotope tagging is lost after averaging
- collapse averagedPeakLists
- test avg params
- metadata() generic?


## compounds
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- do something about negative H explained fragments by MF?
- SusDat MF support
- MetFrag: auto-include suspect results if suspectListScore is selected?
- SIRIUS: message/progress bar when processing a large batch in singular mode?


## formulas
- customize/document ranking column order? (only do rank for sirius?)
- getFormInfoList(): take care of consensus results like getPrecursorFormScores()

## components
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?


## reporting
- add more options to reportPlots argument of reportHTML()?


## Cleanup
- Reduce non-exported class only methods

## MP

- future MP
    - delayBetweenProc?


# Future

## General

- msPurity integration
- suspect screening: add MS/MS qualifiers
- fillPeaks for CAMERA (and RAMClustR?)
- support fastcluster for compounds clustering/int component clusters?
- algorithmObject() generic: for xset, xsa, rc, ...
- newProject(): fix multi line delete (when possible)
- more withr wrapping? (dev, par)
- improve default plotting for plotInt and cluster plot functions
- newProject()
    - concentration column for anaInfo
    - generate more detailed script with e.g. commented examples of subsetting, extraction etc
- support more of the new SIRIUS functionality
	- newProject(): import Bruker seq file?


## Features

- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- parallel enviPick
- OpenMS MetaboliteAdductDecharger support?
- OpenMS: Support KD grouper?
- suspect screening: tag fGroups with suspect instead of converting fGroups object (and add filter to remove non-hits)
- Integration of mzMine features (package pending...), MS-DIAL and KPIC2, peakonly, SIRIUS?


## MSPeakLists

- DA
    - generateMSPeakListsDA: find precursor masses with larger window
    - tests
        - utils? EICs with export/vdiffr?
        - test MS peak lists deisotoping?
- metadata for Bruker peaklists?


## Formulas

- DBE calculation for SIRIUS?
- OM reporting
- as.data.table: option to average per replicate group?


## Compounds

- do something with sirius fingerprints? --> comparison?
- fix compoundViewer
- add new MF HD scorings and make sure default normalization equals that of MF web
- CFM-ID and MS-FINDER integration
- utility functions to make custom DBs for MetFrag and SIRIUS and support to use them with the latter


## components
- mass defect components
- CliqueMS
- split peak correlation and adduct etc annotation? would allow better non-target integration
- intclust
    - optionally take areas instead of intensities
    - cache results

