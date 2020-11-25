# Release

## general
- test negative subset indices
- refs to OpenBabel?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
    - Support OpenMS vendor conversion? (eg thermo)

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
- provide methods for non-implemented functionality (eg consensus)
- find nice way to re-use docs
- filter() for features/fGroups: support ionized masses for mass filters? or just clarify it doesn't.
- handle/test empty objects
- more descriptive messages what's going on with all the avaraging of MSPeakLists
- remove sets argument for some methods (as.data.table, accessors etc)?
    - if keep, be consistent with all classes
- as.data.table() for formulas: average=T will now produce strange averaged ionized formula, for now simply remove this column.. also give a note in docs? or maybe only remove is not all adducts are equal?
- different name/generic for ionize()? makes less sense for annotation classes
- test DA algorithms
- check if more has to be cached (eg merged results from sets)
- compoundsSetMF sub-class (for settings slot)? or is access via setObjects sufficient? may need to explain anyway for other cases like intclust components
- base set class
    - don't use for fGroupsSet?
- components
    - neutralize masses? otherwise document
        - yay: consistent with other set classes
        - nay: might be a bit strange when looking for adducts etc and components are per set anyway
    - intclust: return componentsSet? if not document somewhere...
    - clearly mention that nontarget is done per set now
    - nontarget-set: plotGraph method? and make sure it's used in reportHTML()
- implement XCMS conversion functions? maybe with given set. Could just ionize() it.
- ionize() for compounds/components? if not remove formulas?
- setThreshold filter() argument, and remove argument from generators?
- handle errors when object has <=1 set
    - groupFeaturesScreening()
    - mergeScreeningSetInfos()
- suspect screening
    - handle suspect and fragment mz values somehow
    - implement TASQ?
    - consensus?
    - support recursive screening? or throw error otherwise
- fix empty MS(MS) peaklists if unavailable during merging sets
- use neutral_formula for annotated similarity calculation for formulas
- NEWS
    - [..., reAverage = FALSE] and implications of filtering when setting it to TRUE

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


## formulas
- customize/document ranking column order? (only do rank for sirius?)
- getFormInfoList(): take care of consensus results like getPrecursorFormScores()

## components
- RC: check spearmans correlation
- NT
    - minimum size argument, combine rows for multiple rGroups?
    - update rt range reported after merging series


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

