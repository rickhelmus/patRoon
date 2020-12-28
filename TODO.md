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
- comparison(): support xcms3? (needs missing support for missing raw data)
- Fix: blank filter with multiple replicate groups (and maybe others?)
- Check: units of plotChord() rt/mz graphs seems off
- update docs and checkChromatograms() for mzWindow --> mzExpWindow
- also subset features if fGroups are?
- peak qualities
    - check for final names and whether to keep them in gInfo
    - include them in as.data.table (also with features=T)
        - (re-)normalization
    - add to EIC tool
- checkFeatures()
    - decide on final function name (perhaps compatible for future component checker)
    - report relevant group/feature info in tables
        - group peak score so it can be used for sorting
        - retention times should be linked to input unit
        - configurable columns?
    - save sessions & make filter
    - mode: group (topMost=1), rGroup (topMostRGroup=1), feature (topMost=NULL, onlyPresent configurable)
    - buttons to reset all groups/features
    - EIC preview dissabled by default? re-rendering table is slow...


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

