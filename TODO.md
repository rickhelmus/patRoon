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
    - GenForm scoring: somehow exclude non MS/MS candidates if MS/MS candidates are present?
        - already fine now with minimum rank?
    - add scorings for SIRIUS/DA
    - remove mustExist?
    - convert to YAML
        - update tests
        - refs?
- interface
    - groupFeaturesScreening() --> screenSuspects()?
    - deprecate old screen functions
    - also convert TASQ?
    - newProject(): update for new interface
    - annotateSuspects() --> annotate() latter is a function (but not generic) from ggplot2 and RAMClustR and method from CAMERA, so probably no
    conflicts
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
    - spec similarity:
        - port from TPs someday
        - clarify that intensity filter happens after filtering precursor (make this optional?)
        - add proper refs
    - update filter() and reportHTML() for new annMSMS cols
- misc
    - util to check if there are suspect results? (ie to replace inherits(...))
    - check for empty names in assertion/preparation functions
    - credits to ES
    - update version number
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
- expand reporting
    - eg marking which candidate corresponds to suspect and include suspect name in EICs
        - mark with different row colour and label?
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects
    - adducts from components?
    - update reportCSV() / reportPDF()?
- update docs & handbook
    - renamed rt/mz columns
    - new plotVenn list functionality
- tests
    - automatic InChIKey/formula calculation from InChIs/SMILES
        - already done implicitly?
    - handling empty results
    - filters
    - more?
    - new plotVenn list functionality


## docs
- improve instructions for MF and SIRIUS installation?
- ref docs and exports for getXCMSnSet


## features
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- suspect screening
    - rename patRoonData::targets?
    - rename groupFeaturesScreening?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- importFeaturesXCMS/importFeaturesXCMS3/importFeatureGroupsXCMS: get rid of anaInfo arg requirement? (or make import func?)
- comparison(): support xcms3? (needs missing support for missing raw data)
- Fix: blank filter with multiple replicate groups (and maybe others?)
- Check: units of plotChord() rt/mz graphs seems off


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

