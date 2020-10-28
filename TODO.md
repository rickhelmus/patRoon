# Release

## general
- test negative subset indices
- refs to OpenBabel?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)

## AutoID

- update version number
- update docs & handbook
    - renamed rt/m columns
- credits to ES
- tests
    - automatic InChIKey/formula calculation from InChIs/SMILES
        - already done implicitly?
    - handling empty results
    - more?
- ID level rules
    - GenForm scoring: somehow exclude non MS/MS candidates if MS/MS candidates are present?
        - already fine now with minimum rank?
        - absolute thresholds for iso/comb scores?
    - add scorings for SIRIUS/DA
- annotation
    - docs
    - only add sim, ranking etc columns if data is available?
    - cache
- screenSuspects()
    - combine screenSuspects() and groupFeaturesScreening()
        - screenSuspects() does both the screening and making new fGroups
        - as.data.table() for fGroupsScreening which adds suspect metadata and ID levels
            - optionally un-collapses suspects per row
        - annotate() method which does what annotateSuspects() now does
            - can later also be defined for fGroups?
            - adds data to screenInfo slot
        - filter() method for minimal rankings, matched MS/MS fragments, no/one/multiple hits etc
            - onlyHits = TRUE by default? or by default if any of the others are given?
        - deprecate features method
        - tag hits, but keep all fGroups?
            - fGroups[, hits = TRUE] ?
            - still keep re-naming functionality? not really possible, but put names in reporting functionality
        - same format for TASQ screening import?
        - throw informative error when screenSuspects()/screenTargets() is called
        - also selectBy for fGroups with multiple hits
- newProject(): create template auto ID rule csv?
- util to check if there are suspect results? (ie to replace inherits(...))
- expand reporting
    - eg marking which candidate corresponds to suspect and include suspect name in EICs
        - mark with different row colour and label?
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects
- check for empty names in assertion/preparation functions
- prepareSuspectList(): export? and always calculate neutralMass?
- annotateSuspects() --> annotate() latter is a function (but not generic) from ggplot2 and RAMClustR and method from CAMERA, so probably no conflicts


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

