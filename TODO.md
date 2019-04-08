# Release

## general
- more withr wrapping? (dev, par)
- test negative subset indices
- ref --> blank
- convertMSFiles()
    - Agilent .d is also a directory?


## docs
- More vignettes
- Reference docs
    - Examples?
- update tutorial
- improve instructions for MF and SIRIUS installation?


## features
- getXcmsSet --> export?
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- tests to verify getXCMSSet...
- rename screenTargets?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- fList --> features?


## MSPeakLists
- isotope tagging is lost after averaging
- most intense analysis as alternative to averaging?
- collapse averagedPeakLists
- test avg params
- metadata() generic?


## compounds
- MetFrag: (buggy) trivial name fetching not needed anymore?
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- add new MF stat and other scorings and make sure default normalization equals that of MF web
- MF: add support for all databases
- Fix MF database section in docs
- timeouts for SIRIUS?
- improve formula scoring


## formulas
- set column order?
- customize/document ranking column order? (only do rank for sirius?)
- how to handle ranking of consensus results?
    - rank by one normalized column per algo? (GenForm: either MS_match or comb_match)
    - or simply don't and mention limitation in doc? (one limitation: filtering)
- reporting incorrect GF adducts


## components
- RC: check spearmans correlation


## reporting
- add more options to reportPlots argument of reportMD()?


## Cleanup
- Reduce non-exported class only methods


# Future

## General

- msPurity integration
- suspect screening: add MS/MS qualifiers, calculate ion masses
- newProject(): generate Rmd?
- fillPeaks for CAMERA (and RAMClustR?)
- support fastcluster for compounds clustering/int component clusters?
- algorithmObject() generic: for xset, xsa, rc, ...

## Features

- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- parallel enviPick
- OpenMS MetaboliteAdductDecharger support?
- OpenMS: Support KD grouper?


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

- do something with sirius fingerprints?
- fix compoundViewer


## components
- mass defect components


## Reporting
- report spectra/tables?

