# Release

## general
- test negative subset indices
- convertMSFiles()
    - Agilent .d is also a directory?
- consistent src file names, split utils?
- improve default plotting for plotInt and cluster plot functions
- multiproc
    - better error reporting
    - log executed command
- get RAMClustR from CRAN?
- set max processes to threads-1? --> decide after bench


## docs
- improve instructions for MF and SIRIUS installation?
- handbook
    - semi-quant advanced section
    - NT components graph plotting


## features
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- tests to verify getXCMSSet
- suspect screening
    - rename screenTargets?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- fList --> features?
- move rGroups from filter() to subset?


## MSPeakLists
- isotope tagging is lost after averaging
- collapse averagedPeakLists
- test avg params
- metadata() generic?
- change subset arg order?


## compounds
- MetFrag: (buggy) trivial name fetching not needed anymore?
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- improve formula scoring
- do something about negative H explained fragments by MF?


## formulas
- customize/document ranking column order? (only do rank for sirius?)


## components
- RC: check spearmans correlation


## reporting
- add more options to reportPlots argument of reportHTML()?


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
- newProject(): fix multi line delete (when possible)
- more withr wrapping? (dev, par)


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
- add new MF HD scorings and make sure default normalization equals that of MF web


## components
- mass defect components


## Reporting
- report spectra/tables?

