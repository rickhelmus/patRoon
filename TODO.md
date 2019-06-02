# Release

## general
- more withr wrapping? (dev, par)
- test negative subset indices
- ref --> blank
- convertMSFiles()
    - Agilent .d is also a directory?
- consistent src file names, split utils?
- newProject(): don't necessarily add MSPeakLists when DA formulas are choosen
- improve default plotting for plotInt and cluster plot functions
- better error reporting for multiproc commands
- get RAMClustR from CRAN
- ref -> blank


## docs
- improve instructions for MF and SIRIUS installation?
- handbook
    - semi-quant advanced section


## features
- getXcmsSet --> export?
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- tests to verify getXCMSSet...
- suspect screening
    - rename screenTargets?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- fList --> features?
- move rGroups from filter() to subset?
- as.data.table() and others: optionally report areas instead of intensities


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
- plotHeatMap() / plotSilhouettes() for comp clusters?
- do something about negative H explained fragments by MF?


## formulas
- customize/document ranking column order? (only do rank for sirius?)
- minExplainedFragPeaks -> minExplainedPeaks?


## components
- RC: check spearmans correlation


## reporting
- add more options to reportPlots argument of reportMD()?
- reportMD() --> reportHTML()


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
- add new MF HD scorings and make sure default normalization equals that of MF web


## components
- mass defect components


## Reporting
- report spectra/tables?

