# Release

## general
- test negative subset indices
- convertMSFiles()
    - Agilent .d is also a directory?
- runWithoutCache? runWithCacheMode()? shortcut to withr::with_options(patRoon.cache.mode=...)


## docs
- improve instructions for MF and SIRIUS installation?


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
- NT: minimum size argument, combine rows for multiple rGroups?


## reporting
- add more options to reportPlots argument of reportHTML()?


## Cleanup
- Reduce non-exported class only methods

## MP
- new interface
    - batch mode
    - logging --> return stdout/stderr from run()
    - printOutput, eg for SIRIUS... or not?
    - showProgress
    - waitTimeout/delayBetweenProc?
- changes
    - logging done in executeMultiProc()
    - caching done in executeMultiProc()
    - configure "engine" through options() --> "classic" and "future"
        - let executeMultiProc() handle general things like caching and logging
        - rename "classic"?
    - somehow let user set future.apply scheduling and chunk.size options, probably through options()
    - logging:
        - enabled through an option
        - implement in engines
        - somehow disable for featuresOptimizerOpenMS
- default GenForm to classic?
- docs
    - write vignette or handbook chapter
    - update sirius-args.R (mention maxProcAmount)
    - update handbook etc for maxProcAmount
    - patRoon.logPath, multiProcMethod
    - update ref docs for executeMultiProcess
    - update docs and NEWS for SIRBatchSize --> splitBatches
    - mention that reportHTML never uses future method
- tests
- limitation: local databases need to be available for MF on the localhost and nodes, both for caching and checking score types. Any workarounds? Otherwise document that location for them should be fixed.
    - force that file should be available locally so it can be used for hashing (with makeFileHash?)/scoring
        --> implement hashing
- remove pre44 SIRIUS support: doesn't work with future procs at the moment, and doesn't seem useful to keep
- remove/update PWizBatches?
- move main future handler to separate function?


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

