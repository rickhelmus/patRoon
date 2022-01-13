# patRoon 2.0.1

- Fixed: `newProject()` would not add suspect annotation to the output script if the example suspect list or sets workflows were chosen.
- Fixed: Default optimization range for KPIC2 `min_width` was incorrect (PR #31, thanks to @@coltonlloyd)
- `installPatRoon()` improvements in determining what is already installed
- Fixed: Group qualities/scores were not transferred to new `featureGroups` objects after calling `screenSuspects()` or `unset()`
- Checking of MS file extensions (e.g. for `generateAnalysisInfo()`) is now case insensitive (see issue #34)


# patRoon 2.0

This release adds a significant amount of new functionality and changes. Please see the updated Handbook and sections below for more information.

Users of previous `patRoon` versions should inform themselves with the important changes highlighted in the next section. Furthermore, it is highly recommended to remove any cached data, i.e. by running `clearCache("all")` or manually removing the `cache.sqlite` file from your project directory.


## Important changes

- Features
    - XMCS(3): Renamed argument `exportedData` to `loadRawData`.
    - The `mzWindow` and `EICMzWindow` arguments were renamed to `mzExpWindow` / `EICMzExpWindow` and are now with slightly different meaning (please see the reference manual).
    - OpenMS: `minFWHM`/`maxFWHM` defaults lowered for `findFeatures` and feature optimization.
- Annotation
    - ggplot2 support for several plotting functions (i.e. `useGGPlot2` argument) is removed (not often used and a maintenance burden).
    - the `precursor` argument to the `plotSpectrum()`, `annotatedSpectrum()` and `plotScores()` methods for `formulas` now expects the neutral formula instead of the ionized formula. This change was done for consistency with compound annotations and sets workflows.
    - The way of obtaining candidate formulae from different features in the same group (i.e. _feature formula_ consensus) was changed.
        - Fixes were applied to improve thresholding with `featThreshold`.
        - A second and new threshold, `featThresholdAnn`, only takes annotated features into account.
        - The default of `featThreshold` is now `0`, for `featThresholdAnn` it is the same as the previous default for `featThreshold`.
        - Candidate results: renamed `analysis` column to `analysis_from` and added `analyses` column that lists all analyses from which the consensus was made.
        - if multiple annotations are available for a single MS/MS peak (eg due to differences between feature annotations) then only the annotation with lowest m/z deviation is kept (and a warning is emitted).
        - Scores of annotated fragments from different features are now averaged.
    - The storage classes and interface for formula and compound annotation was harmonized
        - The `formulas` and `compounds` classes now derive from the new `featureAnnotations` class. Most of the functionality common to formulas/compounds are defined for this class.
        - Storage of formula annotation results mostly follow the format that was already used for compound annotations.
        - The `maxFormulas`/`maxFragFormulas` argument for `as.data.table()` were removed, as these don't make much sense with the new format.
        - The `elements` filter now applies to neutral formulae for both formula and compound annotations (`fragElements` still applies to the ionized fragment formula).
    - Formula annotation with Bruker now require `MSPeakLists`. Since all algorithms now require peak lists, `generateFormulas` now has a mandatory `MSPeakLists` argument (similar to `generateCompounds`).
    - Formula candidates (formula and compound annotations) are now reported in the `ion_formula` (ionized) and `neutral_formula` (neutral) columns. Similarly, the `formula_mz` column was renamed to `ion_formula_mz`.
- Suspect screening
    - The methodology to match m/z values of suspects and features was changed. This was mainly for consistency and compatibility with sets workflows. Please see the updated section on suspect screening in the Handbook.


## Major new functionality

### Transformation product screening

The most important new functionality in `patRoon 2.0` are transformation product (TP) screening workflows. This release adds functionality to predict TPs (with `BioTransformer` or metabolic logic) or search TPs in `PubChem` or custom databases. Furthermore,  other data such as MS/MS similarity or feature classification data can be used to relate parent/TP features. Other TP screening functionality includes TP prioritization and automatic generation of TP compound database for `MetFrag` annotation. The workflows follow the classical design of `patRoon`, where flexible workflows can be executed with a combination of established algorithms and new functionality. For more information, please see the dedicated chapter about TP screening in the Handbook.

### Sets workflows

Another major change in this release is the addition of sets workflows. These workflows are typically used to simultaneously process positive and negative ionization data. Advantages of sets workflows include simplification of data processing, combining positive and negative data to improve e.g. feature annotations and easily comparing features across polarities. A sets workflow requires minimal changes to a 'classical workflow', and most of the additional work needed to process both polarities is done automatically behind the scenes. For more information, please see the dedicated chapter about sets workflows in the Handbook.

### Features

The following new feature detection/grouping algorithms were integrated: `SIRIUS`, `KPIC2` and `SAFD`. Furthermore, integration with `MetaClean` was implemented for the calculation of peak qualities and machine learning based classification of pass/fail peaks. In addition, the peak qualities are used to calculate peak scores, which can be used for quick assessment and prioritization.

### Data curation

Interactive curation of feature data with `checkChromatograms()` was replaced with `checkFeatures()`, which is much faster, is better suitable for larger datasets, customizable and has an improved user interface. Furthermore, this tool can be used for training/assessing `MetaClean` models. Similarly, `checkComponents()` is a function that allows interactive curation of component data.

The `delete()` generic function allows straightforward deletion of parts of the workflow data, such as features, components and annotations. Furthermore, this function makes it easy to implement customized filters.

### Adduct annotation

The algorithms of `OpenMS` (`MetaboliteAdductDecharger`) and `cliqueMS` were integrated for additional ways to detect adducts/isotopes through componentization. Furthermore, the new `selectIons()` method uses these annotations to prioritize features (e.g. by only retaining those with preferable adducts). In addition, this function stores the adduct annotations for the remaining feature groups, which can then be automatically used for e.g. formula and compound annotation.

## Other new functionality

- `newProject`
    - Updated for new functionality such as sets and TP workflows and adduct annotation.
    - Completely re-designed code generation to improve extensibility. The generated code will have a slightly different layout and some parameter defaults were changed.
- Features
    - `as.data.table()`
        - intensity normalization (`normFunc` argument)
        - customized averaging (`averageFunc` argument)
        - calculation of Fold-changes (`FCParams` argument)
        - report peak qualities/scores (`qualities` argument)
    - New `plotVolcano()` method function to plot fold changes.
    - `topMostByRGroup/EICTopMostByRGroup` arguments for plotting/reporting EIC data of only the top most intense replicate(s).
    - `reportHTML()` now only plots the EIC of the most intense feature of each replicate group (i.e. `EICTopMostByRGroup=TRUE` and `EICTopMost=1`).
    - `XCMS3`
        - `loadRawData` argument for feature grouping and `comparison()`
        - `...` argument for `findFeaturesXCMS3`
        - `preGroupParam` to specify grouping parameters used prior to RT alignment (suggested by Ricardo Cunha)
    - The internal `XCMS` feature (group) objects are synchronized as much as possible when feature data is changed.
    - Feature groups: print feature counts with `show()` and `filter()` methods.
    - `OpenMS` feature finding: `useFFMIntensities` argument to speed up intensity loading (experimental).
    - `reportHTML()` now reports general feature information in a separate tab.
    - Feature groups: `results` argument to `[` (subset) and `filter()` to quickly synchronize feature groups between objects (e.g. to quickly remove feature groups without annotation results).
- Annotation
    - The methodology of `plotSpectrum()` to automatically calculate the space necessary for formula annotation texts and candidate structures was improved. Annotation texts are now automatically resized if there is insufficient space, and the maximum size and resolution for candidate structures can be controlled with the `maxMolSize`/`molRes` parameters.
    - `filter()` method for `MSPeakLists`: `minMSMSPeaks` filter to only retain MSMS peak lists with a minimum number of peaks.
    - `filter()` method for `MSPeakLists`: `annotatedBy` filter to only keep peaks with formula/compound annotations.
- Suspect screening
    - The `screenSuspects()` method now supports the `amend` argument, which allows combining results of different `screenSuspects()` calls (see the Handbook for details).
- Components
    - A new algorithm, `specclust`, which generates components based on hierarchically clustering feature groups with high MS/MS similarities.


## Minor changes

- Features
    - `groupFeatures`: renamed the `feat` argument to `obj`.
    - Improved performance for some feature group filters.
    - `reportHTML()`: EICs are shared between tabs to avoid duplicated plotting
    - The `features` object embedded in `featureGroups` objects is now synchronized, and any features not present in any group are removed accordingly. This reduces memory usage and indirectly causes `reportCSV()` to only report features still present.
    - `plotInt()`: now has `xnames` and `showLegend` arguments to adjust plotting.
- Annotation
    - The `[` (subset) and `filter()` methods for `MSPeakLists` now only re-average peak lists if the new `reAverage` argument is set to `TRUE` (default `FALSE`). This change was mainly done as (1) the effects are usually minor and (2) re-averaging invalidates any formula/compound annotations done prior to filtering.
    - `filter()` method for `MSPeakLists`: the `withMSMS` filter is now applied after all other filters.
    - `MetFrag`: the raw unprocessed annotation formulas are now additionally stored in the `fragInfo` tables (`formula_MF` column).
    - `MetFrag`: the precursor ion m/z is now taken from peak list data instead of the feature group to improve annotation.
    - `MetFrag`: the `useSmiles` parameter is now set to `true` as this seems to improve results sometimes.
    - `as.data.table()` method for `formulas`: if `average=TRUE` then all column data that cannot be reasonably averaged are excluded.
    - `annotatedPeakList()`: also add annotation columns for missing results (for consistency).
    - Compound MS/MS annotations do not include peak intensities anymore (already stored in MS peak lists).
    - The `minMaxNormalization` argument to the `consensus()` method for `compounds` was removed (unused).
    - `filter()` for `formulas`/`compounds`: if algorithm consensus results are filtered with `scoreLimits`, and a score term exists multiple times for a candidate, only one of the terms needs to fall within the specified limits for the candidate to be kept (was all).
    - `plotSpectrum()` for `compounds`: `plotStruct` is now defaulted to `FALSE`.
    - `MSPeakLists` data now store an unique identifier for each mass peak in the `ID` column. These IDs are used by e.g. formula/compound annotations, and stored in the `PLID` column in their `fragInfo` data. This replaces the `PLIndex` column in `fragInfo` data, which was only row based, and therefore invalidated in case peak lists were filtered afterwards.
- Adducts
    - `GenFormAdducts()` and `MetFragAdducts()` now additionally return adducts in generic format and use cached data for efficiency. 
    - `err` argument to `as.character()` to control if an error or `NA` should returned if conversion fails.
    - `as.adduct()` now removes any whitespace and performs stricter format checks to make conversion more robust.
    - Standardized GenForm/MetFrag element addition/subtraction data to improve consistency for conversions (eg NH4 --> H4N).
    - Conversion from/to adduct formats of OpenMS (`MetaboliteAdductDecharger`) and `cliqueMS`.
    - `calculateIonFormula()` and `calculateNeutralFormula()` now Hill sort their result
    - The embedded GenForm code was updated to the latest version.
- Suspect screening
    - `as.data.table()`: Suspect screening specific columns are now prefixed with `susp_`.
    - The `suspFormRank` and `suspCompRank` suspect annotation data columns were renamed to `formRank`/`suspCompRank` (the previous change made prefixing unnecessary).
    - Several updates for Bruker TASQ import.
    - `logPath` argument for `annotateSuspects()` to specify the file path for log files are disable logging completely.
    - suspect names are now trimmed to 150 characters to avoid logging issues on e.g. Windows
- Components
    - Intensity clusters now use `fastcluster` for hierarchical clustering.
    - Changed column `rt` to `ret` for consistency.
    - `show()`: show unique feature group counts.
    - `filter()`: allow negative `rtIncrement` values.
    - `nontarget`: replaced `extraOpts` argument with `...`.
    - `nontarget`: store links as character string indices instead of numeric indices.
    - `RAMClustR`: moved position of `ionization` argument to improve consistency.
    - The 'reduced components' mechanism, where a components object was returned without any algorithm specific data (using the `componentsReduced` class) when filtering/subsetting components, was removed. This system was quite unintuitive and imposed unnecessary limitations. Instead, functions that cannot work after component data is changed (e.g. those specific to intensity clustering) will throw an error if needed.
    - The objects returned from `intclust` components are now derived from a general `componentsClust` class, which is shared with `specclust` components. The common functionality for both algorithms is implemented for this class.
- Misc
    - `show()` methods now print class inheritance tree
    - The `progressr` package is not used anymore, thus, it is not necessary to set up progress bars with `future` based multiprocessing.
    - `newProject`: Moved order of componentization step (now before annotation & suspect screening).
    - Plots of chromatograms, spectra etc that are without data now reflect this in the generated plot.


## Fixes

- Features
    - Blank filter: don't subtract blanks from each other
    - Fixed: when `xlim`/`ylim` was used with `plotChroms` then peaks were not always correctly filled
    - `retMin` argument to `plot()` method for `featureGroupsComparison` wasn't properly used/defaulted.
- Annotations
    - Fixed: `plotSpectrum()` if `xlim` is set and this yields no data then an empty plot is shown.
    - Fixed: `plotSpectrum()` automatic `ylim` determination was incorrect if only one peak is shown.
    - Fixed: consensus from feature formulas possibly could have fragment m/zs not in group MS/MS peak lists.
    - Fixed: consensus from feature formulas possibly could have fragment m/zs that deviated from those in in group MS/MS peaklists.
    - Fixed: formula algorithm consensus wrongly ranked candidates not ubiquitously present in all algorithms.
    - Fixed: the `scoreLimits` filter for formulas could ignore results not obtained with MS/MS data.
    - Fixed: MetFrag was using a wrong/inconsistent cache name.
    - Fixed: `as.data.table(compounds, fragments=TRUE)` returned empty results for candidates without fragment annotations.
    - Fixed: `topX` arguments for the `MSPeakLists` method for `filter()` would re-order peak lists, thereby invaliding any annotations.
    - Fixed: conversion of adducts with multiple additions/subtractions to GenForm/MetFrag format failed.
    - Fixed: Hill ordering: H wasn't sorted alphabetically if no C is present.
    - Several fixes were applied to improve handling of `SIRIUS` 'adduct fragments'.
    - formula/compound annotation consensus ranking is now properly scaled.
    - Fixed: `generateMSPeakListsDAFMF()` potentially used wrong DA compound data in case features were filtered.
- Suspect screening
    - Fixed: `numericIDLevel()` now properly handles `NA` values.
    - `importFeatureGroupsBrukerTASQ()`: Improved handling of absent analyses in imported results files.
    - Fixed: Automatic _m/z_ calculation for suspects
        - Improperly handled electron masses for adducts involving element subtract (e.g. `[M-H]-`), resulting in ~1.5 mDa deviations
        - Adduct conversion didn't handle multiple molecules (e.g. `[2M+H]+`) and multiple charges (e.g. `[M+2H]2+`)
- Components
    - `RAMClustR`: ensure that columns are the right type if all values are NA.
    - `CAMERA`: correctly handle cases when `minSize` filter results in zero components.
    - `plotGraph()`: improve error handling with empty objects.
- Misc
    - Future multiprocessing: make sure that logs are created even when an error occurs.
    - Classic multiprocessing: intermediate results are cached again.
    - Fixed: parallelization issues with cached data (thanks to https://blog.r-hub.io/2021/03/13/rsqlite-parallel/)
    - `newProject()`: correctly handle DIA with Bruker MS peak lists.


# patRoon 1.2.1

- Fixed: XCMS feature grouping didn't work when the `peakgroups` alignment method was used (fixes issue #22)
- Fixed: (harmless) `mapply` warning was shown with `newProject()`
- `newProject()`: don't show _Remove_ button in analyses select screen when the script option is selected, as this will not work properly.
- `IPO`: add default limits for OpenMS `traceTermOutliers`
- `IPO` optimization fix: integer parameters are properly rounded
- Fixed: `generateFeatureOptPSet("xcms3", method="matchedFilter")` would return a parameter set with `step` instead of `binSize` (issue #23)
- Fixed: `newProject()` would generate an ID levels configuration file even when no suspect list was selected.
- MCS calculation: handle `NULL` values that may occasionally be returned by `rcdk::get.mcs`
- Fixed: intensity filter failed if previous filters lead to zero feature groups.
- Fixed: `reportHTML()` annotation table was paged.
- Fixed: Check final path lengths of log files and truncate where necessary (reported by Corey Griffith)
- Fixed: in some cases the checking of `analysisInfo` validity may result in an error (reported by Tiago Sobreira)
- Fixed: `convertMSFiles()` error with `dirs=TRUE` (reported by Tiago Sobreira)
- Small updates/fixes for `installPatRoon()`
- Fixed: `screenSuspects()` did not take `onlyHits` into account for caching
- `screenSuspects()`: The original suspect name is stored in the `name_orig` column
- `XCMS3` feature group optimization: `binSize` and `minFraction` values were rounded while they shouldn't (issue #27)


# patRoon 1.2.0

This releases focuses on a significantly changed suspect screening interface, which brings several utilities to assist suspect annotation, prioritization, mixing of suspect and full NTA workflows and other general improvements.

**IMPORTANT**: The suspect screening interface has changed significantly. Please read the documentation (`?screenSuspects` and the handbook) for more details. If you want to quickly update your code without using any new functionality:

Change your existing code, e.g.

```r
scr <- screenSuspects(fGroups, suspectList, ...)
fGroupsScr <- groupFeaturesScreening(fGroups, scr)
```

to

```r
fGroupsScr <- screenSuspects(fGroups, suspectList, ..., onlyHits = TRUE)
```

**Major changes**

* New suspect screening interface
    * By default, feature groups without suspect hit are _not_ removed (unless `onlyHits=TRUE`). This allows straightforward mixing of suspect and full non-target workflows.
    * The feature groups are _not_ renamed tot the suspect name anymore. If you quickly want to assess which suspects were found, use the `screenInfo()` or `as.data.table()` methods.
    * Subsetting of suspsect screening results can be done with the `suspects` argument to `[`, e.g. `fGroupsScr[, suspects = "carbamazepine"]`
    * A new method, `annotateSuspects()`, allows combining the annotation workflow data (peak lists, formulas, compounds) to perform a detailed annotation for the suspects found during the workflow. This method calculates properties such as
        * Rankings: where is the suspect formula/compound ranked in the candidates from the workflow data
        * Annotation similarity: how well does the MS/MS spectrum of the suspect matches with formula/compound annotations.
        * An _estimation_ of identification levels to assist in quickly evaluating how well the suspect was annotated. The rules for identification levels are fully configurable.
    * A dedicated `filter()` method for suspect screening results, which allows you to easily prioritize data, for instance, by selecting minimum annotation ranks and similarities, identification levels and automatically choosing the best match in case multiple suspects are assigned to one feature (and vice versa).
    * A dedicated `as.data.table()` method and reporting functionality for suspect screening results to quickly inspect their annotation data.
    * Please refer to the updated suspect screening sections in the handbook and `?screenSuspects` and `?annotateSuspects` for  more information.
* Changes to suspect lists
    * Whenever possible, suspect information such as formulae, neutral masses, InChIKeys etc will be calculated for the input suspect list (obtainable afterwards with `screenInfo()`).
    * The suspect names will be checked to be file compatible, and automatically adjusted if necessary.
    * If MS/MS fragments are known for a suspect (formula or `m/z`), these can be included in the suspect list to improve suspect annotation.
    * The old suspect screening support for `features` objects was removed. The same and much more functionality can be obtained by the workflow for feature groups.
* The `reportCSV()` function was simplified and uses `as.data.table()` to generate the CSV data. This should give more consistent results.
* The `individualMoNAScore` MetFrag scoring is now enabled by default.

Other changes

* `reportHTML()` now allows toggling visibility for the columns shown in the feature annotation table.
* The `plotVenn()` method for `featureGroups` now allows to compare combinations of multiple replicate groups with each other. See `?plotVenn` for more information.
* Fix: locating `SIRIUS` binary on `macOS` did not work properly
* Fix: timeout warning for `GenForm` resulted in an error (https://github.com/rickhelmus/patRoon/issues/18)
* Fix: plotting structures resulted in a Java error on the RStudio Docker image (https://github.com/rickhelmus/patRoon/issues/18)

            
# patRoon 1.1

* **IMPORTANT**: The `plotEIC()`, `groups()` and `plotSpec()` methods were renamed to `plotChroms()`, `groupTable()` and `plotSpectrum()`. This was done to avoid name clashes with `XCMS` and `CAMERA`. The old functions still work (with a warning), but please update your scripts as these will be removed in the future.
* **IMPORTANT**: Major changes to the parallelization functionality
    * `patRoon` now supports an additional method to perform parallelization for tools such as `MetFrag`, `SIRIUS` etc. The main purpose of this method is to allow you to perform such calculations on external computer clusters. Please see the updated parallelization section in the handbook for more details.
    * The `logPath` and `maxProcAmount` arguments to functions such `generateFormulas`, `generateCompounds` etc were removed. These should now solely be configured through package options (see `?patRoon`).
    * The `patRoon.maxProcAmount` package option was renamed to `patRoon.MP.maxProcs`.
* Changes related to `SIRIUS`
    * **IMPORTANT:** Support for SIRIUS 4.5.0. Please update to this version since these changes break support for older versions.
    * Fix: SIRIUS formula calculation with `calculateFeatures=TRUE` would try to calculate formulas for features even if not present (eg after being removed by subsetting or filtering steps).
    * The `SIRBatchSize` argument to formula and compound generation functions was renamed to `splitBatches` and its meaning has slightly changed. Please see the reference manual (e.g. `?generateFormulas`) for more details.
* Changes related to MetFrag
    * Paths to local database files for MetFrag are now normalized, which makes handling of relative paths more reliable.
    * Changes in the specified local MetFrag database files are now inspected to improve caching.
    * Consistency: `generateCompoundsMetfrag` was renamed to `generateCompoundsMetFrag`.
* Optimized loading of spectra and EIC data.
* New utility functions
    * `withOpt()` to temporarily change (`patRoon`) package options.
    * `printPackageOpts()`: display current package options of `patRoon`.
* Finding features with `OpenMS`: potentially large temporary files are removed when possible to avoid clogging up disk space (especially relevant on some Linux systems where `/tmp` is small).
* Several packages such as `XCMS` are not attached by default, which significantly speeds up loading `patRoon` (e.g. with `library()`).
* The `compoundViewer()` function was marked as defunct, as it hasn;t been working for some time and its functionality is largely replaced by `reportHTML()`.
* `generateComponentsNontarget()`: update homolog statistics for merged series.
* `checkChromatograms()`: fix error when `fGroups` has only one replicate group
* `convertMSFiles()`: If `algorithm="pwiz"` and vendor centroiding is used then any extra filters are now correctly put after the `peakPicking` filter.
* `getXCMSnExp()` is now properly exported and documented.

# patRoon 1.0.4
* Small compatibility fixes for macOS
* Updated support for latest PubChemLite
* The `annoTypeCount` score for annotated compounds with PubChemLite is now not normalized by default anymore when reporting results.
* `reportHTML()` now correctly handles relative paths while opening the final report in a browser.


# patRoon 1.0.3
* `componentsNT`: include algorithm data returned by `nontarget::homol.search` in `homol` slot (suggested by Vittorio Albergamo)
* several `convertMSFiles()` fixes (issue #14)
    * prevent error when no input files are found
    * only allow one input/output format (didn't properly work before)
    * recognize that Waters files are directories
    * `cwt` option is now available for conversion with ProteoWizard
* minor fixes for subsetting XCMS `features` objects
* Fixed: `generateCompoundsMetFrag()`: compound names could be sometimes be interpreted as dates (reported by Corey Griffith)
* Fixed: on very rare cases empty peaklists could be present after averaging.
* Fixed: `SIRIUS` annotation didn't use set adduct but used default instead
* `SIRIUS` results are better handled if choosen adduct is not `[M+H]+` or `[M+H]+`
* More fixes for loading `data.table` objects properly from cache.
* RStudio Docker image: see the updated installation instructions in the handbook (thanks to Thanh Wang for some last minute fixes!)


# patRoon 1.0.2

* Fixed: avoid errors when SIRIUS returns zero results (reported by Vittorio Albergamo)
* Fixed: `plotGraph()` didn't properly handle components without linked series (reported by Vittorio Albergamo)
* Keep signal to noise data when importing/exporting XCMS data (`sn` column) (suggested by Ricardo Cunha)
* Reversed argument order of `exportedData`/`verbose` to `getXCMSSet()` functions to avoid ambiguities
* Automated tests for importing/exporting XCMS(3) data + small fixes for surfaced bugs
* `generateComponentsNontarget()`: allow wider _m/z_ deviation for proper linkage of different series (controlled by `absMzDevLink` argument).
* Fixed: `addAllDAEICs()` sometimes used wrong names for EICs
* Improved handling of empty feature groups objects when reporting
* Fixed: `reportPDF()` may report formula annotated spectra of results not present in input `featureGroups`
* Fixed: Loading `data.table` data from cache now calls `data.table::setalloccol()` to ensure proper behavior if `data.table::set()` is called on cached data.
* Fixed: plotSpec() for `compounds` with `useGGPlot2=TRUE` would try to plot formulas for non-annotated peaks (resulting in many extra diagonal lines)  
* Fixed: some functions involved in caching plot data for HTML reports sometimes returned invalid data.
* Fixed: EICs plotted by `reportPDF()` where not properly placed in a grid (as specified by `EICGrid` argument)
* Small tweaks/fixes for `reportHTML()`
    * now displays subscripted neutral formulae
    * Fixed: x axis title of EICs in annotation tab was cut-off
    * Fixed: The rt vs mz plot in the summary page now uses minutes for retention times if `retMin=TRUE`
* Updates for SIRIUS 4.4.29


# patRoon 1.0.1

* Perform neutral mass calculation for suspect screening with OpenBabel to avoid some possible Java/RCDK bugs on Linux.


# patRoon 1.0

## June 2020

* Fixed: `newProject()` didn't show polarity selection if only a compound identification algorithm was selected.
* Updated external dependency versions in installer script. 
* Fixed: `groupFeaturesXCMS3()` didn't properly cache results.
* `MSPeakLists`: results for averaged peak lists are now the same order as the input feature groups
* Fixed: XCMS(3) feature group import used wrong variable name internally (reported by Ricardo Cunha)

## May 2020

* **IMPORTANT** Major changes were made related to `SIRIUS` support 
    * Multiple features can now be annotated at once by `SIRIUS` (configurable with new `SIRBatchSize` function argument). This dramatically improves overal calculation times (thanks to Markus Fleischauer for pointing out this possibility!).
    * `generateFormulasSirius()` and `generateCompoundsSirius()` are now properly capitalized to `generateFormulasSIRIUS()` and `generateCompoundsSIRIUS()`
    * Support for `SIRIUS` 4.4.
    * If all features are annotated at once then `SIRIUS` output is directly shown on the console.
    * The amount of cores used by `SIRIUS` can be specified with the `cores` function arguments.
    * More extra commandline options can be given to `SIRIUS`
* Fixed: `groupNames()`, `analyses()` and similar methods sometimes returned `NULL` instead of an empty `character` vector for empty objects.
* `plotHeatMap()` with `interactive=TRUE`: switch from now removed `d3heatmap` package to `heatmaply`
* Fixed: `reportHTML()` didn't split PubChem URLs when multiple identifiers were reported.
* `PWizBatchSize` argument for `convertMSFiles()`


## April 2020

* `extraOptsRT`/`extraOptsGroup` arguments for OpenMS feature grouping to allow custom command line options.
* `importFeatureGroupsBrukerTASQ`
    * now correctly takes retention times of suspects into account when creating feature groups.
    * retention times / _m/z_ values are now averaged over grouped suspects.
* The `plot()` method for `featureGroups` now allows drawing legends when `colourBy="fGroups"` and sets `colourBy="none"` by default, both for consistency with `plotEIC()`.
* All documentation is now available as PDF files on the website (https://rickhelmus.github.io/patRoon/)
* `newProject()` now uses XCMS3 algorithms instead of the older XCMS interface.
* Fixed: features in objects generated by `xcms` (not `xcms3`) could not be subset with zero analyses (which resulted in errors by e.g. `unique()` and `reportHTML()`). Reported by Corey Griffith.


## March 2020

* Fixed: Normalization of scorings for formulae/compounds potentially uses wrong data after subsetting/filtering of `formulas`/`compounds` objects
* Suspect screening
    * Fixed: Errors/Warnings of missing data in suspect list were not shown if using cached data
    * If a value is missing in any of the columns from the suspect list used for automatic ion mass calculation (e.g. SMILES, formula, ...) then data from another suitable column is tried.
    * Fixed: invalid neutral mass calculation for suspects with charged SMILES/InChIs
    * Default adduct can be specified in `newProject()` dialog
* Small compatibility fix for feature finding with OpenMS 2.5 (reported by Thanh Wang)
* RAMClustR is now supported from CRAN and no need to install github package anymore
* pubchemlite identifiers are now URL linked in HTML reports
* related CIDs are now reported for PubChemLite results.
* MetFrag compound generation: removed `addTrivialNames` option as it never worked very well.
* `reportHTML()`: only components with reported feature groups are now reported.


## February 2020

* Several small improvements and fixes for TASQ import
* Suspect screening:
    * now also support chemical formula for automatic m/z calculation
    * more robust loading of suspect lists (e.g. skip suspects with missing/invalid data)


## January 2020

* Ignore user specified scorings for local databases such as CompTox that are not actually present in the DB. This makes it easier to use e.g. different DB versions with differing scorings.
* Add scorings from wastewater and smoking metadata comptox MetFrag databases
* Windows install script now install latest (March2019) CompTox
* Updates for latest PubChemLite relaease (Jan2020)
* Suspect screening now doesn't require pre-calculated ion `m/z` values. Instead, suspect lists can contain SMILES, InChI or neutral mass values which are used for automatic ion `m/z` calculation. See `?screenSuspects` for more details.


## December 2019

* Added missing score terms for latest CompTox MetFrag database
* labels parameter for formulas/compounds methods of `consensus()`
* Fixed: colour assignment for scores plotting of merged formulae/compound results might be incorrect (reported by Emma Schymanski)
* Fixed: analysis table in `newProject()` UI only showed partial amount of rows.
* Fixed: don't print normalized instead of absolute design parameters when only one parameter is optimized in DoEs (fixes issue #10)


## November 2019

* **IMPORTANT** The `addFormulaScoring()` function now uses a different algorithm to calculate formula scores for compound candidates. The score is now based on the actual formula ranking in the provided `formulas` object, and is fixed between _zero_ (no match) and _one_ (best match).
* Formula feature consensus:
    * All scorings are now averaged, including those that are not fragment specific (e.g. precursor m/z error)
    * This also improves ranking in certain specific cases
* Vectorized plotting of MS spectra to make it potentially faster
* Added PubChemLite support for MetFrag


## October 2019

* Fixed: `convertMSFiles` correctly checks if input exists
* Specific optimizations after benchmarking results:
    * `maxProcAmount` (i.e. number of parallel processes) now defaults to amount of physical cores instead of total number of CPU threads.
    * Decreased `batchSize` to `8` for GenForm formula calculation.
* `plot()` for `featureGroups` can now highlight unique/shared features across replicates (suggested by V Albergamo)
* Linking of homologous series:
    * Improved info descriptions for `plotGraph()`
    * Series are now properly unlinked when merging (was too greedy)
    * Better algorithm to detect conflicting series
    * Fixed bug when updating removed links
* `concs` option for `generateAnalysisInfo()` to set concentration data


## September 2019

* Labels for objects in a `featureGroupsComparison` can be customized (useful for e.g. plotting)
* Caching and progress bar support for suspect screening
* Updated/Fixed JDK installation for installation script
* Fixed missing pipe operator import (`%>%`)


## August 2019

* `topMost` argument for GenForm formula calculation.
* Added XCMS3 support for finding and grouping features, importing/exisiting data and parameter optimization (i.e. mostly on-par with classic XCMS support).
* Changed compound result column name from InChi to InChI


## June 2019

* **IMPORTANT** Several things are renamed for clarity/consistency
    * The column to specify replicate groups for blank subtraction in the analysis information is re-named from `ref` to `blank`. Similarly, the `refs` argument to `generateAnalysisInfo()` is now called `blanks`.
    * `reportMD()` is renamed to `reportHTML()`
    *  `filter()` method for `formulas`: `minExplainedFragPeaks` is now called `minExplainedPeaks`
    * `screenTargets` and its `targets` parameter have been renamed to `screenSuspects()` / `suspects`
* Fixed incorrect selection after feature table (or other interactive tables) have been manually re-ordered (reported by Thanh Wang)
* `groups()` and `as.data.table()` methods for `featureGroups`: optionally consider feature areas instead of peak intensities.
* `plotSilhouettes()` method for `compoundsCluster`
* Added `rGroups` argument to subset operator for `featureGroups` to subset by replicate groups (equivalent to `rGroups` argument to `filter()`).
* Improved logging of output from CLI tools (e.g. OpenMS, MetFrag, SIRUS, ...) 

## May 2019

* Formula updates
    * `GenForm` formula calculation with `MSMode="both"` (the default): instead of repeating calculations with and without MS/MS data and combining the data, it now simply does either of the two depending on MS/MS data availability. The old behavior turned out to be redundant, hence, calculation is now a bit faster.
    * `GenForm` now perform _precursor isolation_ to cleanup MS1 data prior to formula calculation. During this step any mass peaks that are unlikely part of the isotopic pattern of the feature are removed, which otherwise would penalize the isotopic scoring. The result is that isotopic scoring is dramatically improved now. This filter step is part of new filter functionality for `MSPeakLists`, see `?MSPeakLists` and `?generateFormulas` for more information.
    * When formula consensus are made from multiple features the scorings and mass errors are now averaged (instead of taking the values from the best ranked feature).
    * Improved ranking of candidates from a consensus of multiple formula objects (see `?formulas`).
* Consensus for compounds are now similarly ranked as formulas.
* More consistent minimum abundance arguments for `consensus()` (`absMinAbundance` and `relMinAbundance`)
* `MetFrag`: for-ident database and new statistical scores are now supported
* `as.data.table()` / `as.data.frame()` for `featureGroups` now optionally reports regression information, which may be useful for quantitative purposes. This replaces the (defunct) `regression()` method and limited support from `screenTargets()`.
* `plotGraph()` method to visually inspect linked homologous series.

## April 2019

* Misc small tweaks and fixes for `newProject()` (e.g. loading of example data).
* Improved graphical output of various common plotting functions.
* Updated tutorial vignette and added handbook


## March 2019
* `reportMD()`: most time consuming plots are now cached. Hence, re-reporting should be signficiantly faster now.
* Updates to MS data file conversion:
    * `convertMSFiles()` now (optionally) takes analysis information (`anaInfo`) for file input.
    * `convertMSFiles()` now supports Bruker DataAnalysis as conversion algorithm (replaces now deprecated `exportDAFiles()` function).
    * `MSFileFormats()` function to list supported input conversion formats.
    * `generateAnalysisInfo()` now recognizes more file formats. This is mainly useful so its output can be used with `convertMSFiles()`.
    * `convertMSFiles()` now has the `centroid` argument to more easily perform centroiding.
* Updates to `newProject()`:
    * The analyses selector recognizes more data file formats. This way you can select analyses that have not been converted yet.
    * Data pre-treatment options now include more sophisticated file conversion options (_e.g._ using ProteoWizard). This and the new analysis selector functionality ensures that data files in all major vendor formats do not have to be converted prior to generating a script.
    * Re-organized tabs to mirror non-target workflow.
    * Suspect screening support.
    * Improved layout of output script.
* `withMSMS` filter for MS peak lists.
* Timeout for formula calculation with GenForm to avoid excessive calculation times.
* `importFeatures()` generic function
* Reporting functions renamed arguments related to compounds reporting (e.g. compoundTopMost to compound**s**TopMost)


## February 2019
* Compound scorings are now normalized towards their original min/max values. This ensures that the `score` column of MetFrag results stays correct.
* plotScores(): Option to only report scorings that have been used for ranking
* as.data.table()/as.data.frame() method for compounds: optionally normalize scores.
* `reportPDF()`/`reportMD()` now report only 5 top most candidate compounds by default (controlled by `compoundsTopMost` argument).
* metadata for MS peak lists
* `plotSpec()` now displays subscripted formulae
* **IMPORTANT** Several major changes were made to the `filter()` methods for `features` and `featureGroups`. Please carefully read the updated documentation for these methods! (i.e. `` ?`filter,features-method` `` and `` ?`filter,featureGroups-method` ``).
    * Most argument have been renamed for consistency, simplicity and clarity.
    * The order when multiple filters are specified to the `featureGroups` method was adjusted, notably to improve reliability of blank filtration. Again, please see `` ?`filter,featureGroups-method` ``.
    * The following new filters were added:
        * mass defect range (`mzDefectRange` argument)
        * maximum relative standard deviation (RSD) of intensities between replicates (`maxReplicateIntRSD` argument)
        * minimum number of features within analyses (`absMinFeatures` and `relMinFeatures` arguments).
        * pre-intensity filters (`preAbsMinIntensity` and `preRelMinIntensity` arguments)
        * most existing filters now accept both relative and absolute values.
    * The script generation functionality of `newScript()` has been updated and supports more filter types.
    * The `repetitions` argument is not needed anymore for the new algorithm and has been removed.
    * `Inf` values now should be used to specify no maximum for range filters (was `-1`).
* Fixed: GenForm now always uses Hill sorting.
* `annotatedPeakList()` method for `formulas` and `compounds`. Also used by `reportMD` for improved annotation peak tables.
* Tweaked default mzR MS peak lists settings (halved `maxRtMSWidth` and `precursorMzWindow`)
* Fixed: Make sure that MetFrag web doesn't try to set unsupported database
* **IMPORTANT** Several changes were made to improve clarity and consensistency for arguments that specify retention/mz windows or allowable deviations.
    * Functions with changed argument names: `generateComponentsNontarget`, `generateComponentsRAMClustR`, `generateCompoundsSirius`, `generateFormulasGenForm`, `generateFormulasSirius`, `generateMSPeakListsDA`, `generateMSPeakListsMzR`, `importFeatureGroupsBrukerPA`
    * The `maxRtMSWidth` argument to `generateMSPeakListsDA`, `generateMSPeakListsMzR` (now `maxMSRtWindow`) now specifies a retention time window (\emph{i.e.} +/- retention time feature) instead of total retention width around a feature. Hence, _current input values should be halved_.
* CAMERA and RAMClustR components: both now have `minSize` and `relMinReplicates` (replaces `ubiquitous` for CAMERA) arguments. Note that their defaults may filter out (feature groups from) components. See their documentation for more info.
* Changed capitalisation of MetFrag CL location option from `patRoon.path.metFragCL` to `patRoon.path.MetFragCL`. The old name still works for backward compatability.
* Documented usage of the CompTox database with MetFrag. See `?generateCompounds`.
* Default normalization of MetFrag scorings now follows MetFrag web behaviour.
* `topMostFormulas` argument for SIRIUS compound generation.
* Fixed GenForm ranking in case both MS and MS/MS formulae are present.
* `reportPDF()`/`reportMD()` now report only 5 top most candidate formulae by default (controlled by `formulasTopMost` argument).
* Added `verifyDependencies()` function to let the user verify if external tools can be found.
* The meaning of the `dirs` argument to `convertMSFiles()` was slightly changed: if `TRUE` (the default) the input can either be paths to analyses files or to directories containing the analyses files.
* More effective locating ProteoWizard binaries by using the Windows registry.
* Nicer default graphics for `featureGroups` method for `plot()`.
* `reportMD()`: Don't plot Chord if <3 (non-empty) replicate groups are available. 
* All `filter()` methods now support negation by `negate` argument.


## January 2019
* minSize and ubiquitous arguments for CAMERA component generation. See ?generateComponentsCamera.
* Various tweaks for plotEIC() and plotSpec() methods
* Various small additions to newProject()
* `reportMD()`: added table with annotated fragments for compounds/formulas
* `consensus()` updates
    * `consensus()` methods now support extracting unique data. This also replaces the `unique()` method that was defined for `featureGroupsComparison`.
    * `comparison()` now automatically determines object names from algorithm (consistency with `consensus()` method for other objects).
    * Fixed: coverage calculation for consensus formulas now correctly based on precursor overlap (was overlap of precursor+fragment).    
* `plotVenn()` and `plotUpSet()` methods to compare different compounds or formulas objects.
* `filter()` method for components.
* DataAnalysis formula generation: fixed neutral formula calculation if `MSMode="msms"`, now needs `adduct` argument.
* Neutral loss filter for compounds.
* **IMPORTANT** Adduct specification is now simplified and more generic by usage of a new `adduct` class. This means that `generateCompounds()` and `generateFormulas()` now expect slightly differing arguments. Please see their manual pages.
* Workaround for homologous series generation with nontarget (see https://github.com/blosloos/nontarget/issues/6)
* Improvements to terminate background commandline processes when e.g. R is terminated. 
* `clearCache()` now supports removal of caches via regular expressions.
* Added/Improved `topMost` and `extraOpts` arguments for SIRIUS formula/compound generation.
* Annotated fragments from SIRIUS compounds now correctly contain charged molecular form.
* `filter()` method for compounds now support generic scoring filtering and on elements of precursor and fragment formulae.
* **IMPORTANT** Several changes were made to the MetFrag compound generation interface in order to simplify it and make it more generic. See `?generateCompounds` for more details (notably the Scorings section).
* More MS peak list updates
    * Precursor peaks are now flagged in MS peak list data and `plotSpec()`
    * Prune MS peak lists (not MS/MS) if no precursor could be determined (enabled by default, see `pruneMissingPrecursorMS` option in `?generateMSPeakLists`).
    * Better retain precursor peaks after filtering steps: only intensity thresholds may remove precursors (always for MS data, optional for MS/MS with `retainPrecursorMSMS` function arguments, see `?MSPeakLists` and `?generateMSPeakLists`).
* All major workflow classes now have `algorithm()` and `as.data.table()/as.data.frame()` methods. The latter replaces and enhances the `makeTable()` (`formulas` class) and `groupTable()` (`featureGroups` class) methods.


## December 2018
* Moved OpenMS XML writing code from `R` to `C++`: significantly reduces time required for grouping large amount of features.
* Several updates for functionality that uses Bruker DataAnalyses
    * Improved verification and consistency for handling processed data from DataAnalysis
    * Automatic saving & closing of analyses processed with DataAnalysis. Files are now generally closed by default to limit the resources needed by DataAnalysis. 
    * `revertDAAnalyses()` function: brings back set of Bruker analyses to their unprocessed state.
    * Minimum intensity arguments for Bruker DataAnalysis MS peak lists.
    * Slightly different `doFMF` behaviour for DataAnalysis feature finding.
* Several important updates were made to fomula calculation functionality.
    * The interface has been simplified as the functionality from the `formula` and `formulaConsensus` classes are now merged: there is no need to call `consensus()` anymore after `generateFormulas()`.
    * Formulae can now directly be calculated for feature groups by using group averaged MS peak lists (by setting `calculateFeatures=FALSE`). This can greatly speed up calulcation, especially with many analyses.
    * The new `filter()` and `as.data.table()`/`as.data.frame` methods bring new functionalities related to filtering, extracting data and performing several processing steps commonly performed for organic matter (OM) characterization.
    * Other updates on formulas
        * length now returns number of unique precursor formulas (was total number of results)
        * Fixed: Reported fragment formulas from SIRIUS were incorrectly assumed to be charged. Charged fragment formulas are now calculated manually (the neutral form is stored in the `frag_neutral_formula` column). This ensures correct comparison when a consensus is made.
        * `reportCSV()` now splits formulas for each feature group in separate CSV files (similar to `compounds` reporting).
        * Fixed: `reportPDF()` now actually includes formula annotations in annotated compound spectra when formulas are specified.
        * New oc argument when using GenForm: if enabled only organic formulae are accepted (i.e. with at least one carbon atom). Also incorporated a small fix for the actual GenForm binary to make this option work (https://sourceforge.net/p/genform/tickets/1/).
        * Fixed: coverage calculation of formulae across features treated formulae calculated only from MS data separately.
        * GenForm now also includes precursor annotation from MS/MS data.
* `file` argument for `clearCache()`
* Updates on MS peak lists
    * More consistent naming for algorithm specific MS peak list generators (i.e. `generateMSPeakListsX` where X is the algo).
    * Additional MS peak lists are generated by averaging the lists of features within a feature group.
    * `generateCompounds()` and plotting functionality now uses averaged group peak lists instead of peak list of most intense analysis.
    * `plotSpec()` method for MSPeakLists: plot (non-annotated) MS and MS/MS spectra.
    * Minimum intensity filter option that is applied after averaging.
    * Now uses "hclust" method for averaging by default, which now uses the [fastcluster] package.


## November 2018
* Default value for `maxRtMSWidth` argument used for peak list generation.
* Fixed: `maxRtMSWidth` argument for mzR peak list generation had no effect.
* Preliminary EPA DSSTox support (via LocalCSV).
* Added `addAllDAEICs()` function.
* Renamed `mzWidth` argument of `addDAEIC()` to `mzWindow`.
* Normalization of compound scores: normalization method can now be set and specified scorings can be excluded.
* Store/report IUPACName (as compoundName) from MetFrag PubChem data.
* Renamed trivialName to compoundName for compound tables.
* `convertMSFiles`: changed interface with more options, parallelization and ProteoWizard support.
* Automatic optimization of parameters necessary for feature finding and grouping. Heavily based on the IPO R package. See the 'feature-optimization' manual page.
* **IMPORTANT** `getXcmsSet()` is renamed to `getXCMSSet()`
* verbose option for `findFeatures()` / `groupFeatures()`
* Changed `nintersects` default for plotUpSet so that all intersections are plotted by default.
* plotChord() now properly stops if nothing overlaps.
* replicateGroupSubtract() now removes replicate groups that were subtracted.
* Fixed: replicateGroupSubtract() now correctly takes maximum mean intensity for threshold determination when multiple rGroups are specified.
* Fixed: Wrong compound clusters plotted in reportMD().
* Fixed: Added timeout after restarting failed command (e.g. MetFrag CL) to prevent rare error "The requested operation cannot be performed on a file with a user-mapped section open".


## October 2018
* OpenMS `features` class objects now store number of isotopes found for each feature.
* **IMPORTANT** Added all relevant options of FeatureFinderMetabo as function arguments to findFeaturesOpenMS() and renamed/reordered current options for more conistent style. Please check ?findFeatures for the updated function arguments!
* openReport option for reportMD(). If TRUE the generated report will be opened with a web browser.
* reportPlots option for reportMD() which collapses reportFGroups, reportChord and reportFormulaSpectra and adds control to plot Venn and UpSet diagrams.
* plotUpSet() methods to compare feature groups by UpSet plots. See e.g. http://caleydo.org/tools/upset/
* filter() method for features.
* EICs now loaded via faster C++ code thats uses mzR instead of XCMS
* Moved feature intensity loading code for OpenMS features to C++. This results in much faster feature finding.


## September 2018
* Removed filterBy methods: these are now deprecated with new subset operators and groupNames()/analyses() methods. Example: `fGroups <- fGroups[, groupNames(compounds)]`
* subset/extraction operators ("[", "[[" and "$") for features, featureGroups, MSPeakLists, formulas, formulaConsensus, compounds, compoundsCluster and components classes.
* analyses() and groupNames() generics to get analyses and feature group names of the data within an object.
* "[" method for featureGroups: empty feature groups now always dropped, drop argument now ignored.
* reportMD(): The layout to show compounds, formulas and components is now done with DataTables (DT package). This change allows faster initial loading of results. Furthermore, several small tweaks were done to improve general design.
* plotSpec() (compounds method): remove unused normalizeScores flag
* plotSpec() (compounds method): plotting of embedded structure now optional (plotStruct argument)
* plotSpec() (compounds method): automatic calculation of necessary extra height to plot labels/structure


## Augustus 2018
* The XML code required to load feature (group) data generated by OpenMS is now moved to a C++ interface that uses [Rcpp] and [pugixml]. This results in a significant reduction of required processing time. In addition, files are now processed in chunks, allowing even very large feature sets (>10000) without clogging up system memory.
* Improved general numeric comparisons, resulting in e.g. improved EIC generation.
* Tweaked OpenMS feature intensity loading: now takes intensity from data point closest to retention time instead of max intensity from datapoints in the search window. Furthermore, the search window for datapoints was reduced and made configurable.


## July 2018
* getMCS() method for compounds
* plotStructure() method for compounds will draw MCS when mutiple indices are specified


## June 2018
* Added removeRefAnalyses argument to filter() (featureGroups method) to easily remove e.g. analyses that are used as blanks after blank subtraction.
* Added filterBy() method which removes any feature groups from a featureGroups object of which no results are present in a specified object. Methods are defined for MSPeakLists, formulaConsenus, compounds and components. This method replaces some of the functionality of the filter() method for featureGroups (formConsensus and compounds arguments).
* Added mz and chromatographic peak width range options to filter() method for feature groups.
* Moved intensity clustering code (makeHCluster) to new component type (see componentsIntClust class documentation).


## May 2018
* Added compound clustering (see makeHCluster method for compounds). This is an useful tool to get an overview of all the candidate chemical structures after compound identification. The clustering will reduce data complexity. Furthermore, maximum common sucstructures (MCS) can be calculated and plotted for each cluster to get a quick impression of the different structures of candidates.
* Added function arguments checks using [checkmate]. This guards all exported functions and methods from wrong user input arguments. If any problems are found (e.g. a wrong data type or range was specified) then the user is informed about this and what kind of input is to be expected.
* Added workaround (removed latex dependencies added automatically by `kableExtra` package) that may cause memory leakage when `reportMD()` is called repeatedly.


## April 2018
* Added unit tests (using [testthat]) and fixed several small bugs that were revealed in the process.
* Continuous integration (CI) with running tests on [CircleCI] (Linux builds) and [AppVeyor] (Windows builds) and testing coverage on [Codecov]. Docker images with patRoon and all its dependencies are automatically pushed on [Docker Hub][DH].
* Many small bug fixes.




[Rcpp]: http://www.rcpp.org/
[pugixml]: https://pugixml.org/
[checkmate]: https://github.com/mllg/checkmate
[testthat]: https://github.com/r-lib/testthat
[CircleCI]: https://circleci.com/gh/rickhelmus/patRoon
[AppVeyor]: https://ci.appveyor.com/project/rickhelmus/patroon/branch/master
[Codecov]: https://codecov.io/gh/rickhelmus/patRoon
[DH]: https://hub.docker.com/r/patroonorg/patroon/
[fastcluster]: https://cran.r-project.org/web/packages/fastcluster/index.html
