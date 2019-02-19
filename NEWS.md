# patRoon 0.1.0.9000

## February 2019
* Compound scorings are now normalized towards their original min/max values. This ensures that the `score` column of MetFrag results stays correct.
* plotScores(): Option to only report scorings that have been used for ranking
* as.data.table()/as.data.frame() method for compounds: optionally normalize scores.
* reportPDF()/reportMD() now report only 5 top most candidate compounds by default (controlled by compoundTopMost argument).
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
