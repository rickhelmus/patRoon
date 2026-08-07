# Report feature group data (legacy interface)

Functionality to report data produced by most workflow steps such as
features, feature groups, calculated chemical formulae and tentatively
identified compounds. This is the legacy interface, for the updated
interface see
[reporting](https://rickhelmus.github.io/patRoon/reference/reporting.md).

## Usage

``` r
reportCSV(
  fGroups,
  path = "report",
  reportFeatures = FALSE,
  formulas = NULL,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
    "annotHitCount", "libMatch"),
  compsCluster = NULL,
  components = NULL,
  retMin = TRUE,
  clearPath = FALSE
)

reportPDF(
  fGroups,
  path = "report",
  reportFGroups = TRUE,
  formulas = NULL,
  formulasTopMost = 5,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  reportFormulaSpectra = TRUE,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
    "annotHitCount", "libMatch"),
  compoundsOnlyUsedScorings = TRUE,
  compoundsTopMost = 5,
  compsCluster = NULL,
  components = NULL,
  MSPeakLists = NULL,
  retMin = TRUE,
  EICGrid = c(2, 1),
  EICParams = getDefEICParams(window = 20, topMost = 1, topMostByReplicate = TRUE),
  clearPath = FALSE
)

# S4 method for class 'featureGroups'
reportCSV(
  fGroups,
  path = "report",
  reportFeatures = FALSE,
  formulas = NULL,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
    "annotHitCount", "libMatch"),
  compsCluster = NULL,
  components = NULL,
  retMin = TRUE,
  clearPath = FALSE
)

# S4 method for class 'featureGroups'
reportPDF(
  fGroups,
  path = "report",
  reportFGroups = TRUE,
  formulas = NULL,
  formulasTopMost = 5,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  reportFormulaSpectra = TRUE,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
    "annotHitCount", "libMatch"),
  compoundsOnlyUsedScorings = TRUE,
  compoundsTopMost = 5,
  compsCluster = NULL,
  components = NULL,
  MSPeakLists = NULL,
  retMin = TRUE,
  EICGrid = c(2, 1),
  EICParams = getDefEICParams(),
  clearPath = FALSE
)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that should be used for reporting data.

- path:

  The destination file path for files generated during reporting. Will
  be generated if needed.

- reportFeatures:

  If set to `TRUE` then for each analysis a `.csv` file will be
  generated with information about its detected features.

- formulas, compounds, compsCluster, components:

  Further objects
  ([`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md),
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md),
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md),
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md))
  that should be reported. Specify `NULL` to skip reporting a particular
  object.

- compoundsNormalizeScores, formulasNormalizeScores:

  A `character` that specifies how normalization of annotation scorings
  occurs. Either `"none"` (no normalization), `"max"` (normalize to max
  value) or `"minmax"` (perform min-max normalization). Note that
  normalization of negative scores (e.g. output by `SIRIUS`) is always
  performed as min-max. Furthermore, currently normalization for
  `compounds` takes the original min/max scoring values into account
  when candidates were generated. Thus, for `compounds` scoring,
  normalization is not affected when candidate results were removed
  after they were generated (*e.g.* by use of `filter`).

- compoundsExclNormScores, formulasExclNormScores:

  A `character` vector specifying any compound scoring names that should
  *not* be normalized. Set to `NULL` to normalize all scorings. Note
  that whether any normalization occurs is set by the
  `compoundsExclNormScores,formulasExclNormScores` argument.

  For `compounds`: By default `score` and `individualMoNAScore` are set
  to mimic the behavior of the `MetFrag` web interface.

- retMin:

  If `TRUE` then report retention times in minutes (otherwise seconds).

- clearPath:

  If `TRUE` then the destination path will be (recursively) removed
  prior to reporting.

- reportFGroups:

  If `TRUE` then feature group data will be reported.

- formulasTopMost, compoundsTopMost:

  Only this amount of top ranked candidate formulae/compounds are
  reported. Lower values may significantly speed up reporting. Set to
  `NULL` to ignore.

- reportFormulaSpectra:

  If `TRUE` then explained MS/MS spectra (if available) for candidate
  formulae will be reported. Specifying `formulas` and setting this
  argument to `FALSE` still allows further annotation of compound MS/MS
  spectra.

- compoundsOnlyUsedScorings:

  If `TRUE` then only scorings are plotted that actually have been used
  to rank data (see the `scoreTypes` argument to
  [`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md)
  for more details).

- MSPeakLists:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that is *mandatory* when spectra for formulae and/or compounds
  will be reported.

- EICGrid:

  An integer vector in the form `c(columns, rows)` that is used to
  determine the plotting grid when reporting EICs in PDF files.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

## Details

These functions are usually called at the very end of the workflow. It
is used to report various data on features and feature groups. In
addition, these functions may be used for reporting formulae and/or
compounds that were generated for the specified feature groups. Data can
be reported in tabular form (*i.e.* `.csv` files) by `reportCSV` or
graphically by `reportPDF` and `reportHTML`. The latter functions will
plot for instance chromatograms and annotated mass spectra, which are
useful to get a graphical overview of results.

All functions have a wide variety of arguments that influence the
reporting process. Nevertheless, most parameters are optional and only
required to be given for fine tuning. In addition, only those objects
(*e.g.* formulae, compounds, clustering) that are desired to be reported
need to be specified.

`reportCSV` generates tabular data (*i.e.* `.csv` files) for given data
to be reported. This may also be useful to allow import by other tools
for post processing.

`reportPDF` will report graphical data (*e.g.* chromatograms and mass
spectra) within PDF files. Compared to `reportHTML` this function may be
faster and yield smaller report files, however, its functionality is a
bit more basic and generated data is more 'scattered' around.

## Note

Any formulae and compounds for feature groups which are not present
within `fGroups` (*i.e.* because it has been subset afterwards) will not
be reported.

The `topMost`, `topMostByReplicate` and `onlyPresent` [EIC
parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
may be ignored, *e.g.*, when generating overview plots.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by the report functions to process HRMS (or IMS-HRMS)
data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## See also

reporting
