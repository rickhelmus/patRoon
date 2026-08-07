# Report workflow data

Functionality to report data produced by most workflow steps such as
features, feature groups, formula and compound annotations, and TPs.

## Usage

``` r
report(
  fGroups,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  compsCluster = NULL,
  components = NULL,
  TPs = NULL,
  settingsFile = system.file("report", "settings.yml", package = "patRoon"),
  path = NULL,
  EICParams = getDefEICParams(topMost = 1, topMostByReplicate = TRUE),
  EIMParams = getDefEIMParams(topMost = 1, topMostByReplicate = TRUE),
  specSimParams = getDefSpecSimParams(),
  clearPath = FALSE,
  openReport = TRUE,
  parallel = FALSE,
  overrideSettings = list()
)

# S4 method for class 'featureGroups'
report(
  fGroups,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  compsCluster = NULL,
  components = NULL,
  TPs = NULL,
  settingsFile = system.file("report", "settings.yml", package = "patRoon"),
  path = NULL,
  EICParams = getDefEICParams(topMost = 1, topMostByReplicate = TRUE),
  EIMParams = getDefEIMParams(topMost = 1, topMostByReplicate = TRUE),
  specSimParams = getDefSpecSimParams(),
  clearPath = FALSE,
  openReport = TRUE,
  parallel = FALSE,
  overrideSettings = list()
)

genReportSettingsFile(out = "report.yml", baseFrom = NULL)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that should be used for reporting data.

- MSPeakLists, formulas, compounds, compsCluster, components, TPs:

  Further objects
  ([`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md),
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md),
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md),
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md),
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md),
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md))
  that should be reported. Specify `NULL` to skip reporting a particular
  object. Note that `MSPeakLists` must be set if either `formulas` or
  `compounds` is set.

- settingsFile:

  The path to the report settings file used for report configuration
  (see `Report settings`).

- path:

  The destination file path for files generated during reporting. Will
  be generated if needed. If `path=NULL` then the destination path is
  taken from the report settings (see below).

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- EIMParams:

  A named `list` with parameters used for extracted ion mobilogram (EIM)
  creation. See the [EIM
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- clearPath:

  If `TRUE` then the report destination path will be (recursively)
  removed prior to reporting.

- openReport:

  If set to `TRUE` then the output report file will be opened with the
  system browser.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

  **NOTE**: parallelization is disabled by default, as it may slow down
  reporting on some systems (*e.g.* Windows) and under some
  circumstances. It is best to experiment with this setting to see if it
  speeds up report generation for your system and data.

- overrideSettings:

  A `list` with settings that override those from the report settings
  file. Example: `overrideSettings=list(compounds=list(topMost=25))`.

- out:

  The output file path.

- baseFrom:

  An existing report file to which the report settings should be based
  from. This is primarily used to update old settings files: the output
  settings file will be based on the old settings and amended with any
  missing.

## Details

The reporting functionality is typically used at the very end of the
workflow. It is used to overview the data generated during the workflow,
such as features, their annotations and TP screening results.

`report` reports all workflow data in an interactive HTML file. The
reports include both tabular data (*e.g.* retention times, annotation
properties, screening results) and various plots (*e.g.* chromatograms,
(annotated) mass spectra and many more). This function uses
functionality from other R packages, such as
[rmarkdown](https://CRAN.R-project.org/package=rmarkdown),
[knitr](https://CRAN.R-project.org/package=knitr) and
[bslib](https://CRAN.R-project.org/package=bslib).

The `genReportSettingsFile` function generates a new template `YAML`
file to configure report settings (see the next section).

## Note

No data will be reported for feature groups in any of the reported
objects (`formulas`, `compounds` etc) which are *not* present in the
input
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object (`fGroups`).

The `topMost`, `topMostByReplicate` and `onlyPresent` [EIC
parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
may be ignored, *e.g.*, when generating overview plots.

## Report settings

The report generation can be customized with a variety of settings that
are read from a `YAML` file. This is especially useful if you want to
change more advanced settings or want to add or remove the parts that
are reported. The report settings file is specified through the
`settingsFile` argument. If not specified then default settings will be
used. To ease creation of a new template settings file, the
`genReportSettingsFile` function can be used.

The following settings are currently available:

- General

  - `version`: version of the settings file.

  - `format`: the report format. Currently this can only be `"html"`.

  - `path`: the destination path (ignored if the `path` argument is
    specified).

  - `keepUnusedPlots`: the number of days that unused plot files are
    kept (see `Plot file caching`).

  - `selfContained`: If `true` then the output `report.html` embeds all
    graphics and script dependencies. Otherwise these files are read
    from the `report_files/` directory. Self-contained reports are
    generally smaller and easily shared, since only the `report.html`
    needs to be copied. However, they are slower to generate and plots
    cannot be cached.

  - `noDate`: Set to `true` to omit the date from the report. Mainly
    used for internal purposes.

- `summary`: defines the plots on the summary page: `chord`, `venn`
  and/or `upset`.

- `features`

  - `retMin`: if `true` then retention times are reported in minutes.

  - `chromatograms`

    - `large`: inclusion of large chromatograms (used in feature group
      table and TP parent chromatogram view).

    - `small`: inclusion of small chromatograms (feature group table).

    - `features`: inclusion of chromatograms for individual features
      (features view). Set to `all` to also include plots for analyses
      in which a feature was not found (or removed afterwards).

    - `intMax`: Method to determine the maximum intensity plot range:
      `eic` or `feature`. Sets the `intMax` argument to `plotChroms`.

  - `mobilograms`

    - `large`, `small`, `features`: inclusion of mobilogram plots, see
      `chromatograms` above.

  - `intensityPlots`: inclusion of intensity trend plots.

  - `aggregateConcs`, `aggregateTox`: function name used for
    concentration and toxicity aggregation, *e.g.* `mean`.

- `MSPeakLists`

  - `spectra`: inclusion of MS and MS/MS spectra (not annotated).

- `formulas`

  - `include`: whether formula results are reported (formula view). If
    `false` then the input `formulas` object is still used to amend
    *e.g.* compound annotated spectra.

  - `normalizeScores`, `exclNormScores`: controls score normalization
    and which score fields to exclude from normalization; these are
    forwarded to *e.g.* `plotScores`.

  - `topMost`: only report this number of top ranked candidates. This
    number can be lowered to speed-up report generation.

- `compounds`

  - `normalizeScores`, `exclNormScores`, `topMost`: same as `formulas`,
    see above.

  - `onlyUsedScorings`: if `true` only scorings used by the current
    dataset are considered when normalizing or reporting compound
    scores.

- `TPs`

  - `graphs`: inclusion of TP hierarchy graphs (generated with
    [`plotGraph`](https://rickhelmus.github.io/patRoon/reference/generics.md)).

  - `graphStructuresMax`: maximum number of structures to plot in
    hierarchy graphs (sets `structuresMax` argument of
    [`plotGraph`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)).

- `internalStandards`

  - `graph`: inclusion of internal standard network plot
    ([`plotGraph`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md)).

## Plot file caching

When a new report is generated the plot files are stored inside the
`report_files` sub-directory inside the destination path of the report.
The plot files are kept so they can be reused to speed-up re-creation of
reports (*e.g.* with different report settings). After the report is
generated, any unused plot files are removed unless they were recently
created (controlled by the `keepUnusedPlots` setting, see previous
section). The `clearPath` argument can be used to completely remove any
old files.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `report` to process HRMS (or IMS-HRMS) data. Please
see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## References

Creating MetFrag landing page URLs based on code from
[MetFamily](https://github.com/Treutler/MetFamily) R package.\
\
Xie Y (2014). “knitr: A Comprehensive Tool for Reproducible Research in
R.” In Stodden V, Leisch F, Peng RD (eds.), *Implementing Reproducible
Computational Research*. Chapman and Hall/CRC. ISBN 978-1466561595.\
\
Xie Y (2015). *Dynamic Documents with R and knitr*, 2nd edition. Chapman
and Hall/CRC, Boca Raton, Florida. ISBN 978-1498716963,
<https://yihui.org/knitr/>.\
\
Xie Y (2025). *knitr: A General-Purpose Package for Dynamic Report
Generation in R*. R package version 1.51, <https://yihui.org/knitr/>.
