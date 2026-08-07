# Plotting of grouped features

Various plotting functions for feature group data.

## Usage

``` r
plotMobilograms(obj, ...)

# S4 method for class 'featureGroups,missing'
plot(
  x,
  groupBy = NULL,
  onlyUnique = FALSE,
  retMin = FALSE,
  showLegend = TRUE,
  IMS = "maybe",
  col = NULL,
  pch = NULL,
  ...
)

# S4 method for class 'featureGroups'
plotInt(
  obj,
  average = FALSE,
  averageFunc = mean,
  areas = FALSE,
  normalized = FALSE,
  xBy = NULL,
  xNames = TRUE,
  groupBy = "fGroups",
  regression = FALSE,
  showLegend = FALSE,
  IMS = "maybe",
  pch = 20,
  type = if (regression) "p" else "b",
  lty = 3,
  xlim = NULL,
  ylim = NULL,
  col = NULL,
  plotArgs = NULL,
  linesArgs = NULL
)

# S4 method for class 'featureGroups'
plotChord(
  obj,
  addSelfLinks = FALSE,
  addRetMzPlots = TRUE,
  aggregate = FALSE,
  groupBy = NULL,
  addIntraOuterGroupLinks = FALSE,
  IMS = "maybe",
  ...
)

# S4 method for class 'featureGroups'
plotChroms(
  obj,
  analysis = analyses(obj),
  groupName = names(obj),
  retMin = FALSE,
  showPeakArea = FALSE,
  showFGroupRect = TRUE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  annotate = c("none", "ret", "mz", "mob"),
  intMax = "eic",
  EICParams = getDefEICParams(),
  showProgress = FALSE,
  IMS = "maybe",
  xlim = NULL,
  ylim = NULL,
  EICs = NULL,
  ...
)

# S4 method for class 'featureGroups'
plotChroms3D(
  obj,
  analysis = analyses(obj),
  groupName = names(obj),
  dim3 = "mz",
  retMin = FALSE,
  showLimits = TRUE,
  rtWindow = defaultLim("retention", "medium"),
  mzWindow = defaultLim("mz", "medium"),
  mobWindow = defaultLim("mobility", "medium"),
  gridSize = 50,
  title = NULL,
  ...
)

# S4 method for class 'featureGroupsSet'
plotChroms3D(obj, analysis = analyses(obj), ...)

# S4 method for class 'featureGroups'
plotMobilograms(
  obj,
  analysis = analyses(obj),
  groupName = names(obj),
  showPeakArea = FALSE,
  showFGroupRect = TRUE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  annotate = c("none", "ret", "mz", "mob"),
  EIMParams = getDefEIMParams(),
  showProgress = FALSE,
  xlim = NULL,
  ylim = NULL,
  EIMs = NULL,
  ...
)

# S4 method for class 'featureGroups'
plotVenn(obj, which = NULL, aggregate = TRUE, IMS = "maybe", ...)

# S4 method for class 'featureGroups'
plotUpSet(
  obj,
  which = NULL,
  aggregate = TRUE,
  IMS = "maybe",
  nsets = NULL,
  nintersects = NA,
  ...
)

# S4 method for class 'featureGroups'
plotVolcano(
  obj,
  FCParams,
  showLegend = TRUE,
  averageFunc = mean,
  normalized = FALSE,
  IMS = "maybe",
  col = NULL,
  pch = 19,
  ...
)

# S4 method for class 'featureGroups'
plotGraph(obj, onlyPresent = TRUE, width = NULL, height = NULL)

# S4 method for class 'featureGroupsSet'
plotGraph(obj, onlyPresent = TRUE, set, ...)

# S4 method for class 'featureGroups'
plotTICs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S4 method for class 'featureGroups'
plotBPCs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- obj, x:

  `featureGroups` object to be used for plotting.

- ...:

  passed to [`plot`](https://rdrr.io/r/base/plot.html) (`plot`,
  `plotChroms`, `plotMobilograms`, `plotTICs` and `plotBPCs`),
  [`filled.contour`](https://rdrr.io/r/graphics/filled.contour.html)
  (`plotChroms3D`),
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  functions (`plotVenn`), `chordDiagram` (`plotChord`),
  [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html) (`plotUpSet`).

- groupBy:

  Specifies how results are grouped in the plot. Should be a name of a
  column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  table which is used to make analysis groups (*e.g.* `"replicate"`).
  Set to `NULL` for no grouping. The `groupBy` argument can also be set
  to `"fGroups"` to group by feature groups (except `plotChord`).

  For `plotInt`: see the `Data aggregation in intensity plots` section
  for more details.

  For `plotChord`: the grouping is used to generate 'outer groups'.

- onlyUnique:

  If `TRUE` and `groupBy` is set to a column name of the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  then only feature groups that are unique to a replicate are plotted.

- retMin:

  Plot retention time in minutes (instead of seconds).

- showLegend:

  Plot a legend if `TRUE`.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  plotting in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  For `plotChroms`: if the `groupName` argument is set to a subset of
  the feature groups, then the `IMS` argument is forced to `"both"` to
  ensure that the selected feature groups are plotted.

- col:

  Colour(s) used. If `col=NULL` then colours are automatically
  generated.

- pch, type, lty:

  Common plotting parameters passed to *e.g.*
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For `plot`: if `pch=NULL` then values are automatically assigned.

  for `plotInt`: the `type` argument defaults to single points (`"p"`)
  when `regression=TRUE`, and lines+points (`"b"`) otherwise. For
  non-aggregated data, it probably makes more sense to set *e.g.*
  `type="p"` to avoid lines being drawn to points with equal x-values.

- average:

  Controls plot data averaging: see the
  `Data aggregation in intensity plots` section.

- averageFunc, normalized:

  Used for intensity data treatment, see the documentation for the
  [`as.data.table method`](https://rickhelmus.github.io/patRoon/reference/feature-table.md).

- areas:

  Set to `TRUE` to use feature areas instead of peak intensities for
  plotting.

- xBy:

  Controls x-value grouping in the plot: see the
  `Data aggregation in intensity plots` section.

- xNames:

  Plot names (`xNames=TRUE`) or a number sequence (`xNames=FALSE`) on
  the x axis.

- regression:

  If `TRUE` then a regression line is plotted for each group, using the
  x values set by `xBy` (see `Data aggregation in intensity plots`). if
  `showLegend=TRUE` then the R-squared value is show in the legend for
  each group (unless there are multiple feature groups and
  `groupBy != "fGroups"`).

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- plotArgs, linesArgs:

  A `list` with further arguments passed to
  [`plot`](https://rdrr.io/r/base/plot.html) and
  [`lines`](https://rdrr.io/r/graphics/lines.html), respectively.

- addSelfLinks:

  If `TRUE` then 'self-links' are added which represent non-shared data.

- addRetMzPlots:

  Set to `TRUE` to enable *m/z* *vs* retention time scatter plots.

- aggregate:

  Specifies how data should be aggregated prior to comparison. Set to
  `FALSE` to compare analyses, `TRUE` to compare replicates or to the
  name of a column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  to compare by a custom grouping of analyses.

- addIntraOuterGroupLinks:

  If `TRUE` then links will be added within outer groups.

- analysis, groupName:

  `character` vector with the analyses/group names to be considered for
  plotting.

  For `plotChroms` and `plotMobilograms`: Compared to subsetting the
  `featureGroups` object (`obj`) upfront this is slightly faster.
  Furthermore, if `onlyPresent=FALSE` in `EICParams` or `EIMParams`,
  this allows plotting chromatograms for feature groups where none of
  the specified analyses contain the feature (which is impossible
  otherwise since subsetting leads to removal of 'empty' feature
  groups). For IMS workflows: if `IMS!="both"` then the `analysis` and
  `groupName` arguments are adjusted for the remaining data after IMS
  selection.

  For `plotChroms3D`: a single analysis and group name should be
  specified (unless `obj` contains only a single feature, *e.g.* created
  as `fGroups[1, 1]`).

- showPeakArea:

  Set to `TRUE` to display integrated peak ranges by filling (shading)
  their areas.

- showFGroupRect:

  Set to `TRUE` to mark the full integration/intensity range of all
  features within a feature group by drawing a rectangle around it.

- title:

  Character string used for title of the plot. If `NULL` a title will be
  automatically generated.

- annotate:

  Set to `"ret"`, `"mz"` and/or `"mob"` to annotate peaks with the
  retention time, *m/z* and/or ion mobility (if available) feature group
  values, respectively.

- intMax:

  Method used to determine the maximum intensity plot limit. Should be
  `"eic"` (from EIC data) or `"feature"` (from feature data). Ignored if
  the `ylim` parameter is specified.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- showProgress:

  if set to `TRUE` then a text progressbar will be displayed when all
  EICs are being plot. Set to `"none"` to disable any annotation.

- EICs, EIMs:

  Internal parameter for now and should be kept at `NULL` (default).

- dim3:

  The third dimension to plot besides retention time and intensity. Can
  be either `"mz"` or `"mobility"`.

- showLimits:

  If `TRUE`, a rectangle is drawn to indicate the feature limits.

- rtWindow, mzWindow, mobWindow:

  Numeric values specifying the window size around the feature for
  retention time, m/z, and ion mobility respectively. Values `>0` will
  effectively zoom out.

- gridSize:

  The size of the grid for interpolation.

- EIMParams:

  A named `list` with parameters used for extracted ion mobilogram (EIM)
  creation. See the [EIM
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- which:

  A character vector with the selection to compare (*e.g.* replicates,
  as set by the `aggregate` argument). Set to `NULL` to select
  everything.

- nsets, nintersects:

  See [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html). If
  `nsets == NULL` then it will be set to the number of compared items.

- FCParams:

  A parameter list to calculate Fold change data. See `getFCParams` for
  more details.

- onlyPresent:

  Only plot feature groups of internal standards that are still present
  in the `featureGroups` input object (which may be otherwise be removed
  by *e.g.* subsetting or
  [`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md)).

- width, height:

  Passed to `visNetwork`.

- set:

  (**sets workflow**) The set for which data must be plotted.

- retentionRange:

  Range of retention time (in seconds) to collect TIC traces. Should be
  a numeric vector with length of two containing the min/max values. Set
  to NULL to ignore.

- MSLevel:

  Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to
  obtain TIC traces.

## Value

`plotChroms3D` returns the interpolated grid used for plotting
(generated with
[mba.surf](https://finleya.github.io/MBA/reference/mba.surf.html)).

`plotVenn` (invisibly) returns a list with the following fields:

- `gList` the `gList` object that was returned by the utilized
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  function.

- `areas` The total area for each plotted group.

- `intersectionCounts` The number of intersections between groups.

The order for the `areas` and `intersectionCounts` fields is the same as
the parameter order from the used plotting function (see *e.g.*
`draw.pairwise.venn` and `draw.triple.venn`).

`plotGraph` returns the result of `visNetwork`.

## Details

`plot` Generates an *m/z* *vs* retention time plot for all featue
groups. Optionally highlights unique/overlapping presence amongst
replicates.

`plotInt` Generates a line plot with feature intensities.

`plotChord` Generates a chord diagram which can be used to visualize
shared presence of feature groups between analyses or replicate groups.
In addition, analyses/replicates sharing similar properties (*e.g.*
location, age, type) may be grouped to enhance visualization between
these 'outer groups'.

`plotChroms` Plots extracted ion chromatograms (EICs) of feature groups.

`plotChroms3D` generates a 3D plot of chromatographic data for a single
feature group in a single analysis. The plot shows retention time on the
x-axis, m/z or mobility on the y-axis, and intensity as color in a
contour plot. The plot is made with the
[`filled.contour`](https://rdrr.io/r/graphics/filled.contour.html)
function.

`plotMobilograms` Plots extracted ion mobilograms (EIMs) of feature
groups.

`plotVenn` plots a Venn diagram (using
[VennDiagram](https://CRAN.R-project.org/package=VennDiagram)) outlining
unique and shared feature groups between up to five replicates.

`plotUpSet` plots an UpSet diagram (using the
[`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html) function) outlining
unique and shared feature groups between given replicates.

`plotVolcano` Plots Fold change data in a 'Volcano plot'.

`plotGraph` generates an interactive network plot which is used to
explore internal standard (IS) assignments to each feature group. This
requires the availability of IS assignments, see the documentation for
[`normInts`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
for details. The graph is rendered with visNetwork.

`plotTICs` Plots the total ion chromatogram/s (TICs) of the analyses.

`plotTICs` Plots the base peak chromatogram/s (BPCs) of the analyses.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `plotChroms`, `plotMobilograms`, `plotTICs` and
`plotBPCs` to process HRMS (or IMS-HRMS) data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## Sets workflows

The following methods are changed or with new functionality:

- `plotGraph` only plots data per set, and requires the `set` argument
  to be set.

In sets workflows the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
contains an additional `"set"` column, which can be used for arguments
that involve grouping of analyses. For instance, if `groupBy="set"` then
plotting data is grouped per set.

## Data aggregation in intensity plots

the `average`, `xBy` and `groupBy` arguments control how data is
aggregated in intensity plots:

- `average`: controls the averaging of feature intensities prior to
  plotting.

- `xBy`: can map the x value of individual points to analysis metadata.
  For example, exposure time or sample location. Non-numeric values are
  allowed (unless `regression=TRUE`).

- `groupBy`: controls the grouping of points in the plot. Equal groups
  are plotted in sequence so they can be connected with lines and are
  coloured equally. Examples include experiment type or feature groups.

The following values are valid:

- `FALSE` (`average`) or `NULL` (`xBy` and `groupBy`): aggregation is
  disabled.

- `TRUE` (only `average`): results are averaged for each replicate

- a name of a column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md):
  results are aggregated for analyses with the same table column value.

- `"fGroups"` (only `groupBy`): plots are grouped by feature groups.

## References

Gu Z, Gu L, Eils R, Schlesner M, Brors B (2014). “circlize implements
and enhances circular visualization in R.” *Bioinformatics*, **30**,
2811-2812.

Conway JR, Lex A, Gehlenborg N (2017). “UpSetR: an R package for the
visualization of intersecting sets and their properties.”
*Bioinformatics*, **33**(18), 2938-2940.
[doi:10.1093/bioinformatics/btx364](https://doi.org/10.1093/bioinformatics/btx364)
. <http://dx.doi.org/10.1093/bioinformatics/btx364>.\
\
Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H (2014). “UpSet:
Visualization of Intersecting Sets.” *IEEE Transactions on Visualization
and Computer Graphics*, **20**(12), 1983–1992.
[doi:10.1109/tvcg.2014.2346248](https://doi.org/10.1109/tvcg.2014.2346248)
.

## See also

[`featureGroups-class`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md),
[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)

## Author

Rick Helmus \<<r.helmus@uva.nl>\> and Ricardo Cunha \<<cunha@iuta.de>\>
(`plotTICs` and `plotBPCs` functions)
