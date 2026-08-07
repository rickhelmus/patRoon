# Feature group to table conversion

Convert feature group data to a `data.table` (or `data.frame`).

## Usage

``` r
# S4 method for class 'featureGroups'
as.data.table(
  x,
  average = FALSE,
  areas = FALSE,
  features = FALSE,
  qualities = FALSE,
  regression = FALSE,
  regressionBy = NULL,
  averageFunc = mean,
  normalized = FALSE,
  FCParams = NULL,
  concAggrParams = getDefPredAggrParams(),
  toxAggrParams = getDefPredAggrParams(),
  normConcToTox = FALSE,
  anaInfoCols = NULL,
  IMS = "both"
)

# S4 method for class 'featureGroupsScreening'
as.data.table(x, ..., collapseSuspects = ",", onlyHits = FALSE)

# S4 method for class 'featureGroupsScreeningSet'
as.data.table(x, ..., collapseSuspects = ",", onlyHits = FALSE)
```

## Arguments

- x:

  The `featureGroups` like object to be exported as table.

- average:

  Controls the averaging of feature intensities. Averaging also
  influences the calculation of regression parameters. Set to `FALSE` to
  disable averaging, `TRUE` to average per replicate or to the name of a
  column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  to compare by a custom grouping of analyses.

  If `features=TRUE` all numerical properties are averaged and
  non-numericals are collapsed. Each row represents the feature
  aggregated data and the groups used for aggregation (*e.g.*
  replicates) are stored in the `average_group` column. It is also
  possible to average these data per feature group, by setting
  `average="fGroups"`.

- areas:

  If set to `TRUE` then areas are output instead of peak intensities.
  Ignored if `features=TRUE`, as areas of features are always reported.

- features:

  If `TRUE` then feature specific data will be added. Also see the
  `average` argument.

- qualities:

  Adds feature (group) qualities (`qualities="quality"`), scores
  (`qualities="score"`) or both (`qualities="both"`), if this data is
  available (*i.e.* from `calculatePeakQualities`). If `qualities=FALSE`
  then nothing is reported.

- regression, regressionBy:

  Used for regression calculations. See the `Regression calculation`
  section below. Set to `NULL` to ignore.

- averageFunc:

  Function used for averaging. Only used when data is averaged or
  `FCParams != NULL`.

- normalized:

  If `TRUE` then normalized intensities are used (see the
  `Feature intensity normalization` section). If no normalization data
  is available (*e.g.* because `normInts` was not used) then an
  automatic group normalization is performed.

- FCParams:

  A parameter list to calculate fold change data. See `getFCParams` for
  more details. Set to `NULL` to not perform FC calculations. **NOTE**:
  the intensity data is always averaged to calculate fold changes (using
  `averageFunc`) by replicates, irrespective of the `average` argument.

- concAggrParams, toxAggrParams:

  Parameters to aggregate calculated concentrations/toxicities (obtained
  with
  [`calculateConcs`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)/[`calculateTox`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md)).
  See [prediction aggregation
  parameters](https://rickhelmus.github.io/patRoon/reference/pred-aggr-params.md)
  for more information. Set to `NULL` to omit this data.

- normConcToTox:

  Set to `TRUE` to normalize concentrations to toxicities. Only relevant
  if this data is present (see
  [`calculateConcs`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)/[`calculateTox`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md)).

- anaInfoCols:

  A `character` with any additional columns from the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  table. Only supported if `features=TRUE`. If averaging is performed,
  then numeric data is averaged and character data collapsed. Set to
  `NULL` to ignore.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered to be
  returned in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

- ...:

  Passed to the parent `as.data.table` method.

- collapseSuspects:

  If a `character` then any suspects that were matched to the same
  feature group are collapsed to a single row and suspect names are
  separated by the value of `collapseSuspects`. If `NULL` then no
  collapsing occurs, and each suspect match is reported on a single row.
  See the `Suspect collapsing` section below for additional details.

- onlyHits:

  If `TRUE` then only feature groups with suspect hits are reported.

## Details

The `as.data.table` generic function converts most feature group data to
a highly customizable `data.table`. If a classical `data.frame` is
preferred, the `as.data.frame` generic function can be used instead and
accepts the exact same arguments. The methods defined for [suspect
screening
workflows](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
will merge the information from `screenInfo`, such as suspect names and
other properties and annotation data.

## Regression calculation

The `regression` argument controls the calculation of regression
parameters from a regression model calculated with feature intensities
(or areas if `areas=TRUE`). Here, simple linear regression is used,
*i.e.* `y=ax+b` with `a` the slope and `b` the intercept. The value for
`regression` should be the name of a column in the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
table with numerical data to be used for x-values. Alternatively, if
`regression=TRUE` then the `"conc"` column is used. Any `NA` x-values
are ignored, and no regression will be calculated if less than two
(non-NA) x-values are available. The output table will contain
properties such as the slope and correlation coefficient (R-squared).
Furthermore, if `features=TRUE` then x-values will be calculated from
the model and stored in the `x_reg` column.

The `regressionBy` argument can be used to construct separate regression
models for different groups of analysis. It should be set to the name of
a column in the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
table which defines the grouping between samples. If `features=TRUE`
then the grouping is stored in the `regression_group` column of the
output table.

Please see the handbook for examples on how to use the regression
functionality.

## Suspect collapsing

The `as.data.table` method for `featureGroupsScreening` supports an
additional format where each suspect hit is reported on a separate row
(enabled by setting `collapseSuspects=NULL`). In this format the suspect
properties from the `screenInfo` method are merged with each suspect
row. Alternatively, if *suspect collapsing* is enabled (the default)
then the regular `as.data.table` format is used, and amended with the
names and estimated ID levels (if available) of the suspects matched to
a feature group (each separated by the value of the `collapseSuspects`
argument).

Suspect collapsing also influences the reporting of predicted [feature
concentrations](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)
and
[toxicities](https://rickhelmus.github.io/patRoon/reference/pred-tox.md).
In the case that (1) suspects are *not* collapsed in the output table
and (2) predictions are available for a specific suspect hit (*i.e.* if
[`predictRespFactors`](https://rickhelmus.github.io/patRoon/reference/generics.md)
or
[`predictTox`](https://rickhelmus.github.io/patRoon/reference/generics.md)
was called on the feature groups object), then only the suspect specific
data is reported and no aggregation is performed. Hence, this allows you
to obtain specific concentration/toxicity values for each
suspect/feature group pair.

## IMS workflows

If the `IMS` argument is set to `"both"` or `"maybe"` then
`"mobility_collapsed"` and `"CCS_collapsed"` columns will be added that
summarize all mobility/CCS values of the IMS features (or IMS feature
groups) assigned to this IMS precursor. These numbers are currently
rounded to `3` decimals.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
normalization of feature intensities occur per set.

In sets workflows the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
contains an additional `"set"` column, which can be used for arguments
that involve grouping of analyses. For instance, if `regressionBy="set"`
then regression models will be calculated for each set.
