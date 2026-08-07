# Base features class

Holds information for all features present within a set of analysis.

## Usage

``` r
# S4 method for class 'features'
length(x)

# S4 method for class 'features'
show(object)

# S4 method for class 'features'
featureTable(obj)

# S4 method for class 'features'
analysisInfo(obj, df = FALSE)

# S4 method for class 'features'
getFeatureQualityNames(obj, scores = FALSE, totScore = TRUE)

# S4 method for class 'features'
analysisInfo(obj) <- value

# S4 method for class 'features'
analyses(obj)

# S4 method for class 'features'
replicates(obj)

# S4 method for class 'features'
hasIMS(obj)

# S4 method for class 'features'
fromIMS(obj)

# S4 method for class 'features'
as.data.table(x)

# S4 method for class 'features'
filter(
  obj,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  retentionRange = NULL,
  mzRange = NULL,
  mzDefectRange = NULL,
  chromWidthRange = NULL,
  IMSRangeParams = NULL,
  qualityRange = NULL,
  negate = FALSE
)

# S4 method for class 'features,ANY,missing,missing'
x[i, j, ..., ni, reorder = FALSE, drop = TRUE]

# S4 method for class 'features,ANY,missing'
x[[i]]

# S4 method for class 'features'
x$name

# S4 method for class 'features'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'features'
calculatePeakQualities(
  obj,
  weights,
  flatnessFactor,
  featureQualities = NULL,
  EICParams = getDefEICParams(window = 0),
  parallel = TRUE
)

# S4 method for class 'features'
getTICs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'features'
getBPCs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'features'
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

# S4 method for class 'features'
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

# S4 method for class 'featuresSet'
sets(obj)

# S4 method for class 'featuresSet'
show(object)

# S4 method for class 'featuresSet'
as.data.table(x)

# S4 method for class 'featuresSet,ANY,missing,missing'
x[i, ..., sets = NULL, reorder = FALSE, drop = TRUE]

# S4 method for class 'featuresSet'
filter(obj, ..., negate = FALSE, sets = NULL)

# S4 method for class 'featuresSet'
unset(obj, set)

# S4 method for class 'featuresKPIC2'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featuresPiek'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featuresXCMS'
analysisInfo(obj) <- value

# S4 method for class 'featuresXCMS'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featuresXCMS3'
delete(obj, i = NULL, j = NULL, ...)
```

## Arguments

- obj, x, object:

  `features` object to be accessed

- df:

  If `TRUE` then the returned value is a `data.frame`, otherwise a
  `data.table`.

- scores:

  If `TRUE` the score names are returned, otherwise the quality names.

- totScore:

  If `TRUE` (and `scores=TRUE`) then the name of the total score is
  included.

- value:

  A `data.frame` or `data.table` with the new analysis information.

- absMinIntensity, relMinIntensity:

  Minimum absolute/relative intensity for features to be kept. The
  relative intensity is determined from the feature with highest
  intensity (within the same analysis). Set to `0` or `NULL` to skip
  this step.

- retentionRange, mzRange, mzDefectRange, chromWidthRange:

  Range of retention time (in seconds), *m/z*, mass defect (defined as
  the decimal part of *m/z* values) or chromatographic peak width (in
  seconds), respectively. Features outside this range will be removed.
  Should be a numeric vector with length of two containing the min/max
  values. The maximum can be `Inf` to specify no maximum range. Set to
  `NULL` to skip this step.

- IMSRangeParams:

  (**IMS workflow**) A `list` with parameters to be used for filtering
  IMS range data. See
  [`getIMSRangeParams`](https://rickhelmus.github.io/patRoon/reference/getIMSRangeParams.md)
  for details and how to make such a parameter list.

- qualityRange:

  Used to filter features by their peak qualities/scores (see
  [`calculatePeakQualities`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
  Should be a named `list` with min/max ranges for each quality/score to
  be filtered (the
  [`getFeatureQualityNames`](https://rickhelmus.github.io/patRoon/reference/generics.md)
  function can be used to obtain valid names). Example:
  `qualityRange=list(ModalityScore=c(0.3, Inf), SymmetryScore=c(0.5, Inf))`.
  Set to `NULL` to ignore.

- negate:

  If set to `TRUE` then filtering operations are performed in opposite
  manner.

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  analyses by their index or name, respectively (for the order/names see
  [`analyses()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all analyses are selected.\
  \
  For `[[`: should be a scalar value.\
  \
  For `delete`: The data to remove from. `i` are the analyses as numeric
  index, logical or character, `j` the features as numeric index (row)
  of the feature. If either is `NULL` then data for all is removed. `j`
  may also be a function: it will be called for each analysis, with the
  feature table (a `data.table`), the analysis name and any other
  arguments passed as `...` to `delete`. The return value of this
  function specifies the feature indices (rows) to be removed (specified
  as an `integer` or `logical` vector).

- ...:

  For `delete`: passed to the function specified as `j`.

  For `plotTICs` and `plotBPCs`: further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `features` method.

- ni:

  Optional argument. An expression used for subsetting the analyses. The
  [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  is first subset and the remaining rows are used to determine for which
  analyses the results should be kept. The unevaluated `ni` expression
  is used to set the `i` argument of the subset operator of
  [data.table](https://CRAN.R-project.org/package=data.table), which
  therefore brings the advanced subsetting capabilities of
  [data.table](https://CRAN.R-project.org/package=data.table) (see the
  [`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)
  documentation for more details). For instance,
  `fList[replicate == "standard"]` would subset all analyses assigned
  with the replicate `"standard"`.

- reorder:

  If `TRUE` then the order of the analyses is changed to match the order
  of the `i` argument.

  (**sets workflow**) If the `sets` argument is specified (and `i` is
  not) then the order of sets is changed instead.

- drop:

  Ignored.

- name:

  The analysis name (partially matched).

- weights:

  A named `numeric` vector that defines the weight for each score to
  calculate the `totalScore`. The names of the vector follow the score
  names. Unspecified weights are defaulted to `1`. Example:
  `weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)`.

- flatnessFactor:

  Passed to MetaClean as the `flatness.factor` argument to
  [`calculateJaggedness`](https://rdrr.io/pkg/MetaClean/man/calculateJaggedness.html)
  and
  [`calculateModality`](https://rdrr.io/pkg/MetaClean/man/calculateModality.html).

- featureQualities:

  Specifies which feature qualities to calculate. Can be `NULL`
  (default, calculates all qualities), a `character` vector with names
  of qualities to calculate (e.g., `c("FWHM2Base", "Symmetry")`), or a
  `list` of custom quality definitions. See the
  [`featureQualities`](https://rickhelmus.github.io/patRoon/reference/feature-quality.md)
  function for more details.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

- MSLevel:

  Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to
  obtain traces.

- retMin:

  Plot retention time in minutes (instead of seconds).

- title:

  Character string used for title of the plot. If `NULL` a title will be
  automatically generated.

- groupBy:

  Specifies how results are grouped in the plot. Should be a name of a
  column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  table which is used to make analysis groups (*e.g.* `"replicate"`), or
  `"fGroups"` to group by feature group. Set to `NULL` for no grouping.

- showLegend:

  Plot a legend if TRUE.

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- sets:

  (**sets workflow**) For `[` and `filter`: a `character` with name(s)
  of the sets to keep (or remove if `negate=TRUE`).

- set:

  (**sets workflow**) The name of the set.

## Value

`featureTable`: A `list` containing a `data.table` for each analysis
with feature data

`analysisInfo`: The [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
of this `features` object.

`delete` returns the object for which the specified data was removed.

`calculatePeakQualities` returns a modified object amended with peak
qualities and scores.

## Details

This class provides a way to store intensity, retention times, *m/z* and
other data for all features in a set of analyses. The class is `virtual`
and derived objects are created by 'feature finders' such as
`findFeaturesOpenMS`, `findFeaturesXCMS` and `findFeaturesBruker`.

## Methods (by generic)

- `length(features)`: Obtain total number of features.

- `show(features)`: Shows summary information for this object.

- `featureTable(features)`: Get table with feature information

- `analysisInfo(features)`: Get analysis information

- `getFeatureQualityNames(features)`: Returns the present
  chromatographic peak quality and score names for features.

- `analysisInfo(features) <- value`: Modifies analysis information

- `analyses(features)`: returns a `character` vector with the names of
  the analyses for which data is present in this object.

- `replicates(features)`: returns a `character` vector with the names of
  the replicates for which data is present in this object.

- `hasIMS(features)`: Returns `TRUE` if the features object has mobility
  information.

- `fromIMS(features)`: Returns `TRUE` if the features object was
  directly created from IMS data.

- `as.data.table(features)`: Returns all feature data in a table.

- `filter(features)`: Performs common rule based filtering of features.
  Note that this (and much more) functionality is also provided by the
  `filter` method defined for
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).
  However, filtering a `features` object may be useful to avoid grouping
  large amounts of features.

- `x[i`: Subset on analyses.

- `x[[i`: Extract a feature table for an analysis.

- `$`: Extract a feature table for an analysis.

- `delete(features)`: Completely deletes specified features.

- `calculatePeakQualities(features)`: Calculates peak qualities for each
  feature. Please see the
  [`featureQualities`](https://rickhelmus.github.io/patRoon/reference/feature-quality.md)
  function and MetaClean publication (referenced below) for more
  details. For each metric, an additional score is calculated by
  normalizing all feature values (unless the quality metric definition
  has a fixed range) and scale from `0` (worst) to `1` (best). Then, a
  `totalScore` for each feature is calculated by the (weighted) sum of
  all score values.

- `getTICs(features)`: Obtain the total ion chromatogram/s (TICs) of the
  analyses.

- `getBPCs(features)`: Obtain the base peak chromatogram/s (BPCs) of the
  analyses.

- `plotTICs(features)`: Plots the TICs of the analyses.

- `plotBPCs(features)`: Plots the BPCs of the analyses.

## Slots

- `features`:

  List of features per analysis file. Use the `featureTable` method for
  access.

- `analysisInfo`:

  A `data.table` with the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).
  Use the `analysisInfo` method for access.

- `featureQualityNames`:

  Character vector with the names of the chromatographic peak quality
  metrics that are present.

- `hasIMS`:

  A `logical` that is `TRUE` if the features object contain mobility/CCS
  information. Use the `hasIMS` method for access.

- `fromIMS`:

  A `logical` that is `TRUE` if the features object was directly created
  from IMS data (*i.e.* direct mobility assignment workflow). Use the
  `fromIMS` method for access.

## Note

For `calculatePeakQualities`: sometimes MetaClean may return `NA` for
the *Gaussian Similarity* and *Symmetry* metrics, in which case it will
be set to `0`.

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`features`**

    - `featuresSet`

    - `featuresUnset`

    - [`featuresFromFeatGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - [`featuresConsensus`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - `featuresBruker`

    - `featuresEnviPick`

    - `featuresKPIC2`

    - `featuresOpenMS`

    - `featuresPiek`

    - `featuresSAFD`

    - `featuresSIRIUS`

    - `featuresTable`

    - [`featuresBrukerTASQ`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerTASQ.md)

    - `featuresXCMS`

    - `featuresXCMS3`

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `calculatePeakQualities` and TIC/BPC related
functions to process HRMS (or IMS-HRMS) data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## Sets workflows

The `featuresSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `features` and therefore largely follows the
same user interface.

The following methods are specifically defined for sets workflows:

- `sets` Returns the set names for this object.

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`featuresUnset`), which allows it to be used in 'regular'
  workflows. The adduct annotations for the selected set (*e.g.* as
  passed to `makeSet`) are used to convert all feature masses to ionic
  *m/z* values.

The following methods are changed or with new functionality:

- `filter` and the subset operator (`[`) have specific arguments to
  choose/filter by (feature presence in) sets. See the `sets` argument
  description.

- **Important**: the `mzRange`, `mzDefectRange` and `IMSRangeParams`
  filters use neutral feature masses, whereas non-sets workflows use
  *m/z* values. Hence, adjust accordingly to avoid (slightly) different
  results!

## References

Chetnik K, Petrick L, Pandey G (2020). “MetaClean: a machine
learning-based classifier for reduced false positive peak detection in
untargeted LC-MS metabolomics data.” *Metabolomics*, **16**(11).
[doi:10.1007/s11306-020-01738-3](https://doi.org/10.1007/s11306-020-01738-3)
.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)

## Author

Rick Helmus \<<r.helmus@uva.nl>\> and Ricardo Cunha \<<cunha@iuta.de>\>
(`getTICs`, `getBPCs`, `plotTICs` and `plotBPCs` functions)
