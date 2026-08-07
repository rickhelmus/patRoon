# Formula annotations class

Contains data of generated chemical formulae for given feature groups.

## Usage

``` r
# S4 method for class 'formulas'
annotations(obj, features = FALSE)

# S4 method for class 'formulas'
analyses(obj)

# S4 method for class 'formulas'
defaultExclNormScores(obj)

# S4 method for class 'formulas'
show(object)

# S4 method for class 'formulas,ANY,ANY'
x[[i, j]]

# S4 method for class 'formulas'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'formulas'
as.data.table(
  x,
  fGroups = NULL,
  fragments = FALSE,
  countElements = NULL,
  countFragElements = NULL,
  OM = FALSE,
  normalizeScores = "none",
  excludeNormScores = defaultExclNormScores(x),
  average = FALSE
)

# S4 method for class 'formulas'
annotatedPeakList(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  onlyAnnotated = FALSE
)

# S4 method for class 'formulas'
plotSpectrum(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  showLegend = TRUE,
  ...
)

# S4 method for class 'formulas'
plotScores(
  obj,
  index,
  groupName,
  analysis = NULL,
  normalizeScores = "max",
  excludeNormScores = defaultExclNormScores(obj)
)

# S4 method for class 'formulas'
consensus(
  obj,
  ...,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL
)

# S4 method for class 'formulasSet'
show(object)

# S4 method for class 'formulasSet'
delete(obj, i, j, ...)

# S4 method for class 'formulasSet,ANY,missing,missing'
x[i, j, ..., sets = NULL, updateConsensus = FALSE, drop = TRUE]

# S4 method for class 'formulasSet'
filter(obj, ..., sets = NULL, updateConsensus = FALSE, negate = FALSE)

# S4 method for class 'formulasSet'
plotSpectrum(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  showLegend = TRUE,
  perSet = TRUE,
  mirror = TRUE,
  ...
)

# S4 method for class 'formulasSet'
annotatedPeakList(obj, index, groupName, analysis = NULL, MSPeakLists, ...)

# S4 method for class 'formulasSet'
consensus(
  obj,
  ...,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL,
  filterSets = FALSE,
  setThreshold = 0,
  setThresholdAnn = 0,
  setAvgSpecificScores = FALSE
)

# S4 method for class 'formulasSet'
unset(obj, set)

# S4 method for class 'formulasConsensusSet'
unset(obj, set)

# S4 method for class 'formulasSIRIUS'
delete(obj, i = NULL, j = NULL, ...)
```

## Arguments

- obj, x, object:

  The `formulas` object.

- features:

  If `TRUE` returns formula data for features, otherwise for feature
  groups.

- i, j:

  For `[[`: If both `i` and `j` are specified then `i` specifies the
  analysis and `j` the feature group of the feature for which
  annotations should be returned. Otherwise `i` specifies the feature
  group for which group annotations should be returned. `i`/`j` can be
  specified as `integer` index or as a `character` name.

  Otherwise passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- ...:

  For `plotSpectrum`: Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For `delete`: passed to the function specified as `j`.

  For `consensus`: Any further (and unique) `formulas` objects.

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `formulas` method.

- fGroups, fragments, countElements, countFragElements, OM:

  Passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- normalizeScores:

  A `character` that specifies how normalization of annotation scorings
  occurs. Either `"none"` (no normalization), `"max"` (normalize to max
  value) or `"minmax"` (perform min-max normalization). Note that
  normalization of negative scores (e.g. output by `SIRIUS`) is always
  performed as min-max. Furthermore, currently normalization for
  `compounds` takes the original min/max scoring values into account
  when candidates were generated. Thus, for `compounds` scoring,
  normalization is not affected when candidate results were removed
  after they were generated (*e.g.* by use of `filter`).

- excludeNormScores:

  A `character` vector specifying any compound scoring names that should
  *not* be normalized. Set to `NULL` to normalize all scorings. Note
  that whether any normalization occurs is set by the
  `excludeNormScores` argument.

  For `compounds`: By default `score` and `individualMoNAScore` are set
  to mimic the behavior of the `MetFrag` web interface.

- average:

  If set to `TRUE` an 'average formula' is generated for each feature
  group by combining all elements from all candidates and averaging
  their amounts. This obviously leads to non-existing formulae, however,
  this data may be useful to deal with multiple candidate formulae per
  feature group when performing elemental characterization. Setting this
  to `TRUE` disables reporting of most other data.

- index:

  The candidate index (row). For `plotSpectrum` two indices can be
  specified to compare spectra. In this case `groupName` and `analysis`
  (if not `NULL`) should specify values for the spectra to compare.

- groupName:

  The name of the feature group for which a plot should be made. To
  compare spectra, two group names can be specified.

- analysis:

  The name of the analysis for which a plot should be made. If `NULL`
  then data from the feature group averaged peak list is used. When
  comparing spectra, either `NULL` or the analyses for both spectra
  should be specified.

- MSPeakLists:

  The
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object with relevant spectral data.

- onlyAnnotated:

  Set to `TRUE` to filter out any peaks that could not be annotated.

- title:

  The title of the plot. If `NULL` a title will be automatically made.

- normalized:

  Controls intensity normalization. Should be `FALSE` (don't normalize),
  `TRUE` (normalize) or `"multiple"` (only normalizes if multiple
  spectra are plotted).

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- mincex:

  The formula annotation labels are automatically scaled. The `mincex`
  argument forces a minimum `cex` value for readability.

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- showLegend:

  Set to `TRUE` to show a legend.

- absMinAbundance, relMinAbundance:

  Minimum absolute or relative (`0-1`) abundance across objects for a
  result to be kept. For instance, `relMinAbundance=0.5` means that a
  result should be present in at least half of the number of compared
  objects. Set to `NULL` to ignore and keep all results. Limits cannot
  be set when `uniqueFrom` is not `NULL`.

- uniqueFrom:

  Set this argument to only retain formulas that are unique within one
  or more of the objects for which the consensus is made. Selection is
  done by setting the value of `uniqueFrom` to a `logical` (values are
  recycled), `numeric` (select by index) or a `character` (as obtained
  with `algorithm(obj)`). For `logical` and `numeric` values the order
  corresponds to the order of the objects given for the consensus. Set
  to `NULL` to ignore.

- uniqueOuter:

  If `uniqueFrom` is not `NULL` and if `uniqueOuter=TRUE`: only retain
  data that are also unique between objects specified in `uniqueFrom`.

- rankWeights:

  A numeric vector with weights of to calculate the mean ranking score
  for each candidate. The value will be re-cycled if necessary, hence,
  the default value of `1` means equal weights for all considered
  objects.

- labels:

  A `character` with names to use for labelling. If `NULL` labels are
  automatically generated.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`). Note: if `updateConsensus=FALSE` then the
  `setCoverage` column of the annotation results is not updated.

- updateConsensus:

  (**sets workflow**) If `TRUE` then the annonation consensus among set
  results is updated. See the `Sets workflows` section for more details.

- drop:

  Passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- negate:

  Passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- perSet, mirror:

  (**sets workflow**) If `perSet=TRUE` then the set specific mass peaks
  are annotated separately. Furthermore, if `mirror=TRUE` (and there are
  two sets in the object) then a mirror plot is generated.

- filterSets:

  (**sets workflow**) Controls how algorithms concensus abundance
  filters are applied. See the `Sets workflows` section below.

- setThreshold, setThresholdAnn:

  (**sets workflow**) Thresholds used to create the annotation set
  consensus. See
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md).

- setAvgSpecificScores:

  (**sets workflow**) If `TRUE` then set specific annotation scores
  (*e.g.* MS/MS and isotopic pattern match scores) are averaged for the
  set consensus. See
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md).

- set:

  (**sets workflow**) The name of the set.

## Value

`annotations` returns a `list` containing for each feature group (or
feature if `features=TRUE`) a `data.table` with an overview of all
generated formulae and other data such as candidate scoring and MS/MS
fragments.

`consensus` returns a `formulas` object that is produced by merging
results from multiple `formulas` objects.

## Details

`formulas` objects are obtained with
[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md).
This class is derived from the
[`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
class, please see its documentation for more methods and other details.

## Methods (by generic)

- `annotations(formulas)`: Accessor method to obtain generated formulae.

- `analyses(formulas)`: returns a `character` vector with the names of
  the analyses for which data is present in this object.

- `defaultExclNormScores(formulas)`: Returns default scorings that are
  excluded from normalization.

- `show(formulas)`: Show summary information for this object.

- `x[[i`: Extracts a formula table, either for a feature group or for
  features in an analysis.

- `as.data.table(formulas)`: Generates a table with all candidate
  formulae for each feature group and other information such as element
  counts.

- `annotatedPeakList(formulas)`: Returns an MS/MS peak list annotated
  with data from a given candidate formula.

- `plotSpectrum(formulas)`: Plots an annotated spectrum for a given
  candidate formula of a feature or feature group. Two spectra can be
  compared by specifying a two-sized vector for the `index`, `groupName`
  and (if desired) `analysis` arguments.

- `plotScores(formulas)`: Plots a barplot with scoring of a candidate
  formula.

- `consensus(formulas)`: Generates a consensus of results from multiple
  objects. In order to rank the consensus candidates, first each of the
  candidates are scored based on their original ranking (the scores are
  normalized and the highest ranked candidate gets value `1`). The
  (weighted) mean is then calculated for all scorings of each candidate
  to derive the final ranking (if an object lacks the candidate its
  score will be `0`). The original rankings for each object is stored in
  the `rank` columns.

## Slots

- `featureFormulas`:

  A `list` with all generated formulae for each analysis/feature group.
  Use the `annotations` method for access.

- `setThreshold,setThresholdAnn,setAvgSpecificScores`:

  (**sets workflow**) A copy of the equally named arguments that were
  passed when this object was created by
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md).

- `origFGNames`:

  (**sets workflow**) The original (order of) names of the
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that was used to create this object.

- `MS2QuantMeta`:

  (**sets workflow**) A named `list` with for each set the metadata from
  MS2Quant filled in by `predictRespFactors`.

## S4 class hierarchy

- [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)

  - **`formulas`**

    - `formulasConsensus`

    - `formulasSet`

      - `formulasConsensusSet`

    - `formulasUnset`

    - [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)

## Source

Subscripting of formulae for plots generated by `plotSpectrum` is based
on the `chemistry2expression` function from the
[ReSOLUTION](https://github.com/schymane/ReSOLUTION) package.

## Sets workflows

The `formulasSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `formulas` and therefore largely follows the
same user interface.

The following methods are specifically defined for sets workflows:

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`formulasUnset`), which allows it to be used in 'regular'
  workflows. Only the annotation results that are present in the
  specified set are kept (based on the set consensus, see below for
  implications).

The following methods are changed or with new functionality:

- `filter` and the subset operator (`[`) Can be used to select data that
  is only present for selected sets. Depending on the `updateConsenus`,
  both either operate on set consensus or original data (see below for
  implications).

- `annotatedPeakList` Returns a combined annotation table with all sets.

- `plotSpectrum` Is able to highlight set specific mass peaks (`perSet`
  and `mirror` arguments).

- `consensus` Creates the algorithm consensus based on the original
  annotation data (see below for implications). Then, like the sets
  workflow method for
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md),
  a consensus is made for all sets, which can be controlled with the
  `setThreshold` and `setThresholdAnn` arguments. The candidate coverage
  among the different algorithms is calculated for each set (*e.g.*
  `coverage-positive` column) and for all sets (`coverage` column),
  which is based on the presence of a candidate in all the algorithms
  from all sets data. The `consensus` method for sets workflow data
  supports the `filterSets` argument. This controls how the algorithm
  consensus abundance filters (`absMinAbundance`/`relMinAbundance`) are
  applied: if `filterSets=TRUE` then the minimum of all `coverage` set
  specific columns is used to obtain the algorithm abundance. Otherwise
  the overall `coverage` column is used. For instance, consider a
  consensus object to be generated from two objects generated by
  different algorithms (*e.g.* `SIRIUS` and `GenForm`), which both have
  a positive and negative set. Then, if a candidate occurs with both
  algorithms for the positive mode set, but only with the first
  algorithm in the negative mode set, `relMinAbundance=1` will remove
  the candidate if `filterSets=TRUE` (because the minimum relative
  algorithm abundance is `0.5`), while `filterSets=FALSE` will not
  remove the candidate (because based on all sets data the candidate
  occurs in both algorithms).

Two types of annotation data are stored in a `formulasSet` object:

1.  Annotations that are produced from a consensus between set results
    (see `generateFormulas`).

2.  The 'original' annotation data per set, prior to when the set
    consensus was made. This includes candidates that were filtered out
    because of the thresholds set by `setThreshold` and
    `setThresholdAnn`. However, when `filter` or subsetting (`[`)
    operations are performed, the original data is also updated.

In most cases the first data is used. However, in a few cases the
original annotation data is used (as indicated above), for instance, to
re-create the set consensus. It is important to realize that the
original annotation data may have *additional* candidates, and a newly
created set consensus may therefore have 'new' candidates. For instance,
when the object consists of the sets `"positive"` and `"negative"` and
`setThreshold=1` was used to create it, then
`formulas[, sets = "positive", updateConsensus = TRUE]` may now have
additional candidates, *i.e.* those that were not present in the
`"negative"` set and were previously removed due to the consensus
threshold filter.

## See also

The
[`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
base class for more relevant methods and
[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md).
