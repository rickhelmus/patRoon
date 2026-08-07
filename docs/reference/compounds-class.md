# Compound annotations class

Contains data for compound annotations for feature groups.

## Usage

``` r
addFormulaScoring(
  compounds,
  formulas,
  updateScore = FALSE,
  formulaScoreWeight = 1
)

# S4 method for class 'compounds'
defaultExclNormScores(obj)

# S4 method for class 'compounds'
show(object)

# S4 method for class 'compounds'
identifiers(compounds)

# S4 method for class 'compounds'
filter(
  obj,
  minExplainedPeaks = NULL,
  minScore = NULL,
  minFragScore = NULL,
  minFormulaScore = NULL,
  scoreLimits = NULL,
  IMSRangeParams = NULL,
  IMSMatchParams = NULL,
  ...
)

# S4 method for class 'compounds'
addFormulaScoring(
  compounds,
  formulas,
  updateScore = FALSE,
  formulaScoreWeight = 1
)

# S4 method for class 'compounds'
getMCS(obj, index, groupName)

# S4 method for class 'compounds'
plotStructure(obj, index, groupName, width = 500, height = 500)

# S4 method for class 'compounds'
plotScores(
  obj,
  index,
  groupName,
  normalizeScores = "max",
  excludeNormScores = defaultExclNormScores(obj),
  onlyUsed = TRUE
)

# S4 method for class 'compounds'
annotatedPeakList(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  onlyAnnotated = FALSE
)

# S4 method for class 'compounds'
plotSpectrum(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  plotStruct = FALSE,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  showLegend = TRUE,
  maxMolSize = c(0.2, 0.4),
  molRes = c(100, 100),
  ...
)

# S4 method for class 'compounds'
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

# S4 method for class 'compoundsSet'
show(object)

# S4 method for class 'compoundsSet'
delete(obj, i, j, ...)

# S4 method for class 'compoundsSet,ANY,missing,missing'
x[i, j, ..., sets = NULL, updateConsensus = FALSE, drop = TRUE]

# S4 method for class 'compoundsSet'
filter(obj, ..., sets = NULL, updateConsensus = FALSE, negate = FALSE)

# S4 method for class 'compoundsSet'
plotSpectrum(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  plotStruct = FALSE,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  showLegend = TRUE,
  maxMolSize = c(0.2, 0.4),
  molRes = c(100, 100),
  perSet = TRUE,
  mirror = TRUE,
  ...
)

# S4 method for class 'compoundsSet'
addFormulaScoring(
  compounds,
  formulas,
  updateScore = FALSE,
  formulaScoreWeight = 1
)

# S4 method for class 'compoundsSet'
annotatedPeakList(obj, index, groupName, MSPeakLists, formulas = NULL, ...)

# S4 method for class 'compoundsSet'
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

# S4 method for class 'compoundsSet'
unset(obj, set)

# S4 method for class 'compoundsConsensusSet'
unset(obj, set)

# S4 method for class 'compoundsSIRIUS'
delete(obj, i = NULL, j = NULL, ...)
```

## Arguments

- formulas:

  The
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  object that should be used for scoring/annotation. For `plotSpectrum`
  and `annotatedPeakList`: set to `NULL` to ignore.

- updateScore, formulaScoreWeight:

  If `updateScore=TRUE` then the annotation `score` column is updated by
  adding normalized values of the formula score (weighted by
  formulaScoreWeight). Currently, this **only** makes sense for
  annotations performed with `MetFrag`!

- obj, object, compounds, x:

  The `compound` object.

- minExplainedPeaks, scoreLimits:

  Passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- minScore, minFragScore, minFormulaScore:

  Minimum overall score, in-silico fragmentation score and formula
  score, respectively. Set to `NULL` to ignore. The `scoreLimits`
  argument allows for more advanced score filtering.

- IMSRangeParams:

  (**IMS workflow**) A `list` with parameters to be used for filtering
  IMS range data. See
  [`getIMSRangeParams`](https://rickhelmus.github.io/patRoon/reference/getIMSRangeParams.md)
  for details and how to make such a parameter list.

- IMSMatchParams:

  (**IMS workflow**) A `list` with parameters to be used for matching
  IMS data. See
  [`getIMSMatchParams`](https://rickhelmus.github.io/patRoon/reference/getIMSMatchParams.md)
  for details and how to make such a parameter list.

- ...:

  For `plotSpectrum`: Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For `delete`: passed to the function specified as `j`.

  for `filter`: passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

  For `consensus`: any further (and unique) `compounds` objects.

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `compounds` method.

- index:

  The numeric index of the candidate structure.

  For `plotStructure` and `getMCS`: multiple indices (*i.e.* vector with
  length \>=2) should be specified to plot/calculate the most common
  substructure (MCS). Alternatively, `-1` may be specified to select all
  candidates.

  For `plotSpectrum`: two indices can be specified to compare spectra.
  In this case `groupName` should specify values for the spectra to
  compare.

- groupName:

  The name of the feature group for which a plot should be made. To
  compare spectra, two group names can be specified.

- width, height:

  The dimensions (in pixels) of the raster image that should be plotted.

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

- onlyUsed:

  If `TRUE` then only scorings are plotted that actually have been used
  to rank data (see the `scoreTypes` argument to
  [`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md)
  for more details).

- MSPeakLists:

  The
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object with relevant spectral data.

- onlyAnnotated:

  Set to `TRUE` to filter out any peaks that could not be annotated.

- plotStruct:

  If `TRUE` then the candidate structure is drawn in the spectrum.
  Currently not supported when comparing spectra.

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

- maxMolSize:

  Numeric vector of size two with the maximum width/height of the
  candidate structure (relative to the plot size).

- molRes:

  Numeric vector of size two with the resolution of the candidate
  structure (in pixels).

- absMinAbundance, relMinAbundance:

  Minimum absolute or relative (`0-1`) abundance across objects for a
  result to be kept. For instance, `relMinAbundance=0.5` means that a
  result should be present in at least half of the number of compared
  objects. Set to `NULL` to ignore and keep all results. Limits cannot
  be set when `uniqueFrom` is not `NULL`.

- uniqueFrom:

  Set this argument to only retain compounds that are unique within one
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

- i, j, drop:

  Passed to the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`). Note: if `updateConsensus=FALSE` then the
  `setCoverage` column of the annotation results is not updated.

- updateConsensus:

  (**sets workflow**) If `TRUE` then the annonation consensus among set
  results is updated. See the `Sets workflows` section for more details.

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
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).

- setAvgSpecificScores:

  (**sets workflow**) If `TRUE` then set specific annotation scores
  (*e.g.* MS/MS and isotopic pattern match scores) are averaged for the
  set consensus. See
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).

- set:

  (**sets workflow**) The name of the set.

## Value

`addFormulaScoring` returns a `compounds` object updated with formula
scoring.

`getMCS` returns an [rcdk](https://CRAN.R-project.org/package=rcdk)
molecule object (`IAtomContainer`).

`consensus` returns a `compounds` object that is produced by merging
multiple specified `compounds` objects.

## Details

`compounds` objects are obtained from [compound
generators](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).
This class is derived from the
[`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
class, please see its documentation for more methods and other details.

## Methods (by generic)

- `defaultExclNormScores(compounds)`: Returns default scorings that are
  excluded from normalization.

- `show(compounds)`: Show summary information for this object.

- `identifiers(compounds)`: Returns a list containing for each feature
  group a character vector with database identifiers for all candidate
  compounds. The list is named by feature group names, and is typically
  used with the `identifiers` option of
  [`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md).

- `filter(compounds)`: Provides rule based filtering for generated
  compounds. Useful to eliminate unlikely candidates and speed up
  further processing. Also see the
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  method.

- `addFormulaScoring(compounds)`: Adds formula ranking data from a
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  object as an extra compound candidate scoring (`formulaScore` column).
  The formula score for each compound candidate is between `0-1`, where
  *zero* means no match with any formula candidates, and *one* means
  that the compound candidate's formula is the highest ranked.

- `getMCS(compounds)`: Calculates the maximum common substructure (MCS)
  for two or more candidate structures for a feature group. This method
  uses the `get.mcs` function from
  [rcdk](https://CRAN.R-project.org/package=rcdk).

- `plotStructure(compounds)`: Plots a structure of a candidate compound
  using the [rcdk](https://CRAN.R-project.org/package=rcdk) package. If
  multiple candidates are specified (*i.e.* by specifying a `vector` for
  `index`) then the maximum common substructure (MCS) of the selected
  candidates is drawn.

- `plotScores(compounds)`: Plots a barplot with scoring of a candidate
  compound.

- `annotatedPeakList(compounds)`: Returns an MS/MS peak list annotated
  with data from a given candidate compound for a feature group.

- `plotSpectrum(compounds)`: Plots an annotated spectrum for a given
  candidate compound for a feature group. Two spectra can be compared by
  specifying a two-sized vector for the `index` and `groupName`
  arguments.

- `consensus(compounds)`: Generates a consensus of results from multiple
  objects. In order to rank the consensus candidates, first each of the
  candidates are scored based on their original ranking (the scores are
  normalized and the highest ranked candidate gets value `1`). The
  (weighted) mean is then calculated for all scorings of each candidate
  to derive the final ranking (if an object lacks the candidate its
  score will be `0`). The original rankings for each object is stored in
  the `rank` columns.

## Slots

- `MS2QuantMeta`:

  Metadata from MS2Quant filled in by `predictRespFactors`.

  (**sets workflow**) A named `list` with the metadata stored for each
  set.

- `setThreshold,setThresholdAnn,setAvgSpecificScores`:

  (**sets workflow**) A copy of the equally named arguments that were
  passed when this object was created by
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).

- `origFGNames`:

  (**sets workflow**) The original (order of) names of the
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that was used to create this object.

## Note

The values ranges in the `scoreLimits` slot, which are used for
normalization of scores, are based on the *original* scorings when the
compounds were generated (*prior* to employing the `topMost` filter to
[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)).

## IMS workflows

In IMS workflows, reference IMS data to candidates can be assigned with
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_comp.md)
method function. Furthermore, CCS values may be assigned directly to
candidates with
[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
if `database="pubchemlite"`.

This data can be used to prioritize candidates with the `IMSMatchParams`
and `IMSRangeParams` filters.

## S4 class hierarchy

- [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)

  - **`compounds`**

    - `compoundsConsensus`

    - [`compoundsMF`](https://rickhelmus.github.io/patRoon/reference/compoundsMF-class.md)

    - `compoundsSet`

      - `compoundsConsensusSet`

    - `compoundsUnset`

    - [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compoundsSIRIUS-class.md)

## Source

Subscripting of formulae for plots generated by `plotSpectrum` is based
on the `chemistry2expression` function from the
[ReSOLUTION](https://github.com/schymane/ReSOLUTION) package.

## Sets workflows

The `compoundsSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `compounds` and therefore largely follows the
same user interface.

The following methods are specifically defined for sets workflows:

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`compoundsUnset`), which allows it to be used in 'regular'
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
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md),
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
  different algorithms (*e.g.* `SIRIUS` and `MetFrag`), which both have
  a positive and negative set. Then, if a candidate occurs with both
  algorithms for the positive mode set, but only with the first
  algorithm in the negative mode set, `relMinAbundance=1` will remove
  the candidate if `filterSets=TRUE` (because the minimum relative
  algorithm abundance is `0.5`), while `filterSets=FALSE` will not
  remove the candidate (because based on all sets data the candidate
  occurs in both algorithms).

- `addFormulaScoring` Adds the formula scorings to the original data and
  re-creates the annotation set consensus (see below for implications).

Two types of annotation data are stored in a `compoundsSet` object:

1.  Annotations that are produced from a consensus between set results
    (see `generateCompounds`).

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
`compounds[, sets = "positive", updateConsensus = TRUE]` may now have
additional candidates, *i.e.* those that were not present in the
`"negative"` set and were previously removed due to the consensus
threshold filter.

## References

Guha R (2007). “Chemical Informatics Functionality in R.” *Journal of
Statistical Software*, **18**(6).

## See also

The
[`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
base class for more relevant methods and
[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).
