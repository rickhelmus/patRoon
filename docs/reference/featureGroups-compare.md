# Comparing feature groups

Functionality to compare feature groups and make a consensus.

## Usage

``` r
comparison(..., groupAlgo, groupArgs = list(rtalign = FALSE))

# S4 method for class 'featureGroups'
comparison(..., groupAlgo, groupArgs = list(rtalign = FALSE))

# S4 method for class 'featureGroupsComparison'
hasIMS(obj)

# S4 method for class 'featureGroupsComparison,missing'
plot(x, retMin = FALSE, ...)

# S4 method for class 'featureGroupsComparison'
plotVenn(obj, which = NULL, ...)

# S4 method for class 'featureGroupsComparison'
plotUpSet(obj, which = NULL, ...)

# S4 method for class 'featureGroupsComparison'
plotChord(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...)

# S4 method for class 'featureGroupsComparison'
consensus(
  obj,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  verifyAnaInfo = TRUE
)

# S4 method for class 'featureGroupsSet'
comparison(..., groupAlgo, groupArgs = list(rtalign = FALSE))

# S4 method for class 'featureGroupsComparisonSet'
consensus(obj, ...)
```

## Arguments

- ...:

  For `comparison`: `featureGroups` objects that should be compared. If
  the arguments are named (*e.g.* `myGroups = fGroups`) then these are
  used for labelling, otherwise objects are automatically labelled by
  their
  [`algorithm`](https://rickhelmus.github.io/patRoon/reference/generics.md).

  For `plot`, `plotVenn`, `plotChord`: further options passed to `plot`,
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  functions (*e.g.* `draw.pairwise.venn`) and `chordDiagram`
  respectively.

  For `plotUpSet`: any further arguments passed to the `plotUpSet`
  method defined for
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

- groupAlgo:

  The
  [`feature grouping algorithm`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
  that should be used for grouping *pseudo* features (see details).
  Valid values are: `"xcms"`, `xcms3`, `kpic2`, `"openms"` or
  `"greedy"`.

  (**IMS workflow**) For IMS workflows the algorithm selection is
  equally limited to as what is used by
  [`feature grouping`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md),
  *i.e.* only `greedy` is currently supported.

- groupArgs:

  A `list` containing further parameters for
  [`feature grouping`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md).

- x, obj:

  The `featureGroupsComparison` object.

- retMin:

  If `TRUE` retention times are plotted as minutes (seconds otherwise).

- which:

  A character vector specifying one or more labels of compared feature
  groups. For `plotVenn`: if `NULL` then all compared groups are used.

- addSelfLinks:

  If `TRUE` then 'self-links' are added which represent non-shared data.

- addRetMzPlots:

  Set to `TRUE` to enable *m/z* *vs* retention time scatter plots.

- absMinAbundance, relMinAbundance:

  Minimum absolute or relative (`0-1`) abundance across objects for a
  result to be kept. For instance, `relMinAbundance=0.5` means that a
  result should be present in at least half of the number of compared
  objects. Set to `NULL` to ignore and keep all results. Limits cannot
  be set when `uniqueFrom` is not `NULL`.

- uniqueFrom:

  Set this argument to only retain feature groups that are unique within
  one or more of the objects for which the consensus is made. Selection
  is done by setting the value of `uniqueFrom` to a `logical` (values
  are recycled), `numeric` (select by index) or a `character` (as
  obtained with `algorithm(obj)`). For `logical` and `numeric` values
  the order corresponds to the order of the objects given for the
  consensus. Set to `NULL` to ignore.

- uniqueOuter:

  If `uniqueFrom` is not `NULL` and if `uniqueOuter=TRUE`: only retain
  data that are also unique between objects specified in `uniqueFrom`.

- verifyAnaInfo:

  If `FALSE` then the analysis information is not verified to be equal
  for all compared objects. This is mainly only useful when the data is
  the same but stored in different formats (*e.g.* `mzXML`/`mzML`).

## Value

`comparison` returns a
[`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md)
object.

`plotVenn` (invisibly) returns a list with the following fields:

- `gList` the `gList` object that was returned by the utilized
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  function.

- `areas` The total area for each plotted group.

- `intersectionCounts` The number of intersections between groups.

The order for the `areas` and `intersectionCounts` fields is the same as
the parameter order from the used plotting function (see *e.g.*
`draw.pairwise.venn` and `draw.triple.venn`).

`consensus` returns a
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object with a consensus from the compared feature groups.

## Details

Feature groups objects originating from differing feature finding and/or
grouping algorithms (or their parameters) may be compared to assess
their output and generate a consensus.

The `comparison` method generates a
[`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md)
object from given feature groups objects, which in turn may be used for
(visually) comparing presence of feature groups and generating a
consensus. Internally, this function will collapse each feature groups
object to *pseudo* features objects by averaging their retention times,
*m/z* values and intensities, where each original feature groups object
becomes an 'analysis'. All *pseudo* features are then grouped using
[regular feature grouping
algorithms](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
so that a comparison can be made.

`hasIMS` returns `TRUE` if the object has ion mobility information.

`plot` generates an *m/z* *vs* retention time plot.

`plotVenn` plots a Venn diagram outlining unique and shared feature
groups between up to five compared feature groups.

`plotUpSet` plots an UpSet diagram outlining unique and shared feature
groups.

`plotChord` plots a chord diagram to visualize the distribution of
feature groups.

`consensus` combines all compared feature groups and averages their
retention, *m/z* and intensity data. Not yet supported for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
