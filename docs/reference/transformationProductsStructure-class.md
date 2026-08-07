# Base transformation products (TP) class with structure information

Holds information for all TPs for a set of parents, including structural
information.

## Usage

``` r
# S4 method for class 'transformationProductsStructure'
convertToMFDB(TPs, out, includeParents = FALSE)

# S4 method for class 'transformationProductsStructure'
filter(
  obj,
  ...,
  removeParentIsomers = FALSE,
  removeTPIsomers = FALSE,
  removeDuplicates = FALSE,
  minSimilarity = NULL,
  verbose = TRUE,
  negate = FALSE
)

# S4 method for class 'transformationProductsStructure'
plotGraph(
  obj,
  which,
  components = NULL,
  structuresMax = 25,
  prune = TRUE,
  onlyCompletePaths = FALSE,
  width = NULL,
  height = NULL
)

# S4 method for class 'transformationProductsStructure'
plotVenn(obj, ..., commonParents = FALSE, labels = NULL, vennArgs = NULL)

# S4 method for class 'transformationProductsStructure'
plotUpSet(
  obj,
  ...,
  commonParents = FALSE,
  labels = NULL,
  nsets = NULL,
  nintersects = NA,
  upsetArgs = NULL
)

# S4 method for class 'transformationProductsStructure'
consensus(
  obj,
  ...,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  labels = NULL
)
```

## Arguments

- out:

  The file name of the the output `MetFrag` database.

- includeParents:

  Set to `TRUE` to include the parents in the database.

- obj, TPs:

  `transformationProductsStructure` derived object to be accessed

- ...:

  For `filter`: Further argument passed to the base
  [`filter method`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

  For `plotVenn`, `plotUpSet` and `consensus`: further (unique)
  `transformationProductsStructure` objects.

- removeParentIsomers:

  If `TRUE` then TPs with an equal formula as their parent (isomers) are
  removed.

- removeTPIsomers:

  If `TRUE` then all TPs with equal formula as any sibling TPs (isomers)
  are removed. Unlike `removeDuplicates`, *all* TP candidates are
  removed (including the first match). This filter automatically sets
  `removeDuplicates=TRUE` so that TPs are only removed if with different
  structure.

- removeDuplicates:

  If `TRUE` then the TPs of a parent with duplicate structures (SMILES)
  are removed. Such duplicates may occur when different transformation
  pathways yield the same TPs. The first TP candidate with duplicate
  structure will be kept.

- minSimilarity:

  Minimum structure similarity (`0-1`) that a TP should have relative to
  its parent. This data is only available if the `calcSims` argument to
  [`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
  was set to `TRUE`. May be useful under the assumption that parents and
  TPs who have a high structural similarity, also likely have a high
  MS/MS spectral similarity (which can be evaluated after
  componentization with
  [`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md).
  Any values that are `NA` are removed (which only occur when a
  consensus was made from objects that not all have similarity
  information).

- verbose:

  If set to `FALSE` then no text output is shown.

- negate:

  If `TRUE` then filters are performed in opposite manner.

- which:

  Either a `character` or `integer` vector with one or more
  names/indices of the parents to plot.

- components:

  If specified (*i.e.* not `NULL`), a
  [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)
  object that is used for matching the graph with screening results. The
  TPs that were found will be marked. See also the `prune` and
  `onlyCompletePaths` arguments.

- structuresMax:

  An `integer` with the maximum number of structures to plot. Setting a
  maximum is mainly done to avoid long times needed to construct the
  graph.

- prune:

  If `TRUE` and `components` is set, then pathways without *any*
  detected TPs are not shown (pruned). See also the `onlyCompletePaths`
  and `components` arguments.

- onlyCompletePaths:

  If `TRUE` and `components` is set, then only pathways are shown for
  which *all* TPs were detected. See also the `prune` and `components`
  arguments.

- width, height:

  Passed to `visNetwork`.

- commonParents:

  Only consider TPs from parents that are common to all compared
  objects.

- labels:

  A `character` with names to use for labelling. If `NULL` labels are
  automatically generated.

- vennArgs:

  A `list` with further arguments passed to VennDiagram plotting
  functions. Set to `NULL` to ignore.

- nsets, nintersects:

  See [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html). If
  `nsets == NULL` then it will be set to the number of compared items.

- upsetArgs:

  A list with any further arguments to be passed to
  [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html). Set to `NULL` to
  ignore.

- absMinAbundance, relMinAbundance:

  Minimum absolute or relative (`0-1`) abundance across objects for a
  result to be kept. For instance, `relMinAbundance=0.5` means that a
  result should be present in at least half of the number of compared
  objects. Set to `NULL` to ignore and keep all results. Limits cannot
  be set when `uniqueFrom` is not `NULL`.

- uniqueFrom:

  Set this argument to only retain TPs that are unique within one or
  more of the objects for which the consensus is made. Selection is done
  by setting the value of `uniqueFrom` to a `logical` (values are
  recycled), `numeric` (select by index) or a `character` (as obtained
  with `algorithm(obj)`). For `logical` and `numeric` values the order
  corresponds to the order of the objects given for the consensus. Set
  to `NULL` to ignore.

- uniqueOuter:

  If `uniqueFrom` is not `NULL` and if `uniqueOuter=TRUE`: only retain
  data that are also unique between objects specified in `uniqueFrom`.

## Value

`filter` returns a filtered `transformationProductsStructure` object.

`plotGraph` returns the result of `visNetwork`.

`plotVenn` (invisibly) returns a list with the following fields:

- `gList` the `gList` object that was returned by the utilized
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  function.

- `areas` The total area for each plotted group.

- `intersectionCounts` The number of intersections between groups.

The order for the `areas` and `intersectionCounts` fields is the same as
the parameter order from the used plotting function (see *e.g.*
`draw.pairwise.venn` and `draw.triple.venn`).

`consensus` returns a `transformationProductsStructure` object that is
produced by merging results from multiple
`transformationProductsStructure` objects.

## Details

This (virtual) class is derived from the
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
base class, please see its documentation for more details. Objects from
this class are returned by [TP
generators](https://rickhelmus.github.io/patRoon/reference/generateTPs.md).
More specifically, algorithms that works with chemical structures
(*e.g.* `biotransformer`), uses this class to store their results. The
methods defined for this class extend the functionality for the base
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
class.

## Methods (by generic)

- `convertToMFDB(transformationProductsStructure)`: Exports this object
  as a `.csv` file that can be used as a `MetFrag` local database. Any
  duplicate TPs (formed by different pathways or parents) will be merged
  based on their InChIKey.

- `filter(transformationProductsStructure)`: Performs rule-based
  filtering. Useful to simplify and clean-up the data.

- `plotGraph(transformationProductsStructure)`: Plots an interactive
  hierarchy graph of the transformation products. The resulting graph
  can be browsed interactively and allows exploration of the different
  TP formation pathways. Furthermore, results from [TP
  componentization](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
  can be used to match the hierarchy with screening results. The graph
  is rendered with visNetwork.

- `plotVenn(transformationProductsStructure)`: plots a Venn diagram
  (using [VennDiagram](https://CRAN.R-project.org/package=VennDiagram))
  outlining unique and shared candidates of up to five different
  `featureAnnotations` objects.

- `plotUpSet(transformationProductsStructure)`: Plots an UpSet diagram
  (using the [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html)
  function) outlining unique and shared TPs between different
  `transformationProductsStructure` objects.

- `consensus(transformationProductsStructure)`: Generates a consensus
  from different `transformationProductsStructure` objects. Currently
  this removes any hierarchical data, and all TPs are considered to
  originate from the same (original) parent.

## Note

`consensus`: If the `retDir` values differs between matched TPs it will
be set to `0`. If structure similarity data is available (*i.e.*
`calcSims=TRUE` to `generateTPs`) then the mean similarity is
calculated.

## Comparison between objects

The methods that compare different objects (*e.g.* `plotVenn` and
`consensus`) use the InChIKey to match TPs between objects. Moreover,
the parents between objects are matched by their name. Hence, it is
*crucial* that the input parents to
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
(*i.e.* the `parents` argument) are named equally.

## S4 class hierarchy

- [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)

  - **`transformationProductsStructure`**

    - `transformationProductsStructureConsensus`

    - `transformationProductsCTS`

    - [`transformationProductsAnnComp`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnComp-class.md)

    - `transformationProductsBT`

    - `transformationProductsLibrary`

## References

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

The base class
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
for more relevant methods and
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
