# Generate components based on intensity profiles

Generates components based on intensity profiles of feature groups.

## Usage

``` r
generateComponentsIntClust(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsIntClust(
  fGroups,
  method = "complete",
  metric = "euclidean",
  normalized = TRUE,
  average = TRUE,
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1,
  IMS = "maybe"
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which components should be generated.

- ...:

  Any parameters to be passed to the selected component generation
  algorithm.

- method:

  Clustering method that should be applied (passed to
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html)).

- metric:

  Distance metric used to calculate the distance matrix (passed to
  `daisy`).

- normalized, average:

  Passed to
  [`as.data.table`](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
  to perform normalization and averaging of data.

- maxTreeHeight, deepSplit, minModuleSize:

  Arguments used by `cutreeDynamicTree`.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  componentization in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

## Value

The components are stored in objects derived from
[`componentsIntClust`](https://rickhelmus.github.io/patRoon/reference/componentsIntClust-class.md).

## Details

This function uses hierarchical clustering of intensity profiles to
generate components. This function is called when calling
`generateComponents` with `algorithm="intclust"`.

Hierarchical clustering is performed on normalized (and optionally
replicate averaged) intensity data and the resulting dendrogram is
automatically cut with `cutreeDynamicTree`. The distance matrix is
calculated with `daisy` and clustering is performed with
[`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html).
The clustering of the resulting components can be further visualized and
modified using the methods defined for
[`componentsIntClust`](https://rickhelmus.github.io/patRoon/reference/componentsIntClust-class.md).

## IMS workflows

In IMS workflows with post mobility assignment (see
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
it may be necessary to call
[`expandForIMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
when componentization was performed *prior* to mobility assignments, see
its documentation for more details.

If mobilities were already assigned prior to componentization, then the
`IMS` argument selects which feature groups are subjected to
componentization. Data for IMS feature groups that were not considered
(*i.e.* when `IMS` is `FALSE` or `"maybe"`), will be expanded similarly
as is done by
[`expandForIMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md).

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
normalization of feature intensities occur per set.

## References

Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python.” *Journal of Statistical
Software*, **53**(9), 1–18.
[doi:10.18637/jss.v053.i09](https://doi.org/10.18637/jss.v053.i09) .

Schollee JE, Bourgin M, von Gunten U, McArdell CS, Hollender J (2018).
“Non-target screening to trace ozonation transformation products in a
wastewater treatment train including different post-treatments.” *Water
Research*, **142**, 267–278.
[doi:10.1016/j.watres.2018.05.045](https://doi.org/10.1016/j.watres.2018.05.045)
.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
