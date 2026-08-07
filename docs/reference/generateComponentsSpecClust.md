# Generate components based on MS/MS similarity

Generates components based on MS/MS similarity between feature groups.

## Usage

``` r
generateComponentsSpecClust(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsSpecClust(
  fGroups,
  MSPeakLists,
  method = "complete",
  specSimParams = getDefSpecSimParams(),
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

- MSPeakLists:

  The
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object for the given feature groups that should be used for MS
  spectral similarity calculations.

- method:

  Clustering method that should be applied (passed to
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html)).

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

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
[`componentsSpecClust`](https://rickhelmus.github.io/patRoon/reference/componentsSpecClust-class.md).

## Details

This function uses hierarchical clustering of MS/MS spectra to generate
components. This function is called when calling `generateComponents`
with `algorithm="specclust"`.

The similarities are converted to a distance matrix and used as input
for hierarchical clustering, and the resulting dendrogram is
automatically cut with `cutreeDynamicTree`. The clustering is performed
with
[`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html).

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
the spectral similarities for each set are combined as is described for
the
[`spectrumSimilarity`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
method for sets workflows.

## References

Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python.” *Journal of Statistical
Software*, **53**(9), 1–18.
[doi:10.18637/jss.v053.i09](https://doi.org/10.18637/jss.v053.i09) .

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.

## Author

Rick Helmus \<<r.helmus@uva.nl>\> and Bas van de Velde (major
contributions to spectral binning and similarity calculation).
