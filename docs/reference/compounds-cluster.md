# Hierarchical clustering of compounds

Perform hierarchical clustering of structure candidates based on
chemical similarity and obtain overall structural information based on
the maximum common structure (MCS).

## Usage

``` r
makeHCluster(obj, method = "complete", ...)

# S4 method for class 'compounds'
makeHCluster(
  obj,
  method,
  fpType = "extended",
  fpSimMethod = "tanimoto",
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1
)
```

## Arguments

- obj:

  The
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object to be clustered.

- method:

  The clustering method passed to
  [`hclust`](https://rdrr.io/r/stats/hclust.html).

- ...:

  further arguments specified to methods.

- fpType:

  The type of structural fingerprint that should be calculated. See the
  `type` argument of the `get.fingerprint` function of
  [rcdk](https://CRAN.R-project.org/package=rcdk).

- fpSimMethod:

  The method for calculating similarities (i.e. not dissimilarity!). See
  the `method` argument of the `fp.sim.matrix` function of the
  [fingerprint](https://CRAN.R-project.org/package=fingerprint) package.

- maxTreeHeight, deepSplit, minModuleSize:

  Arguments used by `cutreeDynamicTree`.

## Value

`makeHCluster` returns an
[`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md)
object.

## Details

Often many possible chemical structure candidates are found for each
feature group when performing [compound
annotation](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).
Therefore, it may be useful to obtain an overview of their general
structural properties. One strategy is to perform hierarchical
clustering based on their chemical (dis)similarity, for instance, using
the Tanimoto score. The resulting clusters can then be characterized by
evaluating their *maximum common substructure* (MCS).

`makeHCluster` performs hierarchical clustering of all structure
candidates for each feature group within a
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
object. The resulting dendrograms are automatically cut using the
`cutreeDynamicTree` function from the dynamicTreeCut package. The
returned
[`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md)
object can then be used, for instance, for plotting dendrograms and MCS
structures and manually re-cutting specific clusters.

## Source

The methodology applied here has been largely derived from `chemclust.R`
from the metfRag package and the package vignette of
[rcdk](https://CRAN.R-project.org/package=rcdk).

## References

Guha R (2007). “Chemical Informatics Functionality in R.” *Journal of
Statistical Software*, **18**(6).

## See also

compoundsCluster
