# Components based on parent and transformation product (TP) linkage.

This class is derived from
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
and is used to store components that result from linking feature groups
that are (predicted to be) parents with feature groups that (are
predicted to be) transformation products. For more details, see
[`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md).

## Usage

``` r
# S4 method for class 'componentsTPs'
as.data.table(x, candidates = FALSE)

# S4 method for class 'componentsTPs'
filter(
  obj,
  ...,
  retDirMatch = FALSE,
  minSpecSim = NULL,
  minSpecSimPrec = NULL,
  minSpecSimBoth = NULL,
  minTotFragMatches = NULL,
  minTotNLMatches = NULL,
  minFragMatches = NULL,
  minNLMatches = NULL,
  formulas = NULL,
  verbose = TRUE,
  negate = FALSE
)

# S4 method for class 'componentsTPs'
plotGraph(obj, onlyLinked = TRUE, width = NULL, height = NULL)
```

## Arguments

- x, obj:

  A `componentsTPs` object.

- candidates:

  If `TRUE` then candidate specific tables are merged in the result
  table. See
  [`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
  (`Result columns` section) for more details.

- ..., verbose:

  Further arguments passed to the base
  [`filter method`](https://rickhelmus.github.io/patRoon/reference/components-class.md).

- retDirMatch:

  If set to `TRUE`, only keep TPs for which the [retention order
  direction](https://rickhelmus.github.io/patRoon/reference/retDir.md)
  from feature group data matches that of the expected values from TP
  data (`retDir` and `TP_retDir` columns). TPs will never be removed if
  either of the directions is `0` (*i.e.* unknown or not significantly
  different than the parent).

- minSpecSim, minSpecSimPrec, minSpecSimBoth:

  The minimum spectral similarity of a TP compared to its parent
  (`0-1`). The `minSpecSimPrec` and `minSpecSimBoth` apply to binned
  data that is shifted with the `"precursor"` and `"both"` method,
  respectively (see [MS spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  for more details). Set to `NULL` to ignore.

- minTotFragMatches, minTotNLMatches, minFragMatches, minNLMatches:

  Minimum number of (total) parent/TP fragment and neutral loss matches.
  Set to `NULL` to ignore. See the `Result columns` section in
  [`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
  for more details.

- formulas:

  A
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  object. The formula annotation data in this object is to verify if
  elemental additions/subtractions from metabolic logic reactions are
  possible (hence, it only works with data from
  [`generateTPsLogic`](https://rickhelmus.github.io/patRoon/reference/generateTPsLogic.md)).
  To verify elemental additions, only TPs with at least one candidate
  formula that has these elements are kept. Similarly, for elemental
  subtractions, any of the parent candidate formulae must contain the
  subtraction elements. Note that TPs are currently not filtered if
  either the parent or the TP has no formula annotations. Set to `NULL`
  to ignore.

- negate:

  If `TRUE` then filters are applied in opposite manner.

- onlyLinked:

  If `TRUE` then only components with links are plotted.

- width, height:

  Passed to `visNetwork`.

## Value

`filter` returns a filtered `componentsTPs` object.

`plotGraph` returns the result of `visNetwork`.

## Methods (by generic)

- `as.data.table(componentsTPs)`: Returns all component data as a
  `data.table`.

- `filter(componentsTPs)`: Provides various rule based filtering options
  to clean and prioritize TP data.

- `plotGraph(componentsTPs)`: Plots an interactive network graph for
  linked components. Components are linked with each other if one or
  more transformation products overlap. The graph is constructed with
  the igraph package and rendered with visNetwork.

## Slots

- `fromTPs`:

  A `logical` that is `TRUE` when the componentization was performed
  with
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
  data (*i.e.* the TPs argument was not `NULL`).

- `parentsFromScreening`:

  A `logical` that is `TRUE` when the parents were obtained from
  screening data.

- `TPsFromScreening`:

  A `logical` that is `TRUE` when the TPs were obtained from screening
  data.

## Note

The intensity values for components (used by `plotSpectrum`) are set to
a dummy value (1) as no single intensity value exists for this kind of
components.

## S4 class hierarchy

- [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

  - **`componentsTPs`**

## References

Antonov M, Csárdi G, Horvát S, Müller K, Nepusz T, Noom D, Salmon M,
Traag V, Welles BF, Zanini F (2023). “igraph enables fast and robust
network analysis across programming languages.” *arXiv preprint
arXiv:2311.10260*.
[doi:10.48550/arXiv.2311.10260](https://doi.org/10.48550/arXiv.2311.10260)
.\
\
Csárdi G, Nepusz T (2006). “The igraph software package for complex
network research.” *InterJournal*, **Complex Systems**, 1695.
<https://igraph.org>.\
\
Csárdi G, Nepusz T, Traag V, Horvát S, Zanini F, Noom D, Müller K,
Schoch D, Salmon M (2026). *igraph: Network Analysis and Visualization
in R*.
[doi:10.5281/zenodo.7682609](https://doi.org/10.5281/zenodo.7682609) . R
package version 2.3.3, <https://CRAN.R-project.org/package=igraph>.

## See also

[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
for other relevant methods and
[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
