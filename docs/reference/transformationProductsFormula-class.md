# Base transformation products (TP) class with formula information

Holds information for all TPs for a set of parents, including chemical
formulae.

## Usage

``` r
# S4 method for class 'transformationProductsFormula'
plotGraph(
  obj,
  which,
  components = NULL,
  prune = TRUE,
  onlyCompletePaths = FALSE,
  width = NULL,
  height = NULL
)
```

## Arguments

- obj:

  `transformationProductsFormula` derived object to be accessed

- which:

  Either a `character` or `integer` vector with one or more
  names/indices of the parents to plot.

- components:

  If specified (*i.e.* not `NULL`), a
  [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)
  object that is used for matching the graph with screening results. The
  TPs that were found will be marked. See also the `prune` and
  `onlyCompletePaths` arguments.

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

## Value

`plotGraph` returns the result of `visNetwork`.

## Details

This (virtual) class is derived from the
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
base class, please see its documentation for more details. Objects from
this class are returned by [TP
generators](https://rickhelmus.github.io/patRoon/reference/generateTPs.md).
More specifically, algorithms that works with chemical formulae (*e.g.*
`library_formula`), uses this class to store their results. The methods
defined for this class extend the functionality for the base
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
class.

## Methods (by generic)

- `plotGraph(transformationProductsFormula)`: Plots an interactive
  hierarchy graph of the transformation products. The resulting graph
  can be browsed interactively and allows exploration of the different
  TP formation pathways. Furthermore, results from [TP
  componentization](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
  can be used to match the hierarchy with screening results. The graph
  is rendered with visNetwork.

## S4 class hierarchy

- [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)

  - **`transformationProductsFormula`**

    - [`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md)

    - `transformationProductsLibraryFormula`

## See also

The base class
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
for more relevant methods and
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
