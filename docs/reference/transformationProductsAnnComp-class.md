# Transformation products obtained from compound annotations

Class to store results transformation products (TPs) obtained from
compound annotations.

## Usage

``` r
# S4 method for class 'transformationProductsAnnComp'
filter(
  obj,
  ...,
  minFitFormula = 0,
  minFitCompound = 0,
  minSimSusp = 0,
  minFitCompOrSimSusp = c(0, 0),
  minTPScore = 0,
  topMost = NULL,
  verbose = TRUE,
  negate = FALSE
)
```

## Arguments

- obj:

  The `transformationProductsAnnComp` object that should be filtered.

- ...:

  Further arguments passed to the [parent filter
  method](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

- minFitFormula, minFitCompound, minSimSusp, minFitCompOrSimSusp,
  minTPScore:

  Thresholds related to TP scoring. See
  [`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md)
  for more details.

- topMost:

  Only keep this number of top-most TPs (based on `TPScore`) for each
  parent/feature group combination. Set to `NULL` to skip this step.

- verbose:

  If set to `FALSE` then no text output is shown.

- negate:

  If `TRUE` then filters are performed in opposite manner.

## Details

This class is derived from the
[`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)
base class, please see its documentation for more details. Objects from
this class are returned by
[`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md).

## Methods (by generic)

- `filter(transformationProductsAnnComp)`: Performs rule-based
  filtering. Useful to simplify and clean-up the data.

## S4 class hierarchy

- [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

  - **`transformationProductsAnnComp`**

## See also

The base class
[`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)
for more relevant methods and
[`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md)
