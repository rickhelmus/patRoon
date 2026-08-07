# Transformation products obtained from formula annotations

Class to store results transformation products (TPs) obtained from
formula annotations.

## Usage

``` r
# S4 method for class 'transformationProductsAnnForm'
filter(
  obj,
  ...,
  minFitFormula = 0,
  minTPScore = 0,
  topMost = NULL,
  verbose = TRUE,
  negate = FALSE
)
```

## Arguments

- obj:

  The `transformationProductsAnnForm` object that should be filtered.

- ...:

  Further arguments passed to the [parent filter
  method](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

- minFitFormula, minTPScore:

  Thresholds related to TP scoring. See
  [`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md)
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
[`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)
base class, please see its documentation for more details. Objects from
this class are returned by
[`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md).

## Methods (by generic)

- `filter(transformationProductsAnnForm)`: Performs rule-based
  filtering. Useful to simplify and clean-up the data.

## S4 class hierarchy

- [`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)

  - **`transformationProductsAnnForm`**

## See also

The base class
[`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)
for more relevant methods and
[`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md)
