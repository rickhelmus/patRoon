# Conversion between mobility and CCS

Utility functions to convert between mobility and CCS data.

## Usage

``` r
convertMobilityToCCS(mobility, mz, CCSParams, charge = NULL)

convertCCSToMobility(ccs, mz, CCSParams, charge = NULL)
```

## Arguments

- mobility, ccs:

  A `numeric` vector with mobility or CCS values that should be
  converted.

- mz:

  A `numeric` vector with the *m/z* values that map to the input
  mobility or CCS values.

- CCSParams:

  A `list` with parameters for mobility \<–\> CCS conversion. See
  [`getCCSParams`](https://rickhelmus.github.io/patRoon/reference/getCCSParams.md)
  for details and to make such parameter lists. Set to `NULL` to skip
  conversions.

- charge:

  A `numeric` vector with the ion charges that map to the input mobility
  or CCS data. Will be recycled if necessary. If `NULL` then the charge
  configured in `CCSParams` is used.
