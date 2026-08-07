# Parameters to specify a IMS data range

Parameters that define a range of mobility or CCS values.

## Usage

``` r
getIMSRangeParams(param, lower, upper, mzRelative = FALSE)
```

## Arguments

- param, lower, upper, mzRelative:

  Arguments to specify the IMS range parameters, see Details.

## Details

The following parameters are used to define an IMS range:

- `param` Should be `"mobility"` or `"CCS"` to specify a mobility or CCS
  range, respectively.

- `lower`,`upper` The lower and upper range.

- `mzRelative` Set to `TRUE` to specify an IMS range that is normalized
  by *m/z*.

These parameters are passed as a named `list` as the `IMSRangeParams`
argument to functions.

The `getIMSRangeParams` function generates such parameter list with
defaults.
