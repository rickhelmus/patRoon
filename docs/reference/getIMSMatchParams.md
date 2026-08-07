# Parameters for IMS matching

Parameters that define how mobility or CCS values between *e.g.*
features and suspects should be matched.

## Usage

``` r
getIMSMatchParams(param, ...)
```

## Arguments

- param:

  Should be `"mobility"` or `"CCS"` to match by mobility or CCS,
  respectively.

- ...:

  optional named arguments that override defaults.

## Details

The following parameters should be defined:

- `param` Should be `"mobility"` or `"CCS"` to match by mobility or CCS,
  respectively.

- `window`,`relative` The `window` parameter sets the tolerance window
  size used for matching. If `relative=TRUE` then the tolerance is
  relative (`0-1`). The defaults for `window` are (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)):

  - `defaultLim("mobility", "medium")` (`param="mobility"` and
    `relative=FALSE`)

  - `defaultLim("mobility", "medium_rel")` (`param="mobility"` and
    `relative=TRUE`)

  - `defaultLim("CCS", "medium")` (`param="CCS"` and `relative=FALSE`)

  - `defaultLim("CCS", "medium_rel")` (`param="CCS"` and
    `relative=TRUE`)

- `minMatches` The minimum number of mobility/CCS matches for a suspect
  hit. If the number of available reference mobility/CCS values in the
  suspect list is less than `minMatches`, then that number is used as
  threshold. Set to `0` to disable.

These parameters are passed as a named `list` as the `IMSMatchParams`
argument to functions.

The `getIMSMatchParams` function generates such parameter list with
defaults.

## Note

If negation is enabled with suspect filtering and `minMatches>0`, then
the `window` match filter is *not* negated. Negating both would lead to
unexpected results, *i.e.* suspects outside `window` are kept and
increase the number of matched suspects as seen by the `minMatches`
filter.
