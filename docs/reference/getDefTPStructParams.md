# Parameters to handle TP data with structural information

Parameters used by
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
with algorithms that use structural information.

## Usage

``` r
getDefTPStructParams(...)
```

## Arguments

- ...:

  optional named arguments that override defaults.

## Details

The following parameters can be configured:

- `calcLogP` A `character` specifying whether `log P` values should be
  calculated with
  [`rcdk::get.xlogp`](https://rdrr.io/pkg/rcdk/man/get.xlogp.html)
  (`calcLogP="rcdk"`),
  [OpenBabel](https://github.com/openbabel/openbabel)
  (`calcLogP="obabel"`) or not at all (`calcLogP="none"`). The `log P`
  are values of parents and TPs are used for [retention order
  calculation](https://rickhelmus.github.io/patRoon/reference/retDir.md).

- `forceCalcLogP` Force calculation of `Log P` values, even if already
  provided by the TP generation algorithm. This is primarily useful to
  obtain log P values that were consistently calculated with the same
  algorithm, as some algorithms may only partially output these values
  (*e.g.* not for parents).

- `forceCalcRetDir` Force calculation of [retention order
  directions](https://rickhelmus.github.io/patRoon/reference/retDir.md),
  even if already provided by the TP generation algorithm. This is
  primarily intended for re-calculation of library TP data, which may
  have been calculated with different log P values.

- `minLogPDiff` The minimum difference in `log P` values between a
  parent and its TPs to be considered eluting differently. This is used
  for [retention order
  calculation](https://rickhelmus.github.io/patRoon/reference/retDir.md).

- `calcSims` If set to `TRUE` then structural similarities between the
  parent and its TPs are calculated. The [filter
  method](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)
  can be used to threshold structural similarities. This may be useful
  under the assumption that parents and TPs who have a high structural
  similarity, also likely have a high MS/MS spectral similarity (which
  can be evaluated after componentization with
  [`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)).

- `fpType` The type of structural fingerprint that should be calculated.
  See the `type` argument of the `get.fingerprint` function of
  [rcdk](https://CRAN.R-project.org/package=rcdk).

- `fpSimMethod` The method for calculating similarities (i.e. not
  dissimilarity!). See the `method` argument of the `fp.sim.matrix`
  function of the
  [fingerprint](https://CRAN.R-project.org/package=fingerprint) package.

These parameters are passed as a named `list` as the `TPStructParams`
argument to functions.

The `getDefTPStructParams` function generates such parameter list with
defaults.

## References

Guha R (2007). “Chemical Informatics Functionality in R.” *Journal of
Statistical Software*, **18**(6).\
\
OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .
