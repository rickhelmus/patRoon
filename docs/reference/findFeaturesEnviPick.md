# Find features using enviPick

Uses the
[`enviPickwrap`](https://rdrr.io/pkg/enviPick/man/enviPickwrap.html)
function from the enviPick R package to extract features.

## Usage

``` r
findFeaturesEnviPick(analysisInfo, ..., parallel = TRUE, verbose = TRUE)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- ...:

  Further parameters passed to
  [`enviPickwrap`](https://rdrr.io/pkg/enviPick/man/enviPickwrap.html).

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses enviPick to automatically find features. This
function is called when calling `findFeatures` with
`algorithm="envipick"`.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## Note

The analysis files must be in the `mzXML` format.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
