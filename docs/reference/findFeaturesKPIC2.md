# Find features using KPIC2

Uses the [KPIC2](https://github.com/hcji/KPIC2) R package to extract
features.

## Usage

``` r
findFeaturesKPIC2(
  analysisInfo,
  kmeans = TRUE,
  level = 1000,
  ...,
  parallel = TRUE,
  verbose = TRUE
)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- kmeans:

  If `TRUE` then
  [`getPIC.kmeans`](https://rdrr.io/pkg/KPIC/man/getPIC.kmeans.html) is
  used to obtain PICs, otherwise it is
  [`getPIC`](https://rdrr.io/pkg/KPIC/man/getPIC.html).

- level:

  Passed to [`getPIC`](https://rdrr.io/pkg/KPIC/man/getPIC.html) or
  [`getPIC.kmeans`](https://rdrr.io/pkg/KPIC/man/getPIC.kmeans.html)

- ...:

  Further parameters passed to
  [`getPIC`](https://rdrr.io/pkg/KPIC/man/getPIC.html)/[`getPIC.kmeans`](https://rdrr.io/pkg/KPIC/man/getPIC.kmeans.html)

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

This function uses KPIC2 to automatically find features. This function
is called when calling `findFeatures` with `algorithm="kpic2"`.

The MS files should be in the `mzML` or `mzXML` format.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## References

Ji H, Zeng F, Xu Y, Lu H, Zhang Z (2017). “KPIC2: An Effective Framework
for Mass Spectrometry-Based Metabolomics Using Pure Ion Chromatograms.”
*Analytical Chemistry*, **89**(14), 7631–7640.
[doi:10.1021/acs.analchem.7b01547](https://doi.org/10.1021/acs.analchem.7b01547)
.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
