# Find features using SIRIUS

Uses [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) to
find features.

## Usage

``` r
findFeaturesSIRIUS(analysisInfo, verbose = TRUE)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses SIRIUS to automatically find features. This function
is called when calling `findFeatures` with `algorithm="sirius"`.

The features are collected by running the `lcms-align` `SIRIUS` command
for every analysis.

The MS files should be in the `mzML` or `mzXML` format. Furthermore,
this algorithms requires the presence of (data-dependent) MS/MS data.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## Parallelization

`findFeaturesSIRIUS` uses multiprocessing to parallelize computations.
Please see the parallelization section in the handbook for more details
and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

Note that for caching purposes, the analyses files must always exist on
the local host computer, even if it is not participating in
computations.

## References

Duhrkop K, Fleischauer M, Ludwig M, Aksenov AA, Melnik AV, Meusel M,
Dorrestein PC, Rousu J, Bocker S (2019). “SIRIUS 4: a rapid tool for
turning tandem mass spectra into metabolite structure information.”
*Nature Methods*, **16**(4), 299–302.
[doi:10.1038/s41592-019-0344-8](https://doi.org/10.1038/s41592-019-0344-8)
.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
