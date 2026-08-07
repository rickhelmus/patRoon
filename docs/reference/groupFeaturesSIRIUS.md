# Group features using SIRIUS

Uses [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) to
find *and* group features.

## Usage

``` r
groupFeaturesSIRIUS(analysisInfo, verbose = TRUE)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function uses SIRIUS to group features. This function is called
when calling `groupFeatures` with `algorithm="sirius"`.

Finding and grouping features is done by running the `lcms-align`
command on every analyses at once. For this reason, grouping feature
data from other algorithms than `SIRIUS` is not supported.

The MS files should be in the `mzML` or `mzXML` format. Furthermore,
this algorithms requires the presence of (data-dependent) MS/MS data.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## References

Duhrkop K, Fleischauer M, Ludwig M, Aksenov AA, Melnik AV, Meusel M,
Dorrestein PC, Rousu J, Bocker S (2019). “SIRIUS 4: a rapid tool for
turning tandem mass spectra into metabolite structure information.”
*Nature Methods*, **16**(4), 299–302.
[doi:10.1038/s41592-019-0344-8](https://doi.org/10.1038/s41592-019-0344-8)
.

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for more details and other algorithms.
