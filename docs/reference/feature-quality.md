# Feature quality definitions and utilities

Functions to define and work with chromatographic peak quality metrics
used for feature and feature group quality calculations.

## Usage

``` r
featureQualities(qualities = NULL)

featureGroupQualities(qualities = NULL)

featureQualityNames(feat = TRUE, group = TRUE, scores = FALSE, totScore = TRUE)
```

## Arguments

- qualities:

  A character vector specifying which qualities to return. If `NULL`,
  returns all available quality definitions.

- feat:

  If `TRUE` then names specific to features are returned.

- group:

  If `TRUE` then names specific to groups are returned.

- scores:

  If `TRUE` the score names are returned, otherwise the quality names.

- totScore:

  If `TRUE` (and `scores=TRUE`) then the name of the total score is
  included.

## Value

For `featureQualities` and `featureGroupQualities`: A named list
containing quality definitions. Each element contains:

- `func`: The MetaClean function to calculate the quality metric

- `HQ`: "HV" (high values good) or "LV" (low values good)

- `range`: Expected range of values (may be Inf for unbounded ranges)

For `featureQualityNames`: A character vector with quality and/or score
names.

## Details

These functions provide access to quality metrics that are used to
assess the quality of detected features and feature groups by the
[`calculatePeakQualities`](https://rickhelmus.github.io/patRoon/reference/generics.md)
function. The quality metrics are calculated using the MetaClean package
and are useful for filtering out low-quality features before further
analysis.

## Feature qualities

The `featureQualities` function defines quality metrics that are
calculated for individual features. The following quality metrics are
available:

- ApexBoundaryRatio - Ratio between apex and boundary intensities

- FWHM2Base - Full width at half maximum to base ratio

- Jaggedness - Measure of peak smoothness

- Modality - Measure of peak multiplicity

- Symmetry - Measure of peak symmetry (range: -1 to 1)

- GaussianSimilarity - Similarity to Gaussian distribution

- Sharpness - Measure of peak sharpness

- TPASR - Triangle peak area similarity ratio

- ZigZag - Zig-zag index measure

## Feature group qualities

The `featureGroupQualities` function defines quality metrics that are
calculated for feature groups across multiple analyses. The following
quality metrics are available:

- ElutionShift - Measure of retention time consistency across analyses

- RetentionTimeCorrelation - Correlation of retention times across
  analyses

## References

Chetnik K, Petrick L, Pandey G (2020). “MetaClean: a machine
learning-based classifier for reduced false positive peak detection in
untargeted LC-MS metabolomics data.” *Metabolomics*, **16**(11).
[doi:10.1007/s11306-020-01738-3](https://doi.org/10.1007/s11306-020-01738-3)
.

## See also

[`calculatePeakQualities`](https://rickhelmus.github.io/patRoon/reference/generics.md)
