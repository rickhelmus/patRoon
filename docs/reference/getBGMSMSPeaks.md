# Background MS/MS peak detection

Detects background MS/MS peaks by gathering frequently occurring peaks
in MS/MS spectra from blanks.

## Usage

``` r
getBGMSMSPeaks(
  anaInfo,
  replicates = NULL,
  MSLevel = 2,
  retentionRange = NULL,
  mobilityRange = NULL,
  minBPIntensity = 5000,
  avgSpectraParams = getDefAvgPListParams(relMinAbundance = 0.1, topMost = 25),
  avgAnalysesParams = getDefAvgPListParams(relMinAbundance = 0.8, topMost = 25)
)
```

## Arguments

- anaInfo:

  The [analysis
  info](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  object with the analyses to be used for background peak detection.

- replicates:

  A `character` with names of replicates for the analyses used for
  background peak detection. If `NULL` then all the analyses defined in
  `anaInfo` will be used.

- MSLevel:

  The MS level of the spectra to be used for background peak detection.
  This should be `1` or `2`. This function is only tested with
  `MSLevel=2`. And while `MSLevel=1` is possible, it may be rather
  computationally intensive due to the relatively large number of MS
  peaks to iterate through.

- retentionRange, mobilityRange:

  A two-sized `numeric` vector with the minimum and maximum retention
  time and ion mobility to consider, respectively. If `NULL` then the
  full range is used.

- minBPIntensity:

  The minimum basepeak intensity of a spectrum to be considered for
  background peak detection. This is primarily intended to optimize the
  detection procedure.

- avgSpectraParams, avgAnalysesParams:

  A `list` with parameters used for averaging the MS spectra within and
  between analyses, respectively. See
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for more details.

## Value

A `data.table` with the detected background peaks and abundance
statistics. The table can directly be passed to the `removeMZs` argument
of the
[`filter method for MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

## Details

This function iterates through all MS/MS spectra of the given analyses
and collects the most frequently occurring peaks. It first averages all
spectra within the same analyses, and retains those peaks above an
abundance threshold (set by `avgSpectraParams`). The analyses-averaged
spectra are then also averaged, and again peaks with a minimum abundance
are kept (set by `avgAnalysesParams`).

The frequent occurrence of a peak throughout MS/MS spectra, including
DDA spectra different isolation m/z values, is often a sign of
contamination from the analytical system (*e.g.* from the mobile phase,
MS quadrupoles etc). Hence, the analyses used for blank subtraction are
typically just measurements of *e.g.* ultrapure water or another
solvent, and not need to be treated by any extraction procedure.

The output of this function can directly be used to the
[`filter method for MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
to subsequently remove these background peaks, which possibly improves
subsequent feature annotation.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `getBGMSMSPeaks` to process HRMS (or IMS-HRMS) data.
Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

The use of profile *m/z* HRMS data (not IMS-HRMS) is currently not
supported.

## References

Helmus R, Bagdonaite I, de Voogt P, van Bommel MR, Schymanski EL, van
Wezel AP, ter Laak TL (2025). “Comprehensive Mass Spectrometry Workflows
to Systematically Elucidate Transformation Processes of Organic
Micropollutants: A Case Study on the Photodegradation of Four
Pharmaceuticals.” *Environmental Science & Technology*, **59**(7),
3723–3736. ISSN 1520-5851.
[doi:10.1021/acs.est.4c09121](https://doi.org/10.1021/acs.est.4c09121) .
<http://dx.doi.org/10.1021/acs.est.4c09121>.
