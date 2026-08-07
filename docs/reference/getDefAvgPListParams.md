# Parameters for averaging MS peak list data

Create parameter lists for averaging MS peak list data.

## Usage

``` r
getDefAvgPListParams(..., IMS = getLimIMS())
```

## Arguments

- ...:

  Optional named arguments that override defaults.

- IMS:

  A `character` that specifies for which IMS instrument defaults are
  returned. Should be `"bruker"` or `"agilent"`. Defaults to what is
  specified in
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md).

## Value

`getDefAvgPListParams` returns a `list` with the peak list averaging
parameters.

## Details

The parameters set used for averaging peak lists are set by the
`avgFeatParams` and `avgFGroupParams` arguments to
[`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md)
and its related algorithm specific functions. The parameters are
specified as a named `list` with the following values:

- `method`,`clusterMzWindow` The cluster method and window (see
  [clustering
  parameters](https://rickhelmus.github.io/patRoon/reference/cluster-params.md))
  used to average mass spectra. `clusterMzWindow` is defaulted as
  `defaultLim("mz", "medium")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `topMost` Only retain this maximum number of MS peaks when generating
  averaged spectra. Lowering this number may exclude more irrelevant
  (noisy) MS peaks and decrease processing time, whereas higher values
  may avoid excluding lower intense MS peaks that may still be of
  interest.

- `minIntensityPre` MS peaks with intensities below this value will be
  removed (applied prior to selection by `topMost` and averaging).

- `minIntensityPost` MS peaks with intensities below this value will be
  removed (after averaging).

- `minIntensityIMS` MS peaks in spectra of raw IMS frames with
  intensities below this value will be removed (applied prior to any
  other treatment steps).

- `absMinAbundance`,`relMinAbundance` Minimum absolute/relative
  abundance of an MS peak across the spectra that are averaged. If
  `absMinAbundance` exceeds the number of spectra then the threshold is
  automatically lowered to the number of spectra.

- `maxRelCumIntensity` The maximum relative cumulative intensity of a
  peak, calculated in descending order (most intense peaks first). Set
  to `1` to disable.

- `smoothWindowIMS`,`halfWindowIMS`,`maxGapIMS` Parameters used for
  centroiding *m/z* peaks from IMS-HRMS data. See `Centroiding IMS data`
  for more details.

- `withPrecursorMS` For MS data only: ignore any spectra that do not
  contain the precursor peak.

  For IMS data this excludes MS spectra within an IMS frame that do not
  contain the precursor peak, typically due to mobility separation.
  Hence, setting this option performs some crude cleanup of MS spectra,
  even for features for which no mobilities were assigned (*e.g.*
  non-IMS workflows).

- `pruneMissingPrecursorMS` For MS data only: if `TRUE` then peak lists
  without a precursor peak are removed. Note that even when this is set
  to `FALSE`, functionality that relies on MS (not MS/MS) peak lists
  (*e.g.* formulae calculation) will still skip calculation if a
  precursor is not found.

- `retainPrecursorMSMS` For MS/MS data only: if `TRUE` then always
  retain the precursor mass peak even if is not amongst the `topMost`
  peaks. Note that MS precursor mass peaks are always kept. Furthermore,
  note that precursor peaks in both MS and MS/MS data may still be
  removed by intensity thresholds (this is unlike the
  [`filter`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  method function).

The `getDefAvgPListParams` function can be used to generate a default
parameter list. The defaults are (with `IMS="bruker"`):

    list(
      method = "distance_mean",
      clusterMzWindow = 0.005,
      topMost = 50,
      minIntensityPre = 500,
      minIntensityPost = 500,
      minIntensityIMS = 25,
      absMinAbundance = 0,
      relMinAbundance = 0,
      maxRelCumIntensity = 1,
      smoothWindowIMS = 0,
      halfWindowIMS = 2,
      maxGapIMS = 0.005,
      withPrecursorMS = TRUE,
      pruneMissingPrecursorMS = TRUE,
      retainPrecursorMSMS = TRUE
    )

## Centroiding IMS data

With IMS-HRMS data the *m/z* peaks are often not or partially
centroided. The following steps are performed to centroid the data:

1.  Sum up mass spectra within an IMS frame. If the feature has mobility
    data, only spectra within its mobility boundaries are considered.

2.  Use point-distance clustering (see [clustering
    parameters](https://rickhelmus.github.io/patRoon/reference/cluster-params.md))
    with a window defined by `maxGapIMS` to find related mass signals.
    This is primarily meant for non-continuous data, *e.g.* due to
    intensity thresholding. The `maxGapIMS` parameter should be set to a
    value that represents the maximum expected distance between two
    *m/z* datapoints. For some instruments, such as Agilent IMS-QTOF,
    this value may be higher than expected. For that reason, if
    `IMS="agilent"` then the default is set to `0.01`.

3.  Smooth the intensity data using a centered moving average with
    window size `smoothWindowIMS` (set to zero to disable smoothing).

4.  Find local maxima within sliding window with +/- `halfWindowIMS`
    points and eliminate non-centroids. This algorithm is based on the
    `C_localMaxima` function from
    [MALDIquant](https://CRAN.R-project.org/package=MALDIquant).

## References

Gibb S, Strimmer K (2012). “MALDIquant: a versatile R package for the
analysis of mass spectrometry data.” *Bioinformatics*, **28**(17),
2270–2271.
[doi:10.1093/bioinformatics/bts447](https://doi.org/10.1093/bioinformatics/bts447)
.
