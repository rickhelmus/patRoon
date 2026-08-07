# Find features using piek

Uses the `piek` algorithm to find features.

## Usage

``` r
findFeaturesPiek(
  analysisInfo,
  genEICParams = getPiekEICParams(),
  peakParams = getDefPeakParams("chrom", "piek"),
  IMS = FALSE,
  suspects = NULL,
  adduct = NULL,
  assignMethod = "basepeak",
  assignRTWindow = defaultLim("retention", "very_narrow"),
  rtWindowDup = defaultLim("retention", "narrow"),
  mzWindowDup = defaultLim("mz", "medium"),
  mobWindowDup = defaultLim("mobility", "medium"),
  minPeakOverlapDup = 0.25,
  minIntensityIMS = 25,
  EICBatchSize = Inf,
  keepDups = FALSE,
  verbose = TRUE
)

getPiekEICParams(..., IMS = getLimIMS())
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- genEICParams:

  A `list` of parameters for the EIC generation. See the
  `EIC generation parameters` section below. The `getPiekEICParams`
  function is used to generate the parameter list.

- peakParams:

  A `list` of parameters for the peak detection. See
  [`getDefPeakParams`](https://rickhelmus.github.io/patRoon/reference/getDefPeakParams.md)
  for details.

- IMS:

  A `character` that specifies for which IMS instrument defaults are
  returned. Should be `"bruker"` or `"agilent"`. Defaults to what is
  specified in
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md).

- suspects:

  The suspect list to be used for suspect pre-filtering of EIC bins. See
  [suspect
  screening](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
  for details on the suspect list format and `EIC generation parameters`
  to enable suspect filtering.

  **NOTE**: Suspect matching can only be performed by mobilities and not
  CCS values. The
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_susp.md)
  method should be used to convert any CCS data in advance.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. Only needs to be specified if
  `suspects` is set.

- assignMethod:

  Should be `"basepeak"` or `"weighted.mean"`. This parameter sets how
  measured *m/z* or mobilities across the EIC datapoints are handled for
  feature assignment. If `assignMethod="basepeak"`, then the value of
  the base peak (=highest intensity peak) from each EIC datapoint is
  taken. If `assignMethod="weighted.mean"` then the intensity weighted
  mean is calculated of the values that fall within the EIC bin.

- assignRTWindow:

  The retention time window (+/- seconds) used for aggregating EIC
  datapoints to assign feature *m/z* and mobility data, using an
  intensity weighted mean. The maximum window is always bound by the
  feature retention time range. Increasing this number may improve
  accuracy by averaging more points. However, decreasing the window may
  reduce inaccuracies due to inclusion of data from closely eluting
  features (with similar *m/z* and mobility) or noisy data from the
  chromatographic peak extremes. If `assignRTWindow=0` then only the EIC
  datapoint at the feature retention time is used.

  The assignment window is automatically adjusted for the values set for
  `sumWindowMZ` and `sumWindowMob` (see `EIC generation parameters`).

- rtWindowDup, mzWindowDup, mobWindowDup:

  The retention time (seconds), *m/z* and mobility windows used to
  identify duplicate (redundant) features detected in multiple EIC bins.
  These values default to `defaultLim("retention", "very_narrow")`,
  `defaultLim("mz", "medium")` and `defaultLim("mobility", "medium")`,
  respectively (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- minPeakOverlapDup:

  The minimum overlap (fraction between 0 and 1) in retention time
  between two features to be considered a duplicate.

- minIntensityIMS:

  (**IMS workflow**) Raw intensity threshold for IMS data. This is
  primarily intended to speed up raw data processing.

- EICBatchSize:

  The number of EICs to be processed in a single batch. Decreasing this
  number will reduce memory usage, at the cost of speed. Set to `Inf` to
  process all EICs in a single batch.

- keepDups:

  Set to `TRUE` to keep duplicate features and features with
  non-centered *m/z* or mobility values. This is primarily intended for
  debugging, but can be useful to investigate why features are missing
  or optimize tolerance windows for duplicate feature detection.

- verbose:

  If set to `FALSE` then no text output is shown.

- ...:

  Any additional parameters to be set in the returned parameter list.
  These will override the defaults. See the `EIC generation parameters`
  section for details.

## Value

`getPiekEICParams` returns a `list` of parameters for the EIC
generation, which is used to set the `genEICParams` argument to
`findFeaturesPiek`.

## Details

This function uses piek to automatically find features. This function is
called when calling `findFeatures` with `algorithm="piek"`.

The `piek` algorithm extends and improves on the simple and fast feature
detection algorithm introduced by Dietrich C, Wick A, Ternes TA (2021).
“Open‐source feature detection for non‐target LC–MS analytics.” *Rapid
Communications in Mass Spectrometry*, **36**(2). ISSN 1097-0231.
[doi:10.1002/rcm.9206](https://doi.org/10.1002/rcm.9206) .
<http://dx.doi.org/10.1002/rcm.9206>. . This algorithm first forms
extracted ion chromatograms (EICs) and subsequently performs automatic
peak detection to generate features. The `piek` algorithm introduces the
following improvements and changes:

- Support for IMS-HRMS workflows.

- The [msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md)
  interface is used to efficiently form EICs from the raw data. All the
  file formats and types can be used that are supported by
  [msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md).
  This includes IMS data, even if not used for feature detection, which
  allows the use of IMS data directly in non-IMS or [post mobility
  assignment](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)
  workflows.

- The EIC binning approach can be extended with the mobility dimension
  to support [direct mobility
  assignment](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)
  workflows.

- The EIC bins can be filtered with suspect or MS2 data to speed up
  feature detection.

- Several filters are available to eliminate EICs with are likely devoid
  of any signal of interest.

- The original peak detection algorithm was further optimized or can be
  be exchanged with others: see
  [`getDefPeakParams`](https://rickhelmus.github.io/patRoon/reference/getDefPeakParams.md)
  for details.

- Several filters are available to improve the data and reduce
  redundancy:

  - The original redundancy detection, which performs a second feature
    detection with EIC bins that are shifted by 50% width and eliminates
    features with `m/z` values outside the center of any bin, was
    extended for IMS support.

  - Redundant features across bins are eliminated if with close
    retention time, `m/z`, mobility and chromatographic overlap. The
    most intense feature is kept.

  - Data from suspects or MS2 precursors that was used to pre-filter
    EICs, can also be used to filter the final feature list.

- Various small bug fixes and improvements for the original code.

- The output feature tables contain raw intensities/areas and those
  subtracted by the estimated noise level (`intensity`, `intensitySub`,
  `area` and `areaSub` columns, respectively) and the estimated signal
  to noise (`signalToNoise` column).

## IMS workflows

If IMS data is used to resolve features (`IMS=TRUE`), a 'pre-check' is
performed to avoid excessive numbers of two-dimensional bins for EIC
formation and peak detection. These EICs are formed by only considering
the m/z dimension, and subsequently filtered by the parameters described
in the `EIC generation parameters` section. The final EICs for feature
detection are then only formed if they have m/z data that was not
removed during the pre-check.

The `m/z` and mobility data from IMS-HRMS data is typically not or
partially centroided. The feature `m/z` and mobility values are derived
from `m/z` or mobility *versus* intensity profiles. The profiles are
generated for each EIC timepoint, and the value at the maximum intensity
or intensity weighted mean of the profile is used to derive the
intermediate values (configured by `assignMethod`). Several parameters
exist to improve the profile data (see next section).

## EIC generation parameters

The `genEICParams` argument to `findFeaturesPiek` configures the
generation of EICs. The `getPiekEICParams` function should be used to
generate the parameter list.

The following general parameters exist:

- `filter` Controls the pre-filtering of EIC bins with *m/z* data.
  Should be `"none"` (no filtering), `"suspects"` (filter with suspect
  data) or `"ms2"` (filter with data from precursors detected in a
  data-dependent MS/MS experiment).

- `mzRange`,`mzStep` Configures the formation of *m/z* bins. `mzRange`
  is a numeric vector of length two that specifies the min/max *m/z*
  range. `mzStep` specifies the bin widths.

- `retRange` A `numeric` vector of length two that specifies the
  retention time range for the EICs. Data outside this range is
  excluded. Set to `NULL` to use the full range.

- `gapFactor` A `numeric` that configures gap filling for EICs. See
  [`getDefEICParams`](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  for further details.

- `minEICIntensity` The minimum intensity of the highest data point in
  the EIC. Used to filter EICs.

- `minEICAdjTime`,`minEICAdjPoints`,`minEICAdjIntensity` The EIC should
  have at least a continuous signal of `minEICAdjTime` seconds and
  `minEICAdjPoints` data points, where the continuity is defined by data
  points with an intensity of at least `minEICAdjIntensity` high. Set
  `minEICAdjTime` or `minEICAdjPoints` to zero to disable continuity
  checks for time or data points, respectively. Set `minEICAdjIntensity`
  to zero to completely disable continuity checks.

- `topMostEICMZ` Only keep this number of top-most intense EICs. The
  intensity is derived from the data point with the highest intensity in
  the EIC. Set to zero to always select all EICs.

  For IMS workflows, this parameter is *only* used to limit the number
  of EICs resulting from the 'pre-check' in the *m/z* dimension.

The following parameters are specifically used for IMS workflows:

- `filterIMS` Similar to the `filter` parameter, but controls how
  mobility data is used for pre-filtering of EIC bins.

  Different values for `filter` and `filterIMS` can be specified:

  - `filter="none"` and `filterIMS="none"`

  - `filter="suspects"` and `filterIMS="suspects"`

  - `filter="suspects"` and `filterIMS="none"` (only use *m/z*
    filtering)

  - `filter="ms2"` and `filterIMS="ms2"`

  - `filter="ms2"` and `filterIMS="none"`

  Currently only Bruker DDA-PASEF experiments provide the data needed
  for `"ms2"` filtering.

- `mobRange`,`mobStep` Equivalent to `mzRange` and `mzStep`, but for ion
  mobility binning.

- `sumWindowMZ`,`sumWindowMob` The retention time window (+/- s) used to
  sum adjacent datapoints for the determination of intermediate EIC
  *m/z* and mobility values. This data is aggregated to determine the
  final feature values (see also the `assignRTWindow` argument). Set to
  `0` to not sum any adjacent timepoints. Larger values can generally
  improve accuracy for noisy data (*e.g.* from TIMS), but care must be
  taken to stay below the expected minimum chromatographic peak width to
  avoid inclusion of data from other features. Defaults to
  `defaultLim("retention", "very_narrow")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `smoothWindowMZ`,`smoothWindowMob` The window size used to perform
  centered moving average smoothing on intensity data of the *m/z* and
  mobility profiles used to determine intermediate EIC values. Smoothing
  of noisy data (*e.g.* TIMS) is highly recommended to improve accuracy
  and consistency. Set to `0` to disable smoothing.

- `smoothExtMZ`,`smoothExtMob` The `m/z` or mobility window to extend
  the smoothing at the edges of the EIC bin. This is recommended to
  improve smoothing, *e.g.* when the peak profile is only partially
  captured in the bin. Defaults to the bin width, *i.e.* data from an
  adjacent bin on each side is additionally included for smoothing. The
  final smoothed data is only taken from the actual EIC bin. Set to `0`
  to disable extension.

- `saveMZProfiles`,`saveEIMs` Set to `TRUE` to save the *m/z* and
  mobility profiles for each feature. Only the profiles at the feature
  retention time is saved. This can be useful for debugging or parameter
  optimization, but will increase memory usage and processing times.

- `topMostEICMob` Equivalent to `topMostEICMZ`, used to reduce the final
  two-dimensional EIC bins with `m/z` and mobility information.

- `minEICsIMSPreCheck` Only perform the `m/z` pre-check if the number of
  two-dimensional EIC bins is at least `minEICsIMSPreCheck`.

The following parameters are specifically for when suspect data is used
to pre-filter EIC bins:

- `rtWindow`,`mzwindow`,`mobWindow`: The retention time, *m/z* and
  mobility tolerance windows for suspect data. These are used for:

  1.  Pre-filtering of EIC bins with suspect data, *i.e.* larger
      tolerances will lead to more EIC bins being kept. (only applicable
      for `mzWindow` and `mobWindow`).

  2.  Matching the final features to suspect data. `rtWindow=Inf` can be
      used to disable retention time matching.

  Defaults to `defaultLim("retention", "medium")`,
  `defaultLim("mz", "medium")` and `defaultLim("mobility", "medium")`,
  see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md).

- `skipInvalid`,`prefCalcChemProps`,`neutralChemProps` Controls
  preparing the suspect list data. See
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).

The following parameters are specifically for when MS2 data is used to
pre-filter EICs:

- `rtWindow` Eliminates any features without an MS/MS spectrum within
  this retention time window. Set `rtWindow=Inf` to disable this filter.
  Defaults to `defaultLim("retention", "very_narrow")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `mzIsoWindow` The maximum *m/z* window considered for MS/MS precursors
  that were isolated by DDA. These *m/z* isolation windows are used to
  pre-filter EICs and match the final features. Setting `mzIsoWindow` to
  a value lower than typical instrument isolation windows will make
  feature detection more specific, as features need to be more close to
  the triggered DDA precursor `m/z` values. In contrast, larger values
  for `mzIsoWindow` allows to include features that were not
  specifically targeted by DDA, but may still have MS/MS data as their
  *m/z* could still fall within the MS/MS isolation window. The
  effective window used will never exceed the instrumental isolation
  window. Setting `mzIsoWindow=Inf` will always use instrumental
  windows.

  **NOTE**: Sometimes the isolation windows are not exported and cannot
  be deduced automatically (e.g. Agilent data). In that case, the
  `mzIsoWindow` parameter is used as the isolation window and therefore
  needs to be set to a finite value.

- `mobWindow` The mobility tolerance window to match DDA MS/MS
  precursors in IMS workflows. Used for pre-filtering EICs and the final
  features. To match DDA precursor data, the measured mobility range of
  the corresponding MS/MS data is used as the mobility window. This
  window is then adjusted to be at least +/- `mobWindow`. Defaults to
  `defaultLim("mobility", "medium")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md))

- `minTIC` The minimum total ion current (TIC) signal for an MS/MS
  spectrum to be considered. Can be increased to eliminate features with
  low intensity MS/MS data.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `findFeaturesPiek` to process HRMS (or IMS-HRMS)
data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

The use of profile *m/z* HRMS data (not IMS-HRMS) is currently not
supported.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
