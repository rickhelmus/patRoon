# Extracted Ion Chromatogram and Mobilogram parameters

Parameters for creation of extracted ion chromatograms and mobilograms.

## Usage

``` r
getDefEICParams(...)

getDefEIMParams(..., IMS = getLimIMS())
```

## Arguments

- ...:

  optional named arguments that override defaults.

- IMS:

  A `character` that specifies for which IMS instrument defaults are
  returned. Should be `"bruker"` or `"agilent"`. Defaults to what is
  specified in
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md).

## Details

The following parameters exist to configure the creation of extracted
ion chromatograms (EICs) and extracted ion mobilograms (EIMs):

- `window` A value that is subtracted or added to the minimum and
  maximum retention time (EICs) or mobility (EIMs) of the feature. Thus,
  setting this value to `>0` will 'zoom out' on the x-axis of a
  chromatogram or mobilogram. Defaults to
  `defaultLim("retention", "wide")` (EICs) or
  `defaultLim("mobility", "wide")` (EIMs) (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `topMost` Only create EICs/EIMs for this number of top most intense
  features. If `NULL` then EICs/EIMs are created for all features.

- `topMostByReplicate` If set to `TRUE` and `topMost` is set: only
  create EICs/EIMs for the top most features in each replicate. For
  instance, when `topMost=1` and `topMostByReplicate=TRUE`, then only
  the most intense feature of each replicate is considered.

- `onlyPresent` If `TRUE` then EICs/EIMs are created only for analyses
  in which a feature was detected, if `onlyPresent=FALSE` then data is
  generated for **all** analyses. The latter is handy to evaluate if a
  peak was 'missed' during peak detection or removed during *e.g.*
  filtering.

- `mzExpMobWindow` (**IMS workflow**) Additional *m/z* tolerance on top
  of the feature limits. This is for IMS workflows where features were
  detected from centroided LC-MS like data, while EICs/EIMs are
  generated from raw IMS data. In this case the feature *m/z* limits
  were derived from centroided data, which typically has smaller *m/z*
  deviations across scans compared to IMS data. The `mzExpMobWindow`
  parameter sets an additional *m/z* tolerance to specifically handle
  this case. Defaults to `defaultLim("mz", "default")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `minIntensityIMS` (**IMS workflow**) Raw intensity threshold for IMS
  data. This is primarily intended to speed up raw data processing.

if `onlyPresent=FALSE` then the following parameters are also relevant:

- `mzExpWindow`,`mobExpWindow` To create EICs or EIMs for analyses in
  which no feature was found, the *m/z* or mobility value is derived
  from the min/max values of all features in the feature group. The
  value of `mzExpWindow` and `mobExpWindow` further expands this window
  to allow a greater tolerance. Defaults to
  `defaultLim("mz", "very_narrow")` and
  `defaultLim("mobility", "very_narrow")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `setsAdductPos`,`setsAdductNeg` (**sets workflow**) In sets workflows
  the adduct must be known to calculate the ionized *m/z*. If a feature
  is completely absent in a particular set then it follows no adduct
  annotations are available and the value of `setsAdductPos` (positive
  ionization data) or `setsAdductNeg` (negative ionization data) will be
  used instead.

The following additional parameters exist specifically for EICs
(`EICParams`):

- `gapFactor` Bruker TIMS data (and maybe others?) seem to omit zero
  intensity scans, which will lead to time gaps between spectra and
  incorrect EICs. To determine a time gap, the `gapFactor` is multiplied
  with the median of time differences between scans. If a gap is
  detected, then appropriate zero intensity points are added to the EIC.
  Set to `0` to disable this.

The following additional parameters exist specifically for EIMs
(`EIMParams`):

- `maxRTWindow` Maximum retention time window (seconds, +/- feature
  retention time) in which mobilograms are collected and averaged.
  Defaults to `defaultLim("retention", "very_narrow")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `smooth` The smoothing method that is applied to the EIM. Can be
  `"none"` for no smoothing, `"sg"` for Savitzky-Golay (using
  [signal](https://CRAN.R-project.org/package=signal))) or `"ma"` for
  centered moving average (same algorithm as used by
  [`findFeaturesPiek`](https://rickhelmus.github.io/patRoon/reference/findFeaturesPiek.md)).

- `smLength` The smoothing length. If `smooth="sg"` then this is passed
  as the `n` argument to the
  [`signal::sgolayfilt`](https://rdrr.io/pkg/signal/man/sgolayfilt.html)
  function.

- `sgOrder` The smoothing order for Savitzky-Golay. Passed as the `p`
  argument to the
  [`signal::sgolayfilt`](https://rdrr.io/pkg/signal/man/sgolayfilt.html)
  function.

These parameters are passed as a named `list` as the `EICParams` or
`EIMParams` argument to functions that work with EIC or EIM data. The
`getDefEICParams` and `getDefEIMParams` functions generate such
parameter list with defaults.
