# MS spectral similarity calculation parameters

Parameters relevant for calculation of similarities between mass
spectra.

## Usage

``` r
getDefSpecSimParams(...)
```

## Arguments

- ...:

  optional named arguments that override defaults.

## Details

For the calculation of spectral similarities the following parameters
exist:

- `method` The similarity method: either `"cosine"` or `"jaccard"`.

- `removePrecursor` If `TRUE` then precursor peaks (*i.e.* the mass peak
  corresponding to the feature) are removed prior to similarity
  calculation.

- `mzWeight`,`intWeight` Mass and intensity weights used for cosine
  calculation.

- `absMzDev` Maximum absolute *m/z* deviation between mass peaks, used
  for binning spectra. Defaults to `defaultLim("mz", "medium")` (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- `relMinIntensity` The minimum relative intensity for mass peaks
  (`0-1`). Peaks with lower intensities are not considered for
  similarity calculation. The relative intensities are called after the
  precursor peak is removed when `removePrecursor=TRUE`.

- `minPeaks` Only consider spectra that have at least this amount of
  peaks (*after* the spectrum is filtered).

- `shift` If and how shifting is applied prior to similarity
  calculation. Valid options are: `"none"` (no shifting), `"precursor"`
  (all mass peaks of the second spectrum are shifted by the mass
  difference between the precursors of both spectra) or `"both"` (the
  spectra are first binned without shifting, and peaks still unaligned
  are then shifted as is done when `shift="precursor"`).

- `setCombinedMethod` (**sets workflow**) Determines how spectral
  similarities from different sets are combined. Possible values are
  `"mean"`, `"min"` or `"max"`, which calculates the combined value as
  the mean, minimum or maximum value, respectively. `NA` values (*e.g.*
  if a set does not have peak list data to combine) are removed in
  advance.

These parameters are typically passed as a named `list` as the
`specSimParams` argument to functions that do spectral similarity
calculations. The `getDefSpecSimParams` function can be used to generate
such parameter list with defaults.
