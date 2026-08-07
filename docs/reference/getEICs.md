# Obtains extracted ion chromatograms (EICs)

This function generates one or more EIC(s) for given retention time,
*m/z* and optionally mobility ranges.

## Usage

``` r
getEICs(
  analysisInfo,
  ranges,
  gapFactor = 3,
  output = "fill",
  minIntensityIMS = 25
)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- ranges:

  A `list` with for each analysis a `data.frame` with `numeric` columns
  `"retmin"`, `"retmax"`, `"mzmin"`, `"mzmax"` with the lower/upper
  ranges of the retention time and *m/z*. Furthermore, columns
  `"mobmin"` and `"mobmax"` can be added for mobility lower/upper ranges
  in IMS data.

- gapFactor:

  A `numeric` that configures gap filling. See
  [`getDefEICParams`](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  for more details.

- output:

  Should be `"fill"`, `"pad"` or `"raw"`. Internally, EIC data is
  compressed by omitting any zero intensity data points. If
  `output="fill"` then the zero intensity points are re-added to obtain
  continuous chromatograms. If `output="pad"` then zero intensity points
  are only re-added that surround others, which is sufficient for *e.g.*
  plotting. If `output="raw"` then the original compressed data is
  returned.

- minIntensityIMS:

  (**IMS workflow**) Raw intensity threshold for IMS data. This is
  primarily intended to speed up raw data processing.

## Value

A `list` with for each analysis a `list` with EIC data for each of the
rows in `ranges`.

If `output="raw"` then additional columns with *e.g.* mean-averaged and
base peak *m/z* values for each data point are returned. Furthermore,
the `allXValues` attribute is set that can be used to obtain the
original retention time values to reconstruct the original complete
chromatogram.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `getEICs` to process HRMS (or IMS-HRMS) data. Please
see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.
