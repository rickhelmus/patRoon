# Obtains extracted ion chromatograms (EICs)

These methods generate one or more EIC(s). The `data.table` and
`data.frame` methods generate EIC(s) for given retention time, *m/z* and
optionally mobility ranges. The `features` and `featureGroups` methods
generate EICs for all (or selected) features and feature groups,
respectively, whereby the ranges are automatically determined from the
feature data (see the [EIC
parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
for configuration options).

## Usage

``` r
# S4 method for class 'features'
getEICs(
  obj,
  analysis = analyses(obj),
  EICParams = getDefEICParams(),
  output = "fill"
)

# S4 method for class 'featureGroups'
getEICs(
  obj,
  analysis = analyses(obj),
  groupName = names(obj),
  EICParams = getDefEICParams(),
  output = "fill"
)

# S4 method for class 'data.table'
getEICs(obj, ranges, gapFactor = 3, output = "fill", minIntensityIMS = 25)

# S4 method for class 'data.frame'
getEICs(obj, ...)
```

## Arguments

- obj:

  For the `data.table` and `data.frame` methods: a table with [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).
  For the `features` and `featureGroups` methods: the object for which
  EICs should be generated.

- analysis:

  A `character` vector with the analyses for which EICs should be
  generated.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- output:

  Should be `"fill"`, `"pad"` or `"raw"`. Internally, EIC data is
  compressed by omitting any zero intensity data points. If
  `output="fill"` then the zero intensity points are re-added to obtain
  continuous chromatograms. If `output="pad"` then zero intensity points
  are only re-added that surround others, which is sufficient for *e.g.*
  plotting. If `output="raw"` then the original compressed data is
  returned.

- groupName:

  A `character` vector with the names of the feature groups for which
  EICs should be generated.

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

- minIntensityIMS:

  (**IMS workflow**) Raw intensity threshold for IMS data. This is
  primarily intended to speed up raw data processing.

- ...:

  For the `data.frame` method: further arguments passed to the
  `data.table` method.

## Value

A `list` with for each analysis a `list` with EIC data. For the
`data.table` and `data.frame` methods the EICs are ordered according to
the rows in `ranges`. For the `features` and `featureGroups` methods the
EICs are named after the feature IDs and feature group names,
respectively, and analyses without any EIC data are omitted.

If `output="raw"` then additional columns with *e.g.* mean-averaged and
base peak *m/z* values for each data point are returned. Furthermore,
the `allXValues` attribute is set that can be used to obtain the
original retention time values to reconstruct the original complete
chromatogram.

## Functions

- `getEICs(features)`: Generates EICs for all (or selected) features
  (method for `features`).

- `getEICs(featureGroups)`: Generates EICs for all (or selected) feature
  groups (method for `featureGroups`).

- `getEICs(data.table)`: Generates one or more EIC(s) for given
  retention time, *m/z* and optionally mobility ranges (method for
  `data.table`).

- `getEICs(data.frame)`: Wrapper for the `data.table` method (method for
  `data.frame`).

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `getEICs` to process HRMS (or IMS-HRMS) data. Please
see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.
