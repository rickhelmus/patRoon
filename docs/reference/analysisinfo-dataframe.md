# AnalysisInfo data.frame methods

Various parsing and plotting functions for the analysisInfo data.frame.

## Usage

``` r
# S4 method for class 'data.frame'
getTICs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'data.frame'
getBPCs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'data.frame'
plotTICs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S4 method for class 'data.frame'
plotBPCs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S4 method for class 'data.table'
plotTICs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S4 method for class 'data.table'
plotBPCs(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  groupBy = NULL,
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- obj:

  An `analysisInfo` data.frame object as obtained by
  [generateAnalysisInfo](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  function.

- retentionRange:

  Range of retention time (in seconds), m/z, respectively. Should be a
  numeric vector with length of two containing the min/max values. The
  maximum can be Inf to specify no maximum range. Set to NULL to skip
  this step.

- MSLevel:

  Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to
  obtain traces.

- retMin:

  Plot retention time in minutes (instead of seconds).

- title:

  Character string used for title of the plot. If `NULL` a title will be
  automatically generated.

- groupBy:

  Specifies how results are grouped in the plot. Should be a name of a
  column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  table which is used to make analysis groups (*e.g.* `"replicate"`), or
  `"fGroups"` to group by feature group. Set to `NULL` for no grouping.

- showLegend:

  Plot a legend if TRUE.

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- ...:

  Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Functions

- `getTICs(data.frame)`: Obtain the total ion chromatogram/s (TICs) of
  the analyses.

- `getBPCs(data.frame)`: Obtain the base peak chromatogram/s (BPCs) of
  the analyses.

- `plotTICs(data.frame)`: Plots the TICs of the analyses.

- `plotBPCs(data.frame)`: Plots the BPCs of the analyses.

- `plotBPCs(data.table)`: Plots the BPCs of the analyses.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by these functions to process HRMS (or IMS-HRMS) data.
Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## Author

Ricardo Cunha (<cunha@iuta.de>) and Rick Helmus (<r.helmus@uva.nl>)
