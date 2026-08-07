# Find features using Bruker DataAnalysis

Uses the 'Find Molecular Features' (FMF) algorithm of Bruker
DataAnalysis vendor software to find features.

## Usage

``` r
findFeaturesBruker(
  analysisInfo,
  doFMF = "auto",
  startRange = 0,
  endRange = 0,
  save = TRUE,
  close = save,
  verbose = TRUE
)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- doFMF:

  Run the 'Find Molecular Features' algorithm before loading compounds.
  Valid options are: `"auto"` (run FMF automatically if current results
  indicate it is necessary) and `"force"` (run FMF *always*, even if
  cached results exist). Note that checks done if `doFMF="auto"` are
  fairly simplistic, hence set `doFMF="force"` if feature data needs to
  be updated.

- startRange, endRange:

  Start/End retention range (seconds) from which to collect features. A
  0 (zero) for `endRange` marks the end of the analysis.

- close, save:

  If `TRUE` then Bruker files are closed and saved after processing with
  DataAnalysis, respectively. Setting `close=TRUE` prevents that many
  analyses might be opened simultaneously in DataAnalysis, which
  otherwise may use excessive memory or become slow. By default `save`
  is `TRUE` when `close` is `TRUE`, which is likely what you want as
  otherwise any processed data is lost.

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses Bruker to automatically find features. This function
is called when calling `findFeatures` with `algorithm="bruker"`.

The resulting 'compounds' are transferred from DataAnalysis and stored
as features.

This algorithm only works with Bruker data files (`.d` extension) and
requires Bruker DataAnalysis and the RDCOMClient package to be
installed. Furthermore, DataAnalysis combines multiple related masses in
a feature (*e.g.* isotopes, adducts) but does not report the actual
(monoisotopic) mass of the feature. Therefore, it is simply assumed that
the feature mass equals that of the highest intensity mass peak.

## Note

If any errors related to `DCOM` appear it might be necessary to
terminate DataAnalysis (note that DataAnalysis might still be running as
a background process). The `ProcessCleaner` application installed with
DataAnalayis can be used for this.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
