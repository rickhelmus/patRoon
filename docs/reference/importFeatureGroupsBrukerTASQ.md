# Imports feature groups from Bruker TASQ

Imports screening results from Bruker TASQ as feature groups.

## Usage

``` r
importFeatureGroupsBrukerTASQ(
  input,
  analysisInfo,
  clusterRTWindow = defaultLim("retention", "medium")
)
```

## Arguments

- input:

  The file path to an Excel export of the Global results table from
  TASQ, converted to `.csv` format.

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- clusterRTWindow:

  This retention time window (in seconds) is used to group hits across
  analyses together. See also the details section.

## Value

A new `featureGroups` object containing converted screening results from
Bruker TASQ.

## Details

This function imports data from Bruker TASQ. This function is called
when calling `importFeatureGroups` with `type="brukertasq"`.

The feature groups across analyses are formed based on the name of
suspects and their closeness in retention time. The latter is necessary
because TASQ does not necessarily perform checks on retention times and
may therefore assign a suspect to peaks with different retention times
across analyses (or within a single analysis). Hence, suspects with
equal names are hierarchically clustered on their retention times (using
fastcluster) to form the feature groups. The cut-off value for this is
specified by the `clusterRTWindow` argument. The input for this function
is obtained by generating an Excel export of the 'global' results and
subsequently converting the file to `.csv` format.

## Note

This function uses estimated min/max values for retention times and
dummy min/max *m/z* values for conversion to features, since this
information is not (readily) available. Hence, when plotting, for
instance, extracted ion chromatograms (with
[`plotChroms`](https://rickhelmus.github.io/patRoon/reference/generics.md))
the integrated chromatographic peak range shown is incorrect.

This function may use suspect names to base file names used for
reporting, logging etc. Therefore, it is important that these are
file-compatible names.

## References

Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python.” *Journal of Statistical
Software*, **53**(9), 1–18.
[doi:10.18637/jss.v053.i09](https://doi.org/10.18637/jss.v053.i09) .

## See also

[`importFeatureGroups`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroups.md)
for more details and other algorithms.
