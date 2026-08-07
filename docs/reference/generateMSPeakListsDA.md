# Generate peak lists with Bruker DataAnalysis (deprecated)

Uses Bruker DataAnalysis to read the data needed to generate MS peak
lists. This function is now deprecated, please use
[`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md)
instead.

## Usage

``` r
generateMSPeakListsDA(fGroups, ...)

# S4 method for class 'featureGroups'
generateMSPeakListsDA(
  fGroups,
  bgsubtr = TRUE,
  maxMSRTWindow = 5,
  minMSIntensity = 500,
  minMSMSIntensity = 500,
  clear = TRUE,
  close = TRUE,
  save = close,
  MSMSType = "MSMS",
  avgFGroupParams = getDefAvgPListParams()
)

# S4 method for class 'featureGroupsSet'
generateMSPeakListsDA(fGroups, ...)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which MS peak lists should be generated.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- bgsubtr:

  If `TRUE` background will be subtracted using the 'spectral'
  algorithm.

- maxMSRTWindow:

  Maximum chromatographic peak window used for spectrum averaging (in
  seconds, +/- retention time). If `NULL` all spectra from a feature
  will be taken into account. Lower to decrease processing time.

- minMSIntensity, minMSMSIntensity:

  Minimum intensity for peak lists obtained with DataAnalysis. Highly
  recommended to set `>0` as DA tends to report many very low intensity
  peaks.

- clear:

  Remove any existing chromatogram traces/mass spectra prior to making
  new ones.

- close, save:

  If `TRUE` then Bruker files are closed and saved after processing with
  DataAnalysis, respectively. Setting `close=TRUE` prevents that many
  analyses might be opened simultaneously in DataAnalysis, which
  otherwise may use excessive memory or become slow. By default `save`
  is `TRUE` when `close` is `TRUE`, which is likely what you want as
  otherwise any processed data is lost.

- MSMSType:

  The type of MS/MS experiment performed: `"MSMS"` for MRM/AutoMSMS or
  `"BBCID"` for broadband CID.

- avgFGroupParams:

  A `list` with parameters used for averaging of peak lists for feature
  groups. See
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for more details.

## Value

An
[`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
object.

## Details

The MS data should be in the Bruker data format (`.d`). This function
leverages DataAnalysis functionality to support averaging of spectra,
background subtraction and identification of isotopes. In order to
obtain mass spectra TICs will be added in DataAnalysis of the MS and
relevant MS/MS signals.

## Note

The Component column should be active
(Method–\>Parameters–\>Layouts–\>Mass List Layout) in order to add
isotopologue information.

If any errors related to `DCOM` appear it might be necessary to
terminate DataAnalysis (note that DataAnalysis might still be running as
a background process). The `ProcessCleaner` application installed with
DataAnalayis can be used for this.
