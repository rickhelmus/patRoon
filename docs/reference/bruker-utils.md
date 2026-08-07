# Bruker DataAnalysis utilities

Miscellaneous utility functions which interface with Bruker DataAnalysis

## Usage

``` r
showDataAnalysis()

setDAMethod(anaInfo, method, close = TRUE)

revertDAAnalyses(anaInfo, close = TRUE, save = close)

recalibrarateDAFiles(anaInfo, close = TRUE, save = close)

getDACalibrationError(anaInfo)

addDAEIC(
  analysis,
  path,
  mz,
  mzWindow = defaultLim("mz", "medium"),
  ctype = "EIC",
  mtype = "MS",
  polarity = "both",
  bgsubtr = FALSE,
  fragpath = "",
  name = NULL,
  hideDA = TRUE,
  close = FALSE,
  save = close
)

addAllDAEICs(
  fGroups,
  mzWindow = defaultLim("mz", "medium"),
  ctype = "EIC",
  bgsubtr = FALSE,
  name = TRUE,
  onlyPresent = TRUE,
  hideDA = TRUE,
  close = FALSE,
  save = close
)
```

## Arguments

- anaInfo:

  [Analysis info
  table](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)

- method:

  The full path of the DataAnalysis method.

- close, save:

  If `TRUE` then Bruker files are closed and saved after processing with
  DataAnalysis, respectively. Setting `close=TRUE` prevents that many
  analyses might be opened simultaneously in DataAnalysis, which
  otherwise may use excessive memory or become slow. By default `save`
  is `TRUE` when `close` is `TRUE`, which is likely what you want as
  otherwise any processed data is lost.

- analysis:

  Analysis name (without file extension).

- path:

  path of the analysis.

- mz:

  *m/z* (Da) value used for the chromatographic trace (if applicable).

- mzWindow:

  *m/z* window (in Da) used for the chromatographic trace (if
  applicable).

- ctype:

  Type of the chromatographic trace. Valid options are: `"EIC"`
  (extracted ion chromatogram), `"TIC"` (total ion chromatogram, only
  for `addDAEIC`) and `"BPC"` (Base Peak Chromatogram).

- mtype:

  MS filter for chromatographic trace. Valid values are: `"all"`,
  `"MS"`, `"MSMS"`, `"allMSMS"` and `"BBCID"`.

- polarity:

  Polarity filter for chromatographic trace. Valid values: `"both"`,
  `"positive"` and `"negative"`.

- bgsubtr:

  If `TRUE` then background subtraction ('Spectral' algorithm) will be
  performed.

- fragpath:

  Precursor *m/z* used for MS/MS traces (`""` for none).

- name:

  For `addDAEIC`: the name for the chromatographic trace. For
  `addAllEICs`: `TRUE` to automatically set EIC names. Set to `NULL` for
  none.

- hideDA:

  Hides DataAnalysis while adding the chromatographic trace (faster).

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which EICs should be made.

- onlyPresent:

  If `TRUE` then EICs are only generated for analyses where the feature
  was detected.

## Value

`getDACalibrationError` returns a `data.frame` with a column of all
analyses (named `analysis`) and their mass error (named `error`).

## Details

These functions communicate directly with Bruker DataAnalysis to provide
various functionality, such as calibrating and exporting data and adding
chromatographic traces. For this the RDCOMClient package is required to
be installed.

`showDataAnalysis` makes a hidden DataAnalysis window visible again.
Most functions using DataAnalysis will hide the window during processing
for efficiency reasons. If the window remains hidden (*e.g.* because
there was an error) this function can be used to make it visible again.
This function can also be used to start DataAnalysis if it is not
running yet.

`setDAMethod` Sets a given DataAnalysis method (`.m` file) to a set of
analyses. **NOTE**: as a workaround for a bug in DataAnalysis, this
function will save(!), close and re-open any analyses that are already
open prior to setting the new method. The `close` argument only controls
whether the file should be closed after setting the method (files are
always saved).

`revertDAAnalyses` Reverts a given set of analyses to their unprocessed
raw state.

`recalibrarateDAFiles` Performs automatic mass recalibration of a given
set of analyses. The current method settings for each analyses will be
used.

`getDACalibrationError` is used to obtain the standard deviation of the
current mass calibration (in ppm).

`addDAEIC` adds an Extracted Ion Chromatogram (EIC) or other
chromatographic trace to a given analysis which can be used directly
with DataAnalysis.

`addAllDAEICs` adds Extracted Ion Chromatograms (EICs) for all features
within a
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object.

## See also

[`analysis-information`](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
