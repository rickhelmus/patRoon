# Find features using OpenMS

uses the
[FeatureFinderMetabo](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMetabo.html)
TOPP tool (see <http://www.openms.de>) to find features.

## Usage

``` r
findFeaturesOpenMS(
  analysisInfo,
  noiseThrInt = 1000,
  chromSNR = 3,
  chromFWHM = 5,
  mzPPM = defaultLim("mz", "medium_rel"),
  reEstimateMTSD = TRUE,
  traceTermCriterion = "sample_rate",
  traceTermOutliers = 5,
  minSampleRate = 0.5,
  minTraceLength = 3,
  maxTraceLength = -1,
  widthFiltering = "fixed",
  minFWHM = 1,
  maxFWHM = 30,
  traceSNRFiltering = FALSE,
  localRTRange = 10,
  localMZRange = 6.5,
  isotopeFilteringModel = "metabolites (5% RMS)",
  MZScoring13C = FALSE,
  useSmoothedInts = TRUE,
  extraOpts = NULL,
  useFFMIntensities = FALSE,
  verbose = TRUE
)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- noiseThrInt:

  Noise intensity threshold. Sets `algorithm:common:noise_threshold_int`
  option.

- chromSNR:

  Minimum S/N of a mass trace. Sets `algorithm:common:chrom_peak_snr`
  option.

- chromFWHM:

  Expected chromatographic peak width (in seconds). Sets
  `algorithm:common:chrom_fwhm` option.

- mzPPM:

  Allowed mass deviation (ppm) for trace detection. Sets
  `algorithm:mtd:mass_error_ppm`.

- reEstimateMTSD:

  If `TRUE` then enables dynamic re-estimation of m/z variance during
  mass trace collection stage. Sets `algorithm:mtd:reestimate_mt_sd`.

- traceTermCriterion, traceTermOutliers, minSampleRate:

  Termination criterion for the extension of mass traces. See
  [FeatureFinderMetabo](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMetabo.html).
  Sets the `algorithm:mtd:trace_termination_criterion`,
  `algorithm:mtd:trace_termination_outliers` and
  `algorithm:mtd:min_sample_rate` options, respectively.

- minTraceLength, maxTraceLength:

  Minimum/Maximum length of mass trace (seconds). Set negative value for
  maxlength to disable maximum. Sets `algorithm:mtd:min_trace_length`
  and `algorithm:mtd:min_trace_length`, respectively.

- widthFiltering, minFWHM, maxFWHM:

  Enable filtering of unlikely peak widths. See
  [FeatureFinderMetabo](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMetabo.html).
  Sets `algorithm:epd:width_filtering`, `algorithm:epd:min_fwhm` and
  `algorithm:epd:max_fwhm`, respectively.

- traceSNRFiltering:

  If `TRUE` then apply post-filtering by signal-to-noise ratio after
  smoothing. Sets the `algorithm:epd:masstrace_snr_filtering` option.

- localRTRange, localMZRange:

  Retention/MZ range where to look for coeluting/isotopic mass traces.
  Sets the `algorithm:ffm:local_rt_range` and
  `algorithm:ffm:local_mz_range` options, respectively.

- isotopeFilteringModel:

  Remove/score candidate assemblies based on isotope intensities. See
  [FeatureFinderMetabo](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMetabo.html).
  Sets the `algorithm:ffm:isotope_filtering_model` option.

- MZScoring13C:

  Use the 13C isotope as the expected shift for isotope mass traces. See
  [FeatureFinderMetabo](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMetabo.html).
  Sets `algorithm:ffm:mz_scoring_13C`.

- useSmoothedInts:

  If `TRUE` then use LOWESS intensities instead of raw intensities. Sets
  the `algorithm:ffm:use_smoothed_intensities` option.

- extraOpts:

  Named `list` containing extra options that will be passed to
  `FeatureFinderMetabo`. Any options specified here will override any of
  the above. Example:
  `extraOpts=list("-algorithm:common:noise_threshold_int"=1000)`
  (corresponds to setting `noiseThrInt=1000`). Set to `NULL` to ignore.

- useFFMIntensities:

  If `TRUE` then peak intensities are directly loaded from
  `FeatureFinderMetabo` output. Otherwise, intensities are loaded
  afterwards from the input `mzML` files, which is potentially much
  slower, especially with many analyses files. However,
  `useFFMIntensities=TRUE` is still somewhat experimental, may be less
  accurate and requires a recent version of `OpenMS` (\>=2.7).

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses OpenMS to automatically find features. This function
is called when calling `findFeatures` with `algorithm="openms"`.

This functionality has been tested with OpenMS version \>= 2.0. Please
make sure it is installed and configured, e.g. by installing
`patRoonExt` or configuring the path of the binaries with the
`patRoon.path.OpenMS` option or the system PATH variable.

The file format of analyses must be `mzML`.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## Parallelization

`findFeaturesOpenMS` with `useFFMIntensities=FALSE` uses multiprocessing
to parallelize computations. Please see the parallelization section in
the handbook for more details and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

Note that for caching purposes, the analyses files must always exist on
the local host computer, even if it is not participating in
computations.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `findFeaturesOpenMS` with `useFFMIntensities=FALSE`
to process HRMS (or IMS-HRMS) data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## References

Rost HL, Sachsenberg T, Aiche S, Bielow C, Weisser H, Aicheler F,
Andreotti S, Ehrlich H, Gutenbrunner P, Kenar E, Liang X, Nahnsen S,
Nilse L, Pfeuffer J, Rosenberger G, Rurik M, Schmitt U, Veit J, Walzer
M, Wojnar D, Wolski WE, Schilling O, Choudhary JS, Malmstrom L,
Aebersold R, Reinert K, Kohlbacher O (2016). “OpenMS: a flexible
open-source software platform for mass spectrometry data analysis.”
*Nature Methods*, **13**(9), 741–748.
[doi:10.1038/nmeth.3959](https://doi.org/10.1038/nmeth.3959) .\
\
[pugixml](https://pugixml.org/) (via [Rcpp](http://www.rcpp.org/)) is
used to process OpenMS XML output.\
\
Eddelbuettel D (2013). *Seamless R and C++ Integration with Rcpp*.
Springer, New York.
[doi:10.1007/978-1-4614-6868-4](https://doi.org/10.1007/978-1-4614-6868-4)
. ISBN 978-1-4614-6867-7.\
\
Eddelbuettel D, Balamuta J (2018). “Extending R with C++: A Brief
Introduction to Rcpp.” *The American Statistician*, **72**(1), 28-36.
[doi:10.1080/00031305.2017.1375990](https://doi.org/10.1080/00031305.2017.1375990)
.\
\
Eddelbuettel D, François R (2011). “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software*, **40**(8), 1–18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08) .\
\
Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar
I, Bates D, Chambers J (2026). *Rcpp: Seamless R and C++ Integration*. R
package version 1.1.2, <https://www.rcpp.org>.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
