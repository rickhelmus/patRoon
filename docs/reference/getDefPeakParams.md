# Peak detection parameters

Algorithms and parameters for automatic detection of peaks in
chromatograms and mobilograms.

## Usage

``` r
getDefPeakParams(type, algorithm, ...)
```

## Arguments

- type:

  The type of parameter defaults: `"chrom"` for chromatograms and
  `"bruker_ims"` and `"agilent_ims"` for mobilograms coming from Bruker
  and Agilent systems, respectively.

- algorithm:

  The peak detection algorithm: `"openms"`, `"xcms3"`, `"envipick"` or
  `"piek"`.

- ...:

  optional named arguments that override defaults.

## Details

The algorithm and its parameters for peak detection should be in a named
`list` with the format:

`list(algorithm = <algorithm>, param1 = ..., param2 = ..., ...)`

Where `<algorithm>` is the name of the algorithm and `param1`, `param2`
etc are the parameters. The `getDefPeakParams` function generates such
parameter list with the algorithm and default parameters.

The following algorithms are currently supported:

- `"openms"`: uses
  [MRMTransitionGroupPicker](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MRMTransitionGroupPicker.html)
  tool from [OpenMS](http://www.openms.de).

- `"xcms3"`: uses the
  [`xcms::peaksWithCentWave`](https://rdrr.io/pkg/xcms/man/peaksWithCentWave.html)
  function.

- `"envipick"`: uses the
  [`enviPick::mzpick`](https://rdrr.io/pkg/enviPick/man/mzpick.html)
  function.

- `"piek"`: uses the peak detection algorithm from (Dietrich et
  al. 2021) , which was optimized with [OpenMP](https://www.openmp.org/)
  parallelization. See
  [`findFeaturesPiek`](https://rickhelmus.github.io/patRoon/reference/findFeaturesPiek.md)
  for more details.

The parameters are discussed in the next sections.

## Note

The peak detection used by `algorithm="openms"` is different than that
of
[`findFeaturesOpenMS`](https://rickhelmus.github.io/patRoon/reference/findFeaturesOpenMS.md).

The `patRoon.threads` package option sets the number of threads for the
`piek` algorithm.

## General parameters

These parameters are applicable to all algorithms

- `forcePeakWidth` a two-sized `numeric` vector with the minimum and
  maximum width for a peak. Peaks that are more narrow or wide will be
  clamped to this range. This is especially useful for algorithms that
  consider an extensive part of the fronting/tailing noise as part as
  the peak. Set to `c(0, 0)` to disable.

- `relMinIntensity` the minimum intensity threshold for a peak relative
  to the highest peak in the same chromatogram/mobilogram. This is
  *e.g.* useful to exclude noise in mobilograms where normally few peaks
  are expected.

- `calcCentroid` Controls how the peak centroid is calculated, which is
  used for retention time or mobility determination. Valid values are:
  `"algorithm"` (use the centroid as determined by the algorithm),
  `"max"` (use the apex of the peak), `"weighted.mean"` (use the
  intensity weighted mean of all data points in the peak) or
  `"centerOfMass"` (use the center of mass or first statistical moment
  of the peak). The latter two might of interest for assymaterical
  peaks. However, most algorithms, including those not interfaced by
  patRoon, seem to use the peak apex. Hence, `calcCentroid="max"` (or
  `calcCentroid="algorithm"` which is usually the same) seems a good
  default for comparative reasons.

## Parameters for `openms`

The parameters directly map to the command line options for
`MRMTransitionGroupPicker`, please see [its
documentation](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MRMTransitionGroupPicker.html).

- `minPeakWidth` the minimum peak width, sets the `min_peak_width`
  option.

- `backgroundSubtraction` the background subtraction method, sets the
  `-algorithm:background_subtraction` option.

- `SGolayFrameLength` the frame length for Savitzky-Golay smoothing,
  sets the `-algorithm:PeakPickerMRM:sgolay_frame_length` option.

- `SGolayPolyOrder` order of the polynomial, sets the
  `-algorithm:PeakPickerMRM:sgolay_polynomial_order` option.

- `useGauss` set to `TRUE` to use Gaussian smoothing (instead of
  Savitzky-Golay, sets the `-algorithm:PeakPickerMRM:use_gauss` option.

- `gauss_width` the Gaussian width, estimated peak size, sets the
  `-algorithm:PeakPickerMRM:gauss_width` option.

- `SN` signal to noise threshold, sets the
  `-algorithm:PeakPickerMRM:signal_to_noise` option.

- `SNWinLen` `SN` window length, sets the
  `-algorithm:PeakPickerMRM:sn_win_len` option.

- `SNBinCount` `SN` bin count, sets the
  `-algorithm:PeakPickerMRM:sn_bin_count` option.

- `method` peak picking method, sets the
  `-algorithm:PeakPickerMRM:method` option.

- `integrationType` the integration technique, sets the
  `-algorithm:PeakIntegrator:integration_type` option.

- `baselineType` the baseline type, sets the
  `-algorithm:PeakIntegrator:baseline_type` option.

- `fitEMG` if `TRUE` then the EMG model is used for fitting, sets the
  `-algorithm:PeakIntegrator:fit_EMG` option.

## Parameters for `xcms3` and `envipick`

See the documentation for
[`xcms::peaksWithCentWave`](https://rdrr.io/pkg/xcms/man/peaksWithCentWave.html)
and [`enviPick::mzpick`](https://rdrr.io/pkg/enviPick/man/mzpick.html)
for `xcms3` and `envipick`, respectively.

## Parameters for `piek`

- `minIntensity` the minimum intensity of a peak.

- `SN` the signal to noise ratio.

- `peakWidth` two-sized `vector` with the minimum and maximum peak width
  (seconds)

- `RTRange` two-sized `vector` with the minimum and maximum retention
  time range (seconds). Set the 2nd element to `Inf` for no upper limit.

- `maxPeaksPerSignal` upper threshold for consecutive maxima of similar
  size to be regarded as noise.

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

Dietrich C, Wick A, Ternes TA (2021). “Open‐source feature detection for
non‐target LC–MS analytics.” *Rapid Communications in Mass
Spectrometry*, **36**(2). ISSN 1097-0231.
[doi:10.1002/rcm.9206](https://doi.org/10.1002/rcm.9206) .
<http://dx.doi.org/10.1002/rcm.9206>. Benton HP, Want EJ, Ebbels TMD
(2010). “Correction of mass calibration gaps in liquid
chromatography-mass spectrometry metabolomics data.” *BIOINFORMATICS*,
**26**, 2488.\
\
Louail P, Brunius C, Garcia-Aloy M, Kumler W, Storz N, Stanstrup J,
Treutler H, Vangeenderhuysen P, Witting M, Neumann S, Rainer J (2025).
“xcms in Peak Form: Now Anchoring a Complete Metabolomics Data
Preprocessing and Analysis Software Ecosystem.” *Analytical Chemistry*.
[doi:10.1021/acs.analchem.5c04338](https://doi.org/10.1021/acs.analchem.5c04338)
. <https://doi.org/10.1021/acs.analchem.5c04338>.\
\
Smith, C.A., Want, E.J., O'Maille, G., Abagyan,R., Siuzdak, G. (2006).
“XCMS: Processing mass spectrometry data for metabolite profiling using
nonlinear peak alignment, matching and identification.” *Analytical
Chemistry*, **78**, 779–787.\
\
Tautenhahn R, Boettcher C, Neumann S (2008). “Highly sensitive feature
detection for high resolution LC/MS.” *BMC Bioinformatics*, **9**, 504.
