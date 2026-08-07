# Find features using SAFD

Uses [SAFD](https://bitbucket.org/SSamanipour/safd.jl/src/master/) to
obtain features. This functionality is still experimental. Please see
the details below.

## Usage

``` r
findFeaturesSAFD(
  analysisInfo,
  prefCentroid = FALSE,
  mzRange = c(0, 400),
  maxNumbIter = 1000,
  maxTPeakW = 300,
  resolution = 30000,
  minMSW = 0.02,
  RThreshold = 0.75,
  minInt = 2000,
  sigIncThreshold = 5,
  S2N = 2,
  minPeakWS = 3,
  centroidMethod = "RFM",
  centroidDM = 0.005,
  verbose = TRUE
)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- prefCentroid:

  Set to `TRUE` to prefer centroided data over other MS data specified
  in `analysisInfo`.

  **NOTE**: if `prefCentroid=FALSE` but the package option
  patRoon.MS.preferIMS=TRUE (see
  [`msdata`](https://rickhelmus.github.io/patRoon/reference/msdata.md)),
  then centroided data will still be preferred over IMS data.

- mzRange:

  The *m/z* window to be imported.

- maxNumbIter, maxTPeakW, resolution, minMSW, RThreshold, minInt,
  sigIncThreshold, S2N, minPeakWS:

  Parameters directly passed to the `safd_s3D` function.

- centroidMethod, centroidDM:

  Passed to the `safd_s3d_cent` function (`method` and `mdm` arguments,
  respectively).

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses SAFD to automatically find features. This function is
called when calling `findFeatures` with `algorithm="safd"`.

The support for SAFD is still experimental, and its interface might
change in the future.

In order to use SAFD, please make sure that its `Julia` packages are
installed and you have verified that everything works, *e.g.* by running
the test data with `SAFD`.

As of patRoon `3.0`, `findFeaturesSAFD` uses the
[msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md)
interface instead of the MS_Import.jl `Julia` package to read HRMS data.
This means that MS_Import.jl does not need to be installed, and all file
formats supported by
[msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md) are
also supported for `SAFD` feature detection. This includes IMS-HRMS
data, however, in that case IMS resolved spectra are summed and the IMS
dimension is removed to make the data compatible for `SAFD`.

The `SAFD` algorithm was primarily developed to detect features in
profile *m/z* data, but centroided data is also supported. To use
profile data, ensure that the paths are correctly set up in the
[analysisInfo](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).
Furthermore, when using profile data you probably also need to specify
centroided data in the `analysisInfo`, as *e.g.*
[`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md)
currently does not support profile data. If IMS-HRMS data is used it is
treated as profile data, as this data is typically not or partially
centroided (`generateMSPeakLists` supports IMS-HRMS data directly).

## Parallelization

`findFeaturesSAFD` uses multiprocessing to parallelize computations.
Please see the parallelization section in the handbook for more details
and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

Note that for caching purposes, the analyses files must always exist on
the local host computer, even if it is not participating in
computations.

## References

Samanipour S, OBrien JW, Reid MJ, Thomas KV (2019). “Self Adjusting
Algorithm for the Nontargeted Feature Detection of High Resolution Mass
Spectrometry Coupled with Liquid Chromatography Profile Data.”
*Analytical Chemistry*, **91**(16), 10800–10807.
[doi:10.1021/acs.analchem.9b02422](https://doi.org/10.1021/acs.analchem.9b02422)
.

## See also

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
