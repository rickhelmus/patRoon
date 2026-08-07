# Group features using XCMS (new interface)

Uses the new `xcms3` interface from the xcms package to find features.

## Usage

``` r
groupFeaturesXCMS3(feat, ...)

# S4 method for class 'features'
groupFeaturesXCMS3(
  feat,
  rtalign = TRUE,
  loadRawData = TRUE,
  groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$replicate),
  preGroupParam = groupParam,
  retAlignParam = xcms::ObiwarpParam(),
  verbose = TRUE
)

# S4 method for class 'featuresSet'
groupFeaturesXCMS3(
  feat,
  groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$replicate),
  verbose = TRUE
)
```

## Arguments

- feat:

  The
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object with the features to be grouped.

- ...:

  Further parameters passed to the selected grouping algorithm.

- rtalign:

  Set to `TRUE` to enable retention time alignment.

- loadRawData:

  Set to `TRUE` if analyses are available as `mzXML` or `mzML` files.
  Otherwise MS data is not loaded, and some dummy data (*e.g.* file
  paths) is used in the returned object.

- groupParam, retAlignParam:

  parameter object that is directly passed to
  [`xcms::groupChromPeaks`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html)
  and
  [`xcms::adjustRtime`](https://rdrr.io/pkg/xcms/man/adjustRtime.html),
  respectively.

- preGroupParam:

  grouping parameters applied when features are grouped *prior* to
  alignment (only with peak groups alignment).

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function uses XCMS3 to group features. This function is called when
calling `groupFeatures` with `algorithm="xcms3"`.

Grouping of features and alignment of their retention times are
performed with the
[`xcms::groupChromPeaks`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html)
and [`xcms::adjustRtime`](https://rdrr.io/pkg/xcms/man/adjustRtime.html)
functions, respectively. Both of these functions support an extensive
amount of parameters that modify their behavior and may therefore
require optimization.

## Sets workflows

`loadRawData` and arguments related to retention time alignment are
currently not supported for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).

## References

Benton HP, Want EJ, Ebbels TMD (2010). “Correction of mass calibration
gaps in liquid chromatography-mass spectrometry metabolomics data.”
*BIOINFORMATICS*, **26**, 2488.\
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

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for more details and other algorithms.
