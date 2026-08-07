# Find features using XCMS (old interface)

Uses the legacy [`xcmsSet`](https://rdrr.io/pkg/xcms/man/xcmsSet.html)
function from the xcms package to find features.

## Usage

``` r
findFeaturesXCMS(analysisInfo, method = "centWave", ..., verbose = TRUE)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- method:

  The method setting used by XCMS peak finding, see
  [`xcms::findPeaks`](https://rdrr.io/pkg/xcms/man/findPeaks-methods.html)

- ...:

  Further parameters passed to
  [`xcmsSet`](https://rdrr.io/pkg/xcms/man/xcmsSet.html).

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses XCMS to automatically find features. This function is
called when calling `findFeatures` with `algorithm="xcms"`.

This function uses the legacy interface of xcms. It is recommended to
use
[`findFeaturesXCMS3`](https://rickhelmus.github.io/patRoon/reference/findFeaturesXCMS3.md)
instead.

The file format of analyses must be `mzML` or `mzXML`.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

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

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.

[`findFeaturesXCMS3`](https://rickhelmus.github.io/patRoon/reference/findFeaturesXCMS3.md)
