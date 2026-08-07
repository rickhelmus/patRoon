# Imports feature groups from XCMS (old interface)

Imports feature groups from XCMS (old interface)

## Usage

``` r
importFeatureGroupsXCMS(input, analysisInfo)
```

## Arguments

- input:

  An `xcmsSet` object.

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function imports data from XCMS. This function is called when
calling `importFeatureGroups` with `type="xcms"`.

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

[`importFeatureGroups`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroups.md)
for more details and other algorithms.

[`groupFeaturesXCMS`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesXCMS.md)
