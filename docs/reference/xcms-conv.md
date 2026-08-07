# Conversion to XCMS objects

Converts a
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
or
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object to an `xcmsSet` or `XCMSnExp` object.

## Usage

``` r
getXCMSSet(obj, verbose = TRUE, ...)

getXCMSnExp(obj, verbose = TRUE, ...)

# S4 method for class 'features'
getXCMSSet(obj, verbose, loadRawData, IMS = FALSE)

# S4 method for class 'featuresXCMS'
getXCMSSet(obj, verbose = TRUE, ...)

# S4 method for class 'featureGroups'
getXCMSSet(obj, verbose, loadRawData, IMS = FALSE)

# S4 method for class 'featureGroupsXCMS'
getXCMSSet(obj, verbose, loadRawData, ...)

# S4 method for class 'featuresSet'
getXCMSSet(obj, ..., set)

# S4 method for class 'featureGroupsSet'
getXCMSSet(obj, ..., set)

# S4 method for class 'features'
getXCMSnExp(obj, verbose, loadRawData, IMS = FALSE)

# S4 method for class 'featuresXCMS3'
getXCMSnExp(obj, verbose = TRUE, ...)

# S4 method for class 'featureGroups'
getXCMSnExp(obj, verbose, loadRawData, IMS = FALSE)

# S4 method for class 'featureGroupsXCMS3'
getXCMSnExp(obj, verbose, loadRawData, ...)

# S4 method for class 'featuresSet'
getXCMSnExp(obj, ..., set)

# S4 method for class 'featureGroupsSet'
getXCMSnExp(obj, ..., set)
```

## Arguments

- obj:

  The object that should be converted.

- verbose:

  If `FALSE` then no text output is shown.

- ...:

  (**sets workflow**) Further arguments passed to non-sets method.

  Otherwise ignored.

- loadRawData:

  Set to `TRUE` if analyses are available as `mzXML` or `mzML` files.
  Otherwise MS data is not loaded, and some dummy data (*e.g.* file
  paths) is used in the returned object.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  export in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  This should be kept `FALSE` as XCMS export currently does not support
  IMS features.

- set:

  (**sets workflow**) The name of the set to be exported.

## Details

The conversion process will introduce some dummy values for metadata not
present in patRoon objects. If the `features` or `featureGroups` object
was generated with XCMS, then no conversion is performed and the
original XCMS object will be returned, if possible. Conversion may still
occur *e.g.* due to the application of some subsetting or filtering
steps or the re-ordering of analyses.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md),
[`unset`](https://rickhelmus.github.io/patRoon/reference/generics.md) is
used to convert the feature (group) data before the object is exported.

## References

reference Benton HP, Want EJ, Ebbels TMD (2010). “Correction of mass
calibration gaps in liquid chromatography-mass spectrometry metabolomics
data.” *BIOINFORMATICS*, **26**, 2488.\
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
