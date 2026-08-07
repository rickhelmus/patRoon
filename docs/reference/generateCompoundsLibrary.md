# Compound annotation with an MS library

Uses a MS library loaded by
[`loadMSLibrary`](https://rickhelmus.github.io/patRoon/reference/loadMSLibrary.md)
for compound annotation.

## Usage

``` r
generateCompoundsLibrary(fGroups, ...)

# S4 method for class 'featureGroups'
generateCompoundsLibrary(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  MSLibrary,
  minSim = 0.75,
  minAnnSim = minSim,
  absMzDev = defaultLim("mz", "narrow"),
  adduct = NULL,
  checkIons = "adduct",
  spectrumType = "MS2",
  specSimParamsLib = specSimParams,
  minIMSSpecSim = 0
)

# S4 method for class 'featureGroupsSet'
generateCompoundsLibrary(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  MSLibrary,
  minSim = 0.75,
  minAnnSim = minSim,
  absMzDev = defaultLim("mz", "narrow"),
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0,
  setAvgSpecificScores = FALSE
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object which should be annotated. This should be the same or a subset
  of the object that was used to create the specified `MSPeakLists`. In
  the case of a subset only the remaining feature groups in the subset
  are considered.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- MSPeakLists:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- MSLibrary:

  The
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
  object that should be used to find candidates.

- minSim:

  The minimum spectral similarity for candidate records.

- minAnnSim:

  The minimum spectral similarity of a record for it to be used to find
  annotations (see the `Details` section).

- absMzDev:

  The maximum absolute *m/z* deviation between the feature group and
  library record *m/z* values for candidate selection.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- checkIons:

  A `character` that excludes library records with different adduct
  (`checkIons="adduct"`) or MS ionization polarity
  (`checkIons="polarity"`). If `checkIons="none"` then these filters are
  not applied.

- spectrumType:

  A `character` vector which limits library records to the given
  spectrum types (`Spectrum_type` field, *e.g.* `"MS2"`). Set to `NULL`
  to allow all spectrum types.

- specSimParamsLib:

  Like `specSimParams`, but these parameters are used for the
  pre-treatment of library spectra (only the `removePrecursor`,
  `relMinIntensity` and `minPeaks` parameters are used).

- minIMSSpecSim:

  (**IMS workflow**) If the spectrum similarity of an IMS feature group
  compared to its IMS precursor (see
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
  is at least this value, then the IMS feature group will not be
  subjected to the annotation algorithm and all feature annotation
  properties will be copied from its precursor. This assumes that
  feature annotation is primarily influenced by the MS/MS spectrum, and
  can be used to speed up the feature annotation process. All scorings,
  annotation similarities etc. are copied from the IMS precursor. The
  fragment annotations are also copied (`fragInfo` result column),
  however, these are adjusted based on the peak list data of the IMS
  feature group.

- setThreshold:

  (**sets workflow**) Minimum abundance for a candidate among all sets
  (`0-1`). For instance, a value of `1` means that the candidate needs
  to be present in all the set data.

- setThresholdAnn:

  (**sets workflow**) As `setThreshold`, but only taking into account
  the set data that contain annotations for the feature group of the
  candidate.

- setAvgSpecificScores:

  (**sets workflow**) If `TRUE` then set specific scorings (*e.g.* MS/MS
  match) are also averaged.

## Details

This function uses MS library spectra to generate compound candidates.
This function is called when calling `generateCompounds` with
`algorithm="library"`.

This method matches measured MS/MS data (peak lists) with those from an
MS library to find candidate structures. Hence, only feature groups with
MS/MS peak list data are annotated.

The library is searched for candidates with the following criteria:

1.  Only records with ion *m/z* (`PrecursorMZ`), SMILES, InChI, InChIKey
    and `formula` data are considered.

2.  Depending on the value of the `checkIons` argument, records with
    different adduct (`Precursor_type`) or polarity (`Ion_mode`) may be
    ignored.

3.  The *m/z* values of the candidate and feature group should match
    (tolerance set by `absMzDev` argument).

4.  The spectral similarity should not be lower than the value defined
    for the `minSim` argument.

5.  If multiple candidates with the same first-block InChIKey are found
    then only the candidate with the best spectral match is kept.

If the library contains annotations these will be added to the matched
MS/MS peaks. However, since the candidate selected from criterion \#5
above may not contain all the annotation data available from the MS
library, annotations from other records are also considered (controlled
by the `minAnnSim` argument). If this leads to different annotations for
the same mass peak then only the most abundant annotation is kept.

## Note

The `score`, `libMatch` and `annSim` output columns are all equal and
resemble the spectral similarity between the experimental and library
spectra.

## See also

[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
for more details and other algorithms.

[`loadMSLibrary`](https://rickhelmus.github.io/patRoon/reference/loadMSLibrary.md)
to obtain MS library data and the methods for
[`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
to treat the data before using it for annotation.
