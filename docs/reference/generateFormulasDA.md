# Generate formula with Bruker DataAnalysis

Uses Bruker DataAnalysis to generate chemical formulae.

## Usage

``` r
generateFormulasDA(fGroups, ...)

# S4 method for class 'featureGroups'
generateFormulasDA(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  precursorMzSearchWindow = defaultLim("mz", "narrow"),
  MSMode = "both",
  adduct = NULL,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = defaultLim("mz", "narrow"),
  save = TRUE,
  close = save
)

# S4 method for class 'featureGroupsSet'
generateFormulasDA(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  precursorMzSearchWindow = defaultLim("mz", "narrow"),
  MSMode = "both",
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
  object for which formulae should be generated. This should be the same
  or a subset of the object that was used to create the specified
  `MSPeakLists`. In the case of a subset only the remaining feature
  groups in the subset are considered.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- MSPeakLists:

  An
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- precursorMzSearchWindow:

  Search window for *m/z* values (+/- the feature *m/z*) used to find
  back feature data of precursor/parent ions from MS/MS spectra (this
  data is not readily available from `SmartFormula3D` results).

- MSMode:

  Whether formulae should be generated only from MS data (`"ms"`), MS/MS
  data (`"msms"`) or both (`"both"`). Selecting "both" will calculate
  formulae from MS data and MS/MS data and combines the results
  (duplicated formulae are removed). This is useful when poor MS/MS data
  would exclude proper candidates.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- featThreshold:

  If `calculateFeatures=TRUE`: minimum presence (`0-1`) of a formula in
  all features before it is considered as a candidate for a feature
  group. For instance, `featThreshold=0.75` dictates that a formula
  should be present in at least 75% of the features inside a feature
  group.

- featThresholdAnn:

  As `featThreshold`, but only considers features with annotations. For
  instance, `featThresholdAnn=0.75` dictates that a formula should be
  present in at least 75% of the features with annotations inside a
  feature group.

- absAlignMzDev:

  When the group formula annotation consensus is made from feature
  annotations, the *m/z* values of annotated MS/MS fragments may
  slightly deviate from those of the corresponding group MS/MS peak
  list. The `absAlignMzDev` argument specifies the maximum *m/z* window
  used to re-align the mass peaks.

- close, save:

  If `TRUE` then Bruker files are closed and saved after processing with
  DataAnalysis, respectively. Setting `close=TRUE` prevents that many
  analyses might be opened simultaneously in DataAnalysis, which
  otherwise may use excessive memory or become slow. By default `save`
  is `TRUE` when `close` is `TRUE`, which is likely what you want as
  otherwise any processed data is lost.

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

## Value

A
[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
object containing all generated formulae.

## Details

This function uses bruker to generate formula candidates. This function
is called when calling `generateFormulas` with `algorithm="bruker"`.

This method supports scoring based on overlap between measured and
theoretical isotopic patterns (both MS and MS/MS data) and the presence
of 'fitting' MS/MS fragments. The method will iterate through all
features (or "Compounds" in DataAnalysis terms) and call `SmartFormula`
(and `SmartFormula3D` if MS/MS data is available) to generate all
formulae. Parameters affecting formula calculation have to be set in
advance within the DataAnalysis method for each analysis (*e.g.* by
[`setDAMethod`](https://rickhelmus.github.io/patRoon/reference/bruker-utils.md)).

This method requires that features were obtained with
[`findFeaturesBruker`](https://rickhelmus.github.io/patRoon/reference/findFeaturesBruker.md).
It is recommended, but not mandatory, that the
[`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
are also generated by DataAnalysis.

Calculation of formulae with DataAnalysis always occurs with the
'feature approach' (see `Candidate assignment` in
[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)).

## Note

If any errors related to `DCOM` appear it might be necessary to
terminate DataAnalysis (note that DataAnalysis might still be running as
a background process). The `ProcessCleaner` application installed with
DataAnalayis can be used for this.

## See also

[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
for more details and other algorithms.
