# Filtering of grouped features

Basic rule based filtering of feature groups.

## Usage

``` r
replicateSubtract(fGroups, replicates, threshold = 0)

# S4 method for class 'featureGroups'
filter(
  obj,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  preAbsMinIntensity = NULL,
  preRelMinIntensity = NULL,
  absMinMaxIntensity = NULL,
  relMinMaxIntensity = NULL,
  absMinAnalyses = NULL,
  relMinAnalyses = NULL,
  absMinReplicates = NULL,
  relMinReplicates = NULL,
  absMinFeatures = NULL,
  relMinFeatures = NULL,
  absMinReplicateAbundance = NULL,
  relMinReplicateAbundance = NULL,
  absMinConc = NULL,
  relMinConc = NULL,
  absMaxTox = NULL,
  relMaxTox = NULL,
  absMinConcTox = NULL,
  relMinConcTox = NULL,
  maxReplicateIntRSD = NULL,
  maxReplicateIntRSDPres = NULL,
  blankThreshold = NULL,
  retentionRange = NULL,
  mzRange = NULL,
  mzDefectRange = NULL,
  chromWidthRange = NULL,
  featQualityRange = NULL,
  groupQualityRange = NULL,
  replicates = NULL,
  IMS = NULL,
  withIMSPrecursor = FALSE,
  IMSRangeParams = NULL,
  results = NULL,
  removeBlanks = FALSE,
  removeISTDs = FALSE,
  checkFeaturesSession = NULL,
  predAggrParams = getDefPredAggrParams(),
  removeNA = FALSE,
  negate = FALSE,
  applyIMS = "both"
)

# S4 method for class 'featureGroupsSet'
filter(
  obj,
  ...,
  negate = FALSE,
  applyIMS = "both",
  sets = NULL,
  absMinSets = NULL,
  relMinSets = NULL
)

# S4 method for class 'featureGroups'
replicateSubtract(fGroups, replicates, threshold = 0)
```

## Arguments

- fGroups, obj:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object to which the filter is applied.

- replicates:

  A character vector of replicates that should be kept (`filter`) or
  subtracted from (`replicateSubtract`).

- threshold:

  Minimum relative threshold (compared to mean intensity of replicate
  being subtracted) for a feature group to be *not* removed. When `0` a
  feature group is always removed when present in the given replicates.

- absMinIntensity, relMinIntensity:

  Minimum absolute/relative intensity for features to be kept. The
  relative intensity is determined from the feature with highest
  intensity (of all features from all groups). Set to `0` or `NULL` to
  skip this step.

- preAbsMinIntensity, preRelMinIntensity:

  As `absMinIntensity`/`relMinIntensity`, but applied *before* any other
  filters. This is typically used to speed-up subsequent filter steps.
  However, care must be taken that a sufficiently low value is chosen
  that is not expected to affect subsequent filtering steps. See below
  why this may be important.

- absMinMaxIntensity, relMinMaxIntensity:

  Feature groups are only kept if at least one feature in the group has
  an intensity above this absolute/relative threshold.

- absMinAnalyses, relMinAnalyses:

  Feature groups are only kept when they contain data for at least this
  (absolute or relative) amount of analyses. Set to `NULL` to ignore.

- absMinReplicates, relMinReplicates:

  Feature groups are only kept when they contain data for at least this
  (absolute or relative) amount of replicates. Set to `NULL` to ignore.

- absMinFeatures, relMinFeatures:

  Analyses are only kept when they contain at least this (absolute or
  relative) amount of features. Set to `NULL` to ignore.

- absMinReplicateAbundance, relMinReplicateAbundance:

  Minimum absolute/relative abundance that a grouped feature should be
  present within a replicate. If this minimum is not met all features
  within the replicate are removed. Set to `NULL` to skip this step.

- absMinConc, relMinConc:

  The minimum absolute/relative predicted concentration (calculated by
  [`calculateConcs`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md))
  assigned to a feature. The toxicities are first aggregated prior to
  filtering, as controlled by the `predAggrParams` argument. Also see
  the `removeNA` argument.

- absMaxTox, relMaxTox:

  The maximum absolute/relative predicted toxicity (LC50) (calculated by
  [`calculateTox`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md))
  assigned to a feature group. The concentrations are first aggregated
  prior to filtering, as controlled by the `predAggrParams` argument.
  Also see the `removeNA` argument.

- absMinConcTox, relMinConcTox:

  Like `absMinConc`/`relMinConc`, but instead considers the ratio
  between feature concentrations and the toxicity of the feature group.
  For instance, `absMinConcTox=0.1` means that the calculated
  concentration of a feature should be at least `10%` of its toxicity.

- maxReplicateIntRSD, maxReplicateIntRSDPres:

  Maximum relative standard deviation (RSD) of intensity values for
  features within a replicate. If the RSD is above this value all
  features within the replicate are removed. The `maxReplicateIntRSD`
  includes absent features in the RSD calculation (*i.e.* zero
  intensities), whereas these are ignored with `maxReplicateIntRSDPres`.
  Set to `NULL` to ignore.

- blankThreshold:

  Feature groups that are also present in blank analyses (see [analysis
  info](https://rickhelmus.github.io/patRoon/reference/analysis-information.md))
  are filtered out unless their relative intensity is above this
  threshold. For instance, a value of `5` means that only features with
  an intensity five times higher than that of the blank are kept. The
  relative intensity values between blanks and non-blanks are determined
  from the mean of all non-zero blank intensities. Set to `NULL` to skip
  this step.

- retentionRange, mzRange, mzDefectRange, chromWidthRange:

  Range of retention time (in seconds), *m/z*, mass defect (defined as
  the decimal part of *m/z* values) or chromatographic peak width (in
  seconds), respectively. Features outside this range will be removed.
  Should be a numeric vector with length of two containing the min/max
  values. The maximum can be `Inf` to specify no maximum range. Set to
  `NULL` to skip this step.

- featQualityRange:

  Used to filter features by their peak qualities/scores (see
  [`calculatePeakQualities`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
  Should be a named `list` with min/max ranges for each quality/score to
  be filtered (the
  [`getFeatureQualityNames`](https://rickhelmus.github.io/patRoon/reference/generics.md)
  function can be used to obtain valid names). Example:
  `featQualityRange=list(ModalityScore=c(0.3, Inf), SymmetryScore=c(0.5, Inf))`.
  Set to `NULL` to ignore.

- groupQualityRange:

  Like `featQualityRange`, but filters on group specific or averaged
  qualities/scores.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered to be
  kept in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  Set to `NULL` to ignore.

- withIMSPrecursor:

  (**IMS workflow**) only keep IMS feature groups with IMS precursors,
  *i.e.* remove all orphans. Unaffected by `negate=TRUE`.

- IMSRangeParams:

  (**IMS workflow**) A `list` with parameters to be used for filtering
  IMS range data. See
  [`getIMSRangeParams`](https://rickhelmus.github.io/patRoon/reference/getIMSRangeParams.md)
  for details and how to make such a parameter list.

- results:

  Only keep feature groups that have results in the object specified by
  `results`. Valid classes are
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  (*e.g.* formula/compound annotations) and
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md).
  Can also be a `list` with multiple objects: in this case a feature
  group is kept if it has a result in *any* of the objects. Set to
  `NULL` to ignore.

- removeBlanks:

  Set to `TRUE` to remove all analyses that belong to replicates that
  are specified as a blank in the
  [analysis-information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).
  This is useful to simplify the analyses in the specified
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object after blank subtraction. When both `blankThreshold` and this
  argument are set, blank subtraction is performed prior to removing any
  analyses.

- removeISTDs:

  If `TRUE` then all feature groups marked as internal standard (IS) are
  removed. This requires IS assignments done by
  [`normInts`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md),
  see its documentation for more details.

- checkFeaturesSession:

  If set then features and/or feature groups are removed that were
  selected for removal (see
  [check-GUI](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)).
  The session files are typically generated with the
  [`checkFeatures`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  and
  [`predictCheckFeaturesSession`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  functions. The value of `checkFeaturesSession` should either by a path
  to the session file or `TRUE`, in which case the default session file
  name is used. If `negate=TRUE` then all non-selected features/feature
  groups are removed instead.

- predAggrParams:

  Parameters to aggregate calculated concentrations/toxicities (obtained
  with
  [`calculateConcs`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)/[`calculateTox`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md))
  prior to filtering data. See [prediction aggregation
  parameters](https://rickhelmus.github.io/patRoon/reference/pred-aggr-params.md)
  for more information.

- removeNA:

  Set to `TRUE` to remove `NA` values. Currently only applicable to the
  concentration and toxicity filters.

- negate:

  If set to `TRUE` then filtering operations are performed in opposite
  manner.

- applyIMS:

  (**IMS workflow**) whether the filters are only applied to IMS
  precursors (`applyIMS=FALSE`), only to IMS features (`applyIMS=TRUE`)
  or to both (`applyIMS="both"`). Other feature groups will always be
  kept. The `negate` option does not affect `applyIMS`.

- ...:

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  method.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`).

- absMinSets, relMinSets:

  (**sets workflow**) Feature groups are only kept when they contain
  data for at least this (absolute or relative) amount of sets. Set to
  `NULL` to ignore.

## Value

A filtered
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object. Feature groups that are filtered away have their intensity set
to zero. In case a feature group is not present in any of the analyses
anymore it will be removed completely.

## Details

`filter` performs common rule based filtering of feature groups such as
blank subtraction, minimum intensity and minimum replicate abundance.
Removing of features occurs by zeroing their intensity values.
Furthermore, feature groups that are left completely empty (*i.e.* all
intensities are zero) will be automatically removed.

`replicateSubtract` removes feature groups present in a given set of
replicates (unless intensities are above a given threshold). The
replicates that are subtracted will be removed.

## Sets workflows

The following methods are changed or with new functionality:

- `filter` has specific arguments to filter by (feature presence in)
  sets. See the argument descriptions.

- **Important**: the `mzRange`, `mzDefectRange` and `IMSRangeParams`
  filters use neutral feature masses, whereas non-sets workflows use
  *m/z* values. Hence, adjust accordingly to avoid (slightly) different
  results!

## Filter order

When multiple arguments are specified to `filter`, multiple filters are
applied in sequence. Since some of these filters may affect each other,
choosing their order correctly may be important for effective data
filtering. For instance, when an intensity filter removes features from
blank analyses, a subsequent blank filter may not adequately perform
blank subtraction. Similarly, when intensity and blank filters are
executed after the replicate abundance filter it may be necessary to
ensure minimum replicate abundance again as the intensity and blank
filters may have removed some features within a replicate.

With this in mind, filters (if specified) occur in the following order:

1.  Features/feature groups selected for removal by the session
    specified by `checkFeaturesSession`.

2.  Pre-Intensity filters (*i.e.* `preAbsMinIntensity` and
    `preRelMinIntensity`).

3.  Chromatography and mass filters (*i.e* `retentionRange`, `mzRange`,
    `mzDefectRange`, `chromWidthRange`, `featQualityRange` and
    `groupQualityRange`).

4.  Replicate abundance filters (*i.e.* `absMinReplicateAbundance`,
    `relMinReplicateAbundance` and
    `maxReplicateIntRSD`/`maxReplicateIntRSDPres`).

5.  Blank filter (*i.e.* blankThreshold).

6.  Intensity filters (*i.e.* `absMinIntensity` and `relMinIntensity`).

7.  Replicate abundance filters (2nd time, only if previous filters
    affected results).

8.  Minimum-maximum intensity filters (*i.e.* `absMinMaxIntensity` and
    `relMinMaxIntensity`).

9.  General abundance filters (*i.e.* `absMinAnalyses`,
    `relMinAnalyses`, `absMinReplicates`, `relMinReplicates`,
    `absMinFeatures`, `relMinFeatures`), `absMinConc`, `relMinConc`,
    `absMaxTox` and `relMaxTox`.

10. Replicate filter (*i.e.* `replicates`), results filter (*i.e.*
    `results`) and blank analyses / internal standard removal (*i.e.*
    `removeBlanks=TRUE` / `removeISTDs=TRUE`).

If another filtering order is desired then `filter` should be called
multiple times with only one filter argument at a time.

## See also

[`featureGroups-class`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
and
[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
