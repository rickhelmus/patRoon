# Base class for grouped features.

This class holds all the information for grouped features.

## Usage

``` r
groupTable(object, ...)

groupFeatIndex(fGroups)

groupInfo(fGroups)

unique(x, incomparables = FALSE, ...)

overlap(fGroups, which = NULL, aggregate = TRUE, exclusive = FALSE, ...)

selectIons(fGroups, components, prefAdduct, ...)

groupQualities(fGroups)

groupScores(fGroups)

internalStandards(fGroups)

internalStandardAssignments(fGroups, ...)

normInts(
  fGroups,
  featNorm = "none",
  groupNorm = FALSE,
  normFunc = max,
  standards = NULL,
  ISTDRTWindow = 120,
  ISTDMZWindow = 300,
  minISTDs = 3,
  ...
)

concentrations(fGroups, ...)

toxicities(fGroups, ...)

updateGroups(fGroups, ...)

# S4 method for class 'featureGroups'
names(x)

# S4 method for class 'featureGroups'
analyses(obj)

# S4 method for class 'featureGroups'
replicates(obj)

# S4 method for class 'featureGroups'
groupNames(obj)

# S4 method for class 'featureGroups'
length(x)

# S4 method for class 'featureGroups'
hasIMS(obj)

# S4 method for class 'featureGroups'
fromIMS(obj)

# S4 method for class 'featureGroups'
show(object)

# S4 method for class 'featureGroups'
groupTable(object, areas = FALSE, normalized = FALSE)

# S4 method for class 'featureGroups'
analysisInfo(obj, df = FALSE)

# S4 method for class 'featureGroups'
analysisInfo(obj) <- value

# S4 method for class 'featureGroups'
groupInfo(fGroups)

# S4 method for class 'featureGroups'
featureTable(obj)

# S4 method for class 'featureGroups'
getFeatures(obj)

# S4 method for class 'featureGroups'
groupFeatIndex(fGroups)

# S4 method for class 'featureGroups'
groupQualities(fGroups)

# S4 method for class 'featureGroups'
groupScores(fGroups)

# S4 method for class 'featureGroups'
getFeatureQualityNames(
  obj,
  feat = TRUE,
  group = TRUE,
  scores = FALSE,
  totScore = TRUE
)

# S4 method for class 'featureGroups'
annotations(obj)

# S4 method for class 'featureGroups'
internalStandards(fGroups)

# S4 method for class 'featureGroups'
internalStandardAssignments(fGroups)

# S4 method for class 'featureGroups'
adducts(obj)

# S4 method for class 'featureGroups'
adducts(obj) <- value

# S4 method for class 'featureGroups'
concentrations(fGroups)

# S4 method for class 'featureGroups'
toxicities(fGroups)

# S4 method for class 'featureGroups,ANY,ANY,missing'
x[i, j, ..., ni, replicates, IMS, results, reorder = FALSE, drop = TRUE]

# S4 method for class 'featureGroups,ANY,ANY'
x[[i, j]]

# S4 method for class 'featureGroups'
x$name

# S4 method for class 'featureGroups'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featureGroups'
export(obj, type, out, IMS = FALSE)

# S4 method for class 'featureGroups'
unique(
  x,
  incomparables = FALSE,
  which,
  aggregate = TRUE,
  relativeTo = NULL,
  outer = FALSE
)

# S4 method for class 'featureGroups'
overlap(fGroups, which, aggregate, exclusive)

# S4 method for class 'featureGroups'
calculatePeakQualities(
  obj,
  weights,
  flatnessFactor,
  featureQualities = NULL,
  featureGroupQualities = NULL,
  avgFunc = mean,
  EICParams = getDefEICParams(window = 0),
  parallel = TRUE
)

# S4 method for class 'featureGroups'
selectIons(
  fGroups,
  components,
  prefAdduct,
  onlyMonoIso = TRUE,
  chargeMismatch = "adduct"
)

# S4 method for class 'featureGroups'
normInts(
  fGroups,
  featNorm = "none",
  groupNorm = FALSE,
  normFunc = max,
  standards = NULL,
  ISTDRTWindow = 120,
  ISTDMZWindow = 300,
  minISTDs = 3,
  ...
)

# S4 method for class 'featureGroups'
getTICs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'featureGroups'
getBPCs(obj, retentionRange = NULL, MSLevel = 1)

# S4 method for class 'featureGroups'
updateGroups(
  fGroups,
  what = c("ret", "mz", "mobility", "CCS"),
  intWeight = FALSE
)

# S4 method for class 'featureGroupsSet'
sets(obj)

# S4 method for class 'featureGroupsSet'
internalStandardAssignments(fGroups, set = NULL)

# S4 method for class 'featureGroupsSet'
adducts(obj, set, ...)

# S4 method for class 'featureGroupsSet'
adducts(obj, set, reGroup = TRUE) <- value

# S4 method for class 'featureGroupsSet'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featureGroupsSet'
show(object)

# S4 method for class 'featureGroupsSet'
featureTable(obj)

# S4 method for class 'featureGroupsSet,ANY,ANY,missing'
x[i, j, ..., sets = NULL, reorder = FALSE, drop = TRUE]

# S4 method for class 'featureGroupsSet'
export(obj, type, out, ..., set)

# S4 method for class 'featureGroupsSet'
selectIons(fGroups, components, prefAdduct, ...)

# S4 method for class 'featureGroupsSet'
normInts(
  fGroups,
  featNorm = "none",
  groupNorm = FALSE,
  normFunc = max,
  standards = NULL,
  ISTDRTWindow = 120,
  ISTDMZWindow = 300,
  minISTDs = 3,
  ...
)

# S4 method for class 'featureGroupsSet'
unset(obj, set)

# S4 method for class 'featureGroupsKPIC2'
delete(obj, ...)

# S4 method for class 'featureGroupsXCMS'
analysisInfo(obj) <- value

# S4 method for class 'featureGroupsXCMS'
delete(obj, ...)

# S4 method for class 'featureGroupsXCMS3'
delete(obj, ...)
```

## Arguments

- ...:

  For the `"["` operator: ignored.

  For `delete`: passed to the function specified as `j`.

  For `normInts`: passed to
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
  if `featNorm="istd"`.

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `featureGroups` method.

- fGroups, obj, x, object:

  `featureGroups` object to be accessed.

- incomparables:

  Not used. Included for compatibility with the generic.

- which:

  A character vector with the selection to compare (*e.g.* replicates,
  as set by the `aggregate` argument).

  For `overlap`: can also be a `NULL` to compare all elements.

- aggregate:

  Specifies how data should be aggregated prior to comparison. Set to
  `FALSE` to compare analyses, `TRUE` to compare replicates or to the
  name of a column in the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  to compare by a custom grouping of analyses.

- exclusive:

  If `TRUE` then all feature groups are removed that are not unique to
  the given replicates.

- components:

  The `components` object that was generated for the given
  `featureGroups` object. Obviously, the components must be created with
  algorithms that support adduct/isotope annotations, such as those from
  RAMClustR and cliqueMS.

- prefAdduct:

  The 'preferred adduct' (see method description). This is often
  `"[M+H]+"` or `"[M-H]-"`.

- featNorm:

  The method applied for feature normalization: `"istd"`, `"tic"`,
  `"conc"` or `"none"`. See the `Feature intensity normalization`
  section for details.

- groupNorm:

  If `TRUE` then group normalization is performed. See the
  `Feature intensity normalization` section for details.

- normFunc:

  A `function` to combine data for normalization. See the
  `Feature intensity normalization` section for details.

- standards:

  A `data.table` (or `data.frame`) with all internal standards. Should
  follow the format of a [suspect
  list](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
  Only used if `featNorm="istd"`. See the
  `Feature intensity normalization` section for details.

  (**sets workflow**) Can also be a `list` with internal standard lists.

  See the `suspects` argument to
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
  for more details.

- ISTDRTWindow, ISTDMZWindow:

  The retention time and *m/z* windows for IS selection. Only used if
  `featNorm="istd"`. See the `Feature intensity normalization` section
  for details.

- minISTDs:

  The minimum number of IS that should be assigned to each feature (if
  possible). Only used if `featNorm="istd"`. See the
  `Feature intensity normalization` section for details.

- areas:

  If set to `TRUE` then areas are considered instead of peak
  intensities.

- normalized:

  If `TRUE` then normalized intensity data is used (see the
  `Feature intensity normalization` section.

- df:

  If `TRUE` then the returned value is a `data.frame`, otherwise a
  `data.table`.

- value:

  For `analysisInfo<-`: A `data.frame` or `data.table` with updated
  analysis information.

  For `adducts<-`: A `character` with adduct annotations assigned to
  each feature group. The length should equal the number of feature
  groups. Can be named with feature group names to customize the
  assignment order.

- feat, group:

  If `TRUE` then names specific to features and feature groups are
  returned, respectively.

- scores:

  If `TRUE` the score names are returned, otherwise the quality names.

- totScore:

  If `TRUE` (and `scores=TRUE`) then the name of the total score is
  included.

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  analyses/feature groups by their index or name, respectively (for the
  order/names see `analyses()/names()`).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all analyses/feature groups are
  selected.\
  \
  For `[[`: should be a scalar value. If `j` is not specified, `i`
  selects by feature groups instead.\
  \
  For `delete`: The data to remove from. `i` are the analyses as numeric
  index, logical or character, `j` the feature groups as numeric index,
  logical or character. If either is `NULL` then data for all is
  removed. `j` may also be a function: it will be called for each
  feature group, with a vector of the group intensities, the group name
  and any other arguments passed as `...` to `delete`. The return value
  of this function specifies the analyses of the features in the group
  to be removed (same format as `i`).

- ni:

  Optional argument. An expression used for subsetting the analyses. The
  [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  is first subset and the remaining rows are used to determine for which
  analyses the results should be kept. The unevaluated `ni` expression
  is used to set the `i` argument of the subset operator of
  [data.table](https://CRAN.R-project.org/package=data.table), which
  therefore brings the advanced subsetting capabilities of
  [data.table](https://CRAN.R-project.org/package=data.table) (see the
  [`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)
  documentation for more details). For instance,
  `fGroups[replicate == "standard"]` would subset all analyses assigned
  with the replicate `"standard"`.

- replicates:

  An optional `character` vector: if specified only keep results for the
  given replicates (equivalent to the `replicates` argument to
  [`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md)).

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered to be
  kept (`"["`) or exported (`export`) in IMS workflows. The following
  options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  For `"["`: Leave unassigned to perform no IMS selection. For `export`:
  keep `FALSE` as the export formats currently do not support IMS data.

- results:

  Optional argument. If specified only feature groups with results in
  the specified object are kept. The class of `results` should be
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  or
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md).
  Multiple objects can be specified in a `list`: in this case a feature
  group is kept if it has a result in *any* of the objects (equivalent
  to the `results` argument to
  [`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md)).

- reorder:

  If `TRUE` then the order of the analyses is changed to match the order
  of the `i` argument.

  (**sets workflow**) If the `sets` argument is specified (and `i` is
  not) then the order of sets is changed instead.

- drop:

  Ignored.

- name:

  The feature group name (partially matched).

- type:

  The export type: `"brukerpa"` (Bruker ProfileAnalysis), `"brukertasq"`
  (Bruker TASQ) or `"mzmine"` (MZmine).

- out:

  The destination file for the exported data.

- relativeTo:

  A character vector of groupings that should be used for unique
  comparison. The groupings (*e.g.* replicates) are configured by the
  `aggregate` argument. If `NULL` then all groupings are used for
  comparison. Data specified in `which` are ignored.

- outer:

  If `TRUE` then only feature groups are kept which do not overlap among
  those present in the groupings specified by the `which` parameter.

- weights:

  A named `numeric` vector that defines the weight for each score to
  calculate the `totalScore`. The names of the vector follow the score
  names. Unspecified weights are defaulted to `1`. Example:
  `weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)`.

- flatnessFactor:

  Passed to MetaClean as the `flatness.factor` argument to
  [`calculateJaggedness`](https://rdrr.io/pkg/MetaClean/man/calculateJaggedness.html)
  and
  [`calculateModality`](https://rdrr.io/pkg/MetaClean/man/calculateModality.html).

- featureQualities:

  Specifies which feature qualities to calculate. Can be `NULL`
  (default, calculates all qualities), a `character` vector with names
  of qualities to calculate (e.g., `c("FWHM2Base", "Symmetry")`), or a
  `list` of custom quality definitions. See the
  [`featureQualities`](https://rickhelmus.github.io/patRoon/reference/feature-quality.md)
  function for more details.

- featureGroupQualities:

  Analogous to `featureQualities` for feature groups metrics. See the
  [`featureGroupQualities`](https://rickhelmus.github.io/patRoon/reference/feature-quality.md)
  function for more details.

- avgFunc:

  The function used to average the peak qualities and scores for each
  feature group.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

- onlyMonoIso:

  Set to `TRUE` to only keep feature groups that were annotated as
  monoisotopic. If multiple annotations are available, any monoisotopic
  annotation is sufficient for the feature group to be kept. Feature
  groups are never removed by this setting if no isotope annotations are
  available.

- chargeMismatch:

  Specifies how to deal with a mismatch in charge between adduct and
  isotope annotations. Valid values are: `"adduct"` (ignore isotope
  annotation), `"isotope"` (ignore adduct annotation), `"none"` (ignore
  both annotations) and `"ignore"` (don't check for charge mismatches).
  *Important*: when `OpenMS` is used to find features, it already
  removes any detected non-monoisotopic features by default. Hence, in
  such case setting `chargeMismatch="adduct"` is more appropriate.

- retentionRange:

  Range of retention time (in seconds) to collect TIC traces. Should be
  a numeric vector with length of two containing the min/max values. Set
  to NULL to ignore.

- MSLevel:

  Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to
  obtain TIC traces.

- what:

  A `character` vector specifying which group-wise values to update.
  Valid values are `"ret"` (retention time, in seconds), `"mz"` and
  `"mobility"`. At least one must be specified. If `"mobility"` is
  specified an no IMS information is present, it will be ignored.

- intWeight:

  If `TRUE`, calculate intensity-weighted means per group; otherwise use
  simple arithmetic means.

- set:

  (**sets workflow**) The name of the set.

- reGroup:

  (**sets workflow**) Set to `TRUE` to re-group the features after the
  adduct annotations are changed. See the `Sets workflow` section for
  more details.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep.

## Value

`delete` returns the object for which the specified data was removed.

`calculatePeakQualities` returns a modified object amended with peak
qualities and scores.

`selectIons` returns a `featureGroups` object with only the selected
feature groups and amended with adduct annotations.

`normInts` returns a `featureGroups` object, amended with data in the
`ISTDs` and `ISTDAssignments` slots if `featNorm="istd"`.

## Details

The `featureGroup` class is the workhorse of patRoon: almost all
functionality operate on its instantiated objects. The class holds all
information from grouped features (obtained from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)).
This class itself is `virtual`, hence, objects are not created directly
from it. Instead, 'feature groupers' such as
[`groupFeaturesXCMS`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesXCMS.md)
return a `featureGroups` derived object after performing the actual
grouping of features across analyses.

## Methods (by generic)

- `names(featureGroups)`: Obtain feature group names.

- `analyses(featureGroups)`: returns a `character` vector with the names
  of the analyses for which data is present in this object.

- `replicates(featureGroups)`: returns a `character` vector with the
  names of the replicates for which data is present in this object.

- `groupNames(featureGroups)`: Same as `names`. Provided for consistency
  to other classes.

- `length(featureGroups)`: Obtain number of feature groups.

- `hasIMS(featureGroups)`: Returns `TRUE` if the feature groups object
  has ion mobility information.

- `fromIMS(featureGroups)`: Returns `TRUE` if the features were directly
  generated from IMS data.

- `show(featureGroups)`: Shows summary information for this object.

- `groupTable(featureGroups)`: Accessor for `groups` slot.

- `analysisInfo(featureGroups)`: Obtain analysisInfo (see analysisInfo
  slot in
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)).

- `analysisInfo(featureGroups) <- value`: Modifies the [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  of this `features` object. This is primarily intended to change or add
  analysis metadata columns or can be used to re-order analysis. The
  removal or addition of analyses and changes to the `"analysis"` column
  are not supported. This function performs several internal updates
  after analysis information modifications. Hence, never attempt to
  change the `analysisInfo` slot directly.

- `groupInfo(featureGroups)`: Accessor for `groupInfo` slot.

- `featureTable(featureGroups)`: Obtain feature information (see
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)).

- `getFeatures(featureGroups)`: Accessor for `features` slot.

- `groupFeatIndex(featureGroups)`: Accessor for `ftindex` slot.

- `groupQualities(featureGroups)`: Accessor for `groupQualities` slot.

- `groupScores(featureGroups)`: Accessor for `groupScores` slot.

- `getFeatureQualityNames(featureGroups)`: Returns feature quality names
  that were calculated for this object.

- `annotations(featureGroups)`: Accessor for `annotations` slot.

- `internalStandards(featureGroups)`: Accessor for `ISTDs` slot.

- `internalStandardAssignments(featureGroups)`: Accessor for
  `ISTDAssignments` slot.

- `adducts(featureGroups)`: Returns a named `character` with adduct
  annotations assigned to each feature group (if available).

- `adducts(featureGroups) <- value`: Sets adduct annotations for feature
  groups.

- `concentrations(featureGroups)`: Accessor for `concentrations` slot.

- `toxicities(featureGroups)`: Accessor for `toxicities` slot.

- `x[i`: Subset on analyses/feature groups.

- `x[[i`: Extract intensity values.

- `$`: Extract intensity values for a feature group.

- `delete(featureGroups)`: Completely deletes specified feature groups.

- `export(featureGroups)`: Exports feature groups to a `.csv` file that
  is readable to Bruker ProfileAnalysis (a 'bucket table'), Bruker TASQ
  (an analyte database) or that is suitable as input for the
  `Targeted peak detection` functionality of
  [MZmine](http://mzmine.github.io/).

- `unique(featureGroups)`: Obtain a subset with unique feature groups
  present in one or more analyses, replicates etc.

- `overlap(featureGroups)`: Obtain a subset with feature groups that
  overlap between a set of specified replicate(s).

- `calculatePeakQualities(featureGroups)`: Calculates peak and group
  qualities for all features and feature groups. The peak qualities (and
  scores) are calculated with the [features method of this
  function](https://rickhelmus.github.io/patRoon/reference/features-class.md),
  and subsequently averaged per feature group. Group metrics are then
  calculated and scored and scaled by normalizing qualities among all
  groups and scaling them from `0` (worst) to `1` (best). The
  `totalScore` for each group is then calculated as the weighted sum
  from all feature (group) scores. The
  [`getMCTrainData`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  and
  [`predictCheckFeaturesSession`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  functions can be used to train and apply Pass/Fail ML models from
  MetaClean.

- `selectIons(featureGroups)`: uses
  [componentization](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
  results to select feature groups with preferred adduct ion and/or
  isotope annotation. Typically, this means that only feature groups are
  kept if they are (de-)protonated adducts and are monoisotopic. The
  adduct annotation assignments for the selected feature groups are
  copied from the components to the `annotations` slot. If the adduct
  for a feature group is unknown, its annotation is defaulted to the
  'preferred' adduct, and hence, the feature group will never be
  removed. Furthermore, if a component does not contain an annotation
  with the preferred adduct, the most intense feature group is selected
  instead. Similarly, if no isotope annotation is available, the feature
  group is assumed to be monoisotopic and thus not removed. An important
  advantage of `selectIons` is that it may considerably simplify your
  dataset. Furthermore, the adduct assignments allow formula/compound
  annotation steps later in the workflow to improve their annotation
  accuracy. On the other hand, it is important the componentization
  results are reliable. Hence, it is highly recommended that, prior to
  calling `selectIons`, the settings to
  [`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
  are optimized and its results are reviewed with
  [`checkComponents`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md).
  Finally, the `adducts<-` method can be used to manually correct adduct
  assignments afterwards if necessary.

- `normInts(featureGroups)`: Provides various methods to normalizes
  feature intensities for each sample analysis or of all features within
  a feature group. See the `Feature intensity normalization` section
  below.

- `getTICs(featureGroups)`: Obtain the total ion chromatogram/s (TICs)
  of the analyses.

- `getBPCs(featureGroups)`: Obtain the base peak chromatogram/s (BPCs)
  of the analyses.

- `updateGroups(featureGroups)`: Recalculate group information from
  feature data.

## Slots

- `groups`:

  Matrix (`data.table`) with intensities for each feature group
  (columns) per analysis (rows). Access with `groups` method.

- `features`:

  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  class associated with this object. Access with`featureTable` methods.

- `groupInfo`:

  `data.table` with retention time (`ret` column, in seconds) and *m/z*
  (`mz` column) for each feature group. Access with `groupInfo` method.

- `ftindex`:

  Matrix (`data.table`) with feature indices for each feature group
  (columns) per analysis (rows). Each index corresponds to the row
  within the feature table of the analysis (see
  [`featureTable`](https://rickhelmus.github.io/patRoon/reference/generics.md)).

- `groupQualities,groupScores`:

  A `data.table` with qualities/scores for each feature group (see the
  `calculatePeakQualities` method).

- `annotations`:

  A `data.table` with adduct annotations for each group (see the
  `selectIons` method).

- `ISTDs`:

  A `data.table` with screening results for internal standards (filled
  in by the `normInts` method).

- `ISTDAssignments`:

  A `list`, where each item is named by a feature group and consists of
  a vector with feature group names of the internal standards assigned
  to it (filled in by the `normInts` method).

- `concentrations,toxicities`:

  A `data.table` with predicted concentrations/toxicities for each
  feature group. Assigned by the
  [`calculateConcs`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)/[`calculateTox`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md)
  methods. Use the `concentratrions`/`toxicities` methods for access.

- `groupAlgo,groupArgs,groupVerbose`:

  (**sets workflow**) Grouping parameters that were used when this
  object was created. Used by `adducts<-` and `selectIons` when these
  methods perform a re-grouping of features.

- `annotations,ISTDAssignments`:

  (**sets workflow**) As the `featureGroups` slots, but contains the
  data per set.

- `annotationsChanged`:

  Set internally by `adducts()<-` and applied as soon as `reGroup=TRUE`.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `calculatePeakQualities` to process HRMS (or
IMS-HRMS) data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`featureGroups`**

    - `featureGroupsSet`

      - [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

    - `featureGroupsUnset`

    - [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

      - [`featureGroupsSetScreeningUnset`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

    - `featureGroupsBruker`

    - [`featureGroupsConsensus`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - `featureGroupsEnviMass`

    - `featureGroupsGreedy`

    - `featureGroupsIMS`

    - `featureGroupsKPIC2`

    - `featureGroupsOpenMS`

    - `featureGroupsSIRIUS`

    - `featureGroupsTable`

    - [`featureGroupsBrukerTASQ`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerTASQ.md)

    - `featureGroupsXCMS`

    - `featureGroupsXCMS3`

## Feature intensity normalization

The `normInts` method performs normalization of feature intensities (and
areas). These values are amended in the `features` slot, while the
original intensities/areas are kept. To use the normalized intensities
set `normalized=TRUE` to methods such as
[`plotInt`](https://rickhelmus.github.io/patRoon/reference/generics.md),
[`generateComponentsIntClust`](https://rickhelmus.github.io/patRoon/reference/generateComponentsIntClust.md)
and
[`as.data.table`](https://rickhelmus.github.io/patRoon/reference/feature-table.md).
Please see the `normalized` argument documentation for these methods for
more details.

The `normInts` method supports several methods to normalize
intensities/areas of features within the same analysis. Most methods are
influenced by the *normalization concentration* (`norm_conc` in the
[analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md))
set for each sample analysis. For `NA` or zero values the output will be
zero. If `norm_conc` is completely absent from the analysis information
or all values are `NA`, the normalization concentration is defaulted to
one.

The different normalization methods are:

1.  `featNorm="istd"` Uses *internal standards* (IS) for normalization.
    The IS are screened internally by the
    [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
    function. Hence, the IS specified by the `standards` argument should
    follow the format of a [suspect
    list](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
    Note that labelled elements in IS formulae should be specified with
    the [rcdk](https://CRAN.R-project.org/package=rcdk) format, *e.g.*
    `"[13]C"` for 13C, `"[2]H"` for a deuterium etc. Example IS lists
    are provided with the patRoonData and patRoonDataIMS packages.

    The assignment of IS to features is automatically performed, using
    the following criteria:

    1.  Only analyses are considered with a defined normalization
        concentration.

    2.  The IS must be detected in all of the analyses in which the
        feature was detected.

    3.  The retention time and *m/z* are reasonably close
        (`ISTDRTWindow`/`ISTDMZWindow` arguments). However, additional
        IS candidates outside these windows will be chosen if the number
        of candidates is less than the `minISTDs` argument. In this case
        the next close(st) candidate(s) will be chosen.

    Normalization of features within the same feature group always occur
    with the same IS. If multiple IS are assigned to a feature then
    normalization occurs with the combined intensity (area), which is
    calculated with the function defined by the `normFunc` argument. The
    (combined) IS intensity is then normalized by the normalization
    concentration, and finally used for feature normalization.

2.  `featNorm="tic"` Uses the Total Ion Current (TIC) to normalize
    intensities. The TIC is calculated by combining all intensities with
    the function defined by the `normFunc` argument. For this reason,
    you may need to take care to perform normalization before *e.g.*
    suspect screening or other prioritization techniques. The TIC
    normalized intensities are finally divided by the normalization
    concentration.

3.  `featNorm="conc"` Simply divides all intensities (areas) with the
    normalization concentration defined for the sample.

4.  `featNorm="none"` Performs no normalization. The raw intensity
    values are simply copied. This is mainly useful if you only want to
    do group normalization (described below).

The meaning of the normalization concentration differs for each method:
for `"istd"` it resembles the IS concentration of a sample analysis,
whereas for `"tic"` and `"conc"` it is used to normalize different
sample amounts (*e.g.* injection volume).

If `groupNorm=TRUE` then feature intensities (areas) will be normalized
by the combined values for its feature group (again, combination occurs
with `normFunc`). This *group normalization* always occurs *after*
aforementioned normalization methods. Group normalization was the only
method with patRoon `<2.1`, and still occurs automatically if `normInts`
was not called when a method is executed that requests normalized data.

In IMS workflows with post mobility assignment (see
[assignMobilities](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)),
any IMS features are excluded for the assignment of internal standards
(`featNorm="istd"`) or calculation of TICs (`featNorm="tic"`).
Furthermore, the normalized intensities and areas for IMS features are
copied from their IMS precursors.

## Sets workflows

The `featureGroupsSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `featureGroups` and therefore largely follows
the same user interface.

The following methods are specifically defined for sets workflows:

- `sets` Returns the set names for this object.

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`featureGroupsUnset`), which allows it to be used in 'regular'
  workflows. The adduct annotations for the selected set are used to
  convert all feature (group) masses to ionic *m/z* values. The
  annotations persist in the converted object.

The following methods are changed or with new functionality:

- `adducts`, `adducts<-` require the `set` argument. The order of the
  data that is returned/changed follows that of the `annotations` slot.
  Furthermore, `adducts<-` will perform a re-grouping of features when
  its `reGroup` parameter is set to `TRUE`. The implications for this
  are discussed below. Note that no adducts are changed *until*
  `reGroup=TRUE`.

- the subset operator (`[`) has specific arguments to choose (feature
  presence in) sets. See the argument descriptions.

- `export` Only allows to export data from one set. The `unset` method
  is used prior to exporting the data.

- `overlap` and `unique` allow to handle data per set. See the `sets`
  argument description.

- `selectIons` Will perform a re-grouping of features. The implications
  of this are discussed below.

- `normInts` Performs normalization for each set *independently*.

A re-grouping of features occurs if `selectIons` is called or
`adducts<-` is used with `reGroup=TRUE`. Afterwards, it is very likely
that feature group names are changed. Since data generated later in the
workflow (*e.g.* annotation steps) rely on feature group names, these
objects are **not valid** anymore, and **must** be re-generated.

## References

Chetnik K, Petrick L, Pandey G (2020). “MetaClean: a machine
learning-based classifier for reduced false positive peak detection in
untargeted LC-MS metabolomics data.” *Metabolomics*, **16**(11).
[doi:10.1007/s11306-020-01738-3](https://doi.org/10.1007/s11306-020-01738-3)
.

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for generating feature groups,
[feature-filtering](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md),
[feature-table](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
and
[feature-plotting](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md)
for more advanced `featureGroups` methods.

## Author

Rick Helmus \<<r.helmus@uva.nl>\> and Ricardo Cunha \<<cunha@iuta.de>\>
(`getTICs` and `getBPCs` functions)
