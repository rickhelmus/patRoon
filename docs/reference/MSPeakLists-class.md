# Class containing MS Peak Lists

Contains all MS (and MS/MS where available) peak lists for a
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object.

## Usage

``` r
peakLists(obj, ...)

averagedPeakLists(obj, ...)

spectrumSimilarity(obj, ...)

spectrumSimilarityIMS(obj, ...)

# S4 method for class 'MSPeakLists'
peakLists(obj)

# S4 method for class 'MSPeakLists'
averagedPeakLists(obj)

# S4 method for class 'MSPeakLists'
analyses(obj)

# S4 method for class 'MSPeakLists'
groupNames(obj)

# S4 method for class 'MSPeakLists'
length(x)

# S4 method for class 'MSPeakLists'
show(object)

# S4 method for class 'MSPeakLists,ANY,ANY,missing'
x[i, j, ..., reAverage = FALSE, drop = TRUE]

# S4 method for class 'MSPeakLists,ANY,ANY'
x[[i, j]]

# S4 method for class 'MSPeakLists'
x$name

# S4 method for class 'MSPeakLists'
as.data.table(x, fGroups = NULL, averaged = TRUE)

# S4 method for class 'MSPeakLists'
delete(obj, i = NULL, j = NULL, k = NULL, reAverage = FALSE, ...)

# S4 method for class 'MSPeakLists'
filter(
  obj,
  MSLevel = 1:2,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  topMostPeaks = NULL,
  minPeaks = NULL,
  maxMZOverPrec = NULL,
  absMinAbundanceFeat = NULL,
  relMinAbundanceFeat = NULL,
  absMinAbundanceFGroup = NULL,
  relMinAbundanceFGroup = NULL,
  maxRelCumIntensity = NULL,
  isolatePrec = NULL,
  deIsotope = FALSE,
  removeMZs = NULL,
  withMSMS = FALSE,
  annotatedBy = NULL,
  retainPrecursor = TRUE,
  mzWindow = defaultLim("mz", "medium"),
  reAverage = FALSE,
  negate = FALSE
)

# S4 method for class 'MSPeakLists'
plotSpectrum(
  obj,
  groupName,
  analysis = NULL,
  MSLevel = 1,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  xlim = NULL,
  ylim = NULL,
  showLegend = TRUE,
  ...
)

# S4 method for class 'MSPeakLists'
spectrumSimilarity(
  obj,
  groupName1,
  groupName2 = NULL,
  analysis1 = NULL,
  analysis2 = NULL,
  MSLevel = 1,
  specSimParams = getDefSpecSimParams(),
  NAToZero = FALSE,
  drop = TRUE
)

# S4 method for class 'MSPeakLists'
spectrumSimilarityIMS(obj, fGroups, doFGroups = TRUE, warn = TRUE, ...)

# S4 method for class 'MSPeakListsSet'
analysisInfo(obj, df = FALSE)

# S4 method for class 'MSPeakListsSet'
show(object)

# S4 method for class 'MSPeakListsSet,ANY,ANY,missing'
x[i, j, ..., reAverage = FALSE, sets = NULL, drop = TRUE]

# S4 method for class 'MSPeakListsSet'
as.data.table(x, fGroups = NULL, averaged = TRUE)

# S4 method for class 'MSPeakListsSet'
delete(obj, i = NULL, j = NULL, k = NULL, reAverage = FALSE, ...)

# S4 method for class 'MSPeakListsSet'
filter(
  obj,
  ...,
  removeMZs = NULL,
  withMSMS = FALSE,
  annotatedBy = NULL,
  retainPrecursor = TRUE,
  mzWindow = defaultLim("mz", "medium"),
  reAverage = FALSE,
  negate = FALSE,
  sets = NULL
)

# S4 method for class 'MSPeakListsSet'
plotSpectrum(
  obj,
  groupName,
  analysis = NULL,
  MSLevel = 1,
  title = NULL,
  normalized = "multiple",
  specSimParams = getDefSpecSimParams(),
  xlim = NULL,
  ylim = NULL,
  perSet = TRUE,
  mirror = TRUE,
  ...
)

# S4 method for class 'MSPeakListsSet'
spectrumSimilarity(
  obj,
  groupName1,
  groupName2 = NULL,
  analysis1 = NULL,
  analysis2 = NULL,
  MSLevel = 1,
  specSimParams = getDefSpecSimParams(),
  NAToZero = FALSE,
  drop = TRUE
)

# S4 method for class 'MSPeakListsSet'
unset(obj, set)

getDefIsolatePrecParams(...)
```

## Arguments

- obj, x, object:

  The `MSPeakLists` object to access.

- ...:

  For the `"["` operator: ignored.

  For `delete`: passed to the function specified as `j`.

  For `plotSpectrum`: passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For `spectrumSimilarityIMS`: passed to `spectrumSimilarity`

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `MSPeakLists` method.

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  analyses/feature groups by their index or name, respectively (for the
  order/names see `analyses()/groupNames()`).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all analyses/feature groups are
  selected.\
  \
  For `[[`: should be a scalar value. If `j` is not specified, `i`
  selects by feature groups instead.\
  \
  For `delete`: The data to remove from. `i` are the feature groups as
  numeric index, logical or character, `j` the MS peaks as numeric
  indices (rows). If either is `NULL` then data for all is removed. `j`
  may also be a function: it will be called for each feature group, with
  the peak list table (a `data.table`), feature group name, analysis
  (`NULL` if averaged function specifies the peak list indices (rows) to
  be removed (specified as an `integer` or `logical`

- reAverage:

  Set to `TRUE` to regenerate group averaged MS peak lists. **NOTE** it
  is very important that any annotation data relying on MS peak lists
  (formulae/compounds) are regenerated afterwards! Otherwise it is
  likely that *e.g.* plotting methods will use wrong MS/MS data.

- drop:

  If set to `TRUE` and if the comparison is made between two spectra
  then [`drop`](https://rdrr.io/r/base/drop.html) is used to reduce the
  `matrix` return value to a `numeric` vector.

- name:

  The feature group name (partially matched).

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that was used to generate this object. If not `NULL` it is used
  to add feature group information (retention and *m/z* values).

- averaged:

  If `TRUE` then feature group averaged peak list data is used.

- k:

  A vector with analyses (`character` with names or `integer` with
  indices) for which the data should be deleted. If `k!=NULL` then
  deletions will *not* occur on group averaged peak lists. Otherwise, if
  `k=NULL` then deletion occurs on *both* group averaged and analysis
  specific peak lists.

- MSLevel:

  The MS level for which data is plotted or filtered: `1` for regular
  MS, `2` for MSMS.

  For `filter`: can also be `1:2` to specify both.

- absMinIntensity, relMinIntensity:

  Absolute/relative intensity threshold for peaks. Set to `NULL` for
  none.

- topMostPeaks:

  Only consider this number of most intense peaks. Set to `NULL` to
  consider all.

- minPeaks:

  If the number of peaks in an MS/MS peak list (**excluding** the
  precursor peak) is lower than this it will be completely removed. Set
  to `NULL` to ignore.

- maxMZOverPrec:

  Any mass peaks with an m/z higher than this value (relative to the
  precursor) will be removed. Set to `NULL` to ignore.

- absMinAbundanceFeat, relMinAbundanceFeat:

  The minimum absolute/relative abundance for a mass peak across spectra
  that are averaged for a feature. Setting `reAverage` determines if
  feature group peak lists are also filtered:

  - `reAverage=FALSE` then this filter is also applied to feature group
    data, using the the mean averaged peak abundance from the peak lists
    in the group.

  - `reAverage=TRUE` then this filter is *not* applied to the
    (regenerated) feature group data.

  In most cases `reAverage=TRUE` makes more sense to avoid inconsistent
  filtering approaches between feature and feature group data.

  Set to `NULL` to ignore.

- absMinAbundanceFGroup, relMinAbundanceFGroup:

  The minimum absolute/relative abundance of a mass peak across spectra
  that are averaged for a feature group. Set to `NULL` to ignore.

- maxRelCumIntensity:

  The maximum relative cumulative intensity of a peak, calculated in
  descending order (most intense peaks first). For instance, a value of
  `0.95` means that only the most intense peaks that together account
  for \<=95% of the total intensity are retained. Set to `NULL` to
  ignore.

- isolatePrec:

  If not `NULL` then value should be a `list` with parameters used for
  isolating the precursor and its isotopes in MS peak lists (see
  `Isolating precursor data`). Alternatively, `TRUE` to apply the filter
  with default settings (as given with `getDefIsolatePrecParams`).

- deIsotope:

  Remove any isotopic peaks in peak lists. This may improve data
  processing steps which do not assume the presence of isotopic peaks
  (e.g. MetFrag for MS/MS). Note that `generateMSPeakLists` does not
  (yet) support flagging of isotopes.

- removeMZs:

  A set of m/z values to be removed from the peak lists. This is
  typically used to remove background peaks. The m/z values should be
  specified by either be a `numeric` vector or a `data.frame` with an
  `mz` column. The latter is returned by the
  [`getBGMSMSPeaks`](https://rickhelmus.github.io/patRoon/reference/getBGMSMSPeaks.md)
  function, which attempts to automatically detect background peaks.

  (**sets workflow**) Should be a `list` that specifies the m/z values
  to be removed (in above mentioned format) for each set. The order
  should match that of the sets in the `MSPeakLists` object.

- withMSMS:

  If set to `TRUE` then only results will be retained for which MS/MS
  data is available. if `negate=TRUE` then only results *without* MS/MS
  data will be retained.

- annotatedBy:

  Either a
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  or
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object, or a `list` with both. Any MS/MS peaks that are *not*
  annotated by any of the candidates in the specified objects are
  removed.

  **NOTE**: the `annotatedBy` filter currently only supports filtering
  peak of feature groups (and not of features). Hence, this filter
  cannot be combined with `reAverage=TRUE`. Furthermore, if peak lists
  are re-averaged after application of this filter, any filtered results
  will be *undone*.

- retainPrecursor:

  If `TRUE` then precursor peaks will never be filtered out from MS/MS
  peak lists (note that precursors are never removed from MS peak
  lists). The `negate` argument does not affect this setting.

- mzWindow:

  The m/z window used to find peaks to be removed from the `removeMZs`
  filter.

- negate:

  If `TRUE` then filters are applied in opposite manner.

- groupName:

  The name of the feature group for which a plot should be made. To
  compare spectra, two group names can be specified.

- analysis:

  The name of the analysis for which a plot should be made. If `NULL`
  then data from the feature group averaged peak list is used. When
  comparing spectra, either `NULL` or the analyses for both spectra
  should be specified.

- title:

  The title of the plot. If `NULL` a title will be automatically made.

- normalized:

  Controls intensity normalization. Should be `FALSE` (don't normalize),
  `TRUE` (normalize) or `"multiple"` (only normalizes if multiple
  spectra are plotted).

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- showLegend:

  Set to `TRUE` to show a legend.

- groupName1, groupName2:

  The names of the feature groups for which the comparison should be
  made. If both arguments are specified then a comparison is made with
  the spectra specified by `groupName1` *vs* those specified by
  `groupName2`. The length of either can be `>1` to generate a
  comparison matrix. Alternatively, if `groupName2` is `NULL` then all
  the spectra specified in `groupName1` will be compared with eachother,
  *i.e.* resulting in a square similarity matrix.

- analysis1, analysis2:

  The name of the analysis (analyses) for the comparison. If `NULL` then
  data from the feature group averaged peak list is used. Otherwise,
  should be the same length as `groupName1`/`groupName2`.

- NAToZero:

  Set to `TRUE` to convert `NA` similarities (*i.e.* when no similarity
  could be calculated) to zero values.

- doFGroups:

  Set to `TRUE` to compare spectra of feature groups, `FALSE` to compare
  spectra of features.

- warn:

  Set to `TRUE` to show a warning when no relevant feature group data is
  found.

- df:

  If `TRUE` then a `data.frame` is returned, otherwise a `data.table` is
  returned.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`).

- perSet, mirror:

  (**sets workflow**) If `perSet=TRUE` then the set specific mass peaks
  are annotated separately. Furthermore, if `mirror=TRUE` (and there are
  two sets in the object) then a mirror plot is generated.

- set:

  (**sets workflow**) The name of the set.

## Value

`peakLists` returns a nested list containing MS (and MS/MS where
available) peak lists per feature group and per analysis. The format is:
`[[analysis]][[featureGroupName]][[MSType]][[PeakLists]]` where `MSType`
is either `"MS"` or `"MSMS"` and `PeakLists` a `data.table` containing
all *m/z* values (`mz` column) and their intensities (`intensity`
column). In addition, the peak list tables may contain a `cmp` column
which contains an unique alphabetical identifier to which isotopic
cluster (or "compound") a mass belongs (only supported by MS peak lists
generated by Bruker tools at the moment).

`averagedPeakLists` returns a nested list of feature group averaged peak
lists in a similar format as `peakLists`.

`delete` returns the object for which the specified data was removed.

`spectrumSimilarityIMS` returns a `data.table` with spectral
similarities for each IMS precursor and feature pair.

## Details

Objects for this class are returned by
[`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md).

The `getDefIsolatePrecParams` is used to create a parameter list for
isolating the precursor and its isotopes (see
`Isolating precursor data`).

## Methods (by generic)

- `peakLists(MSPeakLists)`: Accessor method to obtain the MS peak lists.

- `averagedPeakLists(MSPeakLists)`: Accessor method to obtain the
  feature group averaged MS peak lists.

- `analyses(MSPeakLists)`: returns a `character` vector with the names
  of the analyses for which data is present in this object.

- `groupNames(MSPeakLists)`: returns a `character` vector with the names
  of the feature groups for which data is present in this object.

- `length(MSPeakLists)`: Obtain total number of *m/z* values.

- `show(MSPeakLists)`: Shows summary information for this object.

- `x[i`: Subset on analyses/feature groups.

- `x[[i`: Extract a list with MS and MS/MS (if available) peak lists. If
  the second argument (`j`) is not specified the averaged peak lists for
  the group specified by the first argument (`i`) will be returned.

- `$`: Extract group averaged MS peaklists for a feature group.

- `as.data.table(MSPeakLists)`: Returns all MS peak list data in a
  table.

- `delete(MSPeakLists)`: Completely deletes specified peaks from MS peak
  lists.

- `filter(MSPeakLists)`: provides post filtering of generated MS peak
  lists, which may further enhance quality of subsequent workflow steps
  (*e.g.* formulae calculation and compounds identification) and/or
  speed up these processes. The filters are applied to peak lists for
  each feature and feature group. The feature group peak lists are *not*
  re-averaged by default (see the `reAverage` argument). *not* filtered
  afterwards.

- `plotSpectrum(MSPeakLists)`: Plots a spectrum using MS or MS/MS peak
  lists for a given feature group. Two spectra can be compared when two
  feature groups are specified.

- `spectrumSimilarity(MSPeakLists)`: Calculates the spectral similarity
  between two or more spectra.

- `spectrumSimilarityIMS(MSPeakLists)`: Calculates the spectral
  similarity between spectra from IMS features (or feature groups) and
  their IMS precursors in post mobility workflows (see
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)).

## Slots

- `peakLists`:

  Contains a list of all MS (and MS/MS) peak lists. Use the `peakLists`
  method for access.

- `metadata`:

  Metadata for all spectra used to generate peak lists. Follows the
  format of the `peakLists` slot.

- `averagedPeakLists`:

  A `list` with averaged MS (and MS/MS) peak lists for each feature
  group.

- `avgPeakListArgs`:

  A `list` with arguments used to generate feature group averaged
  MS(/MS) peak lists.

- `origFGNames`:

  A `character` with the original input feature group names.

- `analysisInfo`:

  (**sets workflow**) [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).
  Use the `analysisInfo` method for access.

## Isolating precursor data

Formula calculation typically relies on evaluating the measured isotopic
pattern from the precursor to score candidates. Some algorithms
(currently only `GenForm`) penalize candidates if mass peaks are present
in MS1 spectra that do not contribute to the isotopic pattern. Since
these spectra are typically very 'noisy' due to background and
co-eluting ions, an additional filtering step may be recommended prior
to formula calculation. During this precursor isolation step all mass
peaks are removed that are (1) not the precursor and (2) not likely to
be an isotopologue of the precursor. To determine potential isotopic
peaks the following parameters are used:

- `maxIsotopes` The maximum number of isotopes to consider. For
  instance, a value of `5` means that `M+0` (*i.e.* the monoisotopic
  peak) till `M+5` is considered. All mass peaks outside this range are
  removed.

- `mzDefectRange` A two-sized `vector` specifying the minimum (can be
  negative) and maximum *m/z* defect deviation compared to the precursor
  *m/z* defect. When chlorinated, brominated or other compounds with
  strong *m/z* defect in their isotopologues are to be considered a
  higher range may be desired. On the other hand, for natural compounds
  this range may be tightened. Note that the search range is propegated
  with increasing distance from the precursor, *e.g.* the search range
  is doubled for `M+2`, tripled for `M+3` etc.

- `intRange` A two-sized `vector` specifying the minimum and maximum
  relative intensity range compared to the precursor. For instance,
  `c(0.001, 2)` removes all peaks that have an intensity below 0.1% or
  above 200% of that of the precursor.

- `z` The `z` value (*i.e.* absolute charge) to be considerd. For
  instance, a value of `2` would look for `M+0.5`, `M+1` etc. Note that
  the `mzDefectRange` is adjusted accordingly (*e.g.* halved if `z=2`).

- `maxGap` The maximum number of missing adjacent isotopic peaks
  ('gaps'). If the (rounded) *m/z* difference to the previous peak
  exceeds this value then this and all next peaks will be removed.
  Similar to `z`, the maximum gap is automatically adjusted for
  `charge`.

These parameters should be in a `list` that is passed to the
`isolatePrec` argument to `filter`. The default values can be obtained
with the `getDefIsolatePrecParams` function:

`maxIsotopes=5`; `mzDefectRange=c(-0.01, 0.01)`; `intRange=c(0.001, 2)`;
`z=1`; `maxGap=2`

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`MSPeakLists`**

    - `MSPeakListsSet`

    - `MSPeakListsUnset`

## Source

`spectrumSimilarity`: The principles of spectral binning and cosine
similarity calculations were loosely was based on the code from
`SpectrumSimilarity()` function of OrgMassSpecR.

## Sets workflows

The `MSPeakListsSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `MSPeakLists` and therefore largely follows
the same user interface.

The following methods are specifically defined for sets workflows:

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`MSPeakListsUnset`), which allows it to be used in 'regular'
  workflows. Only the MS peaks that are present in the specified set are
  kept.

- `analysisInfo` Returns the [analysis
  info](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  for this object.

The following methods are changed or with new functionality:

- `filter` and the subset operator (`[`) Can be used to select data that
  is only present for selected sets (`sets` argument).

- The `filter` method is applied for each set individually, and
  afterwards the results are combined again (see
  [`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md)).
  Note that this has important implications for *e.g.* intensity filters
  (`absMinIntensity`/`relMinIntensity`), `topMostPeaks` and `minPeaks`.
  Furthermore, when the `annotatedBy` filter is applied, each set
  specific MS peak list is filtered by the annotation results from only
  that set. Finally, the `removeMZs` filter should be set for each set
  separately.

- `plotSpectrum` Is able to highlight set specific mass peaks (`perSet`
  and `mirror` arguments).

- `spectrumSimilarity` First calculates similarities for each spectral
  pair per set (*e.g.* all positive mode spectra are compared and then
  all negative mode spectra are compared). This data is then combined
  into an overall similarity value. How this combination is performed
  depends on the `setCombineMethod` field of the
  [`specSimParams`](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  argument.

## Author

For `spectrumSimilarity`: major contributions by Bas van de Velde for
spectral binning and similarity calculation.
