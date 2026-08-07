# Generation of MS Peak Lists

Functionality to convert MS and MS/MS spectra into MS peak lists.

## Usage

``` r
generateMSPeakLists(fGroups, ...)

# S4 method for class 'featureGroups'
generateMSPeakLists(
  fGroups,
  maxMSRTWindow = defaultLim("retention", "narrow"),
  fixedIsolationWidth = FALSE,
  topMost = NULL,
  avgFeatParams = getDefAvgPListParams(),
  avgFGroupParams = getDefAvgPListParams()
)

# S4 method for class 'featureGroupsSet'
generateMSPeakLists(fGroups, ...)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which MS peak lists should be generated.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- maxMSRTWindow:

  Maximum chromatographic peak window used for spectrum averaging (in
  seconds, +/- retention time). If `NULL` all spectra from a feature
  will be taken into account. Lower to decrease processing time.

- fixedIsolationWidth:

  Configures how MS/MS spectra are selected for a feature:

  - `NA`: select all spectra unconditionally, *i.e.* ignoring any
    precursor information of spectra.

  - `FALSE`: select the spectra that were recorded for the feature m/z,
    using the instrumental precursor isolation window (*e.g.* quadrupole
    m/z width) as the selection tolerance.

  - A numeric value: select the spectra that were recorded for the
    feature m/z within this tolerance window (+/- precursor m/z).

  - A two-sized `numeric` vector: select the spectra that were recorded
    for the feature m/z within the specified minimum and maximum
    tolerance window (relative to the precursor m/z).

  If no isolation was applied to record MS/MS data (*e.g.*
  data-independent MS/MS), then all MS/MS spectra will be always be
  selected.

  **NOTE**: Sometimes the isolation windows are not exported and cannot
  be deduced automatically (e.g. Agilent data). In that case,
  `fixedIsolationWidth` needs to be set to a numeric or `NA` value.

- topMost:

  Only extract MS peak lists from a maximum of `topMost` features with
  highest intensity. If `NULL` all features will be used.

- avgFeatParams:

  Parameters used for averaging MS peak lists of individual features.
  Analogous to `avgFGroupParams`.

- avgFGroupParams:

  A `list` with parameters used for averaging of peak lists for feature
  groups. See
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for more details.

## Value

An
[`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
object.

## Details

Data processing steps that use mass spectral data (*e.g.* [formula
generation](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
and [compound
generation](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md))
typically use 'MS peak lists', which are tables storing the *m/z*,
intensity and other data from the raw mass spectra. The
`generateMSPeakLists` function generates MS and MS/MS peak lists first
for all features (or a subset, if the `topMost` argument is set). During
this step multiple spectra over the feature elution profile are
averaged. Subsequently, peak lists will be generated for each feature
group by averaging peak lists of the features within the group. The data
processing steps that uses peak lists will either use the data from
individual features or from group averaged peak lists. For instance, the
former may be used by formulae calculation, while compound
identification and plotting functionality typically uses group averaged
peak lists.

## Deprecated algorithms

Prior to patRoon 3.0 the `generateMSPeakLists` function was a wrapper to
the now deprecated functions
[`generateMSPeakListsMzR`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakListsMzR.md),
[`generateMSPeakListsDA`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakListsDA.md)
and
[`generateMSPeakListsDAFMF`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakListsDAFMF.md).
These functions are now unmaintained and may not (fully) work anymore.
The current interface is much faster and provides most of the previous
functionality (and additional). However, please provide feedback if you
feel any functionality is missing.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `generateMSPeakLists` to process HRMS (or IMS-HRMS)
data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

The use of profile *m/z* HRMS data (not IMS-HRMS) is currently not
supported.

## Sets workflows

With a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md),
the feature group averaged peak lists are made per set. This is
important, because for averaging peak lists cannot be mixed, for
instance, when different ionization modes were used to generate the
sets. The group averaged peaklists are then simply combined and labelled
in the final peak lists. However, annotation and most other
functionality typically uses only the (uncombined and) set specific peak
lists, as most algorithms cannot work with mixed peak lists.
