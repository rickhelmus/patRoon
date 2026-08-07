# Generate components of transformation products

Generates components by linking feature groups of transformation
products and their parents.

## Usage

``` r
generateComponentsTPs(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsTPs(
  fGroups,
  fGroupsTPs = fGroups,
  TPs = NULL,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  ignoreParents = FALSE,
  minRTDiff = 20,
  specSimParams = getDefSpecSimParams(),
  IMS = "maybe"
)

# S4 method for class 'featureGroupsSet'
generateComponentsTPs(
  fGroups,
  fGroupsTPs = fGroups,
  TPs = NULL,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  ignoreParents = FALSE,
  minRTDiff = 20,
  specSimParams = getDefSpecSimParams(),
  IMS = "maybe"
)
```

## Arguments

- fGroups:

  The input
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  for componentization. See `fGroupsTPs`.

- ...:

  Further arguments specified to the methods.

- fGroupsTPs:

  A
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object containing the feature groups that are expected to be
  transformation products. If a distinction between parents and TPs is
  not yet known, `fGroupsTPs` should equal the `fGroups` argument.
  Otherwise, `fGroups` should only contain the parent feature groups,
  and both `fGroups` and `fGroupsTPs` *must* be a subset of the same
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object.

- TPs:

  A
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
  object. Set to `NULL` to perform linking without this data.

- MSPeakLists, formulas, compounds:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)/[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)/[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object to calculate MS/MS or annotation similarities between parents
  and TPs. If `NULL` then this data is not calculated. For more details
  see the `Linking parents and transformation products` section below.

- ignoreParents:

  If `TRUE` then feature groups present in both `fGroups` and
  `fGroupsTPs` are not considered as TPs.

- minRTDiff:

  Minimum retention time (in seconds) difference between the parent and
  a TP to calculate the [retention order
  direction](https://rickhelmus.github.io/patRoon/reference/retDir.md).

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  componentization in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

## Value

The components are stored in objects derived from
[`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md).

## Details

This function uses transformation product screening to generate
components. This function is called when calling `generateComponents`
with `algorithm="tp"`.

This method typically employs data from [generated transformation
products](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
to find parents and their TPs. However, this data is not necessary, and
components can also be made based on MS/MS similarity and/or other
annotation similarities between the parent and its TPs. For more details
see the `Linking parents and transformation products` section below.

## Note

The `shift` parameter of `specSimParams` is ignored by
`generateComponentsTPs`, since it always calculates similarities with
all supported options.

## Linking parents and transformation products

Each component consists of feature groups that are considered to be
transformation products for one parent (the parent that 'belongs' to the
component can be retrieved with the
[`componentInfo`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
method). The parent feature groups are taken from the `fGroups`
parameter, while the feature groups for TPs are taken from `fGroupsTPs`.
If a feature group occurs in both variables, it may therefore be
considered as both a parent or TP.

If transformation product data is given, *i.e.* the `TPs` argument is
set, then a suspect screening of the parents and/or TPs may need to be
performed in advance to facilitate linkage. This depends on the
algorithm that was used to generate the TPs:

- the parents need to be screened for all algorithms except `logic`
  ([`generateTPsLogic`](https://rickhelmus.github.io/patRoon/reference/generateTPsLogic.md))

- the TPs need to be screened for all algorithms except `ann_form` and
  `ann_comp`
  ([`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md)
  and
  [`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md)).

See
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
to perform the screening and
[`convertToSuspects`](https://rickhelmus.github.io/patRoon/reference/generics.md)
to create the suspect list. To include parents make sure to set
`includeParents=TRUE` when calling `convertToSuspects` or first screen
for the parents and then amend the screening object with TP screening
results by setting `amend=TRUE` to `screenSuspects`. If the the suspect
screening yields multiple TP hits, all will be reported. Similarly, if
the suspect screening contains multiple hits for a parent, a component
is made for each of the parent hits.

In case no transformation product data is provided (`TPs=NULL`), the
componentization algorithm simply assumes that each feature group from
`fGroupsTPs` is a potential TP for every parent feature group in
`fGroups`. For this reason, it is highly recommended to specify which
feature groups are parents/TPs (see the `fGroupsTPs` argument
description above) and *crucial* that the data is post-processed, for
instance by only retaining TPs that have high annotation similarity with
their parents (see the
[`filter`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)
method for
[`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)).

A typical way to distinguish which feature groups are parents or TPs
from two different (groups of) samples is by calculating Fold Changes
(see the
[`as.data.table`](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
method for feature groups and
[`plotVolcano`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
Of course, other statistical techniques from R are also suitable.

## Result columns

During componentization, several characteristics are calculated which
may be useful for post-processing. These can be obtained with *e.g.*
[`as.data.table`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)
or
[`componentTable`](https://rickhelmus.github.io/patRoon/reference/components-class.md).
The properties are either reported for each feature group in a
component, or for each candidate of a feature group in a component (only
if `TPs` was set).

The following properties may be reported for each feature group:

- `specSimilarity`: the MS/MS spectral similarity between the feature
  groups of the TP and its parent (`0-1`).

- `specSimilarityPrec`,`specSimilarityBoth`: as `specSimilarity`, but
  calculated with binned data using the `"precursor"` and `"both"`
  method, respectively (see [MS spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  for more details).

- `totalFragmentMatches` The total number of MS/MS fragment annotations
  that overlap between *all* feature annotation candidates for the TP
  feature group and the feature annotations specifically for the parent
  (based on the assigned fragment formula). If both the `formulas` and
  `compounds` arguments are specified then the annotation data is pooled
  prior to calculation. Each unique match is only counted once.

- `totalNeutralLossMatches` As `totalFragmentMatches`, but counting
  overlapping neutral loss formulae.

- `retDir`,`TP_retDir` The [retention order
  direction](https://rickhelmus.github.io/patRoon/reference/retDir.md)
  derived from the feature groups (`retDir`) or the (expected) value
  from TP data (`TP_retDir`).

- `retDiff`,`mzDiff`, The retention time and *m/z* difference between
  the parent and TP.

The candidate specific properties are stored inside the `candidates`
column in component tables, and can be obtained with `as.data.table` by
setting `candidates=TRUE`. The following properties may be present:

- `fragmentMatches`,`neutralLossMatches` As `totalFragmentMatches` and
  `totalNeutralLossMatches`, but only considering the feature
  annotations specifically for this candidate.

- `formulaDiff` The formula difference between the parent and TP (if
  formula data is available).

- `TPScore`,`annSim`,`fitFormula`,`fitCompound`,`simSusps`: TP scoring
  properties, see
  [`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md)
  and
  [`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md).

## IMS workflows

In IMS workflows with post mobility assignment (see
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
it may be necessary to call
[`expandForIMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
when componentization was performed *prior* to mobility assignments, see
its documentation for more details.

If mobilities were already assigned prior to componentization, then the
`IMS` argument selects which feature groups are subjected to
componentization. Data for IMS feature groups that were not considered
(*i.e.* when `IMS` is `FALSE` or `"maybe"`), will be expanded similarly
as is done by
[`expandForIMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md).

**NOTE**: IMS expansion by
[`expandForIMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
only expands results for TP candidates, *i.e.* no new components from
parents assigned to IMS feature groups will be added.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the component tables are amended with extra information such as
overall/specific set spectrum similarities. As sets data is mixed,
transformation products are able to be linked with a parent, even if
they were not measured in the same set.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
