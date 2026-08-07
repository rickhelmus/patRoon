# Class for suspect screened feature groups.

This class derives from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
and adds suspect screening information.

## Usage

``` r
screenInfo(obj)

# S4 method for class 'featureGroupsScreening'
screenInfo(obj)

# S4 method for class 'featureGroupsScreening'
show(object)

# S4 method for class 'featureGroupsScreening,ANY,ANY,missing'
x[i, j, ..., suspects = NULL, reorder = FALSE, drop = TRUE]

# S4 method for class 'featureGroupsScreening'
delete(obj, i = NULL, j = NULL, k = NULL, ...)

# S4 method for class 'featureGroupsScreening'
filter(
  obj,
  ...,
  onlyHits = NULL,
  IMSMatchParams = NULL,
  selectHitsBy = NULL,
  selectBestFGroups = FALSE,
  maxLevel = NULL,
  maxFormRank = NULL,
  maxCompRank = NULL,
  minAnnSimForm = NULL,
  minAnnSimComp = NULL,
  minAnnSimBoth = NULL,
  absMinFragMatches = NULL,
  relMinFragMatches = NULL,
  minRF = NULL,
  maxLC50 = NULL,
  negate = FALSE,
  applyIMS = "both"
)

# S4 method for class 'featureGroupsScreeningSet'
screenInfo(obj)

# S4 method for class 'featureGroupsScreeningSet'
show(object)

# S4 method for class 'featureGroupsScreeningSet,ANY,ANY,missing'
x[i, j, ..., suspects = NULL, sets = NULL, reorder = FALSE, drop = TRUE]

# S4 method for class 'featureGroupsScreeningSet'
delete(obj, i = NULL, j = NULL, k = NULL, ...)

# S4 method for class 'featureGroupsScreeningSet'
filter(
  obj,
  ...,
  onlyHits = NULL,
  IMSMatchParams = NULL,
  selectHitsBy = NULL,
  selectBestFGroups = FALSE,
  maxLevel = NULL,
  maxFormRank = NULL,
  maxCompRank = NULL,
  minAnnSimForm = NULL,
  minAnnSimComp = NULL,
  minAnnSimBoth = NULL,
  absMinFragMatches = NULL,
  relMinFragMatches = NULL,
  minRF = NULL,
  maxLC50 = NULL,
  negate = FALSE,
  applyIMS = "both"
)

# S4 method for class 'featureGroupsScreeningSet'
unset(obj, set)
```

## Arguments

- obj, object, x:

  The `featureGroupsScreening` object.

- i, j, reorder:

  See
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

- ...:

  Further arguments passed to the base method.

- suspects:

  An optional `character` vector with suspect names. If specified, only
  `featureGroups` will be kept that are assigned to these suspects.

- drop:

  Ignored.

- k:

  The `k` argument is used to delete screening results (instead of
  features) and should be:

  - a `character` vector with suspect names that should be removed

  - a `function` that is called with the screening info table and should
    return a `logical` vector for each suspect row to be removed

  - `NA` to remove all screening results, which is especially useful
    when paired with the `j` argument, *i.e.* to remove all screening
    results for a particular set of feature groups.

  - `NULL` to not touch screening results and only perform deletion as
    the
    [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
    method.

  Setting both `i` and `k` is currently not supported.

- onlyHits:

  If

  - `negate=FALSE` and `onlyHits=TRUE` then all feature groups without
    suspect hits will be removed. Otherwise nothing will be done.

  - `negate=TRUE` then `onlyHits=TRUE` will select feature groups
    without suspect hits, `onlyHits=FALSE` will only retain feature
    groups with suspect matches and this filter is ignored if
    `onlyHits=NULL`.

- IMSMatchParams:

  (**IMS workflow**) A `list` with parameters to be used for matching
  IMS data. See
  [`getIMSMatchParams`](https://rickhelmus.github.io/patRoon/reference/getIMSMatchParams.md)
  for details and how to make such a parameter list.

- selectHitsBy:

  Should be `"intensity"` or `"level"`. For cases where the same suspect
  is matched to multiple feature groups, only the suspect to the feature
  group with highest mean intensity (`selectHitsBy="intensity"`) or best
  identification level (`selectHitsBy="level"`) is kept. In case of ties
  only the first hit is kept. Set to `NULL` to ignore this filter. If
  `negate=TRUE` then only those hits with lowest mean intensity/poorest
  identification level are kept.

- selectBestFGroups:

  If `TRUE` then for any cases where a single feature group is matched
  to several suspects only the suspect assigned to the feature group
  with best identification score is kept. In case of ties only the first
  is kept.

- maxLevel, maxFormRank, maxCompRank, minAnnSimForm, minAnnSimComp,
  minAnnSimBoth:

  Filter suspects by maximum identification level (*e.g.* `"3a"`),
  formula/compound rank or with minimum formula/compound/combined
  annotation similarity. Set to `NULL` to ignore.

- absMinFragMatches, relMinFragMatches:

  Only retain suspects with this minimum number MS/MS matches with the
  fragments specified in the suspect list (*i.e.*
  `fragments_mz`/`fragments_formula`). `relMinFragMatches` sets the
  minimum that is relative (`0-1`) to the maximum number of MS/MS
  fragments specified in the `fragments_*` columns of the suspect list.
  Set to `NULL` to ignore.

- minRF:

  Filter suspect hits by the given minimum predicted response factor (as
  calculated by
  [`predictRespFactors`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
  Set to `NULL` to ignore.

- maxLC50:

  Filter suspect hits by the given maximum toxicity (LC50) (as
  calculated by
  [`predictTox`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
  Set to `NULL` to ignore.

- negate:

  If set to `TRUE` then filtering operations are performed in opposite
  manner.

- applyIMS:

  (**IMS workflow**) whether the filters are only applied to IMS
  precursors (`applyIMS=FALSE`), only to IMS features (`applyIMS=TRUE`)
  or to both (`applyIMS="both"`). Other feature groups will always be
  kept. The `negate` option does not affect `applyIMS`.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`).

- set:

  (**sets workflow**) The name of the set.

## Value

`delete` returns the object for which the specified data was removed.

`filter` returns a filtered `featureGroupsScreening` object.

## Methods (by generic)

- `screenInfo(featureGroupsScreening)`: Returns a table with screening
  information (see `screenInfo` slot).

- `show(featureGroupsScreening)`: Shows summary information for this
  object.

- `x[i`: Subset on analyses, feature groups and/or suspects.

- `delete(featureGroupsScreening)`: Completely deletes specified feature
  groups or screening results.

- `filter(featureGroupsScreening)`: Performs rule based filtering. This
  method builds on the comprehensive filter functionality from the base
  [`filter,featureGroups-method`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md).
  It adds several filters to select *e.g.* the best ranked suspects or
  those with a minimum estimated identification level. **NOTE**: most
  filters *only* affect suspect hits, not feature groups. Set
  `onlyHits=TRUE` to subsequently remove any feature groups that lost
  any suspect matches due to these filter steps.

## Slots

- `screenInfo`:

  A (`data.table`) with results from suspect screening. This table will
  be amended with ID confidence data when
  [`estimateIDConfidence`](https://rickhelmus.github.io/patRoon/reference/id-conf.md)
  is run.

- `MS2QuantMeta`:

  Metadata from MS2Quant filled in by `predictRespFactors`.

  (**sets workflow**) A named `list` with the metadata stored for each
  set.

## Note

`filter` removes suspect hits with `NA` values when any of the filters
related to minimum or maximum values are applied (unless `negate=TRUE`).

## S4 class hierarchy

- [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

  - **`featureGroupsScreening`**

    - `featureGroupsSetScreeningUnset`

## Sets workflows

The `featureGroupsScreeningSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `featureGroupsScreening` and therefore
largely follows the same user interface.

The following methods are specifically defined for sets workflows:

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`featureGroupsScreeningUnset`), which allows it to be used in
  'regular' workflows. Only the screening results present in the
  specified set are kept.

The following methods are changed or with new functionality:

- `estimateIDConfidence` See the `Sets workflows` section in the
  documentation for
  [`estimateIDConfidence`](https://rickhelmus.github.io/patRoon/reference/id-conf.md).

- `filter` All filters related to estimated identification levels and
  formula/compound rankings are applied to the overall set data (see
  above). All others are applied to set specific data: in this case
  candidates are only removed if none of the set data confirms to the
  filter.

This class derives also from
[`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).
Please see its documentation for more relevant details with sets
workflows.

Note that the `formRank` and `compRank` columns are *not* updated when
the data is subset.

## See also

[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
