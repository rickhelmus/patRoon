# Group features using greedy algorithm

Group features using a greedy algorithm that maximizes group scores
based on retention time, m/z, mobility, and intensity similarities.

## Usage

``` r
groupFeaturesGreedy(feat, ...)

# S4 method for class 'features'
groupFeaturesGreedy(
  feat,
  rtalign = FALSE,
  rtWindow = defaultLim("retention", "medium"),
  mzWindow = defaultLim("mz", "medium"),
  mobWindow = defaultLim("mobility", "medium"),
  scoreWeights = c(retention = 1, mz = 1, mobility = 1, intensity = 1),
  verbose = TRUE
)
```

## Arguments

- feat:

  The
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object with the features to be grouped.

- ...:

  Further parameters passed to the selected grouping algorithm.

- rtalign:

  Not yet supported. Provided for consistency with other grouping
  methods.

- rtWindow, mzWindow, mobWindow:

  Numeric tolerances for retention time (seconds), *m/z*, and mobility,
  respectively. The scoring terms are normalized to these values.
  Defaults to `defaultLim("retention", "medium")`,
  `defaultLim("mz", "medium")`, and `defaultLim("mobility", "medium")`,
  respectively (see
  [limits](https://rickhelmus.github.io/patRoon/reference/limits.md)).

- scoreWeights:

  Numeric vector specifying the scoring weights. Should contain the
  following named elements: `"retention"`, `"mz"`, `"mobility"`, and
  `"intensity"`. Missing elements are defaulted to one.

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function uses greedy to group features. This function is called
when calling `groupFeatures` with `algorithm="greedy"`.

The `greedy` algorithm is a simple feature grouping algorithm that can
work with both HRMS and IMS-HRMS data. The algorithm groups features by
iteratively building the best possible groups. Features are processed in
order of decreasing intensity. For each feature, candidate groups are
formed from all other (ungrouped) features within the specified
retention time, *m/z* and mobility windows. Each candidate group only
contains a maximum of one feature per analysis. The candidates are then
scored and the group with the lowest overall variations in retention
time, *m/z*, mobility and replicate intensity is then selected. This
process is repeated until all features have been assigned to a group.
The weights for each of the scoring terms can be configured.

## Note

Any links between IMS precursors and IMS features are removed. This can
occur *e.g.* when `greedy` is used to generate a [feature
consensus](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)
from a [post mobility
assignment](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)
workflow.

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for more details and other algorithms.
