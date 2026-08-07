# Componentization of adducts, isotopes etc. with OpenMS

Uses the
[MetaboliteAdductDecharger](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MetaboliteAdductDecharger.html)
utility (see <http://www.openms.de>) to generate components.

## Usage

``` r
generateComponentsOpenMS(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsOpenMS(
  fGroups,
  ionization = NULL,
  chargeMin = 1,
  chargeMax = 1,
  chargeSpan = 3,
  qTry = "heuristic",
  potentialAdducts = NULL,
  minRTOverlap = 0.66,
  retWindow = defaultLim("retention", "very_narrow"),
  absMzDev = defaultLim("mz", "medium"),
  minSize = 2,
  relMinAdductAbundance = 0.75,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
  prefAdducts = c("[M+H]+", "[M-H]-"),
  extraOpts = NULL
)

# S4 method for class 'featureGroupsSet'
generateComponentsOpenMS(
  fGroups,
  ionization = NULL,
  chargeMin = 1,
  chargeMax = 1,
  chargeSpan = 3,
  qTry = "heuristic",
  potentialAdducts = NULL,
  ...
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which components should be generated.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- ionization:

  Which ionization polarity was used to generate the data: should be
  `"positive"` or `"negative"`. If the `featureGroups` object has adduct
  annotations, and `ionization=NULL`, the ionization will be detected
  automatically.

  (**sets workflow**) This parameter is not supported for sets
  workflows, as the ionization will always be detected automatically.

- chargeMin, chargeMax:

  The minimum/maximum charge to consider. Corresponds to the
  `algorithm:MetaboliteFeatureDeconvolution:charge_min`/`algorithm:MetaboliteFeatureDeconvolution:charge_min`
  options.

- chargeSpan:

  The maximum charge span for a single analyte. Corresponds to
  `algorithm:MetaboliteFeatureDeconvolution:charge_span_max`.

- qTry:

  Sets how charges are determined. Corresponds to
  `algorithm:MetaboliteFeatureDeconvolution:q_try`. Valid options are
  `"heuristic"` and `"all"` (the `"feature"` option from `OpenMS` is
  currently not supported).

- potentialAdducts:

  The adducts to consider. Should be a `numeric` vector with
  probabilities for each adduct, *e.g.*
  `potentialAdducts=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2)`. Note that the
  sum of probabilities should always be `1`. Furthermore, note that
  additions of multiple adducts should be controlled by the
  `chargeMin`/`chargeMax` arguments (and *not* with `potentialAdducts`),
  *e.g.* if `chargeMax=2` then both `[M+H]+` and `[2M+H]2+` may be
  considered. Please see the
  `algorithm:MetaboliteFeatureDeconvolution:potential_adducts` option of
  [MetaboliteAdductDecharger](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MetaboliteAdductDecharger.html)
  for more details. If `NULL` then the a default is chosen with
  [`defaultOpenMSAdducts`](https://rickhelmus.github.io/patRoon/reference/defaultOpenMSAdducts.md)
  (which is *not* the same as `OpenMS`).

  (**sets workflow**) Should be a `list` where each entry specifies the
  potential adducts for a set. Should either be named with the sets
  names or follow the same order as `sets(fGroups)`. Example:
  `potentialAdducts=list(positive=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2), negative=c("[M-H]-" = 0.8, "[M-H2O-H]-" = 0.2))`

- minRTOverlap, retWindow:

  Sets feature retention tolerances when grouping features. Sets the
  `algorithm:MetaboliteFeatureDeconvolution:min_rt_overlap` and
  `"algorithm:MetaboliteFeatureDeconvolution:retention_max_diff"`
  options.

- absMzDev:

  Maximum absolute *m/z* deviation. Sets the
  `algorithm:MetaboliteFeatureDeconvolution:mass_max_diff` option

- minSize:

  The minimum size of a component. Smaller components than this size
  will be removed. See note below.

- relMinAdductAbundance:

  The minimum relative abundance (`0-1`) that an adduct should be
  assigned to features within the same feature group. See the
  `Feature components` section for more details.

- adductConflictsUsePref:

  If set to `TRUE`, and not all adduct assigments to the features within
  a feature group are equal and at least one of those adducts is a
  preferential adduct (`prefAdducts` argument), then only the features
  with (the lowest ranked) preferential adduct are considered. In all
  other cases or when `adductConflictsUsePref=FALSE` only features with
  the most frequently assigned adduct is considered. See the
  `Feature components` section for more details.

- NMConflicts:

  The strategies to employ when not all neutral masses within a
  component are equal. Valid options are: `"preferential"`,
  `"mostAbundant"` and `"mostIntense"`. Multiple strategies are
  possible, and will be executed in the given order until one succeeds.
  See the `Feature components` section for more details.

- prefAdducts:

  A `character` vector with one or more *preferential adducts*. See the
  `Feature components` section for more details.

- extraOpts:

  Named character vector with extra command line parameters directly
  passed to `MetaboliteAdductDecharger`. Set to `NULL` to ignore.

## Value

A
[`componentsFeatures`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
derived object.

## Details

This function uses OpenMS to generate components. This function is
called when calling `generateComponents` with `algorithm="openms"`.

Features that show highly similar chromatographic elution profiles are
grouped, and subsequently annotated with their adducts.

## Feature components

The returned components are based on so called *feature components*.
Unlike other algorithms, components are first made on a feature level
(per analysis), instead of for complete feature groups. In the final
step the feature components are converted to 'regular' components by
employing a consensus approach with the following steps:

1.  If an adduct assigned to a feature only occurs as a minority
    compared to other adduct assigments within the same feature group,
    it is considered as an outlier and removed accordingly (controlled
    by the `relMinAdductAbundance` argument).

2.  For features within a feature group, only keep their adduct
    assignment if it occurs as the most frequent or is preferential
    (controlled by `adductConflictsUsePref` and `prefAdducts`
    arguments).

3.  Components are made by combining the feature groups for which at
    least one of their features are jointly present in the same feature
    component.

4.  Conflicts of neutral mass assignments within a component (*i.e.* not
    all are the same) are dealt with. Firstly, all feature groups with
    an unknown neutral mass are split in another component. Then, if
    conflicts still occur, the feature groups with similar neutral mass
    (determined by `absMzDev` argument) are grouped. Depending on the
    `NMConflicts` argument, the group with one or more preferential
    adduct(s) or that is the largest or most intense is selected,
    whereas others are removed from the component. In case multiple
    groups contain preferential adducts, and `>1` preferential adducts
    are available, the group with the adduct that matches first in
    `prefAdducts` 'wins'. In case of ties, one of the next strategies in
    `NMConflicts` is tried.

5.  If a feature group occurs in multiple components it will be removed
    completely.

6.  the `minSize` filter is applied.

## IMS workflows

The componentization algorithm is not aware of the IMS dimension. For
this reason, no IMS feature groups will be considered for
componentization, and direct IMS workflows (see
[`assignMobilitities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
are currently not supported.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the componentization is first performed for each set independently. The
resulting components are then all combined in a
[`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
object. Note that the components themselves are never merged. The
components are renamed to include the set name from which they were
generated (*e.g.* `"CMP1"` becomes `"CMP1-positive"`).

## Parallelization

generateComponentsOpenMS uses multiprocessing to parallelize
computations. Please see the parallelization section in the handbook for
more details and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

## References

Bielow C, Ruzek S, Huber CG, Reinert K (2010). “Optimal Decharging and
Clustering of Charge Ladders Generated in ESI-MS.” *Journal of Proteome
Research*, **9**(5), 2688–2695.
[doi:10.1021/pr100177k](https://doi.org/10.1021/pr100177k) .

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
