# Componentization of adducts, isotopes etc. with cliqueMS

Uses [cliqueMS](https://github.com/osenan/cliqueMS) to generate
components using the
[`cliqueMS::getCliques`](https://rdrr.io/pkg/cliqueMS/man/getCliques.html)
function.

## Usage

``` r
generateComponentsCliqueMS(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsCliqueMS(
  fGroups,
  ionization = NULL,
  maxCharge = 1,
  maxGrade = 2,
  ppm = 10,
  adductInfo = NULL,
  absMzDev = defaultLim("mz", "medium"),
  minSize = 2,
  relMinAdductAbundance = 0.75,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
  prefAdducts = c("[M+H]+", "[M-H]-"),
  extraOptsCli = NULL,
  extraOptsIso = NULL,
  extraOptsAnn = NULL,
  parallel = TRUE
)

# S4 method for class 'featureGroupsSet'
generateComponentsCliqueMS(fGroups, ionization = NULL, ...)
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

- maxCharge, maxGrade, ppm:

  Arguments passed to
  [`cliqueMS::getIsotopes`](https://rdrr.io/pkg/cliqueMS/man/getIsotopes.html)
  and/or
  [`cliqueMS::getAnnotation`](https://rdrr.io/pkg/cliqueMS/man/getAnnotation.html).

- adductInfo:

  Sets the `adinfo` argument to
  [`cliqueMS::getAnnotation`](https://rdrr.io/pkg/cliqueMS/man/getAnnotation.html).
  If `NULL` then the default adduct information from cliqueMS is used
  (*i.e.* the `positive.adinfo`/`negative.adinfo` package datasets).

- absMzDev:

  Maximum absolute *m/z* deviation.

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

- extraOptsCli, extraOptsIso, extraOptsAnn:

  Named `list` with further arguments to be passed to
  [`cliqueMS::getCliques`](https://rdrr.io/pkg/cliqueMS/man/getCliques.html),
  [`cliqueMS::getIsotopes`](https://rdrr.io/pkg/cliqueMS/man/getIsotopes.html)
  and
  [`cliqueMS::getAnnotation`](https://rdrr.io/pkg/cliqueMS/man/getAnnotation.html),
  respectively. Set to `NULL` to ignore.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

## Value

A
[`componentsFeatures`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
derived object.

## Details

This function uses cliqueMS to generate components. This function is
called when calling `generateComponents` with `algorithm="cliquems"`.

The grouping of features in each component ('clique') is based on high
similarity of chromatographic elution profiles. All features in each
component are then annotated with the
[`cliqueMS::getIsotopes`](https://rdrr.io/pkg/cliqueMS/man/getIsotopes.html)
and
[`cliqueMS::getAnnotation`](https://rdrr.io/pkg/cliqueMS/man/getAnnotation.html)
functions.

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

## References

Senan O, Aguilar-Mogas A, Navarro M, Capellades J, Noon L, Burks D,
Yanes O, Guimera R, Sales-Pardo M (2019). “CliqueMS: a computational
tool for annotating in-source metabolite ions from LC-MS untargeted
metabolomics data based on a coelution similarity network.”
*Bioinformatics*, **35**(20), 4089–4097.
[doi:10.1093/bioinformatics/btz207](https://doi.org/10.1093/bioinformatics/btz207)
.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
