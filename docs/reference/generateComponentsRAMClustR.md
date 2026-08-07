# Componentization of adducts, isotopes etc. with RAMClustR

Uses [RAMClustR](https://github.com/cbroeckl/RAMClustR) to generate
components from feature groups which follow similar chromatographic
retention profiles and annotate their relationships (*e.g.* adducts and
isotopes).

## Usage

``` r
generateComponentsRAMClustR(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsRAMClustR(
  fGroups,
  ionization = NULL,
  st = NULL,
  sr = NULL,
  maxt = 12,
  hmax = 0.3,
  normalize = "TIC",
  absMzDev = defaultLim("mz", "narrow"),
  relMzDev = defaultLim("mz", "narrow_rel"),
  minSize = 2,
  relMinReplicates = 0.5,
  RCExperimentVals = list(design = list(platform = "LC-MS"), instrument = list(ionization
    = ionization, MSlevs = 1)),
  extraOptsRC = NULL,
  extraOptsFM = NULL
)

# S4 method for class 'featureGroupsSet'
generateComponentsRAMClustR(fGroups, ionization = NULL, ...)
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

- st, sr, maxt, hmax, normalize:

  Arguments to tune the behaviour of feature group clustering. See their
  documentation from
  [`ramclustR`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html). When
  `st` is `NULL` it will be automatically calculated as the half of the
  median for all chromatographic peak widths.

- absMzDev:

  Maximum absolute *m/z* deviation. Sets the `mzabs.error` argument to
  [`do.findmain`](https://rdrr.io/pkg/RAMClustR/man/do.findmain.html)

- relMzDev:

  Maximum relative mass deviation (ppm). Sets the `ppm.error` argument
  to
  [`do.findmain`](https://rdrr.io/pkg/RAMClustR/man/do.findmain.html).

- minSize:

  The minimum size of a component. Smaller components than this size
  will be removed. See note below. Sets the `minModuleSize` argument to
  [`ramclustR`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html).

- relMinReplicates:

  Feature groups within a component are only kept when they contain data
  for at least this (relative) amount of replicate analyses. For
  instance, `0.5` means that at least half of the replicates should
  contain data for a particular feature group in a component. In this
  calculation replicates that are fully absent within a component are
  not taken in to account. See note below.

- RCExperimentVals:

  A named `list` containing two more `list`s: `design` and `instrument`.
  These are used to construct the `ExpDes` argument passed to
  [`ramclustR`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html).

- extraOptsRC, extraOptsFM:

  Named `list` with further arguments to be passed to
  [`ramclustR`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html) and
  [`do.findmain`](https://rdrr.io/pkg/RAMClustR/man/do.findmain.html).
  Set to `NULL` to ignore.

## Value

A
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
(derived) object containing all generated components.

## Details

This function uses RAMClustR to generate components. This function is
called when calling `generateComponents` with `algorithm="ramclustr"`.

This method uses the
[`ramclustR`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html)
functions for generating the components, whereas
[`do.findmain`](https://rdrr.io/pkg/RAMClustR/man/do.findmain.html) is
used for annotation.

## Note

The default value for `relMinReplicates` results in extra filtering,
hence, the final results may be different than what the algorithm
normally would return.

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

Broeckling, Heuberger CD;, Prince AL;, Ingelsson JA;, Prenni E;, E. J
(2013). “Assigning precursor-product ion relationships in indiscriminant
MS/MS data from non-targeted metabolite profiling studies.” *Analytical
Chemistry*, **9**, 33-43.\
\
Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE (2014).
“RAMClust: A Novel Feature Clustering Method Enables
Spectral-Matching-Based Annotation for Metabolomics Data.” *Analytical
Chemistry*, **86 (14)**, 6812–6817.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
