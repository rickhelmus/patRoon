# Componentization of adducts, isotopes etc. with CAMERA

Interfaces with
[CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html)
to generate components from known adducts, isotopes and in-source
fragments.

## Usage

``` r
generateComponentsCAMERA(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsCAMERA(
  fGroups,
  ionization = NULL,
  onlyIsotopes = FALSE,
  minSize = 2,
  relMinReplicates = 0.5,
  extraOpts = NULL
)

# S4 method for class 'featureGroupsSet'
generateComponentsCAMERA(fGroups, ionization = NULL, ...)
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

- onlyIsotopes:

  Logical value. If `TRUE` only isotopes are considered when generating
  components (faster). Corresponds to `quick` argument of
  [`CAMERA::annotate`](https://rdrr.io/pkg/CAMERA/man/annotate-methods.html).

- minSize:

  The minimum size of a component. Smaller components than this size
  will be removed. See note below.

- relMinReplicates:

  Feature groups within a component are only kept when they contain data
  for at least this (relative) amount of replicate analyses. For
  instance, `0.5` means that at least half of the replicates should
  contain data for a particular feature group in a component. In this
  calculation replicates that are fully absent within a component are
  not taken in to account. See note below.

- extraOpts:

  Named character vector with extra arguments directly passed to
  [`CAMERA::annotate`](https://rdrr.io/pkg/CAMERA/man/annotate-methods.html).
  Set to `NULL` to ignore.

## Value

A
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
(derived) object containing all generated components.

## Details

This function uses CAMERA to generate components. This function is
called when calling `generateComponents` with `algorithm="camera"`.

The specified `featureGroups` object is automatically converted to an
`xcmsSet` object using
[`getXCMSSet`](https://rickhelmus.github.io/patRoon/reference/xcms-conv.md).

## Note

The default value for `minSize` and `relMinReplicates` results in extra
filtering, hence, the final results may be different than what the
algorithm normally would return.

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

Kuhl C, Tautenhahn R, Boettcher C, Larson TR, Neumann S (2012). “CAMERA:
an integrated strategy for compound spectra extraction and annotation of
liquid chromatography/mass spectrometry data sets.” *Analytical
Chemistry*, **84**, 283–289.
<http://pubs.acs.org/doi/abs/10.1021/ac202450g>.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
