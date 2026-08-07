# Initiate sets workflows

Initiate sets workflows from specified feature data.

## Usage

``` r
makeSet(obj, ...)

# S4 method for class 'features'
makeSet(obj, ..., adducts, labels = NULL)

# S4 method for class 'featuresSet'
makeSet(obj, ...)

# S4 method for class 'featureGroups'
makeSet(
  obj,
  ...,
  groupAlgo,
  groupArgs = NULL,
  verbose = TRUE,
  adducts = NULL,
  labels = NULL
)

# S4 method for class 'featureGroupsSet'
makeSet(obj, ...)
```

## Arguments

- obj, ...:

  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  or
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  objects that should be used for the [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).

- adducts:

  The adduct assignments to each set. Should either be a `list` with
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  objects or a `character` vector (*e.g.* `"[M+H]+"`). The order should
  follow that of the objects given to the `obj` and `...` arguments.

  For the
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  method: if `NULL` then adduct annotations are used.

- labels:

  The labels, or *set names*, for each set to be created. The order
  should follow that of the objects given to the `obj` and `...`
  arguments. If `NULL`, then labels are automatically generated from the
  polarity of the specified `adducts` argument (*e.g.* `"positive"`,
  `"negative"`).

- groupAlgo:

  groupAlgo The name of the feature grouping algorithm. See the
  `algorithm` argument of
  [`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
  for details.

- groupArgs:

  A `list` with arguments directly passed to `groupFeatures` (can be
  named). Example: `groupArgs=list(maxAlignMZ=0.002)`.

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

Either a
[`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
object (`features` method) or
[`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object (`featureGroups` method).

## Details

The `makeSet` method function is used to initiate a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
The features from input objects are combined and then *neutralized* by
replacing their *m/z* values by neutral monoisotopic masses. After
neutralization features measured with *e.g.* different ionization
polarities can be grouped since their neutral mass will be the same.

The [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
for this object is updated with all analyses, and a `set` column is
added to designate the set of each analysis. Note that currently, all
analyses names **must** be unique across different sets.

`makeSet` supports two types of input:

1.  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
    objects: `makeSet` combines the input objects into a `featuresSet`
    object, which is then grouped in the 'usual way' with
    [`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md).

2.  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
    objects: In this case the features from the input objects are first
    neutralized and feature groups between sets are then combined with
    `groupFeatures`.

The advantage of the `featureGroups` method is that it preserves any
adduct annotations already present (*e.g.* as set by `selectIons` or
`adducts<-`). Furthermore, this approach allows more advanced workflows
where the input `featureGroups` are first pre-treated with *e.g.* filter
before the sets object is made. On the other hand, the `features` method
is easier, as it doesn't require intermediate feature grouping steps and
is often sufficient since adduct annotations can be made afterwards with
`selectIons`/`adducts<-` and most `filter` operations do not need to be
done per individual set.

The adduct information used for feature neutralization is specified
through the `adducts` argument. Alternatively, when the `featureGroups`
method of `makeSet` is used, then the adduct annotations already present
in the input objects can also by used by setting `adducts=NULL`. The
adduct information is also used to add adduct annotations to the output
of `makeSet`.

## Note

Initiating a sets workflow recursively, *i.e.* with `featuresSet` or
`featureGroupsSet` objects as input, is currently not supported.
