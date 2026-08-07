# Grouping of features

Group equal features across analyses.

## Usage

``` r
groupFeatures(obj, algorithm, ...)

# S4 method for class 'features'
groupFeatures(obj, algorithm, ..., verbose = TRUE)

# S4 method for class 'data.frame'
groupFeatures(obj, algorithm, ..., verbose = TRUE)
```

## Arguments

- obj:

  Either a
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object to be grouped, or a `data.frame` with [analysis
  info](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  to be passed to `groupFeaturesSIRIUS`

- algorithm:

  A `character` that specifies the algorithm to be used: either
  `"openms"`, `"xcms"`, `"xcms3"`, `"kpic2"` or `"greedy"` (`features`
  method), or `"sirius"` (`data.frame` method).

- ...:

  Further parameters passed to the selected grouping algorithm.

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

After [features have been
found](https://rickhelmus.github.io/patRoon/reference/findFeatures.md),
the next step is to align and group them across analyses. This process
is necessary to allow comparison of features between multiple analyses,
which otherwise would be difficult due to small deviations in retention
and mass data. Thus, algorithms of 'feature groupers' are used to
collect features with similar retention and mass data. In addition,
advanced retention time alignment algorithms exist to enhance grouping
of features even with relative large retention time deviations (*e.g.*
possibly observed from analyses collected over a long period). Like
[findFeatures](https://rickhelmus.github.io/patRoon/reference/findFeatures.md),
various algorithms are supported which may have many parameters that can
be fine-tuned. This fine-tuning is likely to be necessary, since optimal
settings often depend on applied methodology and instrumentation.

`groupFeatures` is a generic function that will groupFeatures by one of
the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `groupFeaturesOpenMS` and
`groupFeaturesXCMS3`. While these functions may be called directly,
`groupFeatures` provides a generic interface and is therefore usually
preferred.

The `data.frame` method for `groupFeatures` is a special case that
currently only supports the `"sirius"` algorithm.

## See also

The
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
output class and its methods and the algorithm specific functions:
[`groupFeaturesOpenMS`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesOpenMS.md),
[`groupFeaturesXCMS`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesXCMS.md),
[`groupFeaturesXCMS3`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesXCMS3.md),
[`groupFeaturesKPIC2`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesKPIC2.md),
[`groupFeaturesGreedy`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesGreedy.md),
[`groupFeaturesSIRIUS`](https://rickhelmus.github.io/patRoon/reference/groupFeaturesSIRIUS.md)
