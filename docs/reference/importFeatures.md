# Import features

Generic function to import features produced by other software.

## Usage

``` r
importFeatures(input, type, ...)
```

## Arguments

- input:

  The input object or path that should be imported. See the algorithm
  specific functions for more details.

- type:

  What type of data should be imported: `"xcms"`, `"xcms3"`, `"kpic2"`,
  `"table"`, or `"envimass"`.

- ...:

  Further arguments passed to the selected import algorithm function.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

`importFeatures` is a generic function that will import features by one
of the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `importFeaturesXCMS3` and
`importFeaturesTable`. While these functions may be called directly,
`importFeatures` provides a generic interface and is therefore usually
preferred.

## See also

The
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
output class and its methods and the algorithm specific functions:
[`importFeaturesXCMS`](https://rickhelmus.github.io/patRoon/reference/importFeaturesXCMS.md),
[`importFeaturesXCMS3`](https://rickhelmus.github.io/patRoon/reference/importFeaturesXCMS3.md),
[`importFeaturesKPIC2`](https://rickhelmus.github.io/patRoon/reference/importFeaturesKPIC2.md),
[`importFeaturesTable`](https://rickhelmus.github.io/patRoon/reference/importFeaturesTable.md),
[`importFeaturesEnviMass`](https://rickhelmus.github.io/patRoon/reference/importFeaturesEnviMass.md)

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
to find new features.
