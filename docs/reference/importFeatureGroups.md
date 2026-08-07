# Import feature groups from files

Generic function to import feature groups produced by other software
from files.

## Usage

``` r
importFeatureGroups(input, type, ...)
```

## Arguments

- input:

  The input object or path that should be imported. See the algorithm
  specific functions for more details.

- type:

  What type of data should be imported: `"xcms"`, `"xcms3"`, `"kpic2"`,
  `"table"`, `"brukerpa"` (Bruker ProfileAnalysis), `"brukertasq"`
  (Bruker TASQ) or `"envimass"`.

- ...:

  Further arguments passed to the selected import algorithm function.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

`importFeatureGroups` is a generic function that will import feature
groups from files by one of the supported algorithms. The actual
functionality is provided by algorithm specific functions such as
`importFeatureGroupsXCMS3` and `importFeatureGroupsTable`. While these
functions may be called directly, `importFeatureGroups` provides a
generic interface and is therefore usually preferred.

## See also

The
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
output class and its methods and the algorithm specific functions:
[`importFeatureGroupsXCMS`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsXCMS.md),
[`importFeatureGroupsXCMS3`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsXCMS3.md),
[`importFeatureGroupsKPIC2`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsKPIC2.md),
[`importFeatureGroupsTable`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsTable.md),
[`importFeatureGroupsBrukerPA`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerPA.md),
[`importFeatureGroupsBrukerTASQ`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerTASQ.md),
[`importFeatureGroupsEnviMass`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsEnviMass.md)

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
to group features.
