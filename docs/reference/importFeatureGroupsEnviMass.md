# Imports feature groups from enviMass

Imports a 'profiles' produced by enviMass.

## Usage

``` r
importFeatureGroupsEnviMass(input, feat, positive)
```

## Arguments

- input:

  The path of the enviMass project.

- feat:

  The
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object obtained with
  [`importFeaturesEnviMass`](https://rickhelmus.github.io/patRoon/reference/importFeaturesEnviMass.md).

- positive:

  Whether data from positive (`TRUE`) or negative (`FALSE`) should be
  loaded.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function imports data from enviMass. This function is called when
calling `importFeatureGroups` with `type="envimass"`.

This function *only* imports 'raw' profiles, *not* any results from
further componentization steps performed in enviMass. Furthermore, this
functionality has only been tested with older versions of enviMass.
Finally, please note that this function only supports features imported
by
[`importFeaturesEnviMass`](https://rickhelmus.github.io/patRoon/reference/importFeaturesEnviMass.md)
(obviously, the same project should be used for both importing
functions).

## See also

[`importFeatureGroups`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroups.md)
for more details and other algorithms.
