# Feature groups comparison class

This class is used for comparing different
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
objects.

## Usage

``` r
# S4 method for class 'featureGroupsComparison'
names(x)

# S4 method for class 'featureGroupsComparison'
length(x)

# S4 method for class 'featureGroupsComparison,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'featureGroupsComparison,ANY,missing'
x[[i, j]]

# S4 method for class 'featureGroupsComparison'
x$name
```

## Arguments

- x:

  A `featureGroupsComparison` object.

- i:

  For `[`/`[[`: A numeric or character value which is used to select
  labels by their index or name, respectively (for the order/names see
  [`names()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all labels are selected.\
  \
  For `[[`: should be a scalar value.

- ...:

  Ignored.

- drop, j:

  ignored.

- name:

  The label name (partially matched).

## Details

Objects from this class are returned by
[`comparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md).

## Methods (by generic)

- `names(featureGroupsComparison)`: Obtain the labels that were given to
  each compared feature group.

- `length(featureGroupsComparison)`: Number of feature groups objects
  that were compared.

- `x[i`: Subset on labels that were assigned to compared feature groups.

- `x[[i`: Extract a
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object by its label.

- `$`: Extract a compound table for a feature group.

## Slots

- `fGroupsList`:

  A `list` of
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that were compared

- `comparedFGroups`:

  A *pseudo* `featureGroups` object containing grouped feature groups.
