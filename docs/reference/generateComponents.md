# Grouping feature groups in components

Functionality to automatically group related feature groups (*e.g.*
isotopes, adducts and homologues) to assist and simplify annotation.

## Usage

``` r
generateComponents(fGroups, algorithm, ...)

# S4 method for class 'featureGroups'
generateComponents(fGroups, algorithm, ...)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which components should be generated.

- algorithm:

  A character string describing the algorithm that should be used:
  `"ramclustr"`, `"camera"`, `"nontarget"`, `"intclust"`, `"openms"`,
  `"cliquems"`, `"specclust"`, `"tp"`

- ...:

  Any parameters to be passed to the selected component generation
  algorithm.

## Value

A
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
(derived) object containing all generated components.

## Details

Several algorithms are provided to group feature groups that are related
in some (chemical) way to each other. How feature groups are related
depends on the algorithm: examples include adducts, statistics and
parents/transformation products. The linking of this data is generally
useful for annotation purposes and reducing data complexity.

`generateComponents` is a generic function that will generateComponents
by one of the supported algorithms. The actual functionality is provided
by algorithm specific functions such as `generateComponentsRAMClustR`
and `generateComponentsNontarget`. While these functions may be called
directly, `generateComponents` provides a generic interface and is
therefore usually preferred.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the componentization data is generated differently depending on the used
algorithm. Please see the details in the algorithm specific functions
linked in the `See Also` section.

## See also

The
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
output class and its methods and the algorithm specific functions:
[`generateComponentsRAMClustR`](https://rickhelmus.github.io/patRoon/reference/generateComponentsRAMClustR.md),
[`generateComponentsCAMERA`](https://rickhelmus.github.io/patRoon/reference/generateComponentsCAMERA.md),
[`generateComponentsNontarget`](https://rickhelmus.github.io/patRoon/reference/generateComponentsNontarget.md),
[`generateComponentsIntClust`](https://rickhelmus.github.io/patRoon/reference/generateComponentsIntClust.md),
[`generateComponentsOpenMS`](https://rickhelmus.github.io/patRoon/reference/generateComponentsOpenMS.md),
[`generateComponentsCliqueMS`](https://rickhelmus.github.io/patRoon/reference/generateComponentsCliqueMS.md),
[`generateComponentsSpecClust`](https://rickhelmus.github.io/patRoon/reference/generateComponentsSpecClust.md),
[`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
