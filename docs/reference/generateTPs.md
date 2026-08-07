# Generation of transformation products (TPs)

Functionality to automatically obtain transformation products for a
given set of parent compounds.

## Usage

``` r
generateTPs(algorithm, ...)
```

## Arguments

- algorithm:

  A character string describing the algorithm that should be used:
  `"biotransformer"`, `"logic"`, `"library"`, `"library_formula"`,
  `"cts"`, `"ann_form"`, `"ann_comp"`

- ...:

  Any parameters to be passed to the selected TP generation algorithm.

## Value

A
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
(derived) object containing all generated TPs.

## Details

`generateTPs` is a generic function that will generate transformation
products by one of the supported algorithms. The actual functionality is
provided by algorithm specific functions such as
`generateTPsBioTransformer` and `generateTPsLogic`. While these
functions may be called directly, `generateTPs` provides a generic
interface and is therefore usually preferred.

## See also

The
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
output class and its methods and the algorithm specific functions:
[`generateTPsBioTransformer`](https://rickhelmus.github.io/patRoon/reference/generateTPsBioTransformer.md),
[`generateTPsLogic`](https://rickhelmus.github.io/patRoon/reference/generateTPsLogic.md),
[`generateTPsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibrary.md),
[`generateTPsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibraryFormula.md),
[`generateTPsCTS`](https://rickhelmus.github.io/patRoon/reference/generateTPsCTS.md),
[`generateTPsAnnForm`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnForm.md),
[`generateTPsAnnComp`](https://rickhelmus.github.io/patRoon/reference/generateTPsAnnComp.md)

The derived classes
[`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md),
[`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md),
[`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md)
and
[`transformationProductsAnnComp`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnComp-class.md)
for more specific methods to post-process TP data.
