# Obtain transformation products (TPs) with metabolic logic

Automatically calculate potential transformation products with
*metabolic logic*.

## Usage

``` r
generateTPsLogic(fGroups, minMass = 40, ...)

# S4 method for class 'featureGroups'
generateTPsLogic(fGroups, minMass = 40, adduct = NULL, transformations = NULL)

# S4 method for class 'featureGroupsSet'
generateTPsLogic(fGroups, minMass = 40, transformations = NULL)
```

## Arguments

- fGroups:

  A
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which TPs should be calculated.

- minMass:

  A `numeric` that specifies the minimum mass of calculated TPs. If
  below this mass it will be removed.

- ...:

  Further arguments specified to the methods.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- transformations:

  A `data.frame` with transformation reactions to be used for
  calculating the TPs (see details below). If `NULL`, a default table
  from Schollee *et al.* is used, that can be obtained with
  [`TPLogicTransformations`](https://rickhelmus.github.io/patRoon/reference/TPLogicTransformations.md).

## Value

A
[`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
(derived) object containing all generated TPs.

## Details

This function uses metabolic logic to obtain transformation products.
This function is called when calling `generateTPs` with
`algorithm="logic"`.

With this algorithm TPs are predicted from common (environmental)
chemical reactions, such as hydroxylation, demethylation etc. The
generated TPs result from calculating the mass differences between a
parent feature after it underwent the reaction. While this only results
in little information on chemical properties of the TP, an advantage of
this method is that it does not rely on structural information of the
parent, which may be unknown in a full non-target analysis.

## Transformation reactions

The `transformations` argument specifies custom rules to calculate
transformation products. This should be a `data.frame` with the
following columns:

- `transformation` The name of the chemical transformation

- `add` The elements that are added by this reaction (*e.g.* `"O"`).

- `sub` The elements that are removed by this reaction (*e.g.* `"H2O"`).

- `retDir` The expected [retention order
  direction](https://rickhelmus.github.io/patRoon/reference/retDir.md).

## Source

The algorithms using transformation reactions are directly based on the
work done by Schollee *et al.* (see references).

## References

Schollee JE, Schymanski EL, Avak SE, Loos M, Hollender J (2015).
“Prioritizing Unknown Transformation Products from Biologically-Treated
Wastewater Using High-Resolution Mass Spectrometry, Multivariate
Statistics, and Metabolic Logic.” *Analytical Chemistry*, **87**(24),
12121–12129.
[doi:10.1021/acs.analchem.5b02905](https://doi.org/10.1021/acs.analchem.5b02905)
.

## See also

[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
for more details and other algorithms.
