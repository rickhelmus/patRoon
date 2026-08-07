# Obtain default rules for metabolic logic

This function returns a `data.frame` with the default rules for
metabolic logic, which can be used by
[`generateTPsLogic`](https://rickhelmus.github.io/patRoon/reference/generateTPsLogic.md)
and
[`genFormulaTPLibrary`](https://rickhelmus.github.io/patRoon/reference/genFormulaTPLibrary.md).

## Usage

``` r
TPLogicTransformations()
```

## Value

A `data.frame` with columns describing each transformation rule.

## Source

The table is based on the work done by Schollee *et al.* (see
references).

## References

Schollee JE, Schymanski EL, Avak SE, Loos M, Hollender J (2015).
“Prioritizing Unknown Transformation Products from Biologically-Treated
Wastewater Using High-Resolution Mass Spectrometry, Multivariate
Statistics, and Metabolic Logic.” *Analytical Chemistry*, **87**(24),
12121–12129.
[doi:10.1021/acs.analchem.5b02905](https://doi.org/10.1021/acs.analchem.5b02905)
.
