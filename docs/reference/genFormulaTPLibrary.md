# Automatically generate a transformation product library with formula data.

Functionality to automatically generate a TP library with formula data
from a set of transformation rules, which can be used with
[`generateTPsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibraryFormula.md).
TP calculation will be skipped if the transformation involves
subtraction of elements not present in the parent.

## Usage

``` r
genFormulaTPLibrary(
  parents,
  transformations = NULL,
  minMass = 40,
  generations = 1,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE
)
```

## Arguments

- parents:

  The parents to which the given transformation rules should be used to
  generate the TP library. Should be either a suspect list (see [suspect
  screening](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
  for more information) or the resulting output of
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).

- transformations:

  A `data.frame` with transformation reactions to be used for
  calculating the TPs (see details below). If `NULL`, a default table
  from Schollee *et al.* is used, that can be obtained with
  [`TPLogicTransformations`](https://rickhelmus.github.io/patRoon/reference/TPLogicTransformations.md).

- minMass:

  The minimum mass for a TP to be kept.

- generations:

  An `integer` that specifies the number of transformation generations
  that should be calculated. If `generations>1` then TPs are calculated
  by applying the transformation rules to the TPs generated in the
  previous generation.

- skipInvalid:

  Set to `TRUE` to skip parents without formula information. Otherwise
  an error is thrown.

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the parent
  suspect list. For efficiency reasons it is recommended to set this to
  `TRUE`. See the `Validating and calculating chemical properties`
  section for more details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

## Value

A `data.table` that is suitable for the `TPLibrary` argument to
[`generateTPsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibraryFormula.md).

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

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the parent
suspect list are automatically validated and calculated if
missing/invalid.

The internal validation/calculation process performs the following
steps:

- Validation of SMILES, InChI, InChIKey and formula data (if present).
  Invalid entries will be set to `NA`.

- If `neutralChemProps=TRUE` then chemical data (SMILES, formulae etc.)
  is neutralized by (de-)protonation (using the `–neutralized` option of
  `OpenBabel`). An additional column `molNeutralized` is added to mark
  those molecules that were neutralized. Note that neutralization
  requires either SMILES or InChI data to be available.

- The SMILES and InChI data are used to calculate missing or invalid
  SMILES, InChI, InChIKey and formula data. If `prefCalcChemProps=TRUE`
  then existing InChIKey and formula data is overwritten by calculated
  values whenever possible.

- The chemical formulae which were *not* calculated are verified and
  normalized. This process may be time consuming, and is potentially
  largely avoided by setting `prefCalcChemProps=TRUE`.

- Neutral masses are calculated for missing values
  (`prefCalcChemProps=FALSE`) or whenever possible
  (`prefCalcChemProps=TRUE`).

Note that calculation of formulae for molecules that are isotopically
labelled is currently only supported for deuterium (2H) elements.

This functionality relies heavily on
[OpenBabel](https://github.com/openbabel/openbabel), please make sure it
is installed.

## References

Schollee JE, Schymanski EL, Avak SE, Loos M, Hollender J (2015).
“Prioritizing Unknown Transformation Products from Biologically-Treated
Wastewater Using High-Resolution Mass Spectrometry, Multivariate
Statistics, and Metabolic Logic.” *Analytical Chemistry*, **87**(24),
12121–12129.
[doi:10.1021/acs.analchem.5b02905](https://doi.org/10.1021/acs.analchem.5b02905)
.

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .

## See also

[`generateTPsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibraryFormula.md)
and
[`generateTPsLogic`](https://rickhelmus.github.io/patRoon/reference/generateTPsLogic.md)
