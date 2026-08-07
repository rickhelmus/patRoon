# Obtain transformation products (TPs) from a library with formula data

Automatically obtains transformation products from a library with
formula data.

## Usage

``` r
generateTPsLibraryFormula(
  parents = NULL,
  TPLibrary,
  generations = 1,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  matchParentsBy = "name",
  matchGenerationsBy = "name"
)
```

## Arguments

- parents:

  The parents for which transformation products should be obtained. This
  can be

  - a suspect list (see [suspect
    screening](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
    for more information)

  - the output of
    [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
    in which case the suspects hits are used as parents

  The parents need to have formula information available.

- TPLibrary:

  A `data.frame`. See the details below.

- generations:

  An `integer` that specifies the number of transformation generations.
  TPs for subsequent iterations obtained by repeating the library search
  where the TPs from the previous generation are considered parents.

- skipInvalid:

  If set to `TRUE` then the parents will be skipped (with a warning) for
  which insufficient information (*e.g.* SMILES) is available.

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

- matchParentsBy:

  A `character` that specifies how the input parents are matched with
  the data from the TP library. Valid options are: `"InChIKey"`,
  `"InChIKey1"`, `"InChI"`, `"SMILES"`, `"formula"`, `"name"`. If the
  parent from the TP library is matched with multiple input parents then
  only the first is considered.

- matchGenerationsBy:

  Similar to `matchParentsBy`, but specifies how parents/TPs are matched
  when `generations>1`.

## Value

The TPs are stored in an object derived from the
[`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)
class.

## Details

This function uses a library to obtain transformation products. This
function is called when calling `generateTPs` with
`algorithm="library_formula"`.

This function is similar to
[`generateTPsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibrary.md),
however, it only require formula information of the parent and TPs.

## Note

Unlike
[`generateTPsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibrary.md),
this function defaults the `matchParentsBy` and `matchGenerationsBy`
arguments to `"name"`. While matching by `formula` is also possible, it
is likely that duplicate parent formulae (*i.e.* isomers) are present in
`parents` and/or `TPLibrary`, making matching by formula unsuitable.
However, if you are sure that no duplicate formulae are present, it may
be better to set the matching method to `"formula"`.

## TP libraries

The `TPLibrary` argument is used to specify a custom TP library. This
should be a `data.frame` where each row specifies a TP for a parent,
with the following columns:

- `parent_name` and `TP_name`: The name of the parent/TP.

- `parent_formula` and `TP_formula` The formula of the parent/TP
  structure.

- `retDir` The expected [retention order
  direction](https://rickhelmus.github.io/patRoon/reference/retDir.md).
  (**optional**)

  For `generateTPsLibrary`: If not specified or `forceCalcRetDir=TRUE`
  from
  [`TPStructParams`](https://rickhelmus.github.io/patRoon/reference/getDefTPStructParams.md),
  then the `log P` values below may be used to calculate retention order
  directions.

- `parent_LogP` and `TP_LogP` The `log P` values for the parent/TP.
  (**optional**)

- `logPDiff` The difference between parent and TP `Log P` values.
  Ignored if *both* `parent_LogP` and `TP_LogP` are specified.
  (**optional**)

Other columns are allowed, and will be included in the final object.
Multiple TPs for a single parent are specified by repeating the value
within `parent_` columns.

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

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .

## See also

[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
for more details and other algorithms.

[`generateTPsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateTPsLibrary.md)
to generate TPs from a library that contains structural information.

[`genFormulaTPLibrary`](https://rickhelmus.github.io/patRoon/reference/genFormulaTPLibrary.md)
to automatically generate formula TP libraries.
