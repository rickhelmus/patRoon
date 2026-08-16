# Adduct utilities

Several utility functions to work with adducts.

## Usage

``` r
GenFormAdducts()

MetFragAdducts()

as.adduct(
  x,
  format = "generic",
  isPositive = NULL,
  charge = NULL,
  adductInfo = NULL,
  err = TRUE
)

calculateIonFormula(formula, adduct)

calculateNeutralFormula(formula, adduct)
```

## Arguments

- x:

  The object that should be converted. Should be a `character` string, a
  `numeric` `MetFrag` adduct identifier (`adduct_mode` column obtained
  with `MetFragAdducts`) or an
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (in which case no conversion occurs).

- format:

  A `character` that specifies the source format.

  `"generic"` is an internally used generic format that supports full
  textual conversion. Examples: `"[M+H]+"`, `"[2M+H]+"`, `"[M+3H]3+"`.

  `"sirius"` Is the format used by `SIRIUS`. It is similar to `generic`
  but does not allow multiple charges/molecules. See the SIRIUS manual
  for more details.

  `"genform"` and `"metfrag"` support fixed types of adducts which can
  be obtained with the `GenFormAdducts` and `MetFragAdducts` functions,
  respectively.

  `"openms"` is the format used by the `MetaboliteAdductDecharger` tool.

  `"cliquems"` is the format used by cliqueMS.

  `"nontarget"` is the format used by nontarget/enviPat and requires
  `adductInfo` to be set.

- isPositive:

  A logical that specifies whether the adduct should be positive. Should
  only be set when `format="metfrag"` and `x` is a `numeric` identifier.

- charge:

  The final charge. Only needs to be set when `format="openms"`.

- adductInfo:

  A `data.frame` with adduct info from *e.g.*
  [enviPat::adducts](https://rdrr.io/pkg/enviPat/man/adducts.html). Only
  needs to be set when `format="nontarget"`.

- err:

  If `TRUE` then an error will be thrown if conversion fails, otherwise
  returns without data.

- formula:

  A `character` vector with formulae to convert.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with `as.adduct`).
  Examples: `"[M-H]-"`, `"[M+Na]+"`.

## Details

`GenFormAdducts` returns a table with information on adducts supported
by `GenForm`.

`MetFragAdducts` returns a table with information on adducts supported
by `MetFrag`.

`as.adduct` Converts an object in to an
[`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
object.

`calculateIonFormula` Converts one or more neutral formulae to adduct
ions.

`calculateNeutralFormula` Converts one or more adduct ions to neutral
formulae.

## Examples

``` r
as.adduct("[M+H]+")
as.adduct("[M+H2]2+")
as.adduct("[2M+H]+")
as.adduct("[M-H]-")
as.adduct("+H", format = "genform")
as.adduct(1, isPositive = TRUE, format = "metfrag") # MetFrag adduct ID 1 --> returns [M+H]+

calculateIonFormula("C2H4O", "[M+H]+") # C2H5O
calculateNeutralFormula("C2H5O", "[M+H]+") # C2H4O
```
