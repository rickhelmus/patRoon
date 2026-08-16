# Generic adduct class

Objects from this class are used to specify adduct information in an
algorithm independent way.

## Usage

``` r
adduct(...)

# S4 method for class 'adduct'
show(object)

# S4 method for class 'adduct'
as.character(x, format = "generic", adductInfo = NULL, err = TRUE)
```

## Arguments

- x, object:

  An `adduct` object.

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

- adductInfo:

  A `data.frame` with adduct info from *e.g.*
  [enviPat::adducts](https://rdrr.io/pkg/enviPat/man/adducts.html). Only
  needs to be set when `format="nontarget"`.

- err:

  If `TRUE` then an error will be thrown if conversion fails, otherwise
  returns without data.

- ...:

  Any of `add`, `sub`, `molMult` and/or `charge`. See `Slots`.

## Methods (by generic)

- `show(adduct)`: Shows summary information for this object.

- `as.character(adduct)`: Converts an `adduct` object to a specified
  `character` format.

## Slots

- `add,sub`:

  A `character` with one or more formulas to add/subtract.

- `molMult`:

  How many times the original molecule is present in this molecule
  (*e.g.* for a dimer this would be `2`). Default is `1`.

- `charge`:

  The final charge of the adduct (default `1`).

## See also

[`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)
for easy creation of `adduct` objects and [adduct
utilities](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)
for other adduct functionality.

## Examples

``` r
adduct("H") # [M+H]+
adduct(sub = "H", charge = -1) # [M-H]-
adduct(add = "K", sub = "H2", charge = -1) # [M+K-H2]+
adduct(add = "H3", charge = 3) # [M+H3]3+
adduct(add = "H", molMult = 2) # [2M+H]+

as.character(adduct("H")) # returns "[M+H]+"
```
