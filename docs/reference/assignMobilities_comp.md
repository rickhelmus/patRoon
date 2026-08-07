# Assign IMS data to a `compounds` object.

Assigns ion mobility and CCS values to the candidates in a
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
object.

## Usage

``` r
# S4 method for class 'compounds'
assignMobilities(
  obj,
  fGroups,
  IMS = TRUE,
  from = NULL,
  matchFromBy = "InChIKey1",
  overwrite = FALSE,
  adduct = NULL,
  CCSParams = NULL,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  virtualenv = "patRoon-C3SDB"
)

# S4 method for class 'compoundsSet'
assignMobilities(obj, fGroups, IMS = TRUE, from = NULL, ...)
```

## Arguments

- obj:

  The
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object for which IMS assignments should be performed.

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which the compound candidates (`obj`) were calculated.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  IMS assignments in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

- from, matchFromBy, overwrite, CCSParams, prefCalcChemProps,
  neutralChemProps, virtualenv:

  Passed to the [method for
  suspects](https://rickhelmus.github.io/patRoon/reference/assignMobilities_susp.md)

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

## Details

The `assignMobilities` method for
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
is used to (1) add predicted or known IMS data to annotation candidates
and (2) convert (previously added) mobility \<–\> CCS values.
Internally, the [assignMobilities method for
suspects](https://rickhelmus.github.io/patRoon/reference/assignMobilities_susp.md)
is used to perform these operations, please see its documentation for
more details.

If both adduct specific and non-adduct specific mobility and CCS data is
available (and not `NA`), then the non-adduct specific data is assigned
to the compound candidate. Otherwise, data corresponding to the `adduct`
argument or the adduct assigned to the feature group is taken. The
[`filter`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
method can be used to filter out candidates with deviating IMS data.

## Note

`SIRIUS` does currently not report InChIKey values, hence,
`matchFromBy="InChIKey"` is not supported in this case.

If compound annotation is performed with
[`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md)
with `database="pubchemlite"` then CCS data is already added, provided
the local database has it. If you want to overwrite this data, set
`overwrite=TRUE`.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the calculations are performed and stored independently for each set, as
adducts, charges and *m/z* values typically differ.

The `from` argument may be a `list` that specifies the `from` values for
each set. This is primarily intended when tables are used to set `from`.

## See also

The [assignMobilities method for
suspects](https://rickhelmus.github.io/patRoon/reference/assignMobilities_susp.md).
