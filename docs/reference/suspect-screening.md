# Target and suspect screening

Utilities to screen for analytes with known or suspected identity.

## Usage

``` r
screenSuspects(
  fGroups,
  suspects,
  rtWindow = defaultLim("retention", "medium"),
  mzWindow = defaultLim("mz", "medium"),
  IMSMatchParams = NULL,
  adduct = NULL,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  onlyHits = FALSE,
  ...
)

# S4 method for class 'featureGroups'
screenSuspects(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  IMSMatchParams,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits
)

# S4 method for class 'featureGroupsScreening'
screenSuspects(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  IMSMatchParams,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits,
  amend = FALSE
)

# S4 method for class 'featureGroupsSet'
screenSuspects(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  IMSMatchParams,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits
)

# S4 method for class 'featureGroupsScreeningSet'
screenSuspects(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  IMSMatchParams,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits,
  amend = FALSE
)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that should be screened.

- suspects:

  A `data.frame` with suspect information. See the `Suspect list format`
  section below.

  (**sets workflow**) Can also be a `list` with suspect lists to be used
  for each set (otherwise the same suspect lists is used for all sets).
  The `list` can be named with the sets names to mark which suspect list
  is to be used with which set (*e.g.*
  `suspects=list(positive=suspsPos, negative=suspsNeg)`).

- rtWindow, mzWindow:

  The retention time window (in seconds) and *m/z* window that will be
  used for matching a suspect (+/- feature data).

- IMSMatchParams:

  (**IMS workflow**) A `list` with parameters to be used for matching
  IMS data. See
  [`getIMSMatchParams`](https://rickhelmus.github.io/patRoon/reference/getIMSMatchParams.md)
  for details and how to make such a parameter list.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. May be `NULL`, see
  `Suspect list format` and `Matching of suspect masses` sections below.

- skipInvalid:

  If set to `TRUE` then suspects with invalid data (*e.g.* missing names
  or other missing data) will be ignored with a warning. Similarly, any
  suspects for which mass calculation failed (when no `mz` column is
  present in the suspect list), for instance, due to invalid `SMILES`,
  will be ignored with a warning.

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the suspect
  list. For efficiency reasons it is recommended to set this to `TRUE`.
  See the `Validating and calculating chemical properties` section for
  more details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

- onlyHits:

  If `TRUE` then all feature groups not matched by any of the suspects
  will be removed.

- ...:

  Further arguments specified to the methods.

- amend:

  If `TRUE` then screening results will be *amended* to the original
  object.

## Value

`screenSuspects` returns a
[`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)
object, which is a copy of the input `fGroups` object amended with
additional screening information.

## Details

Besides 'full non-target analysis', where compounds may be identified
with little to no prior knowledge, a common strategy is to screen for
compounds with known or suspected identity. This may be a generally
favorable approach if possible, as it can significantly reduce the load
on data interpretation.

`screenSuspects` is used to perform suspect screening. The input
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object will be screened for suspects by *m/z* values and optionally
retention times. Afterwards, any feature groups not matched may be kept
or removed, depending whether a full non-target analysis is desired.

## Note

`screenSuspects` may use the suspect names to base file names used for
reporting, logging etc. Therefore, it is important that these are
file-compatible names. For this purpose, `screenSuspects` will
automatically try to convert long, non-unique and/or otherwise
incompatible suspect names.

## Suspect list format

the `suspects` argument for `screenSuspects` should be a `data.frame`
with the following mandatory and optional columns:

- `name` The suspect name. Must be file-compatible. (**mandatory**)

- `rt` The retention time (in seconds) for the suspect. If specified the
  suspect will only be matched if its retention matches the experimental
  value (tolerance defined by the `rtWindow` argument). (**optional**)

- `neutralMass`,`formula`,`SMILES`,`InChI` The neutral monoisotopic
  mass, chemical formula, SMILES or InChI for the suspect. (data from
  one of these columns are **mandatory** in case no value from the `mz`
  column is available for a suspect)

- `mz` The ionized *m/z* of the suspect. (**mandatory** unless it can be
  calculated from one of the aforementioned columns)

- `adduct` A `character` that can be converted with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md).
  Can be used to automatically calculate values for the `mz` column.
  (**mandatory** unless data from the `mz` column is available, the
  `adduct` argument is set or `fGroups` has adduct annotations)

- `fragments_mz`,`fragments_formula` One or more MS/MS fragments
  (specified as *m/z* or formulae, respectively). Multiple values can be
  specified by separating them with a semicolon (`;`). This data is used
  by
  [`estimateIDConfidence`](https://rickhelmus.github.io/patRoon/reference/id-conf.md)
  to report detected MS/MS fragments and calculate identification
  levels. (**optional**)

- `mobility`,`CCS` The mobility or CCS value of the suspect. These
  values may be used to filter out suspects, see the `IMSMatchParams`
  argument. Multiple values for a single suspect can be specified by
  separating them with a semicolon(`;`). Adduct specific columns may be
  added by suffixing the adduct to the column name, *e.g.*
  `mobility_[M+H]+` and `CCS_[M-H]-`. (**optional**)

## Matching of suspect masses

How the mass of a suspect is matched with the mass of a feature depends
on the available data:

- If the suspect has data from the `mz` column of the suspect list, then
  this data is matched with the detected feature *m/z*.

- Otherwise, if the suspect has data in the `adduct` column of the
  suspect list, this data is used to calculate its `mz` value, which is
  then used like above.

- In the last case, the neutral mass of the suspect is matched with the
  neutral mass of the feature. Hence, either the `adduct` argument needs
  to be specified, or the `featureGroups` input object must have adduct
  annotations.

## IMS reference assignment

If both adduct specific and non-adduct specific reference values are
available, then non-adduct specific data is chosen (unless `NA`) as
reference for the suspect hit. Otherwise, data is taken from the adduct
specific data corresponding to the adduct assigned to the feature group
(or `adduct` argument). If multiple mobility or CCS values for a suspect
are specified in the suspect list, then the reference value is chosen
which is the closest to that of the feature.

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the suspect
list are automatically validated and calculated if missing/invalid.

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

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md),
`screenSuspects` performs suspect screening for each set separately, and
the screening results are combined afterwards. The `sets` column in the
`screenInfo` data marks in which sets the suspect hit was found.

## References

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .

## See also

`featureGroupsScreening`
