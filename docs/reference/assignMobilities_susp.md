# Assign IMS data to suspects

Adds calculated mobility and/or CCS data to a suspect list.

## Usage

``` r
# S4 method for class 'data.table'
assignMobilities(
  obj,
  from = NULL,
  matchFromBy = "InChIKey1",
  overwrite = FALSE,
  adducts = c("[M+H]+", "[M-H]-", NA),
  predictAdductOnly = TRUE,
  CCSParams = NULL,
  prepareChemProps = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  virtualenv = "patRoon-C3SDB"
)

# S4 method for class 'data.frame'
assignMobilities(obj, ...)
```

## Arguments

- obj:

  The suspect list to which the mobility and/or CCS data should be
  added. Should be a `data.frame` or `data.table`.

- from:

  Specifies from where IMS data is added to the suspect list. This can
  be the following:

  - `"pubchemlite"`: CCS data is matched from predicted values of the
    [PubChemLite](https://pubchemlite.lcsb.uni.lu/) database. **Note**:
    this requires a local copy of the [CCS amended PubChemLite
    database](https://zenodo.org/records/15311000) (see the Handbook for
    more details), which is automatically installed by patRoonExt.

  - `"c3sdb"`: Uses the [C3SDB](https://github.com/dylanhross/c3sdb)
    `Python` package to predict CCS values. This requires a local
    installation of `C3SDB`, *e.g.* performed with
    [`installC3SDB`](https://rickhelmus.github.io/patRoon/reference/installC3SDB.md).

  - A `data.table` or `data.frame` to which IMS data is matched. Should
    contain the column defined by `matchFromBy` and columns storing
    (non-)adduct specific mobility/CCS columns (see Details).

  - `NULL`: No IMS data is added to the suspect list.

  Any `NA` values in `from` are ignored.

- matchFromBy:

  Which column should be used to match the IMS data from `from` and
  suspects. Valid options are `"InChIKey"`, `"InChIKey1"` (first block
  InChIKey), `"InChI"`, `"SMILES"`, `"name"`. However, this also depends
  on which columns are available in either of the data sources.
  `InChIKey1` values are automatically calculated from `InChIKey`s, if
  possible.

  Matching by `InChiKey1` is recommended by default to allow tolerance
  between different datasources. Note that most compound annotation
  algorithms also match by `InChIKey1` and current IMS resolving power
  is generally insufficient to distinguish the different
  stereoisomers/tautomers specified by the full `InChIKey`.

- overwrite:

  Set to `TRUE` to overwrite any existing suspect IMS data with data
  from `from`.

- adducts:

  A `character` with the adduct(s) to consider for assigning mobility
  data to suspects and mobility \<–\> CCS conversions. This may be
  limited by what is available in the data source specified by `from`
  (see
  [`C3SDBAdducts`](https://rickhelmus.github.io/patRoon/reference/C3SDBAdducts.md)
  for `from="c3sdb"`). Inclusion of `NA` in `adducts` allows the use of
  non-adduct specific values (see Details).

  The value for `adducts` is automatically expanded by the adducts
  specified in the `adduct` column of the suspect list. Hence, `adducts`
  can be empty ([`character()`](https://rdrr.io/r/base/character.html))
  if no calculations for other adducts are desired.

- predictAdductOnly:

  If `from="c3sdb"` and `predictAdductOnly=TRUE` then only predictions
  are performed for the adduct specified in the `adduct` column in the
  suspect list (if present).

- CCSParams:

  A `list` with parameters for mobility \<–\> CCS conversion. See
  [`getCCSParams`](https://rickhelmus.github.io/patRoon/reference/getCCSParams.md)
  for details and to make such parameter lists. Set to `NULL` to skip
  conversions.

- prepareChemProps:

  Set to `TRUE` to perform chemical property calculation and validation
  on the suspect list (described below).

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the suspect
  data (if `prepareChemProps=TRUE`) and `from` data (if a table). For
  efficiency reasons it is recommended to set this to `TRUE`. See the
  `Validating and calculating chemical properties` section for more
  details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

- virtualenv:

  The virtual `Python` environment in which `C3SDB` is installed. This
  is passed to
  [`reticulate::use_virtualenv`](https://rstudio.github.io/reticulate/reference/use_python.html).
  Set to `NULL` to skip this and not setup the environment.

- ...:

  Arguments passed to `data.table` method.

## Details

The `assignMobilities` method for suspect lists is used to (1) add IMS
data to suspects from predictions or library data and (2) convert
(previously added) mobility \<–\> CCS values. These steps are controlled
by the `from` and `CCSParams` arguments, respectively.

Mobility and CCS values assigned in the suspect list are either adduct
specific or not. Adduct specific values are preferred, as the 'correct'
value can be automatically selected during suspect screening based on
the adduct assigned to the feature (or passed as the `adduct` argument
to
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)).
The non-adduct specific values are typically used when the corresponding
adduct for the mobility/CCS value is unknown (or not of interest). These
values get precedence over adduct specific values. The adduct specific
values are stored in `mobility_<adduct>` and `CCS_<adduct>` columns,
where `<adduct>` is the adduct name (*e.g.* `[M+H]+`, `[M-H]-`). The
`mobility` and CCS columns store any non-adduct specific values. The
`adducts` argument ultimately defines the use of adduct and non-adduct
specific values.

The mobility \<–\> CCS conversions occur both ways, *i.e.* missing CCS
values will be converted from mobility values and *vice versa*. If
adduct specific values are converted then the charge value used for
these calculations is taken from the corresponding adduct. For
non-adduct specific values the charge is taken from the adduct specified
in suspect list if present, or from the default charge specified in
`CCSParams` otherwise.

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the suspect
data (if `prepareChemProps=TRUE`) and `from` data (if a table) are
automatically validated and calculated if missing/invalid.

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

Schymanski EL, Kondić T, Neumann S, Thiessen PA, Zhang J, Bolton EE
(2021). “Empowering large chemical knowledge bases for exposomics:
PubChemLite meets MetFrag.” *Journal of Cheminformatics*, **13**(1).
ISSN 1758-2946.
[doi:10.1186/s13321-021-00489-0](https://doi.org/10.1186/s13321-021-00489-0)
. <http://dx.doi.org/10.1186/s13321-021-00489-0>.\
\
Elapavalore A, Ross DH, Grouès V, Aurich D, Krinsky AM, Kim S, Thiessen
PA, Zhang J, Dodds JN, Baker ES, Bolton EE, Xu L, Schymanski EL (2025).
“PubChemLite Plus Collision Cross Section (CCS) Values for Enhanced
Interpretation of Nontarget Environmental Data.” *Environmental Science
& Technology Letters*, **12**(2), 166–174. ISSN 2328-8930.
[doi:10.1021/acs.estlett.4c01003](https://doi.org/10.1021/acs.estlett.4c01003)
. <http://dx.doi.org/10.1021/acs.estlett.4c01003>.\
\
Ross DH, Cho JH, Xu L (2020). “Breaking Down Structural Diversity for
Comprehensive Prediction of Ion-Neutral Collision Cross Sections.”
*Analytical Chemistry*, **92**(6), 4548–4557. ISSN 1520-6882.
[doi:10.1021/acs.analchem.9b05772](https://doi.org/10.1021/acs.analchem.9b05772)
. <http://dx.doi.org/10.1021/acs.analchem.9b05772>.
