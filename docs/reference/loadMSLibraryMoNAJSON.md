# Load MS library data from MassBank of North America (MoNA)

This function loads, verifies and curates MS library data from
[MoNA](https://mona.fiehnlab.ucdavis.edu/) `.json` files.

## Usage

``` r
loadMSLibraryMoNAJSON(
  file,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  potAdducts = TRUE,
  potAdductsLib = TRUE,
  absMzDev = defaultLim("mz", "narrow"),
  calcSPLASH = TRUE
)
```

## Source

Guessing adducts from neutral/ionic mass differences was inspired from
[MetFrag](http://ipb-halle.github.io/MetFrag/).

## Arguments

- file:

  A `character` string that specifies the file path to the JSON library.

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the MS library.
  For efficiency reasons it is recommended to set this to `TRUE`. See
  the `Validating and calculating chemical properties` section for more
  details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

- potAdducts, potAdductsLib:

  If and how missing adducts (`Precursor_type` data) are guessed,
  `potAdducts` should be either:

  - `FALSE`: do not perform adduct guessing.

  - `TRUE`: guesses adducts based on a common set of known adducts
    (currently based on
    [`GenFormAdducts`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)
    and
    [`MetFragAdducts`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
    If `potAdductsLib` is `TRUE` then also any adducts specified in the
    library are used.

  - A `list` with
    [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
    objects or `character` vector that can be converted with
    [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md).
    Only the specified adducts will be used for guessing missing values.

- absMzDev:

  The maximum absolute *m/z* deviation when guessing missing adducts.

- calcSPLASH:

  If set to `TRUE` then missing SPLASH values will be calculated (see
  below).

## Value

The loaded data is returned in an
[`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
object.

## Details

This function uses an efficient `C++` JSON loader to load MS library
data. This function is called when calling `loadMSLibrary` with
`algorithm="json"`.

This function uses `C++` with
[Rcpp](https://CRAN.R-project.org/package=Rcpp) and
[rapidjsonr](https://CRAN.R-project.org/package=rapidjsonr) to
efficiently load and parse JSON files from
[MoNA](https://mona.fiehnlab.ucdavis.edu/). An advantage compared to
[`loadMSLibraryMSP`](https://rickhelmus.github.io/patRoon/reference/loadMSLibraryMSP.md)
is that this function supports loading spectral annotations.

The record field names are converted to those used in `.msp` files.

## Automatic curation of library data

Several strategies are applied to automatically verify and improve
library data. This is important, since library records may have
inconsistent or erroneous data, which makes them unsuitable in automated
workflows such as compounds annotation with
[`generateCompoundsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsLibrary.md).

The loaded library data is post-treated as follows:

- The `DB#` field is renamed to `DB_ID` to improve compatibility with R
  column names.

- Synonyms (`Synon` fields) are merged together, mainly to save memory
  usage.

- Inconsistently formatted `NA` data (*e.g.* `"n/a"`, `"N/A"` or empty
  strings) are set to regular R `NA` values.

- The case of record field names are made consistent.

- The `Formula` and `ExactMass` fields are renamed to `formula` and
  `neutralMass`, respectively. This is for consistency with other data
  generated with patRoon.

- `character` field data is trimmed from leading/trailing whitespace.

- Mass data is verified to be properly numeric, and set to `NA`
  otherwise.

- The format of formulae data is made consistent: ionic species (with or
  without square brackets) or converted to a regular formula format.

- Chemical identifiers such as SMILES and formulae are verified and
  missing values are calculated if possible. See below for more details.

- Shortened data in the `Ion_mode` field (P/N) is converted to the long
  format (`POSITIVE`/`NEGATIVE`).

- Many different adduct flavors typically found as `Precursor_type` data
  are converted and normalized to the generic textual format used by
  patRoon (see
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).

- If `potAdducts!=FALSE` then missing or invalid adduct data in
  `Precursor_type` is guessed based on the difference between the
  neutral and ionic mass. If multiple adducts explain the mass
  difference the result is `NA`.

- Missing ion *m/z* data (`PrecursorMZ` field) is calculated from adduct
  data, if possible.

- Missing [SPLASH](https://splash.fiehnlab.ucdavis.edu/) data is
  calculated with the splashR package if `calcSPLASH=TRUE`.

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the MS
library are automatically validated and calculated if missing/invalid.

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

Wohlgemuth G, Mehta SS, Mejia RF, Neumann S, Pedrosa D, Pluskal T,
Schymanski EL, Willighagen EL, Wilson M, Wishart DS, Arita M, Dorrestein
PC, Bandeira N, Wang M, Schulze T, Salek RM, Steinbeck C, Nainala VC,
Mistrik R, Nishioka T, Fiehn O (2016). “SPLASH, a hashed identifier for
mass spectra.” *Nature Biotechnology*, **34**(11), 1099–1101.
[doi:10.1038/nbt.3689](https://doi.org/10.1038/nbt.3689) .\
\
Ruttkies C, Schymanski EL, Wolf S, Hollender J, Neumann S (2016).
“MetFrag relaunched: incorporating strategies beyond in silico
fragmentation.” *Journal of Cheminformatics*, **8**(1).
[doi:10.1186/s13321-016-0115-9](https://doi.org/10.1186/s13321-016-0115-9)
.\
\
Eddelbuettel D (2013). *Seamless R and C++ Integration with Rcpp*.
Springer, New York.
[doi:10.1007/978-1-4614-6868-4](https://doi.org/10.1007/978-1-4614-6868-4)
. ISBN 978-1-4614-6867-7.\
\
Eddelbuettel D, Balamuta J (2018). “Extending R with C++: A Brief
Introduction to Rcpp.” *The American Statistician*, **72**(1), 28-36.
[doi:10.1080/00031305.2017.1375990](https://doi.org/10.1080/00031305.2017.1375990)
.\
\
Eddelbuettel D, François R (2011). “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software*, **40**(8), 1–18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08) .\
\
Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar
I, Bates D, Chambers J (2026). *Rcpp: Seamless R and C++ Integration*. R
package version 1.1.2, <https://www.rcpp.org>.

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .

## See also

[`loadMSLibrary`](https://rickhelmus.github.io/patRoon/reference/loadMSLibrary.md)
for more details and other algorithms.

The
[`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
documentation for various methods to post-process the data and
[`generateCompoundsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsLibrary.md)
for annotation of features with the library data.
