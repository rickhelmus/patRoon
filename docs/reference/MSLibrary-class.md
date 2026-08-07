# Class to store data from a loaded MS library

Stores the spectra and metadata from the records of an MS library.

## Usage

``` r
records(obj)

spectra(obj)

# S4 method for class 'MSLibrary'
records(obj)

# S4 method for class 'MSLibrary'
spectra(obj)

# S4 method for class 'MSLibrary'
length(x)

# S4 method for class 'MSLibrary'
names(x)

# S4 method for class 'MSLibrary'
show(object)

# S4 method for class 'MSLibrary,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'MSLibrary,ANY,missing'
x[[i, j]]

# S4 method for class 'MSLibrary'
x$name

# S4 method for class 'MSLibrary'
as.data.table(x)

# S4 method for class 'MSLibrary'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'MSLibrary'
filter(
  obj,
  properties = NULL,
  massRange = NULL,
  mzRangeSpec = NULL,
  relMinIntensity = NULL,
  topMost = NULL,
  onlyAnnotated = FALSE,
  negate = FALSE
)

# S4 method for class 'MSLibrary'
convertToSuspects(
  obj,
  adduct,
  spectrumType = "MS2",
  avgSpecParams = getDefAvgPListParams(minIntensityPre = 0, minIntensityPost = 2, topMost
    = 10),
  collapse = TRUE,
  suspects = NULL,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE
)

# S4 method for class 'MSLibrary'
export(obj, type = "msp", out)

# S4 method for class 'MSLibrary,MSLibrary'
merge(x, y, ...)
```

## Arguments

- x, obj, object:

  `MSLibrary` object to be accessed.

- i:

  For `[`/`[[`: A numeric or character value which is used to select
  records by their index or name, respectively (for the order/names see
  [`names()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all records are selected.\
  \
  For `[[`: should be a scalar value.

- ...:

  Unused.

- drop, j:

  ignored.

- name:

  The record name (partially matched).

- properties:

  A named `list` with properties to be filtered. Each item in the `list`
  should be named with the name of the property, and should be a vector
  with allowed values. To obtain the possible properties, run *e.g.*
  `names(records)`. Example:
  `properties=list(Instrument_type=c("LC-ESI-QTOF","LC-ESI-TOF"))`. Set
  to `NULL` to ignore.

- massRange:

  Records with a neutral mass outside this range will be removed. Should
  be a two-sized `numeric` vector with the lower and upper mass range.
  Set to `NULL` to ignore.

- mzRangeSpec:

  Similar to the `massRange` argument, but removes any peaks from
  recorded mass spectra outside the given *m/z* range.

- relMinIntensity:

  The minimum relative intensity (`0-1`) of a mass peak to be kept. Set
  to `NULL` to ignore.

- topMost:

  Only keep `topMost` number of mass peaks for each spectrum. This
  filter is applied after others. Set to `NULL` to ignore.

- onlyAnnotated:

  If `TRUE` then only recorded spectra that are formula annotated are
  kept.

- negate:

  If `TRUE` then filters are performed in opposite manner.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Any records with a different adduct (`Precursor_type`) are not
  considered. Alternatively, `adduct` can be set to `NULL` to not filter
  out any records. However, in this case *no* MS/MS fragments will be
  added to the returned suspect list.

- spectrumType:

  A `character` vector which limits library records to the given
  spectrum types (`Spectrum_type` field, *e.g.* `"MS2"`). Set to `NULL`
  to allow all spectrum types.

- avgSpecParams:

  A `list` with parameters used for averaging spectra. See
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for more details.

- collapse:

  Whether records with the same first-block InChIKey should be
  collapsed. See the `Suspect conversion` section for details.

- suspects:

  If not `NULL` then this should be a suspect list (see
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md))
  which will be amended with spectra data. See the `Suspect conversion`
  section for details.

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the input
  suspect list to `convertToSuspects`. For efficiency reasons it is
  recommended to set this to `TRUE`. See the
  `Validating and calculating chemical properties` section for more
  details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

- type:

  The export type. Currently just `"msp"`.

- out:

  The file path to the output library file.

- y:

  The `MSLibrary` to be merged with `x`.

## Value

`delete` returns the object for which the specified data was removed.

`filter` returns a filtered `MSLibrary` object.

`convertToSuspects` return a suspect list (`data.table`), which can be
used with
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).

`merge` returns a merged `MSLibrary` object.

## Details

This class is used by
[`loadMSLibrary`](https://rickhelmus.github.io/patRoon/reference/loadMSLibrary.md)
to store the loaded MS library data.

## Methods (by generic)

- `records(MSLibrary)`: Accessor method for the `records` slot of an
  `MSLibrary` class.

- `spectra(MSLibrary)`: Accessor method for the `spectra` slot of an
  `MSLibrary` class.

- `length(MSLibrary)`: Obtains the total number of records stored.

- `names(MSLibrary)`: Obtains the names of the stored records (`DB_ID`
  field).

- `show(MSLibrary)`: Shows summary information for this object.

- `x[i`: Subset on records.

- `x[[i`: Extracts a spectrum table for a record.

- `$`: Extracts a spectrum table for a record.

- `as.data.table(MSLibrary)`: Converts all the data (spectra and
  metadata) to a single `data.table`.

- `delete(MSLibrary)`: Completely deletes specified full records or
  spectra.

- `filter(MSLibrary)`: Performs rule-based filtering of records and
  spectra. This may be especially to improve annotation with
  [`generateCompoundsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsLibrary.md).

- `convertToSuspects(MSLibrary)`: Converts the MS library data to a
  suspect list, which can be used with
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
  See the `Suspect conversion` section for details.

- `export(MSLibrary)`: Exports the library data to a `.msp` file. The
  export is accelerated by an `C++` interface with
  [Rcpp](https://CRAN.R-project.org/package=Rcpp).

- `merge(x = MSLibrary, y = MSLibrary)`: Merges two `MSLibrary` objects
  (`x` and `y`). The records from `y` that are unique are added to `x`.
  Records that were already in `x` are simply ignored. The
  [SPLASH](https://splash.fiehnlab.ucdavis.edu/) values are used to test
  equality between records, hence, the `calcSPLASH` argument to
  [`loadMSLibrary`](https://rickhelmus.github.io/patRoon/reference/loadMSLibrary.md)
  should be `TRUE`.

## Slots

- `records`:

  A `data.table` with metadata for all records. Use the `records` method
  for access.

- `spectra`:

  A `list` with all (annotated) spectra. Each spectrum is stored in a
  `data.table`. Use the `spectra` method for access.

## Note

`export` does not split any `Synon` data that was merged when the
library was loaded.

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`MSLibrary`**

## Suspect conversion

The `convertToSuspects` method converts MS library data to a suspect
list, which can be used with *e.g.*
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
Furthermore, this function can also amend existing suspect lists with
spectral data.

Conversion occurs in either of the following three methods:

1.  *Direct* (`collapse=FALSE` and `suspects=NULL`): each record is
    considered a suspect, and the resulting suspect list is generated
    directly by converting the records metadata. The `fragments_mz`
    column for each suspect is constructed from the mass peaks of the
    corresponding record.

2.  *Collapse* (`collapse=TRUE` and `suspects=NULL`): All records with
    the same first-block InChIKey are first merged, and their spectra
    are averaged using the parameters from the `avgSpecParams` argument
    (see
    [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)).
    The suspect list is based on the merged records, where the
    `fragments_mz` column is constructed from the averaged spectra. This
    is generally a good default, especially with large MS libraries.

3.  *Amend* (`suspects` is not `NULL`): only those records are
    considered if their first-block InChIKey is present in the suspect
    list. The remaining records and their spectra are then collapsed as
    described for the *Collapse* method, and the `fragments_mz` column
    for each suspect is set from the averaged spectra. If a suspect is
    not present in the library, its `fragments_mz` value will be empty.
    Note that any existing `fragments_mz` data will be overwritten.

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the input
suspect list to `convertToSuspects` are automatically validated and
calculated if missing/invalid.

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
