# Compound annotation with SIRIUS

Uses [SIRIUS](https://bright-giant.com/sirius-features/) for compound
annotation.

## Usage

``` r
generateCompoundsSIRIUS(fGroups, ...)

# S4 method for class 'featureGroups'
generateCompoundsSIRIUS(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  adduct = NULL,
  config = NULL,
  topMost = 100,
  login = "check",
  alwaysLogin = FALSE,
  minIMSSpecSim = 0,
  projectPath = NULL,
  runMode = "execute",
  SIRIUSAPI = NULL,
  verbose = TRUE
)

# S4 method for class 'featureGroupsSet'
generateCompoundsSIRIUS(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  adduct = NULL,
  config = NULL,
  login = "check",
  alwaysLogin = FALSE,
  minIMSSpecSim = 0,
  projectPath = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0,
  setAvgSpecificScores = FALSE
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object which should be annotated. This should be the same or a subset
  of the object that was used to create the specified `MSPeakLists`. In
  the case of a subset only the remaining feature groups in the subset
  are considered.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- MSPeakLists:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- config:

  A
  [`RSirius::JobSubmission`](https://rdrr.io/pkg/RSirius/man/JobSubmission.html)
  configuration object, typically obtained with
  [`getSIRIUSConfig`](https://rickhelmus.github.io/patRoon/reference/getSIRIUSConfig.md).
  If `NULL`, the default `SIRIUS` configuration is used.

- topMost:

  Only keep this number of candidates (per feature group) with highest
  score. Setting this to a high number may result in significant usage
  of CPU/RAM resources for large numbers of candidates.

- login, alwaysLogin:

  Specifies if and how account logging of SIRIUS should be handled:

  `login=FALSE`: no automatic login is performed and the active login
  status is not checked.

  `login="check"`: aborts if no active login is present.

  `login="interactive"`: interactively ask for login (using
  [getPass](https://CRAN.R-project.org/package=getPass)).

  `login=c(username="...", password="...")`: perform the login with the
  given details. For security reasons, please do not enter the details
  directly, but use e.g. environment variables or store/retrieve them
  with the [keyring](https://CRAN.R-project.org/package=keyring)
  package.

  if `alwaysLogin=TRUE` then a login is always performed, otherwise only
  if SIRIUS reports no active login.

  See the [SIRIUS
  website](https://v6.docs.sirius-ms.io/account-and-license/) and
  patRoon handbook for more information.

  **NOTE**: By loggin in you will accept the terms of the Service and
  Privacy Policy of the SIRIUS Webservice.

- minIMSSpecSim:

  (**IMS workflow**) If the spectrum similarity of an IMS feature group
  compared to its IMS precursor (see
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
  is at least this value, then the IMS feature group will not be
  subjected to the annotation algorithm and all feature annotation
  properties will be copied from its precursor. This assumes that
  feature annotation is primarily influenced by the MS/MS spectrum, and
  can be used to speed up the feature annotation process. All scorings,
  annotation similarities etc. are copied from the IMS precursor. The
  fragment annotations are also copied (`fragInfo` result column),
  however, these are adjusted based on the peak list data of the IMS
  feature group.

- runMode, projectPath:

  Whether to execute a `SIRIUS` processing job (`runMode="execute"`) or
  load results from an existing `SIRIUS` project (`runMode"read"`). If
  `runMode="execute"` then `projectPath` can be `NULL` and a temporary
  project will be used, otherwise `projectPath` must point to an
  existing project.

  **NOTE:** if `runMode="execute"` then any existing project at
  `projectPath` will be removed.

  **NOTE:** This is primarily intended for internal purposes, but may be
  of interest to e.g. re-import SIRIUS results.

  (**sets workflow**) `projectPath` should be a `character` specifying
  the paths for each set.

- SIRIUSAPI:

  An `rsirius_api` object for connecting to the `SIRIUS` API. If `NULL`,
  a new connection will be started automatically.

- verbose:

  If `TRUE` then more output is shown in the terminal.

- setThreshold:

  (**sets workflow**) Minimum abundance for a candidate among all sets
  (`0-1`). For instance, a value of `1` means that the candidate needs
  to be present in all the set data.

- setThresholdAnn:

  (**sets workflow**) As `setThreshold`, but only taking into account
  the set data that contain annotations for the feature group of the
  candidate.

- setAvgSpecificScores:

  (**sets workflow**) If `TRUE` then set specific scorings (*e.g.* MS/MS
  match) are also averaged.

## Value

A
[`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compoundsSIRIUS-class.md)
object.

## Details

This function uses SIRIUS to generate compound candidates. This function
is called when calling `generateCompounds` with `algorithm="sirius"`.

Similar to
[`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md),
candidate formulae are generated with SIRIUS. These results are then fed
to `CSI:FingerID` to acquire candidate structures. Candidate formulae
without any assigned structure will be removed (unlike
[`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)).
This method requires the availability of MS/MS data, and feature groups
without will be ignored.

## Running SIRIUS

By default, patRoon tries to connect to a running instance of `SIRIUS`.
This is generally faster and may be useful for debugging by *e.g.*
checking the logs in `SIRIUS`. Otherwise, an attempt will be made to
start `SIRIUS` automatically. The binaries are searched from the
patRoon.path.SIRIUS package option, patRoonExt package or the system
`PATH` environment variable. Any automatically started `SIRIUS`
instances are automatically closed if jobs are finished. By default, a
temporary `SIRIUS` project is made for `SIRIUS` data processing and
removed afterwards. See the `projectPath` to change this.

## SIRIUS 6 functionality

The interface to `SIRIUS 6` is still in development and may be extended
in the future. There is a vast amount of functionality available, which
will require quite some effort to support all. However, the current
functionality in patRoon is mostly equal to what was supported with
previous `SIRIUS` releases. Any feedback on the inclusion of specific
functionality is welcome!

## References

Hoffmann MA, Nothias L, Ludwig M, Fleischauer M, Gentry EC, Witting M,
Dorrestein PC, Dührkop K, Böcker S (2021). “High-confidence structural
annotation of metabolites absent from spectral libraries.” *Nature
Biotechnology*, **40**(3), 411–421. ISSN 1546-1696.
[doi:10.1038/s41587-021-01045-9](https://doi.org/10.1038/s41587-021-01045-9)
. <http://dx.doi.org/10.1038/s41587-021-01045-9>.\
\
Dührkop K, Nothias L, Fleischauer M, Reher R, Ludwig M, Hoffmann MA,
Petras D, Gerwick WH, Rousu J, Dorrestein PC, Böcker S (2020).
“Systematic classification of unknown metabolites using high-resolution
fragmentation mass spectra.” *Nature Biotechnology*, **39**(4), 462–471.
ISSN 1546-1696.
[doi:10.1038/s41587-020-0740-8](https://doi.org/10.1038/s41587-020-0740-8)
. <http://dx.doi.org/10.1038/s41587-020-0740-8>.\
\
Duhrkop K, Fleischauer M, Ludwig M, Aksenov AA, Melnik AV, Meusel M,
Dorrestein PC, Rousu J, Bocker S (2019). “SIRIUS 4: a rapid tool for
turning tandem mass spectra into metabolite structure information.”
*Nature Methods*, **16**(4), 299–302.
[doi:10.1038/s41592-019-0344-8](https://doi.org/10.1038/s41592-019-0344-8)
.\
\
Duhrkop K, Bocker S (2015). “Fragmentation Trees Reloaded.” In Przytycka
TM (ed.), *Research in Computational Molecular Biology*, 65–79. ISBN
978-3-319-16706-0.\
\
Duhrkop K, Shen H, Meusel M, Rousu J, Bocker S (2015). “Searching
molecular structure databases with tandem mass spectra using
CSI:FingerID.” *Proceedings of the National Academy of Sciences*,
**112**(41), 12580–12585.
[doi:10.1073/pnas.1509788112](https://doi.org/10.1073/pnas.1509788112)
.\
\
Bocker S, Letzel MC, Liptak Z, Pervukhin A (2008). “SIRIUS: decomposing
isotope patterns for metabolite identification.” *Bioinformatics*,
**25**(2), 218–224.
[doi:10.1093/bioinformatics/btn603](https://doi.org/10.1093/bioinformatics/btn603)
.

## See also

[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
for more details and other algorithms.
