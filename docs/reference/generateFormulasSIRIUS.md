# Generate formula with SIRIUS

Uses [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) to
generate chemical formulae candidates.

## Usage

``` r
generateFormulasSIRIUS(fGroups, ...)

# S4 method for class 'featureGroups'
generateFormulasSIRIUS(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  relMzDev = defaultLim("mz", "narrow_rel"),
  adduct = NULL,
  projectPath = NULL,
  elements = "CHNOP",
  profile = "qtof",
  database = NULL,
  noise = NULL,
  cores = NULL,
  getFingerprints = FALSE,
  topMost = 100,
  login = FALSE,
  alwaysLogin = FALSE,
  extraOptsGeneral = NULL,
  extraOptsFormula = NULL,
  calculateFeatures = FALSE,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = defaultLim("mz", "narrow"),
  minIMSSpecSim = 0,
  verbose = TRUE,
  splitBatches = FALSE,
  dryRun = FALSE
)

# S4 method for class 'featureGroupsSet'
generateFormulasSIRIUS(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  relMzDev = defaultLim("mz", "narrow_rel"),
  adduct = NULL,
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
  object for which formulae should be generated. This should be the same
  or a subset of the object that was used to create the specified
  `MSPeakLists`. In the case of a subset only the remaining feature
  groups in the subset are considered.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- MSPeakLists:

  An
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- relMzDev:

  Maximum relative deviation between the measured and candidate formula
  *m/z* values (in ppm). Sets the --ppm-max command line option.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- projectPath, dryRun:

  These are mainly for internal purposes. `projectPath` sets the output
  directory for the `SIRIUS` output (a temporary directory if `NULL`).
  If `dryRun` is `TRUE` then no computations are done and only the
  results from `projectPath` are processed.

  (**sets workflow**) `projectPath` should be a `character` specifying
  the paths for each set.

- elements:

  Elements to be considered for formulae calculation. This will heavily
  affects the number of candidates! Always try to work with a minimal
  set by excluding elements you don't expect. The minimum/maximum number
  of elements can also be specified, for example: a value of
  `"C[5]H[10-15]O"` will only consider formulae with up to five carbon
  atoms, between ten and fifteen hydrogen atoms and any amount of oxygen
  atoms. Sets the --elements command line option.

- profile:

  Name of the configuration profile, for example: "qtof", "orbitrap",
  "fticr". Sets the --profile commandline option.

- database:

  If not `NULL`, use a database for retrieval of formula candidates.
  Possible values are: "pubchem", "bio", "kegg", "hmdb". Sets the
  --database commandline option.

- noise:

  Median intensity of the noise (`NULL` ignores this parameter). Sets
  the --noise commandline option.

- cores:

  The number of cores `SIRIUS` will use. If `NULL` then the default of
  all cores will be used.

- getFingerprints:

  Set to `TRUE` to load `SIRIUS-CSI:FingerID` MS/MS fingerprints for the
  formula candidates. This is currently only supported with
  `calculateFeatures=FALSE` to avoid heavy server traffic. The
  fingerprints are stored in the `fingerprints` slot of the returned
  [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)
  object, and are used by the
  [`predictTox`](https://rickhelmus.github.io/patRoon/reference/generics.md)
  and
  [`predictRespFactors`](https://rickhelmus.github.io/patRoon/reference/generics.md)
  methods.

- topMost:

  Only keep this number of candidates (per feature group) with highest
  score. Sets the --candidates command line option.

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
  website](https://boecker-lab.github.io/docs.sirius.github.io/account-and-license/)
  and patRoon handbook for more information.

- extraOptsGeneral, extraOptsFormula:

  a `character` vector with any extra commandline parameters for
  `SIRIUS`. For `SIRIUS` versions `<4.4` there is no distinction between
  general and formula options. Otherwise commandline options specified
  in `extraOptsGeneral` are added prior to the `formula` command, while
  options specified in `extraOptsFormula` are added in afterwards. See
  the `SIRIUS` manual for more details. Set to `NULL` to ignore.

- calculateFeatures:

  If `TRUE` fomulae are first calculated for all features prior to
  feature group assignment (see `Candidate assignment` in
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)).

- featThreshold:

  If `calculateFeatures=TRUE`: minimum presence (`0-1`) of a formula in
  all features before it is considered as a candidate for a feature
  group. For instance, `featThreshold=0.75` dictates that a formula
  should be present in at least 75% of the features inside a feature
  group.

- featThresholdAnn:

  As `featThreshold`, but only considers features with annotations. For
  instance, `featThresholdAnn=0.75` dictates that a formula should be
  present in at least 75% of the features with annotations inside a
  feature group. @param topMost Only keep this number of candidates (per
  feature group) with highest score. Sets the --candidates command line
  option.

- absAlignMzDev:

  When the group formula annotation consensus is made from feature
  annotations, the *m/z* values of annotated MS/MS fragments may
  slightly deviate from those of the corresponding group MS/MS peak
  list. The `absAlignMzDev` argument specifies the maximum *m/z* window
  used to re-align the mass peaks.

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

- verbose:

  If `TRUE` then more output is shown in the terminal.

- splitBatches:

  If `TRUE` then the calculations done by `SIRIUS` will be evenly split
  over multiple `SIRIUS` calls (which may be run in parallel depending
  on the [set package
  options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)).
  If `splitBatches=FALSE` then all feature calculations are performed
  from a single `SIRIUS` execution, which is often the fastest if
  calculations are performed on a single computer.

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
[`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)
object.

## Details

This function uses sirius to generate formula candidates. This function
is called when calling `generateFormulas` with `algorithm="sirius"`.

Similarity of measured and theoretical isotopic patterns will be used
for scoring candidates. Note that `SIRIUS` requires availability of
MS/MS data.

## Note

For annotations performed with `SIRIUS` it is often the fastest to keep
the default `splitBatches=FALSE`. In this case, all `SIRIUS` output will
be printed to the terminal (unless `verbose=FALSE` or
patRoon.MP.method="future"). Furthermore, please note that only
annotations to be performed for the same adduct are grouped in a single
batch execution.

## Parallelization

generateFormulasSIRIUS uses multiprocessing to parallelize computations.
Please see the parallelization section in the handbook for more details
and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

## References

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

[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
for more details and other algorithms.
