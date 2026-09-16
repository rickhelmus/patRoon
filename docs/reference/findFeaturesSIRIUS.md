# Find features using SIRIUS

Uses [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/) to
find features.

## Usage

``` r
findFeaturesSIRIUS(
  analysisInfo,
  noiseIntensity = NULL,
  alignMaxRTDev = NULL,
  minSNR = NULL,
  login = "check",
  alwaysLogin = FALSE,
  projectPath = NULL,
  runMode = "execute",
  SIRIUSAPI = NULL,
  verbose = TRUE
)

importFeaturesSIRIUS(input, analysisInfo, ...)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- noiseIntensity, alignMaxRTDev, minSNR:

  Parameters for the SIRIUS feature finding algorithm. See the SIRIUS
  documentation for details. Set to `NULL` to use the default values.

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

  If set to `FALSE` then no text output is shown.

- input:

  Sets `projectPath`.

- ...:

  Additional arguments passed to `findFeaturesSIRIUS`.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

This function uses SIRIUS to automatically find features. This function
is called when calling `findFeatures` with `algorithm="sirius"`.

The MS files should be in the `mzML` or `mzXML` format.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

`importFeaturesSIRIUS` is a simple wrapper around `findFeaturesSIRIUS`
to import features from an existing SIRIUS project. It will set
`runMode="read"` and `projectPath` to the provided `input` path.

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

[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
for more details and other algorithms.
