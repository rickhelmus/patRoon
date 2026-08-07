# MS data conversion

Conversion of MS analysis files between several open and closed data
formats.

## Usage

``` r
getMSConversionTypes(algorithm, direction)

getMSConversionFormats(algorithm, direction, type = NULL)

convertMSFilesPWiz(
  inFiles,
  outFiles,
  formatTo = "mzML",
  centroid = TRUE,
  IMS = FALSE,
  minIntensity = 0,
  filters = NULL,
  extraOpts = NULL,
  PWizBatchSize = 1
)

convertMSFilesOpenMS(inFiles, outFiles, formatTo = "mzML", extraOpts = NULL)

convertMSFilesBruker(inFiles, outFiles, formatTo = "mzML", centroid = TRUE)

convertMSFilesIMSCollapse(
  inFiles,
  outFiles,
  typeFrom,
  formatTo = "mzML",
  mzRange = NULL,
  mobilityRange = NULL,
  smoothWindow = 0,
  halfWindow = 2,
  maxGap = 0.005,
  clusterMethod = "distance_mean",
  mzWindow = defaultLim("mz", "medium"),
  minIntensityIMS = 0,
  includeMSMS = FALSE,
  ...
)

convertMSFilesTIMSCONVERT(
  inFiles,
  outFiles,
  formatTo = "mzML",
  centroid = TRUE,
  centroidRaw = FALSE,
  IMS = FALSE,
  extraOpts = NULL,
  virtualenv = "patRoon-TIMSCONVERT"
)

convertMSFilesPaths(
  files,
  formatFrom,
  formatTo = "mzML",
  outPath = NULL,
  dirs = TRUE,
  overwrite = FALSE,
  algorithm = "pwiz",
  ...
)

convertMSFiles(
  anaInfo,
  typeFrom = "raw",
  typeTo = "centroid",
  formatFrom,
  formatTo = "mzML",
  overwrite = FALSE,
  algorithm = "pwiz",
  centroidVendor = TRUE,
  ...
)
```

## Arguments

- algorithm:

  Either `"pwiz"` (ProteoWizard), `"openms"`, `"bruker"` (Bruker
  DataAnalysis) , `"imscollapse"` or `"timsconvert"`.

- direction:

  A `character` specifying the direction of conversion. Either `"input"`
  or `"output"`.

- type, typeFrom, typeTo:

  The type of the input or output files. See `getMSConversionTypes` for
  the supported types.

- inFiles, outFiles:

  A `character` vector with input and output files, respectively.
  Lengths and order should be the same.

- centroid:

  Set to `TRUE` to perform centroiding.

  For `convertMSFilesPWiz`: the value may be `"vendor"` to perform
  centroiding with the vendor algorithm or `"cwt"` to use ProteoWizard's
  wavelet algorithm.

- IMS:

  How to handle IMS data.

  For `convertMSFilesPWiz`: if `TRUE` then IMS data is exported and
  spectra for each IMS frame are combined into a single spectrum (using
  the `–combineIonMobilitySpectra` option), which is the format
  supported by patRoon. Set to `NA` to collapse the IMS data by scan
  summing, which mimics 'regular' HRMS data. Set to `FALSE` for non-IMS
  data. **NOTE**: do not set `IMS=FALSE` if the data has IMS data. This
  will result in very large files where MS spectra are not combined by
  frame, which **cannot** be properly read by patRoon.

  For `convertMSFilesTIMSCONVERT`: set to `TRUE` to keep IMS data or
  `FALSE` to exclude IMS data to mimic 'regular' LC-MS data.

- minIntensity:

  The minimum intensity of the mass peaks to be kept. Applying an
  intensity threshold is especially beneficial to reduce export file
  size when there are a lot of zero or very low intensity mass peaks.
  **NOTE** this currently does *not* work well with IMS data.

- filters:

  A `character` vector specifying one or more filters to `msconvert`.
  The elements of the specified vector are directly passed to the
  `--filter` option (see
  [here](http://proteowizard.sourceforge.net/tools/filters.md))

- extraOpts:

  A `character` vector specifying any extra command line parameters
  passed to `msconvert` or `FileConverter`. Set to `NULL` to ignore. For
  options: see
  [FileConverter](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FileConverter.html)
  and
  [msconvert](http://proteowizard.sourceforge.net/tools/msconvert.md).

- PWizBatchSize:

  The number of analyses to process by a single call to `msconvert`.
  Usually a value of one is most efficient. Set to zero to run all
  analyses all at once from a single call.

- mzRange, mobilityRange:

  A two sized vector specifying the m/z and mobility range to be
  exported, respectively. Set to `NULL` to export the full range.

- smoothWindow, halfWindow, maxGap:

  Centroiding parameters: see
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for details. **NOTE**: As described there, `maxGap` may need to be
  increased for Agilent instruments (*e.g.* `0.01`).

- clusterMethod, mzWindow:

  The clustering method and window (see [clustering
  parameters](https://rickhelmus.github.io/patRoon/reference/cluster-params.md))
  used to find and combine MS/MS spectra of precursors with close *m/z*.

- minIntensityIMS:

  The minimum intensity for MS peaks in raw data.

- includeMSMS:

  Set to `TRUE` to include MS/MS spectra in the output. For IMS
  workflows where IMS data is only collapsed to produce compatible data
  files for feature detection, MS/MS data are not needed and can be
  excluded to reduce computational times and file sizes. Setting
  `includeMSMS=TRUE` is primarily intended to perform 'classical LC-MS
  workflows' with IMS data.

- ...:

  For `convertMSFilesIMSCollapse`: further arguments passed to
  [`mzR::writeMSData`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

  For `convertMSFilesPaths` and `convertMSFiles`: further arguments
  passed to algorithm specific conversion functions.

- centroidRaw:

  Only applicable if `IMS=FALSE`. Sets the `mode` parameter of
  `TIMSCONVERT`: `raw` if `centroidRaw=TRUE` or `centroid` if
  `centroidRaw=FALSE`. See
  <https://gtluu.github.io/timsconvert/local.html#notes-on-mode-parameter>
  for more details.

- virtualenv:

  The virtual Python environment in which `TIMSCONVERT` is installed.
  This is passed to
  [`reticulate::use_virtualenv`](https://rstudio.github.io/reticulate/reference/use_python.html),
  which will ensure that the `TIMSCONVERT` command line utility can be
  found by patRoon. Set to `NULL` to skip this step.

- files, dirs:

  The `files` argument should be a `character` vector with input files.
  If `files` contains directories and `dirs=TRUE` then files from these
  directories are also considered.

- formatFrom, formatTo:

  The input or output format. See `getMSConversionFormats` for the
  supported formats.

- outPath:

  A character vector specifying directories that should be used for the
  output. Will be re-cycled if necessary. If `NULL`, output directories
  will be kept the same as the input directories.

- overwrite:

  Should existing destination file be overwritten (`TRUE`) or not
  (`FALSE`)?

- anaInfo:

  An [analysis info
  table](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  that is used to retrieve the input files. The paths set by
  `path_centroid`, `path_profile` and `path_ims` are used to determine
  the output directories. This function automatically determines if and
  how centroiding and IMS conversions should be applied.

- centroidVendor:

  Only for `algorithm="pwiz"`: whether centroiding should be performed
  with vendor algorithms.

## Details

`getMSConversionTypes` returns a `character` with all supported input or
output conversion types for an algorithm.

`getMSConversionFormats` returns a `character` with all supported input
or output conversion formats for an algorithm, optionally filtered by
the given `type`.

`convertMSFilesPWiz` converts and pre-treats HRMS data with the
`msconvert` tool from
[ProteoWizard](http://proteowizard.sourceforge.net/).

`convertMSFilesOpenMS` converts HRMS data with the `FileConvert` tool of
[OpenMS](http://www.openms.de/).

`convertMSFilesBruker` converts and pre-treats Bruker HRMS data with
Bruker DataAnalysis. Note that TIMS data currently is not supported.

`convertMSFilesIMSCollapse` is used to convert IMS data to data that
mimics 'regular' HRMS data by collapsing the IMS dimension. The raw data
interface of patRoon first sums up all spectra within each IMS frame,
performs centroiding and finally exports the resulting data with the
[`mzR::writeMSData`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function. Several thresholds can be set to speed up the conversion
process and reduce noise, but care should be taken that no mass peaks of
interest are lost.

`convertMSFilesTIMSCONVERT` converts and pre-treats TIMS data with
[TIMSCONVERT](https://gtluu.github.io/timsconvert/). The
[`installTIMSCONVERT`](https://rickhelmus.github.io/patRoon/reference/installTIMSCONVERT.md)
function can be used to automatically install `TIMSCONVERT`.

`convertMSFilesPaths` is a wrapper function that simplifies the use of
algorithm specific MS conversion functions, such as
`convertMSFilesPWiz`, and `convertMSFilesTIMSCONVERT`.

`convertMSFiles` is a wrapper function that simplifies the use of
`convertMSFilesPaths`.

## Parallelization

`convertMSFilesPWiz`, `convertMSFilesOpenMS` and
`convertMSFilesTIMSCONVERT` uses multiprocessing to parallelize
computations. Please see the parallelization section in the handbook for
more details and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `convertMSFilesIMSCollapse` to process HRMS (or
IMS-HRMS) data. Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## References

Rost HL, Sachsenberg T, Aiche S, Bielow C, Weisser H, Aicheler F,
Andreotti S, Ehrlich H, Gutenbrunner P, Kenar E, Liang X, Nahnsen S,
Nilse L, Pfeuffer J, Rosenberger G, Rurik M, Schmitt U, Veit J, Walzer
M, Wojnar D, Wolski WE, Schilling O, Choudhary JS, Malmstrom L,
Aebersold R, Reinert K, Kohlbacher O (2016). “OpenMS: a flexible
open-source software platform for mass spectrometry data analysis.”
*Nature Methods*, **13**(9), 741–748.
[doi:10.1038/nmeth.3959](https://doi.org/10.1038/nmeth.3959) .\
\

Chambers MC, Maclean B, Burke R, Amodei D, Ruderman DL, Neumann S, Gatto
L, Fischer B, Pratt B, Egertson J, Hoff K, Kessner D, Tasman N, Shulman
N, Frewen B, Baker TA, Brusniak M, Paulse C, Creasy D, Flashner L, Kani
K, Moulding C, Seymour SL, Nuwaysir LM, Lefebvre B, Kuhlmann F, Roark J,
Rainer P, Detlev S, Hemenway T, Huhmer A, Langridge J, Connolly B,
Chadick T, Holly K, Eckels J, Deutsch EW, Moritz RL, Katz JE, Agus DB,
MacCoss M, Tabb DL, Mallick P (2012). “A cross-platform toolkit for mass
spectrometry and proteomics.” *Nature Biotechnology*, **30**(10),
918–920. [doi:10.1038/nbt.2377](https://doi.org/10.1038/nbt.2377) .\
\

Luu GT, Freitas MA, Lizama-Chamu I, McCaughey CS, Sanchez LM, Wang M
(2022). “TIMSCONVERT: a workflow to convert trapped ion mobility data to
open data formats.” *Bioinformatics*, **38**(16), 4046–4047. ISSN
1367-4811.
[doi:10.1093/bioinformatics/btac419](https://doi.org/10.1093/bioinformatics/btac419)
. <http://dx.doi.org/10.1093/bioinformatics/btac419>.\
\

Chambers, C. M, Maclean, Brendan, Burke, Robert, Amodei, Dario,
Ruderman, L. D, Neumann, Steffen, Gatto, Laurent, Fischer, Bernd, Pratt,
Brian, Egertson, Jarrett, Hoff, Katherine, Kessner, Darren, Tasman,
Natalie, Shulman, Nicholas, Frewen, Barbara, Baker, A. T, Brusniak,
Mi-Youn, Paulse, Christopher, Creasy, David, Flashner, Lisa, Kani, Kian,
Moulding, Chris, Seymour, L. S, Nuwaysir, M. L, Lefebvre, Brent,
Kuhlmann, Frank, Roark, Joe, Rainer, Paape, Detlev, Suckau, Hemenway,
Tina, Huhmer, Andreas, Langridge, James, Connolly, Brian, Chadick, Trey,
Holly, Krisztina, Eckels, Josh, Deutsch, W. E, Moritz, L. R, Katz, E. J,
Agus, B. D, MacCoss, Michael, Tabb, L. D, Mallick, Parag (2012). “A
cross-platform toolkit for mass spectrometry and proteomics.” *Nat
Biotech*, **30**(10), 918–920.
[doi:10.1038/nbt.2377](https://doi.org/10.1038/nbt.2377) .
<http://dx.doi.org/10.1038/nbt.2377>.\
\
Keller A, Eng J, Zhang N, Li X, Aebersold R (2005). “A uniform
proteomics MS/MS analysis platform utilizing open XML file formats.”
*Mol Syst Biol*.\
\
Kessner D, Chambers M, Burke R, Agus D, Mallick P (2008). “ProteoWizard:
open source software for rapid proteomics tools development.”
*Bioinformatics*, **24**(21), 2534–2536.
[doi:10.1093/bioinformatics/btn323](https://doi.org/10.1093/bioinformatics/btn323)
.\
\
Martens L, Chambers M, Sturm M, Kessner D, Levander F, Shofstahl J, Tang
WH, Rompp A, Neumann S, Pizarro AD, Montecchi-Palazzi L, Tasman N,
Coleman M, Reisinger F, Souda P, Hermjakob H, Binz P, Deutsch EW (2010).
“mzML - a Community Standard for Mass Spectrometry Data.” *Mol Cell
Proteomics*.
[doi:10.1074/mcp.R110.000133](https://doi.org/10.1074/mcp.R110.000133)
.\
\
Pedrioli PGA, Eng JK, Hubley R, Vogelzang M, Deutsch EW, Raught B, Pratt
B, Nilsson E, Angeletti RH, Apweiler R, Cheung K, Costello CE, Hermjakob
H, Huang S, Julian RK, Kapp E, McComb ME, Oliver SG, Omenn G, Paton NW,
Simpson R, Smith R, Taylor CF, Zhu W, Aebersold R (2004). “A common open
representation of mass spectrometry data and its application to
proteomics research.” *Nat Biotechnol*, **22**(11), 1459–1466.
[doi:10.1038/nbt1031](https://doi.org/10.1038/nbt1031) .
