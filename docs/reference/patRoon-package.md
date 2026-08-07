# Workflow solutions for mass-spectrometry based non-target analysis.

Provides an easy-to-use interface to a mass spectrometry based
non-target analysis workflow. Various (open-source) tools are combined
which provide algorithms for extraction and grouping of features,
extraction of MS and MS/MS data, automatic formula and compound
annotation and grouping related features to components. In addition,
various tools are provided for e.g. data preparation and cleanup,
plotting results and automatic reporting.

## Package options

The following package options (see
[`options`](https://rdrr.io/r/base/options.html)) can be set:

- `patRoon.cache.mode`: A `character` setting the current caching mode:
  `"save"` and `"load"` will only save/load results to/from the cache,
  `"both"` (default) will do both and `"none"` to completely disable
  caching. This option can be changed anytime, which might be useful,
  for instance, to temporarily disable cached results before running a
  function.

- `patRoon.cache.fileName`: a `character` specifying the name of the
  cache file (default is `cache.sqlite`).

- `patRoon.cache.maxEntries`: a `numeric` specifying the maximum number
  of entries per cache category (default is 100000). When this limit is
  exceeded, the oldest entries are automatically removed.

- `patRoon.MS.backends`,`patRoon.MS.preferIMS`,`patRoon.path.TDFSDK`:
  Options related to the [raw data
  interface](https://rickhelmus.github.io/patRoon/reference/msdata.md).

- `patRoon.threads`: The number of threads to be used for
  parallelization. This is currently only used by the [raw data
  interface](https://rickhelmus.github.io/patRoon/reference/msdata.md)
  and when the `piek` algorithm is used for [peak
  detection](https://rickhelmus.github.io/patRoon/reference/getDefPeakParams.md).

- `patRoon.MP.maxProcs`: The maximum number of processes that should be
  initiated in parallel. A good starting point is the number of physical
  cores, which is the default as detected by
  [`detectCores`](https://rdrr.io/r/parallel/detectCores.html). This
  option is only used when patRoon.MP.method="classic".

- `patRoon.MP.method`: Either `"classic"` or `"future"`. The former is
  the default and uses
  [processx](https://CRAN.R-project.org/package=processx) to execute
  multiple commands in parallel. When `"future"` the `future.apply`
  package is used for parallelization, which is especially useful for
  *e.g.* cluster computing.

- `patRoon.MP.futureSched`: Sets the `future.scheduling` function
  argument for `future_lapply`. Only used if patRoon.MP.method="future".

- `patRoon.MP.logPath`: The path used for logging of output from
  commands executed by multiprocess. Set to `FALSE` to disable logging.

- `patRoon.path.pwiz`: The path in which the `ProteoWizard` binaries are
  installed. If unset an attempt is made to find this directory from the
  Windows registry and PATH environment variable.

- `patRoon.path.GenForm`: The path to the `GenForm` executable. If not
  set (the default) the internal `GenForm` binary is used. Only set if
  you want to override the executable.

- `patRoon.path.MetFragCL`: The complete file path to the MetFrag CL
  `jar` to be used by
  [`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md).
  Example: `"C:/MetFrag2.4.2-CL.jar"`.

- `patRoon.path.MetFragCompTox`: The complete file path to the CompTox
  database `csv` file. See
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
  for more details.

- `patRoon.path.MetFragPubChemLite`: The complete file path to the
  PubChemLite database `csv` file. See
  [`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
  for more details.

- `patRoon.path.SIRIUS`: The directory in which the SIRIUS binaries are
  installed. Used by all functions that interface with `SIRIUS`, such as
  [`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)
  and
  [`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md).
  Example: `"C:/sirius-win64-3.5.1"`. Note that the location of the
  binaries differs for each operating system.

- `patRoon.path.OpenMS`: The path in which the `OpenMS` binaries are
  installed.

- `patRoon.path.obabel`: The path in which the `OpenBabel` binaries are
  installed.

- `patRoon.path.BiotransFormer` The full file path to the
  `biotransformer` `.jar` command line utility. This needs to be set
  when
  [`generateTPsBioTransformer`](https://rickhelmus.github.io/patRoon/reference/generateTPsBioTransformer.md)
  is used. For more details see
  <https://bitbucket.org/djoumbou/biotransformer/src/master>.

- `patRoon.path.limits` A path to a customized [limits YAML
  file](https://rickhelmus.github.io/patRoon/reference/limits.md).

Most external dependencies are provided by patRoonExt or otherwise found
in the system environment PATH variable. However, the `patRoon.path.*`
options should be set if this fails or you want to override the
location. The
[`verifyDependencies`](https://rickhelmus.github.io/patRoon/reference/verifyDependencies.md)
function can be used to assess if dependencies are found.

## See also

Useful links:

- <https://github.com/rickhelmus/patRoon>

- Report bugs at <https://github.com/rickhelmus/patRoon/issues>

## Author

**Maintainer**: Rick Helmus <r.helmus@uva.nl>
([ORCID](https://orcid.org/0000-0001-9401-3133))

Authors:

- Rick Helmus <r.helmus@uva.nl>
  ([ORCID](https://orcid.org/0000-0001-9401-3133))

Other contributors:

- Olaf Brock ([ORCID](https://orcid.org/0000-0003-4727-8459))
  \[contributor\]

- Vittorio Albergamo ([ORCID](https://orcid.org/0000-0002-5347-1362))
  \[contributor\]

- Andrea Brunner ([ORCID](https://orcid.org/0000-0002-2801-1751))
  \[contributor\]

- Emma Schymanski ([ORCID](https://orcid.org/0000-0001-6868-8145))
  \[contributor\]

- Bas van de Velde ([ORCID](https://orcid.org/0000-0003-1292-3251))
  \[contributor\]

- Leon Saal ([ORCID](https://orcid.org/0000-0002-3522-7729))
  \[contributor\]
