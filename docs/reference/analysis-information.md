# Properties of sample analyses

Properties for the sample analyses used in the workflow and utilities to
automatically generate this information.

## Usage

``` r
generateAnalysisInfo(
  fromRaw = NULL,
  fromCentroid = NULL,
  fromProfile = NULL,
  fromIMS = NULL,
  convCentroid = NULL,
  convProfile = NULL,
  convIMS = NULL,
  ...
)

generateAnalysisInfoFromEnviMass(path)
```

## Arguments

- fromRaw, fromCentroid, fromProfile, fromIMS:

  One or more file paths that should be used for finding analyses that
  are stored as raw, centroided, profile or IMS data, respectively (see
  details below). Set to `NULL` to skip file detection for a particular
  file type.

- convCentroid, convProfile, convIMS:

  These arguments specify the [MS file
  conversion](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
  destination paths for centroided, profile and IMS data, respectively.
  These paths are used for those analyses for which no file with a
  particular file type could be found in the directories specified by
  the respective `from*` arguments. Set to `NULL` to not set any
  destination directory. If multiple paths are specified then these will
  be recycled to fill the table rows.

- ...:

  Any other columns that should be added to the analysis information
  table, such as `replicate` and `blank`. The arguments specified by ...
  should be named. Vectors are recycled to the number of rows of the
  table.

- path:

  The path of the enviMass project.

## Value

`generateAnalysisInformation` returns a `data.frame` with automatically
generated analysis information.

## Details

In patRoon a *sample analysis*, or simply *analysis*, refers to a single
MS analysis file (sometimes also called *sample* or *file*). The
*analysis information* summarizes several properties for the analyses,
and is used in various steps throughout the workflow, such as
[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md),
averaging intensities of feature groups and blank subtraction. The
analysis information should be a `data.frame` or `data.table` with a set
of mandatory and optional columns (described below).

`generateAnalysisInfo` is an utility function that automatically
generates analysis information. It scans given directories for analysis
files, and uses this to automatically fill in the `analysis` and
`path_*` columns. This function automatically groups together analyses
that are stored with different file types and formats (see further
details below).

`generateAnalysisInfoFromEnviMass` loads analysis information from an
enviMass project. Note: this funtionality has only been tested with
older versions of enviMass.

## Mandatory analysis information columns

The following columns should be present in the analysis information:

- `path_raw`, `path_centroid`, `path_profile`, `path_ims` Specifies the
  directory path for the raw, centroided, profile and IMS data,
  respectively. See below for more details. At least one column should
  not be empty for each row.

- `analysis` the file name **without** extension and without directory
  path. Must be **unique** across all table rows.

- `replicate` name of the *replicate*. Used to group analyses together
  that are replicates of each other. Thus, the `replicate` column for
  all analyses considered to be belonging to the same replicate should
  have an equal (but unique) value. Used for *e.g.* averaging and
  [`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md).

- `blank` all analyses within this replicate are used by the
  `featureGroups` method of
  [`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md)
  for blank subtraction. Multiple entries can be entered by separation
  with a comma. May be empty (`""`) if no blank subtraction is desired.

## Analysis paths, file types and file formats

Depending on the workflow step, different *file types* for the same
analysis may be required.

- `raw` Specifies the directory to raw HRMS files (*e.g.* `.raw`, `.d`).
  This is used by *e.g.* [conversion of raw MS
  data](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
  and the [OpenTIMS
  backend](https://rickhelmus.github.io/patRoon/reference/msdata.md).

- `centroid` Specifies the directory to centroided and exported HRMS
  files (`.mzML`, `.mzXML`). These files are required by most feature
  finding algorithms.

- `profile` Specifies the directory to exported but not centroided
  (*i.e.* profile) HRMS data files (`.mzML`, `.mzXML`). This is
  currently only used by
  [`findFeaturesSAFD`](https://rickhelmus.github.io/patRoon/reference/findFeaturesSAFD.md).

- `ims` Specifies the directory to exported IMS-HRMS data (`.mzML`).
  This is required in IMS workflows, unless raw IMS-HRMS data is
  directly loaded with the [OpenTIMS
  backend](https://rickhelmus.github.io/patRoon/reference/msdata.md).
  See *e.g.*
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)
  for more details.

Some workflows may require multiple *file formats* for a same *file
type*. In this case, the file formats should be stored within the same
directory specified by the respective `path_*` column. For instance, if
feature finding algorithms from
[OpenMS](https://rickhelmus.github.io/patRoon/reference/findFeaturesOpenMS.md)
and
[enviPick](https://rickhelmus.github.io/patRoon/reference/findFeaturesEnviPick.md)
are mixed then centroided `.mzML` and `.mzXML` files are needed, and
files with both file formats must be stored in the directory specified
by `path_centroid`.

If non-raw data files are not yet present and should be exported by [MS
file
conversion](https://rickhelmus.github.io/patRoon/reference/MSConversion.md),
then `path_centroid`, `path_profile` and `path_ims` should specify the
desired destination paths of the converted files.

## Optional columns and sample metadata

The following columns may need to be present:

- `conc` a numeric value specifying the 'concentration' for the
  analysis. This can be actually any kind of numeric value such as
  exposure time, dilution factor or anything else which may be used to
  form a linear relationship. This is used by the
  [`as.data.table`](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
  method if `regression=TRUE`. As of patRoon version 3.0, any other
  column than `"conc"` can be used by setting its name with the
  `regression` argument.

- `norm_conc` a numeric value specifying the *normalization
  concentration* for the analysis. See the
  `Feature intensity normalization` section in the [featureGroups
  documentation](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md))
  for more details.

Any other columns that are present will be added to the `features` and
`featureGroups` objects as metadata. This metadata can be used *e.g.* in
various plotting and data subsetting functions.
