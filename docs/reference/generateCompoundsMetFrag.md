# Compound annotation with MetFrag

Uses the metfRag package or `MetFrag CL` for compound identification
(see <http://ipb-halle.github.io/MetFrag/>).

## Usage

``` r
generateCompoundsMetFrag(fGroups, ...)

# S4 method for class 'featureGroups'
generateCompoundsMetFrag(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  method = "CL",
  timeout = 300,
  timeoutRetries = 2,
  errorRetries = 2,
  topMost = 100,
  dbRelMzDev = defaultLim("mz", "narrow_rel"),
  fragRelMzDev = defaultLim("mz", "narrow_rel"),
  fragAbsMzDev = defaultLim("mz", "narrow"),
  adduct = NULL,
  database = "pubchemlite",
  extendedPubChem = "auto",
  chemSpiderToken = "",
  scoreTypes = compoundScorings("metfrag", database, onlyDefault = TRUE)$name,
  scoreWeights = 1,
  preProcessingFilters = c("UnconnectedCompoundFilter", "IsotopeFilter"),
  postProcessingFilters = c("InChIKeyFilter"),
  maxCandidatesToStop = 2500,
  identifiers = NULL,
  extraOpts = NULL,
  minIMSSpecSim = 0
)

# S4 method for class 'featureGroupsSet'
generateCompoundsMetFrag(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  method = "CL",
  timeout = 300,
  timeoutRetries = 2,
  errorRetries = 2,
  topMost = 100,
  dbRelMzDev = defaultLim("mz", "narrow_rel"),
  fragRelMzDev = defaultLim("mz", "narrow_rel"),
  fragAbsMzDev = defaultLim("mz", "narrow"),
  adduct = NULL,
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

- method:

  Which method should be used for MetFrag execution: `"CL"` for
  `MetFragCL` and `"R"` for `MetFragR`. The former is usually much
  faster and recommended.

- timeout:

  Maximum time (in seconds) before a metFrag query for a feature group
  is stopped. Also see `timeoutRetries` argument.

- timeoutRetries:

  Maximum number of retries after reaching a timeout before completely
  skipping the metFrag query for a feature group. Also see `timeout`
  argument.

- errorRetries:

  Maximum number of retries after an error occurred. This may be useful
  to handle e.g. connection errors.

- topMost:

  Only keep this number of candidates (per feature group) with highest
  score. Set to `NULL` to always keep all candidates, however, please
  note that this may result in significant usage of CPU/RAM resources
  for large numbers of candidates.

- dbRelMzDev:

  Relative mass deviation (in ppm) for database search. Sets the
  DatabaseSearchRelativeMassDeviation option.

- fragRelMzDev:

  Relative mass deviation (in ppm) for fragment matching. Sets the
  FragmentPeakMatchRelativeMassDeviation option.

- fragAbsMzDev:

  Absolute mass deviation (in Da) for fragment matching. Sets the
  FragmentPeakMatchAbsoluteMassDeviation option.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- database:

  Compound database to use. Valid values are: `"pubchem"`,
  `"chemspider"`, `"for-ident"`, `"comptox"`, `"pubchemlite"`, `"kegg"`,
  `"sdf"`, `"psv"` and `"csv"`. See section below for more information.
  Sets the `MetFragDatabaseType` option.

- extendedPubChem:

  If `database="pubchem"`: whether to use the *extended* database that
  includes information for compound scoring (*i.e.* number of
  patents/PubMed references). Note that downloading candidates from this
  database might take extra time. Valid values are: `FALSE` (never use
  it), `TRUE` (always use it) or `"auto"` (default, use if specified
  scorings demand it).

- chemSpiderToken:

  A character string with the [ChemSpider security
  token](http://www.chemspider.com/AboutServices.aspx) that should be
  set when the ChemSpider database is used. Sets the ChemSpiderToken
  option.

- scoreTypes:

  A character vector defining the scoring types. See the `Scorings`
  section below for more information. Note that both generic and
  `MetFrag` specific names are accepted (*i.e.* `name` and `metfrag`
  columns returned by
  [`compoundScorings`](https://rickhelmus.github.io/patRoon/reference/compoundScorings.md)).
  When a local database is used, the name should match what is given
  there (`e.g` column names when `database=csv`). Note that MetFrag may
  still report other scoring data, however, these are not used for
  ranking. Sets the MetFragScoreTypes option.

- scoreWeights:

  Numeric vector containing weights of the used scoring types. Order is
  the same as set in `scoreTypes`. Values are recycled if necessary.
  Sets the MetFragScoreWeights option.

- preProcessingFilters, postProcessingFilters:

  A character vector defining pre/post filters applied before/after
  fragmentation and scoring (*e.g.* `"UnconnectedCompoundFilter"`,
  `"IsotopeFilter"`, `"ElementExclusionFilter"`). Some methods require
  further options to be set. For all filters and more information refer
  to the `Candidate Filters` section on the [MetFragR
  homepage](http://ipb-halle.github.io/MetFrag/projects/metfragr/). Sets
  the MetFragPreProcessingCandidateFilter and
  `MetFragPostProcessingCandidateFilter` options.

- maxCandidatesToStop:

  If more than this number of candidate structures are found then
  processing will be aborted and no results this feature group will be
  reported. Low values increase the chance of missing data, whereas too
  high values will use too much computer resources and significantly
  slowdown the process. Sets the MaxCandidateLimitToStop option.

- identifiers:

  A `list` containing for each feature group a character vector with
  database identifiers that should be used to find candidates for a
  feature group (the list should be named by feature group names). If
  `NULL` all relevant candidates will be retrieved from the specified
  database. An example usage scenario is to obtain the list of candidate
  identifiers from a
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object obtained with
  [`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md)
  using the
  [`identifiers`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  method. This way, only those candidates will be searched by MetFrag
  that were generated by SIRIUS+CSI:FingerID. Sets the
  PrecursorCompoundIDs option.

- extraOpts:

  A named `list` containing further settings `MetFrag`. See the
  [MetFragR](http://ipb-halle.github.io/MetFrag/projects/metfragr/) and
  [MetFrag CL](http://ipb-halle.github.io/MetFrag/projects/metfragcl/)
  homepages for all available options. Set to `NULL` to ignore.

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

`generateCompoundsMetFrag` returns a
[`compoundsMF`](https://rickhelmus.github.io/patRoon/reference/compoundsMF-class.md)
object.

## Details

This function uses MetFrag to generate compound candidates. This
function is called when calling `generateCompounds` with
`algorithm="metfrag"`.

Several online compound databases such as
[PubChem](https://pubchem.ncbi.nlm.nih.gov/) and
[ChemSpider](http://www.chemspider.com/) may be chosen for retrieval of
candidate structures. This method requires the availability of MS/MS
data, and feature groups without it will be ignored. Many options exist
to score and filter resulting data, and it is highly suggested to
optimize these to improve results. The `MetFrag` options `PeakList`,
`IonizedPrecursorMass` and `ExperimentalRetentionTimeValue` (in minutes)
fields are automatically set from feature data.

## Scorings

`MetFrag` supports *many* different scorings to rank candidates. The
[`compoundScorings`](https://rickhelmus.github.io/patRoon/reference/compoundScorings.md)
function can be used to get an overview: (some columns are omitted)

|  |  |  |
|----|----|----|
| **name** | **metfrag** | **database** |
| score | Score |  |
| fragScore | FragmenterScore |  |
| metFusionScore | OfflineMetFusionScore |  |
| individualMoNAScore | OfflineIndividualMoNAScore |  |
| numberPatents | PubChemNumberPatents | pubchem |
| numberPatents | Patent_Count | pubchemlite |
| pubMedReferences | PubChemNumberPubMedReferences | pubchem |
| pubMedReferences | ChemSpiderNumberPubMedReferences | chemspider |
| pubMedReferences | NUMBER_OF_PUBMED_ARTICLES | comptox |
| pubMedReferences | PubMed_Count | pubchemlite |
| extReferenceCount | ChemSpiderNumberExternalReferences | chemspider |
| dataSourceCount | ChemSpiderDataSourceCount | chemspider |
| referenceCount | ChemSpiderReferenceCount | chemspider |
| RSCCount | ChemSpiderRSCCount | chemspider |
| smartsInclusionScore | SmartsSubstructureInclusionScore |  |
| smartsExclusionScore | SmartsSubstructureExclusionScore |  |
| suspectListScore | SuspectListScore |  |
| retentionTimeScore | RetentionTimeScore |  |
| CPDATCount | CPDAT_COUNT | comptox |
| TOXCASTActive | TOXCAST_PERCENT_ACTIVE | comptox |
| dataSources | DATA_SOURCES | comptox |
| pubChemDataSources | PUBCHEM_DATA_SOURCES | comptox |
| EXPOCASTPredExpo | EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG/KG-BW/DAY | comptox |
| ECOTOX | ECOTOX | comptox |
| NORMANSUSDAT | NORMANSUSDAT | comptox |
| MASSBANKEU | MASSBANKEU | comptox |
| TOX21SL | TOX21SL | comptox |
| TOXCAST | TOXCAST | comptox |
| KEMIMARKET | KEMIMARKET | comptox |
| MZCLOUD | MZCLOUD | comptox |
| pubMedNeuro | PubMedNeuro | comptox |
| CIGARETTES | CIGARETTES | comptox |
| INDOORCT16 | INDOORCT16 | comptox |
| SRM2585DUST | SRM2585DUST | comptox |
| SLTCHEMDB | SLTCHEMDB | comptox |
| THSMOKE | THSMOKE | comptox |
| ITNANTIBIOTIC | ITNANTIBIOTIC | comptox |
| STOFFIDENT | STOFFIDENT | comptox |
| KEMIMARKET_EXPO | KEMIMARKET_EXPO | comptox |
| KEMIMARKET_HAZ | KEMIMARKET_HAZ | comptox |
| REACH2017 | REACH2017 | comptox |
| KEMIWW_WDUIndex | KEMIWW_WDUIndex | comptox |
| KEMIWW_StpSE | KEMIWW_StpSE | comptox |
| KEMIWW_SEHitsOverDL | KEMIWW_SEHitsOverDL | comptox |
| ZINC15PHARMA | ZINC15PHARMA | comptox |
| PFASMASTER | PFASMASTER | comptox |
| peakFingerprintScore | AutomatedPeakFingerprintAnnotationScore |  |
| lossFingerprintScore | AutomatedLossFingerprintAnnotationScore |  |
| agroChemInfo | AgroChemInfo | pubchemlite |
| bioPathway | BioPathway | pubchemlite |
| drugMedicInfo | DrugMedicInfo | pubchemlite |
| foodRelated | FoodRelated | pubchemlite |
| pharmacoInfo | PharmacoInfo | pubchemlite |
| safetyInfo | SafetyInfo | pubchemlite |
| toxicityInfo | ToxicityInfo | pubchemlite |
| knownUse | KnownUse | pubchemlite |
| disorderDisease | DisorderDisease | pubchemlite |
| identification | Identification | pubchemlite |
| annoTypeCount | AnnoTypeCount | pubchemlite |
| annotHitCount | AnnotHitCount | pubchemlite |

In addition, the
[`compoundScorings`](https://rickhelmus.github.io/patRoon/reference/compoundScorings.md)
function is also useful to programmatically generate a set of scorings
to be used for ranking with `MetFrag`. For instance, the following can
be given to the `scoreTypes` argument to use all default scorings for
PubChem:
`compoundScorings("metfrag", "pubchem", onlyDefault=TRUE)$name`.

For all `MetFrag` scoring types refer to the `Candidate Scores` section
on the [MetFragR
homepage](http://ipb-halle.github.io/MetFrag/projects/metfragr/).

## MetFrag databases

When `database="chemspider"` setting the `chemSpiderToken` argument is
mandatory.

If a local database is chosen via `sdf`, `psv`, or `csv` then its file
location should be set with the `LocalDatabasePath` value via the
`extraOpts` argument. For example:
`extraOpts = list(LocalDatabasePath = "C:/myDB.csv")`.

If `database="pubchemlite"` or `database="comptox"` and patRoonExt is
*not* installed then the file location must be specified as above or by
setting the
`patRoon.path.MetFragPubChemLite`/`patRoon.path.MetFragCompTox` option.
See the installation section in the handbook for more details.

If `database="pubchemlite"` and the local file has CCS predictions (see
<https://zenodo.org/records/15311000>), then CCS values will be copied
from the corresponding adduct of the feature groups (or as specified by
the `adduct` argument). These can be converted to mobilities with
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_comp.md)
and used for candidate filtering with
[`filter`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md).

## Parallelization

generateCompoundsMetFrag uses multiprocessing to parallelize
computations. Please see the parallelization section in the handbook for
more details and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

When local database files are used with `generateCompoundsMetFrag`
(*e.g.* when `database` is set to `"pubchemlite"`, `"csv"` etc.) and
patRoon.MP.method="future", then the database file must be present on
all the nodes. When `pubchemlite` or `comptox` is used, the location for
these databases can be configured on the host with the respective
package options (patRoon.path.MetFragPubChemLite and
patRoon.path.MetFragCompTox) or made available by installing the
patRoonExt package. Note that these files must *also* be present on the
local host computer, even if it is not participating in computations.

If the compound database is *not* local, *e.g.* `database="pubchem"`,
then parallelization is disabled to avoid connection errors that
typically occur otherwise.

## References

Ruttkies C, Schymanski EL, Wolf S, Hollender J, Neumann S (2016).
“MetFrag relaunched: incorporating strategies beyond in silico
fragmentation.” *Journal of Cheminformatics*, **8**(1).
[doi:10.1186/s13321-016-0115-9](https://doi.org/10.1186/s13321-016-0115-9)
.\
\
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

## See also

[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md)
for more details and other algorithms.
