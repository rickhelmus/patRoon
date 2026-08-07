# Functionality to predict quantitative data

Functions to predict response factors and feature concentrations from
SMILES and/or `SIRIUS+CSI:FingerID` fingerprints using the MS2Quant
package.

## Usage

``` r
calculateConcs(fGroups, ...)

# S4 method for class 'featureGroups'
calculateConcs(fGroups, featureAnn, areas = FALSE)

# S4 method for class 'featureGroupsSet'
calculateConcs(fGroups, featureAnn, areas = FALSE)

# S4 method for class 'compounds'
predictRespFactors(
  obj,
  fGroups,
  calibrants,
  eluent,
  organicModifier,
  pHAq,
  concUnit = "ugL",
  calibConcUnit = concUnit,
  updateScore = FALSE,
  scoreWeight = 1,
  parallel = TRUE
)

# S4 method for class 'featureGroupsScreening'
predictRespFactors(
  obj,
  calibrants,
  eluent,
  organicModifier,
  pHAq,
  concUnit = "ugL",
  calibConcUnit = concUnit
)

# S4 method for class 'featureGroupsScreening'
calculateConcs(fGroups, featureAnn = NULL, areas = FALSE)

# S4 method for class 'featureGroupsScreeningSet'
predictRespFactors(obj, calibrants, ...)

# S4 method for class 'featureGroupsScreeningSet'
calculateConcs(fGroups, featureAnn = NULL, areas = FALSE)

# S4 method for class 'compoundsSet'
predictRespFactors(obj, fGroups, calibrants, ...)

# S4 method for class 'compoundsSIRIUS'
predictRespFactors(
  obj,
  fGroups,
  calibrants,
  eluent,
  organicModifier,
  pHAq,
  concUnit = "ugL",
  calibConcUnit = concUnit,
  type = "FP"
)

# S4 method for class 'formulasSet'
predictRespFactors(obj, fGroups, calibrants, ...)

# S4 method for class 'formulasSIRIUS'
predictRespFactors(
  obj,
  fGroups,
  calibrants,
  eluent,
  organicModifier,
  pHAq,
  concUnit = "ugL",
  calibConcUnit = concUnit
)

getQuantCalibFromScreening(
  fGroups,
  concs,
  areas = FALSE,
  average = FALSE,
  IMS = "maybe"
)
```

## Arguments

- fGroups:

  For `predictRespFactors` methods for feature annotations: The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which the annotations were performed.

  For `calculateConcs`: The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which concentrations should be calculated.

  For `getQuantCalibFromScreening`: A feature groups object screened for
  the calibrants with
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- featureAnn:

  A
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)
  object (*e.g.*
  [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)
  or
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md))
  which contains response factors. Optional if `calculateConcs` is
  called on suspect screening results (*i.e.*
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)
  method).

- areas:

  Set to `TRUE` to use peak areas instead of peak heights. Note: for
  `calculateConcs` this should follow what is in the `calibrants` table.

- obj:

  The workflow object for which predictions should be performed, *e.g.*
  feature groups with screening results
  ([`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md))
  or compound annotations
  ([`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)).

- calibrants:

  A `data.frame` with calibrants, see the `Calibration` section below.

  (**sets workflow**) Should be a `list` with the calibrants for each
  set.

- eluent:

  A `data.frame` that describes the LC gradient program. Should have a
  column `time` with the retention time in seconds and a column `B` with
  the corresponding percentage of the organic modifier (`0-100`).

- organicModifier:

  The organic modifier of the mobile phase: either `"MeOH"` (methanol)
  or `"MeCN"` (acetonitrile).

- pHAq:

  The pH of the aqueous part of the mobile phase.

- concUnit:

  The concentration unit for calculated concentrations. Can be molar
  based (`"nM"`, `"uM"`, `"mM"`, `"M"`) or mass based (`"ngL"`, `"ugL"`,
  `"mgL"`, `"gL"`). Furthermore, can be prefixed with `"log "` for
  logarithmic concentrations (*e.g.* `"log mM"`).

- calibConcUnit:

  The concentration unit used in the calibrants table. For possible
  values see the `concUnit` argument.

- updateScore, scoreWeight:

  If `updateScore=TRUE` then the annotation `score` column is updated by
  adding normalized values of the response factor (weighted by
  scoreWeight). Currently, this **only** makes sense for annotations
  performed with `MetFrag`!

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

- type:

  Which types of predictions should be performed: should be `"FP"`
  (`SIRIUS+CSI:FingerID` fingerprints), `"SMILES"` or `"both"`. Only
  relevant for
  [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compoundsSIRIUS-class.md)
  method.

- concs:

  A `data.frame` with concentration data. See the `Calibration` section
  below.

- average:

  Set to `TRUE` to average intensity values within replicates.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  generating calibrant data in IMS workflows. The following options are
  valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  Best to keep this `"maybe"`, as calibration typically doesn't support
  IMS filtered data.

## Value

`predictRespFactors` returns an object amended with response factors
(`RF_SMILES`/`LRF_SIRFP` columns).

`calculateConcs` returns a
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
based object amended with concentrations for each feature group
(accessed with the
[`concentrations`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
method).

## Details

The [MS2Quant](https://github.com/kruvelab/MS2Quant) R package predicts
concentrations from SMILES and/or MS/MS fingerprints obtained with
`SIRIUS+CSI:FingerID`. The `predictRespFactors` method functions
interface with this package to calculate response factors, which can
then be used to calculate feature concentrations with the
`calculateConcs` method function.

## Note

The [rcdk](https://CRAN.R-project.org/package=rcdk) package and
[OpenBabel](https://github.com/openbabel/openbabel) tool are used
internally to calculate molecular weights. Please make sure that
`OpenBabel` is installed.

MS2Quant currently *only* supports `M+H` and `M+` adducts when
performing predictions with `SIRIUS:FingerID` fingerprints. Predictions
for candidates with other adducts, including `M-H`, are skipped with a
warning.

## Calibration

The MS2Quant package requires calibration to convert predicted
ionization efficiencies to instrument/method specific response factors.
The calibration data should be specified with the `calibrants` argument
to `predictRespFactors`. This should be a `data.frame` with intensity
observations at different concentrations for a set of calibrants. Each
row specifies one intensity observation at one concentration. The table
should have the following columns:

- `name` The name of the calibrant. Can be freely chosen.

- `SMILES` The SMILES of the calibrant.

- `rt` The retention time of the calibrant (in seconds).

- `intensity` The peak intensity (or area, see the `areas` argument) of
  the calibrant.

- `conc` The concentration of the calibrant (see the `calibConcUnit`
  argument for specifying the unit).

It is recommended to include multiple calibrants (*e.g.* `>=10`) at
multiple concentrations (*e.g.* `>=5`). The latter is achieved by adding
multiple rows for the same calibrant (keeping the `name`/`SMILES`/`rt`
columns constant). It is also possible to follow the column naming used
by MS2Quant (however retention times should still be in seconds!). For
more details and tips see <https://github.com/kruvelab/MS2Quant>.

The `getQuantCalibFromScreening` function can be used to automatically
generate a calibrants table from a feature groups object with suspect
screening results. Here, the idea is to perform a screening with
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
with a suspect list that contain the calibrants, which is then used to
construct the calibrant table. It is highly recommended to add retention
times for the calibrants in the suspect list to ensure the calibrant is
assigned to the correct feature. Furthermore, it is possible to simply
add the calibrants to the 'regular' suspect list in case a suspect
screening was already part of the workflow. The
`getQuantCalibFromScreening` function still requires you to specify
concentration data, which is achieved via the `concs` argument. This
should be a `data.frame` with a column `name` corresponding to the
calibrant name (*i.e.* same as used by `screenSuspects` above) and
columns with concentration data. The latter columns specify the
concentrations of a calibrant in different replicates (as defined in the
[analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)).
The concentration columns should be named after the corresponding
replicate. Only those replicates that should be used for calibration
need to be included. Furthermore, `NA` values can be used if a replicate
should be ignored for a specific calibrant.

## Predicting response factors

The response factors are predicted with the `predictRespFactors` generic
functions, which accepts the following input:

- [Suspect screening
  results](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
  The SMILES data is used to predict response factors for suspect hits.

- Formula annotation data obtained with `"sirius"` algorithm
  ([`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)).
  The predictions are performed for each formula candidate using
  `SIRIUS+CSI:FingerID` fingerprints. For this reason, the
  `getFingerprint` argument must be set to `TRUE` when generating the
  formula data.

- Compound annotation data obtained with the `"sirius"` algorithm
  ([`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md)).
  The predictions are performed for each annotation candidate using its
  SMILES and/or `SIRIUS+CSI:FingerID` fingerprints. The predictions are
  performed on a per formula basis, hence, response factors for isomers
  will be equal.

- Compound annotation data obtained with algorithms other than
  `"sirius"`. The response factors are predicted from SMILES data.

When SMILES data is used then predictions of response factors are
generally more accurate. However, calculations with
`SIRIUS+CSI:FingerID` fingerprints are faster and only require the
formula and MS/MS spectrum, *i.e.* not the full structure. Hence,
calculations with SMILES are mostly useful in suspect screening
workflows, or with high confidence compound annotation data, whereas
MS/MS fingerprints are suitable with unknowns.

For annotation data the calculations are performed for *all* candidates.
This can especially lead to long running calculations when SMILES data
is used. Hence, it is **strongly** recommended to first prioritize the
annotation results, *e.g.* with the `topMost` argument to the [filter
method](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md).

When response factors are predicted from `SIRIUS+CSI:FingerID`
fingerprints then only formula and MS/MS spectra are used, even if
compound annotations are used for input. The major difference is that
with formula annotation input *all* formula candidates for which a
fingerprint could be generated are considered, whereas with compound
annotations only candidate formulae are considered for which also a
structure could be assigned. Hence, the formula annotation input could
be more comprehensive, whereas predictions from structure annotations
could lead to more representative results as only formulae are
considered for which at least one structure could be assigned.

## Assigning concentrations

The `calculateConcs` generic function is used to assign concentrations
for each feature using the response factors discussed in the previous
section. The function takes response factors from suspect screening
results and/or feature annotation data. If multiple response factors
were predicted for the same feature group, for instance when multiple
annotation candidates or suspect hits for this feature group are
present, then a concentrations is assigned for all response factors.
These values can later be easily aggregated with *e.g.* the
[as.data.table](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
function.

In IMS workflows with post mobility assignments (see
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)),
the intensities of *IMS precursors* are used for the calculation of
concentrations for IMS features (as calibration typically does not
consider differences due to mobility filtering). However, the response
factor assigned to IMS features are still used. This may yield small
differences compared to workflows where concentrations are assigned
prior to mobility assignments, as `assignMobilities` simply copies
concentrations from IMS precursors to IMS features.

## References

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .\
\
Guha R (2007). “Chemical Informatics Functionality in R.” *Journal of
Statistical Software*, **18**(6).

Sepman H, Malm L, Peets P, MacLeod M, Martin J, Breitholtz M, Kruve A
(2023). “Bypassing the Identification: MS2Quant for Concentration
Estimations of Chemicals Detected with Nontarget LC-HRMS from MS2 Data.”
*Analytical Chemistry*, **95**(33), 12329–12338.
[doi:10.1021/acs.analchem.3c01744](https://doi.org/10.1021/acs.analchem.3c01744)
. <https://doi.org/10.1021/acs.analchem.3c01744>.

## See also

[Toxicity
prediction](https://rickhelmus.github.io/patRoon/reference/pred-tox.md)
