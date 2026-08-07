# Functionality to predict toxicities

Functions to predict toxicities from SMILES and/or `SIRIUS+CSI:FingerID`
fingerprints using the MS2Tox package.

## Usage

``` r
calculateTox(fGroups, ...)

# S4 method for class 'featureGroups'
calculateTox(fGroups, featureAnn)

# S4 method for class 'featureGroupsSet'
calculateTox(fGroups, featureAnn)

# S4 method for class 'compounds'
predictTox(
  obj,
  LC50Mode = "static",
  concUnit = "ugL",
  updateScore = FALSE,
  scoreWeight = 1,
  parallel = TRUE
)

# S4 method for class 'featureGroupsScreening'
predictTox(obj, LC50Mode = "static", concUnit = "ugL")

# S4 method for class 'featureGroupsScreening'
calculateTox(fGroups, featureAnn = NULL)

# S4 method for class 'featureGroupsScreeningSet'
predictTox(obj, LC50Mode = "static", concUnit = "ugL")

# S4 method for class 'featureGroupsScreeningSet'
calculateTox(fGroups, featureAnn = NULL)

# S4 method for class 'compoundsSet'
predictTox(obj, ...)

# S4 method for class 'compoundsSIRIUS'
predictTox(obj, type = "FP", LC50Mode = "static", concUnit = "ugL")

# S4 method for class 'formulasSet'
predictTox(obj, ...)

# S4 method for class 'formulasSIRIUS'
predictTox(obj, LC50Mode = "static", concUnit = "ugL")
```

## Arguments

- fGroups:

  For `predictTox` methods for feature annotations: The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which the annotations were performed.

  For `calculateTox`: The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which toxicities should be assigned.

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
  which contains toxicities. Optional if `calculateTox` is called on
  suspect screening results (*i.e.*
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)
  method).

- obj:

  The workflow object for which predictions should be performed, *e.g.*
  feature groups with screening results
  ([`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md))
  or compound annotations
  ([`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)).

- LC50Mode:

  The mode used for predictions: should be `"static"` or `"flow"`.

- concUnit:

  The concentration unit for calculated toxicities. Can be molar based
  (`"nM"`, `"uM"`, `"mM"`, `"M"`) or mass based (`"ngL"`, `"ugL"`,
  `"mgL"`, `"gL"`). Furthermore, can be prefixed with `"log "` for
  logarithmic concentrations (*e.g.* `"log mM"`).

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

## Value

`predictTox` returns an object amended with LC 50 values
(`LC50_SMILES`/`LC50_SIRFP` columns).

`calculateTox` returns a
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
based object amended with toxicity values for each feature group
(accessed with the
[`toxicities`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
method).

## Details

The [MS2Tox](https://github.com/kruvelab/MS2Tox) R package predicts
toxicities from SMILES and/or MS/MS fingerprints obtained with
`SIRIUS+CSI:FingerID`. The `predictTox` method functions interface with
this package to predict toxicities, which can then be assigned to
feature groups with the `calculateTox` method function.

## Note

The [rcdk](https://CRAN.R-project.org/package=rcdk) package and
[OpenBabel](https://github.com/openbabel/openbabel) tool are used
internally to calculate molecular weights. Please make sure that
`OpenBabel` is installed.

## Predicting toxicities

The toxicities are predicted with the `predictTox` generic functions,
which accepts the following input:

- [Suspect screening
  results](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).
  The SMILES data is used to predict toxicities for suspect hits.

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
  performed on a per formula basis, hence, toxicities for isomers will
  be equal.

- Compound annotation data obtained with algorithms other than
  `"sirius"`. The toxicities are predicted from SMILES data.

When SMILES data is used then predictions of toxicities are generally
more accurate. However, calculations with `SIRIUS+CSI:FingerID`
fingerprints are faster and only require the formula and MS/MS spectrum,
*i.e.* not the full structure. Hence, calculations with SMILES are
mostly useful in suspect screening workflows, or with high confidence
compound annotation data, whereas MS/MS fingerprints are suitable with
unknowns.

For annotation data the calculations are performed for *all* candidates.
This can especially lead to long running calculations when SMILES data
is used. Hence, it is **strongly** recommended to first prioritize the
annotation results, *e.g.* with the `topMost` argument to the [filter
method](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md).

When toxicities are predicted from `SIRIUS+CSI:FingerID` fingerprints
then only formula and MS/MS spectra are used, even if compound
annotations are used for input. The major difference is that with
formula annotation input *all* formula candidates for which a
fingerprint could be generated are considered, whereas with compound
annotations only candidate formulae are considered for which also a
structure could be assigned. Hence, the formula annotation input could
be more comprehensive, whereas predictions from structure annotations
could lead to more representative results as only formulae are
considered for which at least one structure could be assigned.

## Assigning toxicities

The `calculateTox` generic function is used to assign toxicities for
each feature using the toxicities discussed in the previous section. The
function takes toxicities from suspect screening results and/or feature
annotation data. If multiple toxicities were predicted for the same
feature group, for instance when multiple annotation candidates or
suspect hits for this feature group are present, then a toxicities is
assigned for all toxicities. These values can later be easily aggregated
with *e.g.* the
[as.data.table](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
function.

## References

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .\
\
Guha R (2007). “Chemical Informatics Functionality in R.” *Journal of
Statistical Software*, **18**(6).

Peets P, Wang W, MacLeod M, Breitholtz M, Martin JW, Kruve A (2022).
“MS2Tox Machine Learning Tool for Predicting the Ecotoxicity of
Unidentified Chemicals in Water by Nontarget LC-HRMS.” *Environmental
Science & Technology*, **56**(22), 15508-15517.
[doi:10.1021/acs.est.2c02536](https://doi.org/10.1021/acs.est.2c02536) .
PMID: 36269851, https://doi.org/10.1021/acs.est.2c02536.

## See also

[Concentration
prediction](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)
