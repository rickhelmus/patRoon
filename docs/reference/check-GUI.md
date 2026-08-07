# Interactive GUI utilities to check workflow data

These functions provide interactive utilities to explore and review
workflow data using a shiny graphical user interface (GUI). In addition,
unsatisfactory data (*e.g.* noise identified as a feature and unrelated
feature groups in a component) can easily be selected for removal.

## Usage

``` r
checkFeatures(
  fGroups,
  session = "checked-features.yml",
  EICParams = getDefEICParams(),
  EIMParams = getDefEIMParams(),
  clearSession = FALSE
)

checkComponents(
  components,
  fGroups,
  session = "checked-components.yml",
  EICParams = getDefEICParams(),
  clearSession = FALSE
)

# S4 method for class 'components'
checkComponents(
  components,
  fGroups,
  session = "checked-components.yml",
  EICParams = getDefEICParams(),
  clearSession = FALSE
)

importCheckFeaturesSession(
  sessionIn,
  sessionOut,
  fGroups,
  rtWindow = defaultLim("retention", "narrow"),
  mzWindow = defaultLim("mz", "narrow"),
  mobWindow = defaultLim("mobility", "narrow"),
  overwrite = FALSE
)

# S4 method for class 'featureGroups'
checkFeatures(
  fGroups,
  session = "checked-features.yml",
  EICParams = getDefEICParams(),
  EIMParams = getDefEIMParams(),
  clearSession = FALSE
)

getMCTrainData(fGroups, session)

predictCheckFeaturesSession(fGroups, session, model = NULL, overwrite = FALSE)
```

## Arguments

- fGroups:

  A
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object.

  This should be the 'new' object for `importCheckFeaturesSession` for
  which the session needs to be imported.

- session:

  The session file name.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- EIMParams:

  A named `list` with parameters used for extracted ion mobilogram (EIM)
  creation. See the [EIM
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- clearSession:

  If `TRUE` the session will be completely cleared before starting the
  GUI. This effectively removes all selections for data removal.

- components:

  The
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
  to be checked.

- sessionIn, sessionOut:

  The file names for the input and output sessions.

- rtWindow, mzWindow, mobWindow:

  The retention time (seconds), *m/z* and mobility window (if present)
  used to relate 'old' with 'new' feature groups.

- overwrite:

  Set to `TRUE` to overwrite the output session file if it already
  exists. If `FALSE`, the function will stop with an error message.

- model:

  The model that was created with MetaClean and that should be used to
  predict pass/fail data. If `NULL`, the example model of the
  MetaCleanData package is used.

## Value

A dataframe with the class predictions as well as the associated
probabilities for each EIC as returned by the
[`MetaClean::getPredicitons`](https://rdrr.io/pkg/MetaClean/man/getPredicitons.html)
function. The dataframe has the four columns: EIC, Pred_Class,
Pred_Prob_Pass, Pred_Prob_Fail.

## Details

The data selected for removal is stored in *sessions*. These are `YAML`
files to allow easy external manipulation. The sessions can be used to
restore the selections that were made for data removal when the GUI tool
is executed again. Furthermore, functionality is provided to import and
export sessions. To actually remove the data the
[`filter`](https://rickhelmus.github.io/patRoon/reference/generics.md)
method should be used with the session file as input.

`checkComponents` is used to review components and their feature groups
contained within. A typical use case is to verify that peaks from
features that were annotated as related adducts and/or isotopes are
correctly aligned.

`importCheckFeaturesSession` is used to import a session file that was
generated from a different
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
object. This is useful to avoid re-doing manual interpretation of
chromatographic peaks when, for instance, feature group data is
re-created with different parameters.

`checkFeatures` is used to review chromatographic information for
feature groups. Its main purpose is to assist in reviewing the quality
of detected feature (groups) and easily select unwanted data such as
features with poor peak shapes or noise.

`getMCTrainData` converts a session created by `checkFeatures` to a
`data.frame` that can be used by the MetaClean to train a new model. The
output format is comparable to that from
[`getPeakQualityMetrics`](https://rdrr.io/pkg/MetaClean/man/getPeakQualityMetrics.html).

`predictCheckFeaturesSession` Uses ML data from MetaClean to predict the
quality (Pass/Fail) of feature group data, and converts this to a
session which can be reviewed with `checkFeatures` and used to remove
unwanted feature groups by
[`filter`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md).

## Note

The `topMost` and `topMostByReplicate` [EIC
parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
are ignored.

`checkComponents`: Some componentization algorithms (*e.g.*
[`generateComponentsNontarget`](https://rickhelmus.github.io/patRoon/reference/generateComponentsNontarget.md)
and
[`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md))
may output components where the same feature group in a component is
present multiple times, for instance, when multiple TPs are matched to
the same feature group. If such a feature group is selected for removal,
then *all* of its result in the component will be marked for removal.

`getMCTrainData` only uses session data for selected feature groups.
Selected features for removal are ignored, as this is not supported by
MetaClean.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by these functions to process HRMS (or IMS-HRMS) data.
Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.

## References

Chetnik K, Petrick L, Pandey G (2020). “MetaClean: a machine
learning-based classifier for reduced false positive peak detection in
untargeted LC-MS metabolomics data.” *Metabolomics*, **16**(11).
[doi:10.1007/s11306-020-01738-3](https://doi.org/10.1007/s11306-020-01738-3)
.
