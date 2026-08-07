# Identification confidence estimation

Functions to estimate the identification confidence for suspects and
annotation candidates.

## Usage

``` r
estimateIDConfidence(obj, ...)

# S4 method for class 'formulas'
estimateIDConfidence(
  obj,
  absMzDev = defaultLim("mz", "medium"),
  normalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = NULL
)

# S4 method for class 'compounds'
estimateIDConfidence(
  obj,
  absMzDev = defaultLim("mz", "medium"),
  MSPeakLists = NULL,
  formulas = NULL,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  formulasNormalizeScores = "max",
  compoundsNormalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = NULL
)

# S4 method for class 'featureGroupsScreening'
estimateIDConfidence(
  obj,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  absMzDev = defaultLim("mz", "medium"),
  checkFragments = c("mz", "formula", "compound"),
  formulasNormalizeScores = "max",
  compoundsNormalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = file.path("log", "ident")
)

# S4 method for class 'featureGroupsScreeningSet'
estimateIDConfidence(
  obj,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  absMzDev = defaultLim("mz", "medium"),
  checkFragments = c("mz", "formula", "compound"),
  formulasNormalizeScores = "max",
  compoundsNormalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = file.path("log", "ident")
)

# S4 method for class 'compoundsSet'
estimateIDConfidence(
  obj,
  absMzDev = defaultLim("mz", "medium"),
  MSPeakLists = NULL,
  formulas = NULL,
  formulasNormalizeScores = "max",
  compoundsNormalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = NULL
)

# S4 method for class 'formulasSet'
estimateIDConfidence(
  obj,
  absMzDev = defaultLim("mz", "medium"),
  normalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = NULL
)

numericIDLevel(level)

genIDLevelRulesFile(out, inLevels = NULL, exLevels = NULL)
```

## Arguments

- obj:

  The object for which identification confidence should be estimated.

- ...:

  Method specific arguments.

- absMzDev:

  Maximum absolute *m/z* deviation.

- normalizeScores, compoundsNormalizeScores, formulasNormalizeScores:

  A `character` that specifies how normalization of annotation scorings
  occurs. Either

  `"max"` (normalize to max value) or `"minmax"` (perform min-max
  normalization). Note that normalization of negative scores (e.g.
  output by `SIRIUS`) is always performed as min-max. Furthermore,
  currently normalization for `compounds` takes the original min/max
  scoring values into account when candidates were generated. Thus, for
  `compounds` scoring, normalization is not affected when candidate
  results were removed after they were generated (*e.g.* by use of
  `filter`).

- IDFile:

  A file path to a YAML file with rules used for estimation of
  identification levels. See the `Suspect annotation` section for more
  details. If not specified then a default rules file will be used.

- logPath:

  A directory path to store logging information. If `NULL` then logging
  is disabled. **NOTE**: To avoid slowdowns by logging for potentially
  large number of candidates, logging is disabled for the
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  and
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  methods by default.

- MSPeakLists, formulas, compounds:

  Annotation data
  ([`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md),
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  and
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)).
  All arguments can be `NULL`, but it is recommended to set them if
  possible to allow the most complete estimations.

- specSimParams:

  A named `list` with parameters that influence the calculation of MS
  spectra similarities. See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- checkFragments:

  Which type(s) of MS/MS fragments from workflow data should be checked
  to evaluate the number of suspect fragment matches (*i.e.* from the
  `fragments_mz`/`fragments_formula` columns in the suspect list). Valid
  values are: `"mz"`, `"formula"`, `"compounds"`. The former uses *m/z*
  values in the specified `MSPeakLists` object, whereas the others use
  the formulae that were annotated to MS/MS peaks in the given
  `formulas` or `compounds` objects. Multiple values are possible: in
  this case the maximum number of fragment matches will be reported.

- level:

  The identification level to be converted.

- out:

  The file path to the target file.

- inLevels, exLevels:

  A [regular expression](https://rdrr.io/r/base/regex.html) for the
  identification levels to include or exclude, respectively. For
  instance, `exLevels="4|5"` would exclude level 4 and 5 from the output
  file. Set to `NULL` to ignore.

## Value

`estimateIDConfidence` amends the input object with aforementioned
identification confidence properties.

## Details

The `estimateIDConfidence` methods are used to estimate various
properties to estimate the confidence of identifications assigned to
suspects and feature annotation candidates. These functions are
typically executed after running
[`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md),
[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
and
[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md).
Afterwards, the following columns are added to the result tables
(obtained with *e.g.*
[`screenInfo`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md),
[`annotations`](https://rickhelmus.github.io/patRoon/reference/generics.md)
and
[`as.data.table`](https://rickhelmus.github.io/patRoon/reference/generics.md)):

- `annSim` The *annotation similarity*, defined as the similarity
  between the MS/MS peak list of a feature with (a) only the peaks that
  were annotated and (b) all the peaks. Thus, a value of one means that
  all MS/MS peaks were annotated. The similarity calculation is
  configured with the `specSimParams` argument to
  `estimateIDConfidence`.

- `annSimForm` The annotation similarity specifically for formula
  annotations (equaling the `annSim` column from formula annotations).
  Only calculated for
  [suspects](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)
  and
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md).

- `annSimBoth` The annotation similarity calculated with the combined
  set of annotated MS/MS peaks from formula and compound annotations.
  Only calculated for
  [suspects](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)
  and
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md).

- `estIDLevel` Provides an *estimation* of the identification level,
  roughly following that of (Schymanski et al. 2014) . However, please
  note that this value is only an estimation, and manual interpretation
  is still necessary to assign final identification levels. The
  estimation is done through a set of rules, see the
  `Identification level rules` section below.

In addition, the following columns are specifically added to suspect
screening results:

- `annSimComp` The annotation similarity specifically for compound
  annotations (this equals the `annSim` column in compound annotations.

- `formRank`,`compRank` The rank of the suspect within the
  formula/compound annotation results.

- `maxFrags` The maximum number of MS/MS fragments that can be matched
  for this suspect (based on the `fragments_*` columns from the suspect
  list).

- `maxFragMatches`,`maxFragMatchesRel` The absolute and relative amount
  of experimental MS/MS peaks that were matched from the fragments
  specified in the suspect list. The value for `maxFragMatchesRel` is
  relative to the value for `maxFrags`. The calculation of this column
  is influenced by the `checkFragments` argument to
  `estimateIDConfidence`.

The data for these columns is only calculated if `estimateIDConfidence`
has the required data to do so. For instance, `annSimForm` and
`formRank` are only calculated if the `formulas` argument is set, and
levels for `estIDLevel` will be poor if no compound annotations are
available.

`numericIDLevel` Extracts the numeric part of a given identification
level (*e.g.* `"3a"` becomes `3`).

`genIDLevelRulesFile` Generates a template YAML file that is used to
configure the rules for automatic estimation of identification levels.
This file can then be used as input for `estimateIDConfidence`.

## Identification level rules

The estimation of identification levels is configured through a YAML
file which specifies the rules for each level. The default file is shown
below.

    1:
        suspectFragments: 3
        retention: 12
    2a:
        or:
            - individualMoNAScore:
                min: 0.9
                higherThanNext: .inf
            - libMatch:
                min: 0.9
                higherThanNext: .inf
        rank:
            max: 1
            type: compound
    3a:
        or:
            - individualMoNAScore: 0.7
            - libMatch: 0.7
    3b:
        suspectFragments: 3
    3c:
        annMSMSSim:
            type: compound
            min: 0.7
    4a:
        annMSMSSim:
            type: formula
            min: 0.7
        isoScore:
            min: 0.5
            higherThanNext: 0.2
        rank:
            max: 1
            type: formula
    4b:
        isoScore:
            min: 0.9
            higherThanNext: 0.2
        rank:
            max: 1
            type: formula
    5:
        all: yes

Most of the file should be self-explanatory. Some notes:

- Each rule is either a field of `suspectFragments` (minimum number of
  MS/MS fragments matched from suspect list), `retention` (maximum
  retention deviation from suspect list), `rank` (the maximum annotation
  rank from formula or compound annotations), `all` (this level is
  always matched) or any of the scorings available from the formula or
  compound annotations.

- In case any of the rules could be applied to either formula or
  compound annotations, the annotation type must be specified with the
  `type` field (`formula` or `compound`).

- Identification levels should start with a number and may optionally be
  followed by a alphabetic character. The lowest levels are checked
  first.

- If `relative=yes` then the relative scoring will be used for testing.

- For `suspectFragments`: if the number of fragments from the suspect
  list (`maxFrags` column) is less then the minimum rule value, the
  minimum is adjusted to the number of available fragments.

- The `or` and `and` keywords can be used to combine multiple
  conditions.

- Any conditions that require suspect data (*e.g.* `suspectFragments`)
  are only met with the suspects method for `estimateIDConfidence`
  method.

A template rules file can be generated with the `genIDLevelRulesFile`
function, and this file can subsequently passed to
`estimateIDConfidence`. The file format is highly flexible and
(sub)levels can be added or removed if desired. Note that the default
file is currently only suitable when annotation is performed with
`GenForm` and `MetFrag`, for other algorithms it is crucial to modify
the rules.

## Sets workflows

`estimateIDConfidence` performs its estimations per set. In addition,
the following overall (not set specific) columns are calculated:

- `formRank` and `compRank` based on the ranking of the formula/compound
  in the set consensus data.

- `estIDLevel`: based on the 'best' estimated identification level among
  the sets data (*i.e.* the lowest). In case there is a tie between
  sub-levels (*e.g.* `3a` and `3b`), then the sub-level is stripped
  (*e.g.* `3`).

- Annotation similarities: taken as the maximum value from the data for
  each set.

## References

Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J
(2014). “Identifying Small Molecules via High Resolution Mass
Spectrometry: Communicating Confidence.” *Environmental Science and
Technology*, **48**(4), 2097–2098.
[doi:10.1021/es5002105](https://doi.org/10.1021/es5002105) .\
\
Stein SE, Scott DR (1994). “Optimization and testing of mass spectral
library search algorithms for compound identification.” *Journal of the
American Society for Mass Spectrometry*, **5**(9), 859–866.
[doi:10.1016/1044-0305(94)87009-8](https://doi.org/10.1016/1044-0305%2894%2987009-8)
.

## Author

Rick Helmus \<<r.helmus@uva.nl>\>, Emma Schymanski
\<<emma.schymanski@uni.lu>\> (contributions to identification level
rules), Bas van de Velde (contributions to spectral similarity
calculation).
