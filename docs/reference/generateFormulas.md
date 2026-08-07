# Automatic chemical formula generation

Automatically calculate chemical formulae for all feature groups.

## Usage

``` r
generateFormulas(
  fGroups,
  MSPeakLists,
  algorithm,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  ...
)

# S4 method for class 'featureGroups'
generateFormulas(
  fGroups,
  MSPeakLists,
  algorithm,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  ...
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which formulae should be generated. This should be the same
  or a subset of the object that was used to create the specified
  `MSPeakLists`. In the case of a subset only the remaining feature
  groups in the subset are considered.

- MSPeakLists:

  An
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- algorithm:

  A character string describing the algorithm that should be used:
  `"bruker"`, `"genform"`, `"sirius"`

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- ...:

  Any parameters to be passed to the selected formula generation
  algorithm.

## Value

A
[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
object containing all generated formulae.

## Details

Several algorithms are provided to automatically generate formulae for
given feature groups. All algorithms use the accurate mass of a feature
to back-calculate candidate formulae. Depending on the algorithm and
data availability, other data such as isotopic pattern and MS/MS
fragments may be used to further improve formula assignment and ranking.

`generateFormulas` is a generic function that will generateFormulas by
one of the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `generateFormulasDA` and
`generateFormulasGenForm`. While these functions may be called directly,
`generateFormulas` provides a generic interface and is therefore usually
preferred.

## Candidate assignment

Formula candidate assignment occurs in one of the following ways:

- Candidates are first generated for each feature and then pooled to
  form consensus candidates for the feature group.

- Candidates are directly generated for each feature group by group
  averaged MS peak list data.

With approach (1), scorings and mass errors are averaged and outliers
are removed (controlled by `featThreshold` and `featThresholdAnn`
arguments). Other candidate properties that cannot be averaged are from
the feature from the analysis as specified in the `"analysis"` column of
the results. The second approach only generates candidate formulae once
for every feature group, and is therefore generally much faster.
However, this inherently prevents removal of outliers.

Note that with either approach subsequent workflow steps that use
formula data (*e.g.*
[`addFormulaScoring`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
and
[reporting](https://rickhelmus.github.io/patRoon/reference/reporting.md)
functions) only use formula data that was eventually assigned to feature
groups.

## Scorings

Each algorithm implements their own scoring system. Their names have
been harmonized where possible. An overview is obtained with the
[`formulaScorings`](https://rickhelmus.github.io/patRoon/reference/formulaScorings.md)
function:

|  |  |  |  |  |
|----|----|----|----|----|
| **name** | **genform** | **sirius** | **bruker** | **description** |
| combMatch | comb_match | \- | \- | MS and MS/MS combined match value |
| isoScore | MS_match | isoScore | \- | How well the isotopic pattern matches |
| mSigma | \- | \- | mSigma | Deviation of the isotopic pattern |
| MSMSScore | MSMS_match | treeScore | \- | How well MS/MS data matches |
| score | \- | score | Score | Overall MS formula score |

## Sets workflows

With a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md),
annotation is first performed for each set. This is important, since the
annotation algorithms typically cannot work with data from mixed
ionization modes. The annotation results are then combined to generate a
*sets consensus*:

- The annotation tables for each feature group from the set specific
  data are combined. Rows with overlapping candidates (determined by the
  neutral formula) are merged.

- Set specific data (*e.g.* the ionic formula) is retained by renaming
  their columns with set specific names.

- The MS/MS fragment annotations (`fragInfo` column) from each set are
  combined.

- The scorings for each set are averaged to calculate overall scores. if
  `setAvgSpecificScores=FALSE` then scorings that are considered set
  specific (*e.g.* MS/MS and isotopic pattern match) are *not* averaged.

- The candidates are re-ranked based on their average ranking among the
  set data (if a candidate is absent in a set it is assigned the poorest
  rank in that set).

- The coverage of each candidate among sets is calculated. Depending on
  the `setThreshold` and `setThresholdAnn` arguments, candidates with
  low abundance are removed.

## See also

The
[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
output class and its methods and the algorithm specific functions:
[`generateFormulasDA`](https://rickhelmus.github.io/patRoon/reference/generateFormulasDA.md),
[`generateFormulasGenForm`](https://rickhelmus.github.io/patRoon/reference/generateFormulasGenForm.md),
[`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)

The [GenForm
manual](https://www.researchgate.net/publication/307964728_MOLGEN-MSMS_Software_User_Manual)
(also known as MOLGEN-MSMS).
