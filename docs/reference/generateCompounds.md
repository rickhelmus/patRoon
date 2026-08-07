# Automatic compound annotation

Automatically perform chemical compound annotation for feature groups.

## Usage

``` r
generateCompounds(
  fGroups,
  MSPeakLists,
  algorithm,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  ...
)

# S4 method for class 'featureGroups'
generateCompounds(
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
  object which should be annotated. This should be the same or a subset
  of the object that was used to create the specified `MSPeakLists`. In
  the case of a subset only the remaining feature groups in the subset
  are considered.

- MSPeakLists:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- algorithm:

  A character string describing the algorithm that should be used:
  `"metfrag"`, `"sirius"`, `"library"`

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- ...:

  Any parameters to be passed to the selected compound generation
  algorithm.

## Value

A
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
derived object containing all compound annotations.

## Details

Several algorithms are provided to automatically perform compound
annotation for feature groups. To this end, measured masses for all
feature groups are searched within online database(s) (*e.g.*
[PubChem](https://pubchem.ncbi.nlm.nih.gov/)) to retrieve a list of
potential candidate chemical compounds. Depending on the algorithm and
its parameters, further scoring of candidates is then performed using,
for instance, matching of measured and theoretical isotopic patterns,
presence within other data sources such as patent databases and
similarity of measured and in-silico predicted MS/MS fragments. Note
that this process is often quite time consuming, especially for large
feature group sets. Therefore, this is often one of the last steps
within the workflow and not performed before feature groups have been
prioritized.

`generateCompounds` is a generic function that will generateCompounds by
one of the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `generateCompoundsMetFrag` and
`generateCompoundsSIRIUS`. While these functions may be called directly,
`generateCompounds` provides a generic interface and is therefore
usually preferred.

## Scorings

Each algorithm implements their own scoring system. Their names have
been simplified and harmonized where possible. The
[`compoundScorings`](https://rickhelmus.github.io/patRoon/reference/compoundScorings.md)
function can be used to get an overview of both the algorithm specific
and generic scoring names.

## Sets workflows

With a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md),
annotation is first performed for each set. This is important, since the
annotation algorithms typically cannot work with data from mixed
ionization modes. The annotation results are then combined to generate a
*sets consensus*:

- The annotation tables for each feature group from the set specific
  data are combined. Rows with overlapping candidates (determined by the
  first-block InChIKey) are merged.

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
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
output class and its methods and the algorithm specific functions:
[`generateCompoundsMetFrag`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsMetFrag.md),
[`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md),
[`generateCompoundsLibrary`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsLibrary.md)
