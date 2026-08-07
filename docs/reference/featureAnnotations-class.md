# Base feature annotations class

Holds information for all feature group annotations.

## Usage

``` r
# S4 method for class 'featureAnnotations'
annotations(obj)

# S4 method for class 'featureAnnotations'
groupNames(obj)

# S4 method for class 'featureAnnotations'
length(x)

# S4 method for class 'featureAnnotations,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'featureAnnotations,ANY,missing'
x[[i, j]]

# S4 method for class 'featureAnnotations'
x$name

# S4 method for class 'featureAnnotations'
as.data.table(
  x,
  fGroups = NULL,
  fragments = FALSE,
  countElements = NULL,
  countFragElements = NULL,
  OM = FALSE,
  normalizeScores = "none",
  excludeNormScores = defaultExclNormScores(x)
)

# S4 method for class 'featureAnnotations'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'featureAnnotations'
filter(
  obj,
  minExplainedPeaks = NULL,
  scoreLimits = NULL,
  elements = NULL,
  fragElements = NULL,
  lossElements = NULL,
  fragFormulas = NULL,
  lossFormulas = NULL,
  topMost = NULL,
  OM = FALSE,
  maxLevel = NULL,
  negate = FALSE
)

# S4 method for class 'featureAnnotations'
plotVenn(obj, ..., labels = NULL, vennArgs = NULL)

# S4 method for class 'featureAnnotations'
plotUpSet(
  obj,
  ...,
  labels = NULL,
  nsets = NULL,
  nintersects = NA,
  upsetArgs = NULL
)
```

## Arguments

- obj, x:

  `featureAnnotations` object to be accessed

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  feature groups by their index or name, respectively (for the
  order/names see
  [`groupNames()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all feature groups are selected.\
  \
  For `[[`: should be a scalar value.\
  \
  For `delete`: The data to remove from. `i` are the feature groups as
  numeric index, logical or character, `j` the candidates as numeric
  indices (rows). If either is `NULL` then data for all is removed. `j`
  may also be a function: it will be called for each feature group, with
  the annotation table (a `data.table`), the feature group name and any
  other arguments passed as `...` to `delete`. The return value of this
  function specifies the candidate indices (rows) to be removed
  (specified as an `integer` or `logical` vector).

- ...:

  For the `"["` operator: ignored.

  For `delete`: passed to the function specified as `j`.

  Others: Any further (and unique) `featureAnnotations` objects.

- drop:

  ignored.

- name:

  The feature group name (partially matched).

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that was used to generate this object. If not `NULL` it is used
  to add feature group information (retention and *m/z* values).

- fragments:

  If `TRUE` then information on annotated fragments will be included.
  Automatically set to `TRUE` if `countFragElements` is set.

- countElements, countFragElements:

  A `character` vector with elements that should be counted for each
  candidate's formula. For instance, `c("C", "H")` adds columns for both
  carbon and hydrogen amounts of each formula. Note that the neutral
  formula (`neutral_formula` column) is used to count elements of
  non-fragmented formulae, whereas the charged formula of fragments
  (`ion_formula` column in `fragInfo` data) is used for fragments. Set
  to `NULL` to not count any elements.

- OM:

  For `as.data.table`: if set to `TRUE` several columns with information
  relevant for organic matter (OM) characterization will be added (e.g.
  elemental ratios, classification). This will also make sure that
  `countElements` contains at least C, H, N, O, P and S.

  For `filter`: If `TRUE` then several filters are applied to exclude
  unlikely formula candidates present in organic matter (OM). See Source
  section for details.

- normalizeScores:

  A `character` that specifies how normalization of annotation scorings
  occurs. Either `"none"` (no normalization), `"max"` (normalize to max
  value) or `"minmax"` (perform min-max normalization). Note that
  normalization of negative scores (e.g. output by `SIRIUS`) is always
  performed as min-max. Furthermore, currently normalization for
  `compounds` takes the original min/max scoring values into account
  when candidates were generated. Thus, for `compounds` scoring,
  normalization is not affected when candidate results were removed
  after they were generated (*e.g.* by use of `filter`).

- excludeNormScores:

  A `character` vector specifying any compound scoring names that should
  *not* be normalized. Set to `NULL` to normalize all scorings. Note
  that whether any normalization occurs is set by the
  `excludeNormScores` argument.

  For `compounds`: By default `score` and `individualMoNAScore` are set
  to mimic the behavior of the `MetFrag` web interface.

- minExplainedPeaks:

  Minimum number of explained peaks. Set to `NULL` to ignore.

- scoreLimits:

  Filter results by their scores. Should be a named `list` that contains
  two-sized numeric vectors with the minimum/maximum value of a score
  (use `-Inf`/`Inf` for no limits). The names of each element should
  follow the name column of the table returned by
  [`formulaScorings`](https://rickhelmus.github.io/patRoon/reference/formulaScorings.md)`$name`
  and
  [`compoundScorings`](https://rickhelmus.github.io/patRoon/reference/compoundScorings.md)`()$name`.
  For instance, `scoreLimits=list(numberPatents=c(10, Inf))` specifies
  that `numberPatents` should be at least `10`. Note that a result
  without a specified scoring is never removed. If a score term exists
  multiple times, *i.e.* due to a consensus, then a candidate is kept if
  at least one of the terms falls within the range. Set to `NULL` to
  skip this filter.

- elements:

  Only retain candidate formulae (neutral form) that match a given
  elemental restriction. The format of `elements` is a `character`
  string with elements that should be present where each element is
  followed by a valid amount or a range thereof. If no number is
  specified then `1` is assumed. For instance,
  `elements="C1-10H2-20O0-2P"`, specifies that `1-10`, `2-20`, `0-2` and
  `1` carbon, hydrogen, oxygen and phosphorus atoms should be present,
  respectively. When `length(elements)>1` formulas are tested to follow
  at least one of the given elemental restrictions. For instance,
  `elements=c("P", "S")` specifies that either one phosphorus or one
  sulfur atom should be present. Set to `NULL` to ignore this filter.

- fragElements, lossElements:

  Specifies elemental restrictions for fragment or neutral loss formulae
  (charged form). Candidates are retained if at least one of the
  fragment formulae follow (or not follow if `negate=TRUE`) the given
  restrictions. See `elements` for the used format.

- fragFormulas, lossFormulas:

  A `character` vector with one or more formulae, which are matched to
  MS/MS fragments (`fragFormulas`) or neutral losses (`lossFormulas`).
  Candidates are only kept with at least one match (if both
  `fragFormulas` and `lossFormulas` are set then there must be at least
  one match from both). The input formulae are hill sorted prior to
  matching. For non-exact matching: see `fragElements` and
  `lossElements`. Set to `NULL` to ignore these filters.

- topMost:

  Only keep a maximum of `topMost` candidates with highest score (or
  least highest if `negate=TRUE`). Set to `NULL` to ignore.

- maxLevel:

  Filter by maximum identification level (*e.g.* `"3a"`). Set to `NULL`
  to ignore.

- negate:

  If `TRUE` then filters are applied in opposite manner.

- labels:

  A `character` with names to use for labelling. If `NULL` labels are
  automatically generated.

- vennArgs:

  A `list` with further arguments passed to VennDiagram plotting
  functions. Set to `NULL` to ignore.

- nsets, nintersects:

  See [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html). If
  `nsets == NULL` then it will be set to the number of compared items.

- upsetArgs:

  A list with any further arguments to be passed to
  [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html). Set to `NULL` to
  ignore.

## Value

`as.data.table` returns a `data.table`.

`delete` returns the object for which the specified data was removed.

`filter` returns a filtered `featureAnnotations` object.

`plotVenn` (invisibly) returns a list with the following fields:

- `gList` the `gList` object that was returned by the utilized
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram) plotting
  function.

- `areas` The total area for each plotted group.

- `intersectionCounts` The number of intersections between groups.

The order for the `areas` and `intersectionCounts` fields is the same as
the parameter order from the used plotting function (see *e.g.*
`draw.pairwise.venn` and `draw.triple.venn`).

## Details

This class stores annotation data for feature groups, such as molecular
formulae, SMILES identifiers, compound names etc. The class of objects
that are generated by formula and compound annotation
([`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
and
[`generateCompounds`](https://rickhelmus.github.io/patRoon/reference/generateCompounds.md))
are based on this class.

## Methods (by generic)

- `annotations(featureAnnotations)`: Accessor for the `groupAnnotations`
  slot.

- `groupNames(featureAnnotations)`: returns a `character` vector with
  the names of the feature groups for which data is present in this
  object.

- `length(featureAnnotations)`: Obtain total number of candidates.

- `x[i`: Subset on feature groups.

- `x[[i`: Extracts annotation data for a feature group.

- `$`: Extracts annotation data for a feature group.

- `as.data.table(featureAnnotations)`: Generates a table with all
  annotation data for each feature group and other information such as
  element counts.

- `delete(featureAnnotations)`: Completely deletes specified
  annotations.

- `filter(featureAnnotations)`: Provides rule based filtering for
  feature group annotations. Useful to eliminate unlikely candidates and
  speed up further processing.

- `plotVenn(featureAnnotations)`: plots a Venn diagram (using
  [VennDiagram](https://CRAN.R-project.org/package=VennDiagram))
  outlining unique and shared candidates of up to five different
  `featureAnnotations` objects.

- `plotUpSet(featureAnnotations)`: plots an UpSet diagram (using the
  [`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html) function)
  outlining unique and shared candidates between different
  `featureAnnotations` objects.

## Slots

- `groupAnnotations`:

  A `list` with for each annotated feature group a `data.table` with
  annotation data. Use the `annotations` method for access.

- `scoreTypes`:

  A `character` with all the score types present in this object.

- `scoreRanges`:

  The minimum and maximum score values of all candidates for each
  feature group. Used for normalization.

## Source

Calculation of the aromaticity index (AI) and related double bond
equivalents (DBE_AI) is performed as described in Koch 2015. Formula
classification is performed by the rules described in Abdulla 2013.
Filtering of OM related molecules is performed as described in Koch 2006
and Kujawinski 2006. (see references).

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`featureAnnotations`**

    - [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasConsensus`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasUnset`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)

    - [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsConsensus`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsMF`](https://rickhelmus.github.io/patRoon/reference/compoundsMF-class.md)

      - [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsUnset`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compoundsSIRIUS-class.md)

## References

Koch BP, Dittmar T (2015). “From mass to structure: an aromaticity index
for high-resolution mass data of natural organic matter.” *Rapid
Communications in Mass Spectrometry*, **30**(1), 250–250.
[doi:10.1002/rcm.7433](https://doi.org/10.1002/rcm.7433) .\
\
Abdulla HA, Sleighter RL, Hatcher PG (2013). “Two Dimensional
Correlation Analysis of Fourier Transform Ion Cyclotron Resonance Mass
Spectra of Dissolved Organic Matter: A New Graphical Analysis of
Trends.” *Analytical Chemistry*, **85**(8), 3895–3902.
[doi:10.1021/ac303221j](https://doi.org/10.1021/ac303221j) .\
\
Koch BP, Dittmar T (2006). “From mass to structure: an aromaticity index
for high-resolution mass data of natural organic matter.” *Rapid
Communications in Mass Spectrometry*, **20**(5), 926–932.
[doi:10.1002/rcm.2386](https://doi.org/10.1002/rcm.2386) .\
\
Kujawinski EB, Behn MD (2006). “Automated Analysis of Electrospray
Ionization Fourier Transform Ion Cyclotron Resonance Mass Spectra of
Natural Organic Matter.” *Analytical Chemistry*, **78**(13), 4363–4373.
[doi:10.1021/ac0600306](https://doi.org/10.1021/ac0600306) .

Conway JR, Lex A, Gehlenborg N (2017). “UpSetR: an R package for the
visualization of intersecting sets and their properties.”
*Bioinformatics*, **33**(18), 2938-2940.
[doi:10.1093/bioinformatics/btx364](https://doi.org/10.1093/bioinformatics/btx364)
. <http://dx.doi.org/10.1093/bioinformatics/btx364>.\
\
Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H (2014). “UpSet:
Visualization of Intersecting Sets.” *IEEE Transactions on Visualization
and Computer Graphics*, **20**(12), 1983–1992.
[doi:10.1109/tvcg.2014.2346248](https://doi.org/10.1109/tvcg.2014.2346248)
.

## See also

[`formulas-class`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
and
[`compounds-class`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

The derived
[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
and
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
classes.
