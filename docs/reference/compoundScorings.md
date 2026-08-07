# Scorings terms for compound candidates

Returns an overview of scorings may be applied to rank candidate
compounds.

## Usage

``` r
compoundScorings(
  algorithm = NULL,
  database = NULL,
  includeSuspectLists = TRUE,
  onlyDefault = FALSE,
  includeNoDB = TRUE
)
```

## Arguments

- algorithm:

  The algorithm: `"metfrag"` or `"sirius"`. Set to `NULL` to return all
  scorings.

- database:

  The database for which results should be returned (*e.g.*
  `"pubchem"`). Set to `NULL` to return all scorings.

- includeSuspectLists, onlyDefault, includeNoDB:

  A logical specifying whether scoring terms related to suspect lists,
  default scoring terms and non-database specific scoring terms should
  be included in the output, respectively.

## Value

A `data.frame` with information on which scoring terms are used, what
their algorithm specific name is and other information such as to which
database they apply and short remarks.

## See also

generateCompounds
