# Parameters to aggregate concentrations/toxicity values assigned to feature groups

Parameters that are used by method functions such
[`as.data.table`](https://rickhelmus.github.io/patRoon/reference/feature-table.md)
to aggregate [predicted
concentrations](https://rickhelmus.github.io/patRoon/reference/pred-quant.md)
or
[toxicities](https://rickhelmus.github.io/patRoon/reference/pred-tox.md).

## Usage

``` r
getDefPredAggrParams(all = mean, ...)
```

## Arguments

- all:

  The default aggregation function for all types, *e.g.* `mean`.

- ...:

  optional named arguments that override defaults.

## Details

Multiple concentration or toxicity values may be assigned to a single
feature group. To ease the interpretation and data handling, several
functions aggregate these values prior their use. Aggregation occurs by
the following data:

- The candidate (*i.e.* suspect or annotation candidate). This is mainly
  relevant for sets workflows, where calculations among sets may yield
  different results for the same candidate.

- The prediction type, *e.g.* all values that were obtained from suspect
  or compound annotation data.

- The feature group.

The aggregation of all data first occurs by the same
candidate/type/feature group, then the same type/feature group and
finally for each feature group. This ensures that *e.g.* large numbers
of data points for a prediction type do not bias results.

The `candidateFunc`, `typeFunc` and `groupFunc` parameters specify the
function that should be used to aggregate data. Commonly, functions such
[`mean`](https://rdrr.io/r/base/mean.html),
[`min`](https://rdrr.io/r/base/Extremes.html) or
[`max`](https://rdrr.io/r/base/Extremes.html) can be used here. Note
that the function does not need to handle `NA` values, as these are
removed in advance.

The `preferType` parameters specifies the *preferred* prediction type.
Any values from other prediction types will be ignored unless the
preferred type is not available for a feature group. Valid values are
`"suspect"` (the default), `"compound"` (results from compound
annotation by SMILES), `"SIRIUS_FP"` (results from formula/compound
annotation with `SIRIUS+CSI:FingerID`) or `"none"`.

These parameters should be stored inside a `list`. The
`getDefPredAggrParams` function can be used to generate such parameter
list with defaults.
