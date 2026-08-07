# Fold change calculation

Fold change calculation

## Usage

``` r
getFCParams(replicates, ...)
```

## Arguments

- replicates:

  A `character` vector with the names of the two replicates to be
  compared.

- ...:

  Optional named arguments that override defaults.

## Details

Fold change calculation can be used to easily identify significant
changes between replicates. The calculation process is configured
through a parameter list, which can be constructed with the
`getFCParams` function. The parameter list has the following entries:

- `replicates` the name of the two replicates to compare (taken from the
  `replicates` argument to `getFCParams`).

- `thresholdFC`: the threshold log FC for a feature group to be
  classified as increasing/decreasing.

- `thresholdPV`: the threshold log P for a feature group to be
  significantly different.

- `zeroMethod`,`zeroValue`: how to handle zero values when calculating
  the FC: `add` adds an offset to zero values, `"fixed"` sets zero
  values to a fixed number and `"omit"` removes zero data. The number
  that is added/set by the former two options is defined by `zeroValue`.

- `PVTestFunc`: a function that is used to calculate P values (usually
  using [`t.test`](https://rdrr.io/r/stats/t.test.html)).

- `PVAdjFunc`: a function that is used to adjust P values (usually using
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html))

## See also

[`featureGroups-class`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
and
[feature-plotting](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md)

## Author

The code to calculate and plot Fold change data was created by Bas van
de Velde.
