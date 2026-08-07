# Group features using KPIC2

Uses the the [KPIC2](https://github.com/hcji/KPIC2) R package for
grouping of features.

## Usage

``` r
groupFeaturesKPIC2(feat, ...)

# S4 method for class 'features'
groupFeaturesKPIC2(
  feat,
  rtalign = TRUE,
  loadRawData = TRUE,
  groupArgs = list(tolerance = c(0.005, 12)),
  alignArgs = list(),
  verbose = TRUE
)

# S4 method for class 'featuresSet'
groupFeaturesKPIC2(
  feat,
  groupArgs = list(tolerance = c(0.005, 12)),
  verbose = TRUE
)
```

## Arguments

- feat:

  The
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object with the features to be grouped.

- ...:

  Further parameters passed to the selected grouping algorithm.

- rtalign:

  Set to `TRUE` to enable retention time alignment.

- loadRawData:

  Set to `TRUE` if analyses are available as `mzXML` or `mzML` files.
  Otherwise MS data is not loaded, and some dummy data (*e.g.* file
  paths) is used in the returned object.

- groupArgs, alignArgs:

  Named `character` vector that may contain extra parameters to be used
  by
  [`KPIC::PICset.group`](https://rdrr.io/pkg/KPIC/man/PICset.group.html)
  and
  [`KPIC::PICset.align`](https://rdrr.io/pkg/KPIC/man/PICset.align.html),
  respectively.

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function uses KPIC2 to group features. This function is called when
calling `groupFeatures` with `algorithm="kpic2"`.

Grouping of features and alignment of their retention times are
performed with the
[`KPIC::PICset.group`](https://rdrr.io/pkg/KPIC/man/PICset.group.html)
and
[`KPIC::PICset.align`](https://rdrr.io/pkg/KPIC/man/PICset.align.html)
functions, respectively.

## Sets workflows

`loadRawData` and arguments related to retention time alignment are
currently not supported for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).

## References

Ji H, Zeng F, Xu Y, Lu H, Zhang Z (2017). “KPIC2: An Effective Framework
for Mass Spectrometry-Based Metabolomics Using Pure Ion Chromatograms.”
*Analytical Chemistry*, **89**(14), 7631–7640.
[doi:10.1021/acs.analchem.7b01547](https://doi.org/10.1021/acs.analchem.7b01547)
.

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for more details and other algorithms.
