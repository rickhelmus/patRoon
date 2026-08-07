# Imports feature groups from KPIC2

Imports feature groups from KPIC2

## Usage

``` r
importFeatureGroupsKPIC2(input, analysisInfo)
```

## Arguments

- input:

  A grouped `PIC set` object (*e.g.* as returned by
  [`KPIC::PICset.group`](https://rdrr.io/pkg/KPIC/man/PICset.group.html)).

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function imports data from KPIC2. This function is called when
calling `importFeatureGroups` with `type="kpic2"`.

## References

Ji H, Zeng F, Xu Y, Lu H, Zhang Z (2017). “KPIC2: An Effective Framework
for Mass Spectrometry-Based Metabolomics Using Pure Ion Chromatograms.”
*Analytical Chemistry*, **89**(14), 7631–7640.
[doi:10.1021/acs.analchem.7b01547](https://doi.org/10.1021/acs.analchem.7b01547)
.

## See also

[`importFeatureGroups`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroups.md)
for more details and other algorithms.

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
