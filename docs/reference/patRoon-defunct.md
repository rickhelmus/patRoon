# Defunct functions.

These functions do not work anymore and may be updated in the future.
future.

## Usage

``` r
compoundViewer(fGroups, MSPeakLists, compounds)

# S4 method for class 'featureGroups,MSPeakLists,compounds'
compoundViewer(fGroups, MSPeakLists, compounds)
```

## Arguments

- MSPeakLists:

  A
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object.

- compounds:

  A
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
  object.

## Details

The `compoundViewer` method is used to view compound identification
results. It will display available candidate information such as
scorings and identifiers, MS/MS spectra with explained peaks and
chemical structures.
