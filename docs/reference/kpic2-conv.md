# Conversion to KPIC2 objects

Converts a
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
object to an KPIC object.

## Usage

``` r
getPICSet(obj, ...)

# S4 method for class 'features'
getPICSet(
  obj,
  loadRawData = TRUE,
  IMS = FALSE,
  EICParams = getDefEICParams(window = 0)
)

# S4 method for class 'featuresKPIC2'
getPICSet(obj, ...)
```

## Arguments

- obj:

  The `features` object that should be converted.

- ...:

  Ignored

- loadRawData:

  Set to `TRUE` if analyses are available as `mzXML` or `mzML` files.
  Otherwise MS data is not loaded, and some dummy data (*e.g.* file
  paths) is used in the returned object.

- IMS:

  (**IMS workflow**) Specifies which feature groups are considered for
  export in IMS workflows. The following options are valid:

  - `"both"`: Selects IMS and non-IMS features.

  - `"maybe"`: Selects non-IMS features and IMS features without
    assigned IMS precursor.

  - `FALSE`: Selects only non-IMS features.

  - `TRUE`: Selects only IMS features.

  This should be kept `FALSE` as KPIC2 export currently does not support
  IMS features.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

## Details

The conversion process will introduce some dummy values for metadata not
present in patRoon objects. If the `features` object was generated with
KPIC2 and no post mobility assignment was performed, then no conversion
is performed and the original KPIC2 object will be returned.

## Use of raw HRMS data

The [raw data
interface](https://rickhelmus.github.io/patRoon/reference/msdata.md) of
patRoon is used by `getPICSet` to process HRMS (or IMS-HRMS) data.
Please see [its
documentation](https://rickhelmus.github.io/patRoon/reference/msdata.md)
for more information on the supported formats and available
configuration options.
