# Finding features

Automatically find features.

## Usage

``` r
findFeatures(analysisInfo, algorithm, ..., verbose = TRUE)
```

## Arguments

- analysisInfo:

  A `data.frame` (or `data.table`) with [Analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md).

- algorithm:

  A character string describing the algorithm that should be used:
  `"bruker"`, `"openms"`, `"xcms"`, `"xcms3"`, `"envipick"`, `"sirius"`,
  `"kpic2"`, `"safd"`, `"piek"`

- ...:

  Further parameters passed to the selected feature finding algorithms.

- verbose:

  If set to `FALSE` then no text output is shown.

## Value

An object of a class which is derived from
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

## Details

Several functions exist to collect features (*i.e.* retention and MS
information that represent potential compounds) from a set of analyses.
All 'feature finders' return an object derived from the
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
base class. The next step in a general workflow is to group and align
these features across analyses with
[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md).
Note that some feature finders have a plethora of options which
sometimes may have a large effect on the quality of results. Fine-tuning
parameters is therefore important, and the optimum is largely dependent
upon applied analysis methodology and instrumentation.

`findFeatures` is a generic function that will find features by one of
the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `findFeaturesOpenMS` and
`findFeaturesXCMS`. While these functions may be called directly,
`findFeatures` provides a generic interface and is therefore usually
preferred.

## Note

In most cases it will be necessary to centroid your MS input files. The
only exception is `Bruker`, however, you will still need centroided
`mzXML`/`mzML` files for *e.g.* plotting chromatograms. In this case the
centroided MS files should be stored in the same directory as the raw
`Bruker` `.d` files. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## See also

The
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
output class and its methods and the algorithm specific functions:
[`findFeaturesBruker`](https://rickhelmus.github.io/patRoon/reference/findFeaturesBruker.md),
[`findFeaturesOpenMS`](https://rickhelmus.github.io/patRoon/reference/findFeaturesOpenMS.md),
[`findFeaturesXCMS`](https://rickhelmus.github.io/patRoon/reference/findFeaturesXCMS.md),
[`findFeaturesXCMS3`](https://rickhelmus.github.io/patRoon/reference/findFeaturesXCMS3.md),
[`findFeaturesEnviPick`](https://rickhelmus.github.io/patRoon/reference/findFeaturesEnviPick.md),
[`findFeaturesSIRIUS`](https://rickhelmus.github.io/patRoon/reference/findFeaturesSIRIUS.md),
[`findFeaturesKPIC2`](https://rickhelmus.github.io/patRoon/reference/findFeaturesKPIC2.md),
[`findFeaturesSAFD`](https://rickhelmus.github.io/patRoon/reference/findFeaturesSAFD.md),
[`findFeaturesPiek`](https://rickhelmus.github.io/patRoon/reference/findFeaturesPiek.md)
