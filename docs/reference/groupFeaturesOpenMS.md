# Group features using OpenMS

Group and align features with OpenMS tools

## Usage

``` r
groupFeaturesOpenMS(feat, ...)

# S4 method for class 'features'
groupFeaturesOpenMS(
  feat,
  rtalign = TRUE,
  QT = FALSE,
  maxAlignRT = defaultLim("retention", "wide"),
  maxAlignMZ = defaultLim("mz", "medium"),
  maxGroupRT = defaultLim("retention", "medium"),
  maxGroupMZ = defaultLim("mz", "medium"),
  extraOptsRT = NULL,
  extraOptsGroup = NULL,
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

- QT:

  If enabled, use `FeatureLinkerUnlabeledQT` instead of
  `FeatureLinkerUnlabeled` for feature grouping.

- maxAlignRT, maxAlignMZ:

  Used for retention alignment. Maximum retention time or m/z difference
  (seconds/Dalton) for feature pairing. Sets
  `-algorithm:pairfinder:distance_RT:max_difference` and
  `-algorithm:pairfinder:distance_MZ:max_difference` otpions,
  respectively.

- maxGroupRT, maxGroupMZ:

  as `maxAlignRT` and `maxAlignMZ`, but for grouping of features. Sets
  `-algorithm:distance_RT:max_difference` and
  `-algorithm:distance_MZ:max_difference` options, respectively.

- extraOptsRT, extraOptsGroup:

  Named `list` containing extra options that will be passed to
  `MapAlignerPoseClustering` or
  `FeatureLinkerUnlabeledQT/FeatureLinkerUnlabeled`, respectively. Any
  options specified here will override any of the above. Example:
  `extraOptsGroup=list("-algorithm:distance_RT:max_difference"=12)`
  (corresponds to setting `maxGroupRT=12`). Set to `NULL` to ignore.

- verbose:

  if `FALSE` then no text output will be shown.

## Value

An object of a class which is derived from
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

## Details

This function uses OpenMS to group features. This function is called
when calling `groupFeatures` with `algorithm="openms"`.

Retention times may be aligned by the
[MapAlignerPoseClustering](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MapAlignerPoseClustering.html)
TOPP tool. Grouping is achieved by either the
[FeatureLinkerUnlabeled](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureLinkerUnlabeled.html)
or
[FeatureLinkerUnlabeledQT](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureLinkerUnlabeledQT.html)
TOPP tools.

## References

Rost HL, Sachsenberg T, Aiche S, Bielow C, Weisser H, Aicheler F,
Andreotti S, Ehrlich H, Gutenbrunner P, Kenar E, Liang X, Nahnsen S,
Nilse L, Pfeuffer J, Rosenberger G, Rurik M, Schmitt U, Veit J, Walzer
M, Wojnar D, Wolski WE, Schilling O, Choudhary JS, Malmstrom L,
Aebersold R, Reinert K, Kohlbacher O (2016). “OpenMS: a flexible
open-source software platform for mass spectrometry data analysis.”
*Nature Methods*, **13**(9), 741–748.
[doi:10.1038/nmeth.3959](https://doi.org/10.1038/nmeth.3959) .\
\
[pugixml](https://pugixml.org/) (via [Rcpp](http://www.rcpp.org/)) is
used to process OpenMS XML output.\
\
Eddelbuettel D (2013). *Seamless R and C++ Integration with Rcpp*.
Springer, New York.
[doi:10.1007/978-1-4614-6868-4](https://doi.org/10.1007/978-1-4614-6868-4)
. ISBN 978-1-4614-6867-7.\
\
Eddelbuettel D, Balamuta J (2018). “Extending R with C++: A Brief
Introduction to Rcpp.” *The American Statistician*, **72**(1), 28-36.
[doi:10.1080/00031305.2017.1375990](https://doi.org/10.1080/00031305.2017.1375990)
.\
\
Eddelbuettel D, François R (2011). “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software*, **40**(8), 1–18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08) .\
\
Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar
I, Bates D, Chambers J (2026). *Rcpp: Seamless R and C++ Integration*. R
package version 1.1.2, <https://www.rcpp.org>.

## See also

[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
for more details and other algorithms.
