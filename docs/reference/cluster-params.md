# Clustering parameters

Parameters for clustering data such as mass spectra and mobilograms.

## Details

Different functionality within patRoon uses clustering to group similar
data together, for instance, to average mass spectra. A fast `C++`
backend based on [Rcpp](https://CRAN.R-project.org/package=Rcpp) is used
to perform the clustering.

The clustering can be configured by the `method` and `window` parameter.
The following clustering methods are available:

- `"hclust"`: uses hierarchical clustering to find similar data points
  (using [hclust-cpp](https://github.com/cdalitz/hclust-cpp), which is
  based on the
  [fastcluster](https://CRAN.R-project.org/package=fastcluster)
  package).

- `"distance_point"`: uses a maximum distance between adjacent sorted
  data points to form clusters.

- `"distance_mean"`: uses a maximum distance between the mean of the
  current cluster and the next sorted data point to form clusters.

- `"bin"`: uses a simple binning approach to cluster data points.

The `hclust` method may give more accurate results and was the default
prior to patRoon 3.0, but is more computationally demanding and
generally unsuitable for IMS workflows due to excessive use of RAM. The
`distance_*` methods are now default and suit most cases.

The `window` parameter defines the clustering tolerance. For
`method="hclust"` this corresponds to the cluster height, for
`method="distance_*"` methods this value sets the maximum distance
between compared data and for `method="bin"` it corresponds to the bin
width. Too small windows will prevent clustering close data points
(*e.g.* resulting in split mass peaks in averaged spectra), whereas too
big windows may cluster unrelated data points together (*e.g.* resulting
in mass inaccuracies).

## Source

Averaging of mass spectra was originally based on algorithms from the
[msProcess](https://github.com/zeehio/msProcess) R package (now archived
on CRAN).

## References

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
package version 1.1.2, <https://www.rcpp.org>.\
\
Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python.” *Journal of Statistical
Software*, **53**(9), 1–18.
[doi:10.18637/jss.v053.i09](https://doi.org/10.18637/jss.v053.i09) .\
\
