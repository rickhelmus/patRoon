# Automatically installs C3SDB

Automatically installs the [C3SDB](https://github.com/dylanhross/c3sdb)
`Python` package

## Usage

``` r
installC3SDB(envname = "patRoon-C3SDB", clearEnv = FALSE, ...)
```

## Arguments

- envname:

  The name of the virtual `Python` environment to install `C3SDB` into.
  Passed to
  [reticulate::py_install](https://rstudio.github.io/reticulate/reference/py_install.html).
  Set to `NULL` to not set any virtual `Python` environment.

- clearEnv:

  Set to `TRUE` to remove the virtual environment if it already exists
  (using
  [reticulate::virtualenv_remove](https://rstudio.github.io/reticulate/reference/virtualenv-tools.html)).
  Ignored if `envname` is `NULL`.

- ...:

  Further arguments passed to
  [`reticulate::py_install`](https://rstudio.github.io/reticulate/reference/py_install.html).

## Details

This function uses
[reticulate](https://CRAN.R-project.org/package=reticulate) to install
the [C3SDB](https://github.com/dylanhross/c3sdb) `Python` package in a
virtual environment.

## References

Ross DH, Cho JH, Xu L (2020). “Breaking Down Structural Diversity for
Comprehensive Prediction of Ion-Neutral Collision Cross Sections.”
*Analytical Chemistry*, **92**(6), 4548–4557. ISSN 1520-6882.
[doi:10.1021/acs.analchem.9b05772](https://doi.org/10.1021/acs.analchem.9b05772)
. <http://dx.doi.org/10.1021/acs.analchem.9b05772>.
