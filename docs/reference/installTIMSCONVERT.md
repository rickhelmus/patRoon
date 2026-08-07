# Automatically installs TIMSCONVERT

Automatically installs the
[TIMSCONVERT](https://gtluu.github.io/timsconvert/) `Python` package

## Usage

``` r
installTIMSCONVERT(envname = "patRoon-TIMSCONVERT", clearEnv = FALSE, ...)
```

## Arguments

- envname:

  The name of the virtual Python environment to install `TIMSCONVERT`
  into. Passed to
  [reticulate::py_install](https://rstudio.github.io/reticulate/reference/py_install.html).
  Set to `NULL` to not set any virtual Python environment.

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
the [TIMSCONVERT](https://gtluu.github.io/timsconvert/) Python package
in a virtual environment.

## References

Luu GT, Freitas MA, Lizama-Chamu I, McCaughey CS, Sanchez LM, Wang M
(2022). “TIMSCONVERT: a workflow to convert trapped ion mobility data to
open data formats.” *Bioinformatics*, **38**(16), 4046–4047. ISSN
1367-4811.
[doi:10.1093/bioinformatics/btac419](https://doi.org/10.1093/bioinformatics/btac419)
. <http://dx.doi.org/10.1093/bioinformatics/btac419>.
