# Temporarily changes package options

This function is inspired by
[`withr::with_options`](https://withr.r-lib.org/reference/with_options.html):
it can be used to execute some code where package options are
temporarily changed. This function uses a shortened syntax, especially
when changing options for `patRoon`.

## Usage

``` r
withOpt(code, ..., prefix = "patRoon.")
```

## Arguments

- code:

  The code to be executed.

- ...:

  Named arguments with options to change.

- prefix:

  A `character` that will be used to prefix given option names.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set max parallel processes to five while performing formula calculations
withOpt(MP.maxProcs = 5, {
    formulas <- generateFormulas(fGroups, "genform", ...)
})
} # }
```
