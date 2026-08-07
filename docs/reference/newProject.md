# Easily create new patRoon projects

The `newProject` function is used to quickly generate a processing R
script. This tool allows the user to quickly select the targeted
analyses, workflow steps and configuring some of their common
parameters. This function requires to be run within a
[RStudio](https://www.rstudio.com/) session. The resulting script is
either added to the current open file or to a new file. The [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
will be written to a `.csv` file so that it can easily be modified
afterwards.

## Usage

``` r
newProject(destPath = NULL)
```

## Arguments

- destPath:

  Set destination path value to this value (useful for debugging). Set
  to `NULL` for a default value.
