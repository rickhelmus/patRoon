# Loading of MS library data

Loads, parses, verifies and curates MS library data, *e.g.* obtained
from MassBank.

## Usage

``` r
loadMSLibrary(file, algorithm, ...)
```

## Arguments

- file:

  A `character` string that specifies the the file path to the library.

- algorithm:

  A character string describing the algorithm that should be used:
  `"msp"`, `"json"`

- ...:

  Any parameters to be passed to the selected MS library loading
  algorithm.

## Value

A
[`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
object containing the loaded library data.

## Details

`loadMSLibrary` is a generic function that will loads MS library data by
one of the supported algorithms. The actual functionality is provided by
algorithm specific functions such as `loadMSLibraryMSP` and
`loadMSLibraryMoNAJSON`. While these functions may be called directly,
`loadMSLibrary` provides a generic interface and is therefore usually
preferred.

## See also

The
[`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
output class and its methods and the algorithm specific functions:
[`loadMSLibraryMSP`](https://rickhelmus.github.io/patRoon/reference/loadMSLibraryMSP.md),
[`loadMSLibraryMoNAJSON`](https://rickhelmus.github.io/patRoon/reference/loadMSLibraryMoNAJSON.md)
