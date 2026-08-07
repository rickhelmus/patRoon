# Launch EIC/EIM GUI

Launches a Shiny app for generating and plotting EICs and EIMs.

## Usage

``` r
# S4 method for class 'data.frame'
launchEICGUI(obj, suspects = NULL, adduct = NULL)

# S4 method for class 'features'
launchEICGUI(obj)

# S4 method for class 'featureGroups'
launchEICGUI(obj)
```

## Arguments

- obj:

  For data.frame method: analysis information data.frame. For other
  methods: features/featureGroups object.

- suspects:

  Optional suspect list data.frame (data.frame method only).

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`.

## Value

A Shiny app object.

## Details

These tools were originally developed as internal tools to debug and
explore IMS related algorithms. They were mostly generated with LLMs,
hence, may contain nonsense and bugs. Nevertheless, they may be useful
for users to explore EICs and EIMs interactively. Please use with
caution and report any issues.
