# Obtain transformation products (TPs) from formula annotation candidates

Transforms and prioritizes [formula annotation
candidates](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
to obtain TPs.

## Usage

``` r
generateTPsAnnForm(
  parents,
  formulas,
  minFitFormula = 0.94,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  parallel = TRUE
)
```

## Arguments

- parents:

  The parents for which transformation products should be obtained. This
  can be

  - a suspect list (see [suspect
    screening](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
    for more information)

  - the output of
    [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md)
    in which case the suspects hits are used as parents

  The parents need to have formula information available.

- formulas:

  The
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
  object containing the formula candidates.

- minFitFormula:

  Minimum `fitFormula` (see Details sections) to filter out unlikely
  candidates.

- skipInvalid:

  If set to `TRUE` then the parents will be skipped (with a warning) for
  which insufficient information (*e.g.* SMILES) is available.

- prefCalcChemProps:

  If `TRUE` then calculated chemical properties such as the formula and
  InChIKey are preferred over what is already present in the parent
  suspect list. For efficiency reasons it is recommended to set this to
  `TRUE`. See the `Validating and calculating chemical properties`
  section for more details.

- neutralChemProps:

  If `TRUE` then the neutral form of the molecule is considered to
  calculate SMILES, formulae etc. Enabling this may improve feature
  matching when considering common adducts (*e.g.* `[M+H]+`, `[M-H]-`).
  See the `Validating and calculating chemical properties` section for
  more details.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

## Value

`generateTPsAnnForm` returns an object of the class
[`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md).
Please see its documentation for *e.g.* filtering steps that can be
performed on this object.

## Details

This function uses formula annotations to obtain transformation
products. This function is called when calling `generateTPs` with
`algorithm="ann_form"`.

The `generateTPsAnnForm` function implements the unknown TP screening
from formula candidates approach as described in (Helmus et al. 2025) .
This algorithm does not rely on any known or predicted TPs and is
therefore suitable for 'full non-target' workflows. *All* formula
candidates are considered as potential TPs and are ranked by the
`TP score`: \$\$TP score = fitFormula + annSim\$\$

With:

- `annSim`: the [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md)

- `fitFormula`: the common element count divided by the total element
  count for the formulae of the parent/TP or TP/parent (maximum is
  taken)

To speed up the calculation process, a threshold for `fitFormula` is
applied to rule out unlikely candidates. The default was derived in
(Helmus et al. 2025) .

Unlike most other TP generation algorithms, no additional suspect
screening step is required.

## Note

Setting `parallel=TRUE` may speed up calculations, but is only favorable
for long calculations due to the overhead of setting up multiple R
processes.

## Validating and calculating chemical properties

Chemical properties such as SMILES, InChIKey and formulae in the parent
suspect list are automatically validated and calculated if
missing/invalid.

The internal validation/calculation process performs the following
steps:

- Validation of SMILES, InChI, InChIKey and formula data (if present).
  Invalid entries will be set to `NA`.

- If `neutralChemProps=TRUE` then chemical data (SMILES, formulae etc.)
  is neutralized by (de-)protonation (using the `–neutralized` option of
  `OpenBabel`). An additional column `molNeutralized` is added to mark
  those molecules that were neutralized. Note that neutralization
  requires either SMILES or InChI data to be available.

- The SMILES and InChI data are used to calculate missing or invalid
  SMILES, InChI, InChIKey and formula data. If `prefCalcChemProps=TRUE`
  then existing InChIKey and formula data is overwritten by calculated
  values whenever possible.

- The chemical formulae which were *not* calculated are verified and
  normalized. This process may be time consuming, and is potentially
  largely avoided by setting `prefCalcChemProps=TRUE`.

- Neutral masses are calculated for missing values
  (`prefCalcChemProps=FALSE`) or whenever possible
  (`prefCalcChemProps=TRUE`).

Note that calculation of formulae for molecules that are isotopically
labelled is currently only supported for deuterium (2H) elements.

This functionality relies heavily on
[OpenBabel](https://github.com/openbabel/openbabel), please make sure it
is installed.

## References

OBoyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR
(2011). “Open Babel: An open chemical toolbox.” *Journal of
Cheminformatics*, **3**(1).
[doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33) .

Helmus R, Bagdonaite I, de Voogt P, van Bommel MR, Schymanski EL, van
Wezel AP, ter Laak TL (2025). “Comprehensive Mass Spectrometry Workflows
to Systematically Elucidate Transformation Processes of Organic
Micropollutants: A Case Study on the Photodegradation of Four
Pharmaceuticals.” *Environmental Science & Technology*, **59**(7),
3723–3736. ISSN 1520-5851.
[doi:10.1021/acs.est.4c09121](https://doi.org/10.1021/acs.est.4c09121) .
<http://dx.doi.org/10.1021/acs.est.4c09121>.

## See also

[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
for more details and other algorithms.
