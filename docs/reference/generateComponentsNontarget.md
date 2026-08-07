# Componentization of homologous series with nontarget

Uses [the nontarget R
package](https://cran.r-project.org/web/packages/nontarget/index.html)
to generate components by unsupervised detection of homologous series.

## Usage

``` r
generateComponentsNontarget(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsNontarget(
  fGroups,
  ionization = NULL,
  rtRange = c(-120, 120),
  mzRange = c(5, 120),
  elements = c("C", "H", "O"),
  rtDev = defaultLim("retention", "wide"),
  absMzDev = defaultLim("mz", "narrow"),
  absMzDevLink = defaultLim("mz", "medium"),
  traceHack = all(R.Version()[c("major", "minor")] >= c(3, 4)),
  ...
)

# S4 method for class 'featureGroupsSet'
generateComponentsNontarget(fGroups, ionization = NULL, ...)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which components should be generated.

- ...:

  Any further arguments passed to
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html).\
  \
  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- ionization:

  Which ionization polarity was used to generate the data: should be
  `"positive"` or `"negative"`. If the `featureGroups` object has adduct
  annotations, and `ionization=NULL`, the ionization will be detected
  automatically.

  (**sets workflow**) This parameter is not supported for sets
  workflows, as the ionization will always be detected automatically.

- rtRange:

  A numeric vector containing the minimum and maximum retention time (in
  seconds) between homologues. Series are always considered from low to
  high *m/z*, thus, a negative minimum retention time allows detection
  of homologous series with increasing *m/z* and decreasing retention
  times. These values set the `minrt` and `maxrt` arguments of
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html).

- mzRange:

  A numeric vector specifying the minimum and maximum *m/z* increment of
  a homologous series. Sets the `minmz` and `maxmz` arguments of
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html).

- elements:

  A character vector with elements to be considered for detection of
  repeating units. Sets the `elements` argument of
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html)
  function.

- rtDev:

  Maximum retention time deviation. Sets the `rttol` to
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html).

- absMzDev:

  Maximum absolute *m/z* deviation. Sets the `mztol` argument to
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html)

- absMzDevLink:

  Maximum absolute *m/z* deviation when linking series. This should
  usually be a bit higher than `absMzDev` to ensure proper linkage.

- traceHack:

  Currently
  [`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html)
  does not work with R `>3.3.3`. This flag, which is enabled by default
  on these R versions, implements a (messy) workaround ([more details
  here](https://github.com/blosloos/nontarget/issues/6)).

## Value

The generated comnponents are returned as an object from the
[`componentsNT`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)
class.

## Details

This function uses nontarget to generate components. This function is
called when calling `generateComponents` with `algorithm="nontarget"`.

In the first step the
[`homol.search`](https://rdrr.io/pkg/nontarget/man/homol.search.html)
function is used to detect all homologous series within each replicate
(analyses within each replicate are averaged prior to detection). Then,
homologous series across replicates are merged in case of full overlap
or when merging of partial overlapping series causes no conflicts.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the componentization is first performed for each set independently. The
resulting components are then all combined in a
[`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)
object. Note that the components themselves are never merged. The
components are renamed to include the set name from which they were
generated (*e.g.* `"CMP1"` becomes `"CMP1-positive"`).

The output class supports additional methods such as `plotGraph`.

## IMS workflows

The componentization algorithm is not aware of the IMS dimension. For
this reason, no IMS feature groups will be considered for
componentization, and direct IMS workflows (see
[`assignMobilitities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
are currently not supported.

## References

Loos M, Singer H (2017). “Nontargeted homologue series extraction from
hyphenated high resolution mass spectrometry data.” *Journal of
Cheminformatics*, **9**(1).
[doi:10.1186/s13321-017-0197-z](https://doi.org/10.1186/s13321-017-0197-z)
.\
\
Loos M, Gerber C, Corona F, Hollender J, Singer H (2015). “Accelerated
Isotope Fine Structure Calculation Using Pruned Transition Trees.”
*Analytical Chemistry*, **87**(11), 5738-5744.
<https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00941>.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.
