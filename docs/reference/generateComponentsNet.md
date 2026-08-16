# Network-based componentization

Uses networking to generate components by grouping features based on
similarity of their chromatographic elution profiles.

## Usage

``` r
generateComponentsNet(fGroups, ...)

# S4 method for class 'featureGroups'
generateComponentsNet(
  fGroups,
  ionization = NULL,
  minSize = 2,
  mzWindow = defaultLim("mz", "medium"),
  componSim = "pearson",
  componMinSim = 0.95,
  componMaxP = 0.05,
  componMethod = "community",
  componArgs = list(),
  groupClust = "complete",
  groupClustH = 0.5,
  annotAlgo = "imss",
  annotAdducts = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M-H]-", "[M-H2O-H]-"),
  annotPrefAdducts = c("[M+H]+", "[M-H]-"),
  annotArgs = list()
)

# S4 method for class 'featureGroupsSet'
generateComponentsNet(fGroups, ionization = NULL, ...)
```

## Source

The componentization approach was inspired by CAMERA and cliqueMS.

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which components should be generated.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- ionization:

  Which ionization polarity was used to generate the data: should be
  `"positive"` or `"negative"`. If the `featureGroups` object has adduct
  annotations, and `ionization=NULL`, the ionization will be detected
  automatically.

  (**sets workflow**) This parameter is not supported for sets
  workflows, as the ionization will always be detected automatically.

- minSize:

  The minimum size of a component. Smaller components than this size
  will be removed. See note below.

- componSim:

  The similarity metric to use: `"pearson"` for Pearson correlation or
  `"cosine"` for cosine similarity.

- componMinSim:

  The minimum similarity threshold (`0-1`). Feature pairs with a lower
  similarity are not connected in the network.

- componMaxP:

  The maximum *p*-value threshold (`0-1`) for Pearson correlation. Only
  applicable when `componSim="pearson"`. Feature pairs with a higher
  *p*-value are not connected.

- componMethod:

  The network componentization method: `"community"` for community
  detection with [igraph](https://CRAN.R-project.org/package=igraph),
  `"cliques"` for maximal cliques (using
  [`igraph::max_cliques`](https://r.igraph.org/reference/cliques.html)),
  `"hcs"` for highly connected subgraphs (using
  [`RBGL::highlyConnSG`](https://rdrr.io/pkg/RBGL/man/highlyConnSG.html)),
  or `"hclust"` for hierarchical clustering (using
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html)).

- componArgs:

  A `list` with additional arguments passed to the network
  componentization function.

  For `componMethod="community"`, arguments are passed to the
  [igraph](https://CRAN.R-project.org/package=igraph) clustering
  function. The default function is
  [`igraph::cluster_walktrap`](https://r.igraph.org/reference/cluster_walktrap.html)
  and can be changed by setting the `func` argument.

  For `componMethod="cliques"`: arguments are passed to
  [`igraph::max_cliques`](https://r.igraph.org/reference/cliques.html).

  For `componMethod="hcs"`: arguments are passed to
  [`RBGL::highlyConnSG`](https://rdrr.io/pkg/RBGL/man/highlyConnSG.html).

  For `componMethod="hclust"`: arguments `method` to set the clustering
  method (see
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html))
  (`"complete" by default`) and `h` to set the clustering height
  threshold (`0.95` by default).

- groupClust, groupClustH:

  Clustering method (see
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html))
  and height at which to cut the tree for the consensus clustering step
  across analyses.

- annotAlgo:

  Annotation algorithm: `"imss"` for annotation using
  [`InterpretMSSpectrum::findMAIN`](https://rdrr.io/pkg/InterpretMSSpectrum/man/findMAIN.html)
  or `"nontarget"` for annotation using
  [`nontarget::pattern.search`](https://rdrr.io/pkg/nontarget/man/pattern.search.html)
  and
  [`nontarget::adduct.search`](https://rdrr.io/pkg/nontarget/man/adduct.search.html).

- annotAdducts:

  A `character` vector or `adduct` object vector with adducts to use for
  annotation. At least two adducts are required. Adducts incompatible
  with the ionization mode are automatically removed.

- annotPrefAdducts:

  A `character` vector with preferential adducts. Only relevant if
  `annotAlgo="nontarget"`, see the `nontarget annotation` section.

- annotArgs:

  A `list` with additional arguments passed to the annotation function.
  For `annotAlgo="imss"`, arguments are passed to
  [`InterpretMSSpectrum::findMAIN`](https://rdrr.io/pkg/InterpretMSSpectrum/man/findMAIN.html).
  For `annotAlgo="nontarget"`, the list should contain `pattern` and
  `adduct` elements, which are themselves lists with arguments passed to
  [`nontarget::pattern.search`](https://rdrr.io/pkg/nontarget/man/pattern.search.html)
  and
  [`nontarget::adduct.search`](https://rdrr.io/pkg/nontarget/man/adduct.search.html),
  respectively.

## Value

A
[`componentsNet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)
object (derived from
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)).
The
[`componentTable`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
of each component contains the following columns:

- `group`:

  The feature group name.

- `ret`:

  The retention time.

- `mz`:

  The mass-to-charge ratio.

- `degreeMin`, `degreeMax`, `degreeMean`:

  The minimum, maximum and mean normalized degree of the feature in the
  network across all analyses.

- `corMin`, `corMax`, `corMean`:

  The minimum, maximum and mean correlation (or similarity) of the
  feature with other features in the component across all analyses.

- `intensity`:

  The mean intensity of the feature group across analyses where it was
  detected.

- `intensity_rel`:

  The relative intensity within the component (`0-1`).

When `annotAlgo="imss"`, the following annotation columns are added:

- `isogroup`:

  The isotope group.

- `isonr`:

  The isotope number within the group.

- `charge`:

  The charge.

- `adduct_ion`:

  The assigned adduct.

- `ppm`:

  The mass error in ppm.

- `neutralMass`:

  The calculated neutral mass.

When `annotAlgo="nontarget"`, the output from nontarget's
`pattern.search` and `adduct.search` is significantly restructured to
simplify interpretation. Isotope and adduct information is collapsed to
one row per feature, the "best" grouping is selected in case of
conflicts (see the **nontarget annotation** section), and links to
related features are expressed as feature group names. The following
annotation columns are added:

- `isogroup`:

  The isotope group ID.

- `iso_interaction`:

  The isotope interaction level(s), i.e. the distance from the
  monoisotope. Multiple values are slash-separated.

- `isotope`:

  The isotope label(s) (e.g. `"13C"`). `"mono"` indicates a monoisotope
  peak. Multiple values are slash-separated.

- `iso_mz_tol`:

  The mass tolerance(s) for the isotope assignment. Multiple values are
  slash-separated.

- `charge`:

  The charge level.

- `iso_link`:

  The feature group name(s) of the origin (monoisotope) peak(s) that
  this feature is an isotope of. `NA` for monoisotope peaks or features
  not part of an isotope pattern.

- `addgroup`:

  The adduct group ID, grouping features that share the same neutral
  mass.

- `adduct_ion`:

  The assigned adduct.

- `neutralMass`:

  The calculated neutral mass.

- `add_link`:

  The feature group name(s) of the peak(s) linked to this feature via an
  adduct relationship.

- `add_link_adduct_ion`:

  The adduct(s) of the linked peak(s).

- `add_link_mz_tol`:

  The mass tolerance(s) for the adduct assignment.

The
[`componentInfo`](https://rickhelmus.github.io/patRoon/reference/components-class.md)
table contains the columns `name` (component name), `cmp_ret` (mean
retention time), `cmp_retsd` (retention time standard deviation),
`neutral_mass` (mean neutral mass) and `size` (number of features).

## Details

This function uses network to generate components. This function is
called when calling `generateComponents` with `algorithm="network"`.

Features are first grouped per analysis into feature components by
constructing a network where nodes represent features and edges
represent similarity between their elution profiles. Similarity is
calculated using either Pearson correlation or cosine similarity. The
network is then clustered using one of several methods (see
`componMethod` argument). The resulting feature components are then
merged across analyses using a consensus clustering approach, based on
the pairwise presence similarity across analyses. Finally, components
are annotated using either the InterpretMSSpectrum or nontarget package
to identify isotopes, in-source fragments and adducts.

Compared to other componentization algorithms supported in patRoon,
`generateComponentsNet()`

- does not need rely on XCMS and uses
  [msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md) for
  EIC extraction, which avoids the need of object conversions and is
  generally faster.

- feature based: this is more sensible where samples are not expected to
  be very similar, e.g. in environmental studies.

- uses a different approach to convert feature components across
  analyses into consensus components, which should (hopefully) perform
  better than other feature-based algorithms
  ([`generateComponentsOpenMS`](https://rickhelmus.github.io/patRoon/reference/generateComponentsOpenMS.md)
  and
  [`generateComponentsCliqueMS`](https://rickhelmus.github.io/patRoon/reference/generateComponentsCliqueMS.md))).

- supports different and configurable componentization methods (see
  `componMethod` argument).

- supports isotope annotation that is not based on 13C differences
  thanks to the nontarget package, eg useful for halogenated compounds.

It is recommended to use the
[`plotGraph()`](https://rickhelmus.github.io/patRoon/reference/generics.md)
method to evaluate the feature componentization, e.g. to evaluate
different methods (`componMethod` argument) or relevant parameters.

## Note

The `generateComponentsNet` function is still experimental. Any feedback
is welcome!

## `nontarget` annotation

When `annotAlgo="nontarget"`, annotation is performed using nontarget's
[`pattern.search`](https://rdrr.io/pkg/nontarget/man/pattern.search.html)
and
[`adduct.search`](https://rdrr.io/pkg/nontarget/man/adduct.search.html)
functions.

**Isotope annotation**: `pattern.search` identifies isotope patterns.
Peaks are grouped by isotope patterns and charge levels. When multiple
charge levels are found for the same isotope group, the "best" grouping
is selected using the following priority: (1) groups containing a 13C
isotope, (2) the largest group, or (3) the lowest charge level.

**Adduct annotation**: `adduct.search` identifies adduct relationships
between peaks. Adducts are grouped by their calculated neutral mass. If
there is a conflict in neutral mass assignment, the "best" adduct group
is selected using the following priority: (1) groups containing
preferred adducts (see `annotPrefAdducts`), (2) the largest group, or
(3) the first group encountered.

## IMS workflows

The componentization algorithm is not aware of the IMS dimension. For
this reason, no IMS feature groups will be considered for
componentization, and direct IMS workflows (see
[`assignMobilitities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
are currently not supported.

## Sets workflows

In a [sets
workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
the componentization is first performed for each set independently. The
resulting components are then all combined in a
[`componentsNetSet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)
object. Note that the components themselves are never merged. The
components are renamed to include the set name from which they were
generated (*e.g.* `"CMP1"` becomes `"CMP1-positive"`).

## References

Kuhl C, Tautenhahn R, Boettcher C, Larson TR, Neumann S (2012). “CAMERA:
an integrated strategy for compound spectra extraction and annotation of
liquid chromatography/mass spectrometry data sets.” *Analytical
Chemistry*, **84**, 283–289.
<http://pubs.acs.org/doi/abs/10.1021/ac202450g>.\
\
Senan O, Aguilar-Mogas A, Navarro M, Capellades J, Noon L, Burks D,
Yanes O, Guimera R, Sales-Pardo M (2019). “CliqueMS: a computational
tool for annotating in-source metabolite ions from LC-MS untargeted
metabolomics data based on a coelution similarity network.”
*Bioinformatics*, **35**(20), 4089–4097.
[doi:10.1093/bioinformatics/btz207](https://doi.org/10.1093/bioinformatics/btz207)
.\
\
Loos M, Singer H (2017). “Nontargeted homologue series extraction from
hyphenated high resolution mass spectrometry data.” *Journal of
Cheminformatics*, **9**(1).
[doi:10.1186/s13321-017-0197-z](https://doi.org/10.1186/s13321-017-0197-z)
.\
\
Loos M, Gerber C, Corona F, Hollender J, Singer H (2015). “Accelerated
Isotope Fine Structure Calculation Using Pruned Transition Trees.”
*Analytical Chemistry*, **87**(11), 5738-5744.
<https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00941>.\
\
Jaeger C, Hoffmann F, Schmitt CA, Lisec J (2016). “Automated Annotation
and Evaluation of In-Source Mass Spectra in GC/Atmospheric Pressure
Chemical Ionization-MS-Based Metabolomics.” *Analytical Chemistry*,
**88**(19), 9386–9390. ISSN 1520-6882.
[doi:10.1021/acs.analchem.6b02743](https://doi.org/10.1021/acs.analchem.6b02743)
. <http://dx.doi.org/10.1021/acs.analchem.6b02743>.\
\
Jaeger C, Méret M, Schmitt CA, Lisec J (2017). “Compound annotation in
liquid chromatography/high‐resolution mass spectrometry based
metabolomics: robust adduct ion determination as a prerequisite to
structure prediction in electrospray ionization mass spectra.” *Rapid
Communications in Mass Spectrometry*, **31**(15), 1261–1266. ISSN
1097-0231. [doi:10.1002/rcm.7905](https://doi.org/10.1002/rcm.7905) .
<http://dx.doi.org/10.1002/rcm.7905>.\
\
Antonov M, Csárdi G, Horvát S, Müller K, Nepusz T, Noom D, Salmon M,
Traag V, Welles BF, Zanini F (2023). “igraph enables fast and robust
network analysis across programming languages.” *arXiv preprint
arXiv:2311.10260*.
[doi:10.48550/arXiv.2311.10260](https://doi.org/10.48550/arXiv.2311.10260)
.\
\
Csárdi G, Nepusz T (2006). “The igraph software package for complex
network research.” *InterJournal*, **Complex Systems**, 1695.
<https://igraph.org>.\
\
Csárdi G, Nepusz T, Traag V, Horvát S, Zanini F, Noom D, Müller K,
Schoch D, Salmon M (2026). *igraph: Network Analysis and Visualization
in R*.
[doi:10.5281/zenodo.7682609](https://doi.org/10.5281/zenodo.7682609) . R
package version 2.3.3, <https://CRAN.R-project.org/package=igraph>.\
\
Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python.” *Journal of Statistical
Software*, **53**(9), 1–18.
[doi:10.18637/jss.v053.i09](https://doi.org/10.18637/jss.v053.i09) .

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
for more details and other algorithms.

[`componentsNet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)
