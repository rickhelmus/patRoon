# Generate formula with GenForm

Uses [GenForm](https://sourceforge.net/projects/genform/) to generate
chemical formula candidates.

## Usage

``` r
generateFormulasGenForm(fGroups, ...)

# S4 method for class 'featureGroups'
generateFormulasGenForm(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  relMzDev = defaultLim("mz", "narrow_rel"),
  adduct = NULL,
  elements = "CHNOP",
  hetero = TRUE,
  oc = FALSE,
  thrMS = NULL,
  thrMSMS = NULL,
  thrComb = NULL,
  maxCandidates = Inf,
  extraOpts = NULL,
  calculateFeatures = FALSE,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = defaultLim("mz", "narrow"),
  MSMode = "both",
  isolatePrec = TRUE,
  minIMSSpecSim = 0,
  timeout = 120,
  topMost = 50,
  batchSize = 8
)

# S4 method for class 'featureGroupsSet'
generateFormulasGenForm(
  fGroups,
  MSPeakLists,
  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
  relMzDev = defaultLim("mz", "narrow_rel"),
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0,
  setAvgSpecificScores = FALSE
)
```

## Arguments

- fGroups:

  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which formulae should be generated. This should be the same
  or a subset of the object that was used to create the specified
  `MSPeakLists`. In the case of a subset only the remaining feature
  groups in the subset are considered.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- MSPeakLists:

  An
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  object that was generated for the supplied `fGroups`.

- specSimParams:

  A named `list` with parameters that influence the calculation of the
  [annotation
  similarity](https://rickhelmus.github.io/patRoon/reference/id-conf.md).
  See the [spectral similarity
  parameters](https://rickhelmus.github.io/patRoon/reference/specSimParams.md)
  documentation for more details.

- relMzDev:

  Maximum relative deviation between the measured and candidate formula
  *m/z* values (in ppm). Sets the ppm command line option.

- adduct:

  An
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  object (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md)).
  Examples: `"[M-H]-"`, `"[M+Na]+"`. If the `featureGroups` object has
  adduct annotations then these are used if `adducts=NULL`.

  (**sets workflow**) The `adduct` argument is not supported for sets
  workflows, since the adduct annotations will then always be used.

- elements:

  Elements to be considered for formulae calculation. This will heavily
  affects the number of candidates! Always try to work with a minimal
  set by excluding elements you don't expect. Sets the el command line
  option.

- hetero:

  Only consider formulae with at least one hetero atom. Sets the het
  commandline option.

- oc:

  Only consider organic formulae (*i.e.* with at least one carbon atom).
  Sets the oc commandline option.

- thrMS, thrMSMS, thrComb:

  Sets the thresholds for the `GenForm` MS score (`isoScore`), MS/MS
  score (`MSMSScore`) and combined score (`combMatch`). Sets the
  thms/thmsms/thcomb command line options, respectively. Set to `NULL`
  for no threshold.

- maxCandidates:

  If this number of candidates are found then `GenForm` aborts any
  further formula calculations. The number of candidates is determined
  *after* any formula filters, hence, the properties and 'quality' of
  the candidates is influenced by options such as `oc` and `thrMS`
  arguments. Note that this is different than `topMost`, which selects
  the candidates after `GenForm` finished. Sets the max command line
  option. Set to `0` or `Inf` for no maximum.

- extraOpts:

  An optional character vector with any other command line options that
  will be passed to `GenForm`. See the `GenForm options` section for all
  available command line options.

- calculateFeatures:

  If `TRUE` fomulae are first calculated for all features prior to
  feature group assignment (see `Candidate assignment` in
  [`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)).

- featThreshold:

  If `calculateFeatures=TRUE`: minimum presence (`0-1`) of a formula in
  all features before it is considered as a candidate for a feature
  group. For instance, `featThreshold=0.75` dictates that a formula
  should be present in at least 75% of the features inside a feature
  group.

- featThresholdAnn:

  As `featThreshold`, but only considers features with annotations. For
  instance, `featThresholdAnn=0.75` dictates that a formula should be
  present in at least 75% of the features with annotations inside a
  feature group. @param topMost Only keep this number of candidates (per
  feature group) with highest score.

- absAlignMzDev:

  When the group formula annotation consensus is made from feature
  annotations, the *m/z* values of annotated MS/MS fragments may
  slightly deviate from those of the corresponding group MS/MS peak
  list. The `absAlignMzDev` argument specifies the maximum *m/z* window
  used to re-align the mass peaks.

- MSMode:

  Whether formulae should be generated only from MS data (`"ms"`), MS/MS
  data (`"msms"`) or both (`"both"`). Selecting `"both"` will fall back
  to formula calculation with only MS data in case no MS/MS data is
  available.

- isolatePrec:

  Settings used for isolation of precursor mass peaks and their
  isotopes. This isolation is highly important for accurate isotope
  scoring of candidates, as non-relevant mass peaks will dramatically
  decrease the score. The value of `isolatePrec` should either be a
  `list` with parameters (see the
  [`filter method`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
  for `MSPeakLists` for more details), `TRUE` for default parameters or
  `FALSE` for no isolation (*e.g.* when you already performed isolation
  with the `filter` method). The `z` parameter (charge) is automatically
  deduced from the adduct used for annotation (unless
  `isolatePrec=FALSE`), hence any custom `z` setting is ignored.

- minIMSSpecSim:

  (**IMS workflow**) If the spectrum similarity of an IMS feature group
  compared to its IMS precursor (see
  [`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md))
  is at least this value, then the IMS feature group will not be
  subjected to the annotation algorithm and all feature annotation
  properties will be copied from its precursor. This assumes that
  feature annotation is primarily influenced by the MS/MS spectrum, and
  can be used to speed up the feature annotation process. All scorings,
  annotation similarities etc. are copied from the IMS precursor. The
  fragment annotations are also copied (`fragInfo` result column),
  however, these are adjusted based on the peak list data of the IMS
  feature group.

  This argument does not affect the annotation results for MS-only
  formulae.

- timeout:

  Maximum time (in seconds) that a `GenForm` command is allowed to
  execute. If this time is exceeded a warning is emitted and the command
  is terminated. See the notes section for more information on the need
  of timeouts.

- topMost:

  Only keep this number of candidates (per feature group) with highest
  score.

- batchSize:

  Maximum number of `GenForm` commands that should be run sequentially
  in each parallel process. Combining commands with short runtimes (such
  as `GenForm`) can significantly increase parallel performance. For
  more information see
  [`executeMultiProcess`](https://rickhelmus.github.io/patRoon/reference/executeMultiProcess.md).
  Note that this is ignored if patRoon.MP.method="future".

- setThreshold:

  (**sets workflow**) Minimum abundance for a candidate among all sets
  (`0-1`). For instance, a value of `1` means that the candidate needs
  to be present in all the set data.

- setThresholdAnn:

  (**sets workflow**) As `setThreshold`, but only taking into account
  the set data that contain annotations for the feature group of the
  candidate.

- setAvgSpecificScores:

  (**sets workflow**) If `TRUE` then set specific scorings (*e.g.* MS/MS
  match) are also averaged.

## Value

A
[`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)
object containing all generated formulae.

## Details

This function uses genform to generate formula candidates. This function
is called when calling `generateFormulas` with `algorithm="genform"`.

When MS/MS data is available it will be used to score candidate formulae
by presence of 'fitting' fragments.

## Note

This function always sets the exist and oei `GenForm` command line
options.

Formula calculation with `GenForm` may produce an excessive number of
candidates for high *m/z* values (*e.g.* above 600) and/or many
elemental combinations (set by `elements`). In this scenario formula
calculation may need a very long time. Timeouts are used to avoid
excessive computational times by terminating long running commands (set
by the `timeout` argument).

## GenForm options

Below is a list of options (generated by running `GenForm` without
commandline options) which can be set by the `extraOpts` parameter.

    Formula calculation from MS and MS/MS data as described in
    Meringer et al (2011) MATCH Commun Math Comput Chem 65: 259-290
    Usage: GenForm ms=<filename> [msms=<filename>] [out=<filename>]
            [exist[=mv]] [m=<number>] [ion=-e|+e|-H|+H|+Na] [cha=<number>]
            [ppm=<number>] [msmv=ndp|nsse|nsae] [acc=<number>] [rej=<number>]
            [thms=<number>] [thmsms=<number>] [thcomb=<number>]
            [sort[=ppm|msmv|msmsmv|combmv]] [el=<elements> [oc]] [ff=<fuzzy formula>]
            [vsp[=<even|odd>]] [vsm2mv[=<value>]] [vsm2ap2[=<value>]] [hcf] [kfer[=ex]]
            [wm[=lin|sqrt|log]] [wi[=lin|sqrt|log]] [exp=<number>] [oei]
            [dbeexc=<number>] [ivsm2mv=<number>] [vsm2ap2=<number>]
            [oms[=<filename>]] [omsms[=<filename>]] [oclean[=<filename>]]
            [analyze [loss] [intens]] [dbe] [cm] [pc] [sc] [max]
    Explanation:
            ms      : filename of MS data (*.txt)
            msms    : filename of MS/MS data (*.txt)
            out     : output generated formulas
            exist   : allow only molecular formulas for that at least one
                      structural formula exists;overrides vsp, vsm2mv, vsm2ap2;
                      argument mv enables multiple valencies for P and S
            m       : experimental molecular mass (default: mass of MS basepeak)
            ion     : type of ion measured (default: M+H)
            ppm     : accuracy of measurement in parts per million (default: 5)
            msmv    : MS match value based on normalized dot product, normalized
                      sum of squared or absolute errors (default: nsae)
            acc     : allowed deviation for full acceptance of MS/MS peak in ppm
                      (default: 2)
            rej     : allowed deviation for total rejection of MS/MS peak in ppm
                      (default: 4)
            thms    : threshold for the MS match value
            thmsms  : threshold for the MS/MS match value
            thcomb  : threshold for the combined match value
            sort    : sort generated formulas according to mass deviation in ppm,
                      MS match value, MS/MS match value or combined match value
            el      : used chemical elements (default: CHBrClFINOPSSi)
            oc      : only organic compounds, i.e. with at least one C atom
            ff      : overwrites el and oc and uses fuzzy formula for limits of
                      element multiplicities
            het     : formulas must have at least one hetero atom
            vsp     : valency sum parity (even for graphical formulas)
            vsm2mv  : lower bound for valency sum - 2 * maximum valency
                      (>=0 for graphical formulas)
            vsm2ap2 : lower bound for valency sum - 2 * number of atoms + 2
                      (>=0 for graphical connected formulas)
            hcf     : apply Heuerding-Clerc filter
            kfer    : apply Kind-Fiehn element ratio (extended) ranges
            wm      : m/z weighting for MS/MS match value
            wi      : intensity weighting for MS/MS match value
            exp     : exponent used, when wi is set to log
            oei     : allow odd electron ions for explaining MS/MS peaks
            dbeexc  : excess of double bond equivalent for ions
            ivsm2mv : lower bound for valency sum - 2 * maximum valency
                      for fragment ions
            ivsm2ap2: lower bound for valency sum - 2 * number of atoms + 2
                      for fragment ions
            oms     : write scaled MS peaks to output
            omsms   : write weighted MS/MS peaks to output
            oclean  : write explained MS/MS peaks to output
            analyze : write explanations for MS/MS peaks to output
            loss    : for analyzing MS/MS peaks write losses instead of fragments
            intens  : write intensities of MS/MS peaks to output
            dbe     : write double bond equivalents to output
            cm      : write calculated ion masses to output
            pc      : output match values in percent
            sc      : strip calculated isotope distributions
            noref   : hide the reference information
            max     : maximum number of final candidates (0 is no limit)

## Parallelization

generateFormulasGenForm uses multiprocessing to parallelize
computations. Please see the parallelization section in the handbook for
more details and [patRoon
options](https://rickhelmus.github.io/patRoon/reference/patRoon-package.md)
for configuration options.

When `futures` are used for parallel processing
(`patRoon.MP.method="future"`), calculations with `GenForm` are done
with batch mode disabled (see `batchSize` argument), which generally
limit overall performance.

## References

Meringer M, Reinker S, Zhang J, Muller A (2011). “MS/MS Data Improves
Automated Determination of Molecular Formulas by Mass Spectrometry.”
*MATCH Commun. Math. Comput. Chem.*, **65**(2), 259–290.

## See also

[`generateFormulas`](https://rickhelmus.github.io/patRoon/reference/generateFormulas.md)
for more details and other algorithms.
