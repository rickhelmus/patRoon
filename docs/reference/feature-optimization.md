# Optimization of feature finding and grouping parameters

Automatic optimization of feature finding and grouping parameters
through Design of Experiments (DoE).

## Usage

``` r
optimizeFeatureGrouping(
  features,
  algorithm,
  ...,
  templateParams = list(),
  paramRanges = list(),
  maxIterations = 50,
  maxModelDeviation = 0.1,
  parallel = TRUE
)

generateFGroupsOptPSet(algorithm, ...)

getDefFGroupsOptParamRanges(algorithm)

optimizeFeatureFinding(
  anaInfo,
  algorithm,
  ...,
  templateParams = list(),
  paramRanges = list(),
  isoIdent = if (algorithm == "openms") "OpenMS" else "IPO",
  checkPeakShape = "none",
  CAMERAOpts = list(),
  maxIterations = 50,
  maxModelDeviation = 0.1,
  parallel = TRUE
)

generateFeatureOptPSet(algorithm, ...)

getDefFeaturesOptParamRanges(algorithm, method = "centWave")
```

## Arguments

- features:

  A
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  object with the features that should be used to optimize grouping.

- algorithm:

  The algorithm used for finding or grouping features (see
  [`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
  and
  [`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)).

- ...:

  One or more lists with parameter sets (see below) (for
  `optimizeFeatureFinding` and `optimizeFeatureGrouping`).
  Alternatively, named arguments that set (and possibly override) the
  parameters that should be returned from `generateFeatureOptPSet` or
  `generateFGroupsOptPSet`.

- templateParams:

  Template parameter set (see below).

- paramRanges:

  A list with vectors containing absolute parameter ranges
  (minimum/maximum) that constrain numeric parameters choosen during
  experiments. See the `getDefFeaturesOptParamRanges` and
  `getDefFGroupsOptParamRanges` functions for defaults. Values should be
  `Inf` when no limit should be used.

- maxIterations:

  Maximum number of iterations that may be performed to find optimimum
  values. Used to restrict neededless long optimization procedures. In
  IPO this was fixed to `50`.

- maxModelDeviation:

  See the `Potential suboptimal results by optimization model` section
  below.

- parallel:

  If set to `TRUE` then code is executed in parallel through the
  [future](https://CRAN.R-project.org/package=future) package. Please
  see the parallelization section in the handbook for more details.

- anaInfo:

  [Analysis info
  table](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  (passed to
  [`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)).

- isoIdent:

  Sets the algorithm used to identify isotopes. Valid values are:
  `"IPO"`, `"CAMERA"` and `"OpenMS"`. The latter can only be used when
  OpenMS is used to find features, and is highly recommended in this
  situation.

- checkPeakShape:

  Additional peak shape checking of isotopes. Only used if
  `isoIdent="IPO"`. Valid values: `"none"`, `"borderIntensity"`,
  `"sinusCurve"` or `"normalDistr"`.

- CAMERAOpts:

  A `list` with additional arguments passed to
  [`CAMERA::findIsotopes`](https://rdrr.io/pkg/CAMERA/man/findIsotopes-methods.html)
  when `isoIdent="CAMERA"`.

- method:

  Method used by XCMS to find features (only if `algorithm="xcms"`).

## Value

The `optimizeFeatureFinding` and `optimizeFeatureGrouping` return their
results in a
[`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md)
object.

## Details

Many different parameters exist that may affect the output quality of
feature finding and grouping. To avoid time consuming manual
experimentation, functionality is provided to largely automate the
optimization process. The methodology, which uses design of experiments
(DoE), is based on the excellent [Isotopologue Parameter Optimization
(IPO) R package](https://github.com/rietho/IPO). The functionality of
this package is directly integrated in patRoon. Some functionality was
added or changed, however, the principle algorithm workings are nearly
identical.

Compared to IPO, the following functionality was added or changed:

- The code was made more generic in order to include support for other
  feature finding/grouping algorithms (*e.g.* OpenMS, enviPick, XCMS3).

- The methodology of `FeatureFinderMetabo` (OpenMS) may be used to find
  isotopes.

- The `maxModelDeviation` parameter was added to potentially avoid
  suboptimal results ([issue discussed
  here](https://github.com/rietho/IPO/issues/61)).

- The use of multiple 'parameter sets' (discussed below) which, for
  instance, allow optimizing qualitative paremeters more easily (see
  `examples`).

- More consistent optimization code for feature finding/grouping.

- More consistent output using S4 classes (*i.e.*
  [`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md)
  class).

- Parallelization is performed via the
  [future](https://CRAN.R-project.org/package=future) package instead of
  BiocParallel. If this is enabled (`parallel=TRUE`) then any
  parallelization supported by the feature finding or grouping algorithm
  is disabled.

## Parameter sets

Which parameters should be optimized is determined by a *parameter set*.
A set is defined by a named `list` containing the minimum and maximum
starting range for each parameter that should be tested. For instance,
the set `list(chromFWHM = c(5, 10), mzPPM = c(5, 15))` specifies that
the `chromFWHM` and `mzPPM` parameters (used by OpenMS feature finding)
should be optimized within a range of `5`-`10` and `5`-`15`,
respectively. Note that this range may be increased or decreased after a
DoE iteration in order to find a better optimum. The absolute limits are
controlled by the `paramRanges` function argument.

Multiple parameter sets may be specified (*i.e.* through the ...
function argument). In this situation, the optimization algorithm is
repeated for each set, and the final optimum is determined from the
parameter set with the best response. The `templateParams` function
argument may be useful in this case to define a template for each
parameter set. Actual parameter sets are then constructed by joining
each parameter set with the set specified for `templateParams`. When a
parameter is defined in both a regular and template set, the parameter
in the regular set takes precedence.

Parameters that should not be optimized but still need to be set for the
feature finding/grouping functions should also be defined in a
(template) parameter set. Which parameters should be optimized is
determined whether its value is specified as a vector range or a single
fixed value. For instance, when a set is defined as
`list(chromFWHM = c(5, 10), mzPPM = 5)`, only the `chromFWHM` parameter
is optimized, whereas `mzPPM` is kept constant at `5`.

Using multiple parameter sets with differing fixed values allows
optimization of qualitative values (see examples below).

The parameters specified in parameter sets are directly passed through
the
[`findFeatures`](https://rickhelmus.github.io/patRoon/reference/findFeatures.md)
or
[`groupFeatures`](https://rickhelmus.github.io/patRoon/reference/groupFeatures.md)
functions. Hence, grouping and retention time alignment parameters used
by XCMS should (still) be set through the `groupArgs` and `retcorArgs`
parameters.

**NOTE:** For XCMS3, which normally uses parameter classes for settings
its options, the parameters must be defined in a named list like any
other algorithm. The set parameters are then used passed to the
constructor of the right parameter class object (e.g. `CentWaveParam`,
`ObiwarpParam`). For grouping/alignment sets, these parameters need to
be specified in nested lists called `groupParams` and `retAlignParams`,
respectively (similar to `groupArgs`/`retcorArgs` for
`algorithm="xcms"`). Finally, the underlying XCMS method to be used
should be defined in the parameter set (*i.e.* by setting the `method`
field for feature parameter sets and the `groupMethod` and
`retAlignMethod` for grouping/aligning parameter sets). See the examples
below for more details.

**NOTE:** Similar to IPO, the `peakwidth` and `prefilter` parameters for
XCMS feature finding should be split in two different values:

- The minimum and maximum ranges for `peakwidth` are optimized by
  setting `min_peakwidth` and `max_peakwidth`, respectively.

- The `k` and `I` parameters contained in `prefilter` are split in
  `prefilter` and `value_of_prefilter`, respectively.

*Similary*, for KPIC2, the following parameters should be split:

- the `width` parameter (feature optimization) is optimized by
  specifying the `min_width` and `max_width` parameters.

- the `tolerance` and `weight` parameters (feature grouping
  optimization) are optimized by setting `mz_tolerance`/`rt_tolerance`
  and `mz_weight`/`rt_weight` parameters, respectively.

## Functions

The `optimizeFeatureFinding` and `optimizeFeatureGrouping` are the
functions to be used to optimize parameters for feature finding and
grouping, respectively. These functions are analogous to
`optimizeXcmsSet` and `optimizeRetGroup` from IPO.

The `generateFeatureOptPSet` and `generateFGroupsOptPSet` functions may
be used to generate a parameter set for feature finding and grouping,
respectively. Some algorithm dependent default parameter optimization
ranges will be returned. These functions are analogous to
`getDefaultXcmsSetStartingParams` and `getDefaultRetGroupStartingParams`
from IPO. However, unlike their IPO counterparts, these functions will
not output default fixed values. The `generateFGroupsOptPSet` will only
generate defaults for density grouping if `algorithm="xcms"`.

The `getDefFeaturesOptParamRanges` and `getDefFGroupsOptParamRanges`
return the default absolute optimization parameter ranges for feature
finding and grouping, respectively. These functions are useful if you
want to set the `paramRanges` function argument.

## Potential suboptimal results by optimization model

After each experiment iteration an optimimum parameter set is found by
generating a model containing the tested parameters and their responses.
Sometimes the actual response from the parameters derived from the model
is actually signficantly lower than expected. When the response is lower
than the maximum reponse found during the experiment, the parameters
belonging to this experimental maximum may be choosen instead. The
`maxModelDeviation` argument sets the maximum deviation in response
between the modelled and experimental maxima. The value is relative: `0`
means that experimental values will always be favored when leading to
improved responses, whereas `1` will effectively disable this procedure
(and return to 'regular' IPO behaviour).

## Source

The code and methodology is a direct adaptation from the [IPO R
package](https://github.com/rietho/IPO).

## References

Libiseller G, Dvorzak M, Kleb U, Gander E, Eisenberg T, Madeo F, Neumann
S, Trausinger G, Sinner F, Pieber T, Magnes C (2015). “IPO: a tool for
automated optimization of XCMS parameters.” *BMC Bioinformatics*,
**16**(1).
[doi:10.1186/s12859-015-0562-8](https://doi.org/10.1186/s12859-015-0562-8)
.

## Examples

``` r
# example data from patRoonData package
dataDir <- patRoonData::exampleDataPath()
anaInfo <- generateAnalysisInfo(dataDir)
anaInfo <- anaInfo[1:2, ] # only focus on first two analyses (e.g. training set)

# optimize mzPPM and chromFWHM parameters
ftOpt <- optimizeFeatureFinding(anaInfo, "openms", list(mzPPM = c(5, 10), chromFWHM = c(4, 8)))

# optimize chromFWHM and isotopeFilteringModel (a qualitative parameter)
ftOpt2 <- optimizeFeatureFinding(anaInfo, "openms",
                                 list(isotopeFilteringModel = "metabolites (5% RMS)"),
                                 list(isotopeFilteringModel = "metabolites (2% RMS)"),
                                 templateParams = list(chromFWHM = c(4, 8)))

# perform grouping optimization with optimized features object
fgOpt <- optimizeFeatureGrouping(optimizedObject(ftOpt), "xcms",
                                 list(groupArgs = list(bw = c(22, 28)),
                                      retcorArgs = list(method = "obiwarp")))

# same, but using the XCMS3 interface
fgOpt2 <- optimizeFeatureGrouping(optimizedObject(ftOpt), "xcms3",
                                  list(groupMethod = "density", groupParams = list(bw = c(22, 28)),
                                       retAlignMethod = "obiwarp"))


# plot contour of first parameter set/DoE iteration
plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "contour")

# generate parameter set with some predefined and custom parameters to be
# optimized.
pSet <- generateFeatureOptPSet("openms", chromSNR = c(3, 9),
                               useSmoothedInts = FALSE)
```
