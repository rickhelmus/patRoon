# Class containing optimization results.

Objects from this class contain optimization results resulting from
design of experiment (DoE).

## Usage

``` r
optimizedParameters(object, paramSet = NULL, DoEIteration = NULL)

optimizedObject(object, paramSet = NULL)

scores(object, paramSet = NULL, DoEIteration = NULL)

experimentInfo(object, paramSet, DoEIteration)

# S4 method for class 'optimizationResult'
algorithm(obj)

# S4 method for class 'optimizationResult'
length(x)

# S4 method for class 'optimizationResult'
lengths(x, use.names = FALSE)

# S4 method for class 'optimizationResult'
show(object)

# S4 method for class 'optimizationResult,missing'
plot(
  x,
  paramSet,
  DoEIteration,
  paramsToPlot = NULL,
  maxCols = NULL,
  type = "contour",
  image = TRUE,
  contours = "colors",
  ...
)

# S4 method for class 'optimizationResult'
optimizedParameters(object, paramSet = NULL, DoEIteration = NULL)

# S4 method for class 'optimizationResult'
optimizedObject(object, paramSet = NULL)

# S4 method for class 'optimizationResult'
scores(object, paramSet = NULL, DoEIteration = NULL)

# S4 method for class 'optimizationResult'
experimentInfo(object, paramSet, DoEIteration)
```

## Arguments

- paramSet:

  Numeric index of the parameter set (*i.e.* the first parameter set
  gets index `1`). For some methods optional: if `NULL` the best will be
  selected.

- DoEIteration:

  Numeric index specifying the DoE iteration within the specified
  `paramSet`. For some methods optional: if `NULL` the best will be
  selected.

- obj, x, object:

  An `optimizationResult` object.

- use.names:

  Ignored.

- paramsToPlot:

  Which parameters relations should be plot. If `NULL` all will be plot.
  Alternatively, a `list` containing one or more `character` vectors
  specifying each two parameters that should be plotted. Finally, if
  only one pair should be plotted, can be a `character` vector
  specifying both parameters.

- maxCols:

  Multiple parameter pairs are plotted in a grid. The maximum number of
  columns can be set with this argument. Set to `NULL` for no limit.

- type:

  The type of plots to be generated: `"contour"`, `"image"` or
  `"persp"`. The equally named functions will be called for plotting.

- image:

  Passed to [`contour`](https://rdrr.io/r/graphics/contour.html) (if
  `type="contour"`).

- contours:

  Passed to [`persp`](https://rdrr.io/r/graphics/persp.html) (if
  `type="persp"`).

- ...:

  Further arguments passed to
  [`contour`](https://rdrr.io/r/graphics/contour.html),
  [`image`](https://rdrr.io/r/graphics/image.html) or
  [`persp`](https://rdrr.io/r/graphics/persp.html) (depending on
  `type`).

## Details

Objects from this class are returned by
[`optimizeFeatureFinding`](https://rickhelmus.github.io/patRoon/reference/feature-optimization.md)
and
[`optimizeFeatureGrouping`](https://rickhelmus.github.io/patRoon/reference/feature-optimization.md).

## Methods (by generic)

- `algorithm(optimizationResult)`: Returns the algorithm that was used
  for finding features.

- `length(optimizationResult)`: Obtain total number of experimental
  design iterations performed.

- `lengths(optimizationResult)`: Obtain number of experimental design
  iterations performed for each parameter set.

- `show(optimizationResult)`: Shows summary information for this object.

- `plot(x = optimizationResult, y = missing)`: Generates response plots
  for all or a selected set of parameters.

- `optimizedParameters(optimizationResult)`: Returns parameter set
  yielding optimal results. The `paramSet` and `DoEIteration` arguments
  can be `NULL`.

- `optimizedObject(optimizationResult)`: Returns the object (*i.e.* a
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
  or
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object) that was generated with optimized parameters. The `paramSet`
  argument can be `NULL`.

- `scores(optimizationResult)`: Returns optimization scores. The
  `paramSet` and `DoEIteration` arguments can be `NULL`.

- `experimentInfo(optimizationResult)`: Returns a `list` with
  optimization information from an DoE iteration.

## Slots

- `algorithm`:

  A character specifying the algorithm that was optimized.

- `paramSets`:

  A `list` with detailed results from each parameter set that was
  tested.

- `bestParamSet`:

  Numeric index of the parameter set yielding the best response.

## Examples

``` r
if (FALSE) { # \dontrun{
# ftOpt is an optimization object.

# plot contour of all parameter pairs from the first parameter set/iteration.
plot(ftOpt, paramSet = 1, DoEIteration = 1)

# as above, but only plot two parameter pairs
plot(ftOpt, paramSet = 1, DoEIteration = 1,
     paramsToPlot = list(c("mzPPM", "chromFWHM"), c("chromFWHM", "chromSNR")))

# plot 3d perspective plots
plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "persp")
} # }
```
