# Component class

Contains data for feature groups that are related in some way. These
*components* commonly include adducts, isotopes and homologues.

## Usage

``` r
componentTable(obj)

componentInfo(obj)

findFGroup(obj, fGroup)

expandForIMS(obj, ...)

# S4 method for class 'components'
componentTable(obj)

# S4 method for class 'components'
componentInfo(obj)

# S4 method for class 'components'
groupNames(obj)

# S4 method for class 'components'
length(x)

# S4 method for class 'components'
names(x)

# S4 method for class 'components'
show(object)

# S4 method for class 'components,ANY,ANY,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'components,ANY,ANY'
x[[i, j]]

# S4 method for class 'components'
x$name

# S4 method for class 'components'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'components'
as.data.table(x)

# S4 method for class 'components'
expandForIMS(obj, fGroups)

# S4 method for class 'components'
filter(
  obj,
  size = NULL,
  adducts = NULL,
  isotopes = NULL,
  rtIncrement = NULL,
  mzIncrement = NULL,
  checkComponentsSession = NULL,
  negate = FALSE,
  verbose = TRUE
)

# S4 method for class 'components'
findFGroup(obj, fGroup)

# S4 method for class 'components'
plotSpectrum(obj, index, markFGroup = NULL, xlim = NULL, ylim = NULL, ...)

# S4 method for class 'components'
plotChroms(obj, index, fGroups, EICParams = getDefEICParams(window = 5), ...)

# S4 method for class 'components'
consensus(obj, ...)

# S4 method for class 'componentsCamera'
expandForIMS(obj, ...)

# S4 method for class 'componentsFeatures'
show(object)

# S4 method for class 'componentsCliqueMS'
expandForIMS(obj, ...)

# S4 method for class 'componentsSet'
show(object)

# S4 method for class 'componentsSet,ANY,ANY,missing'
x[i, j, ..., sets = NULL, drop = TRUE]

# S4 method for class 'componentsSet'
filter(obj, ..., negate = FALSE, sets = NULL)

# S4 method for class 'componentsSet'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'componentsSet'
consensus(obj, ...)

# S4 method for class 'componentsSet'
expandForIMS(obj, fGroups)

# S4 method for class 'componentsSet'
unset(obj, set)

# S4 method for class 'componentsNT'
expandForIMS(obj, ...)

# S4 method for class 'componentsNTSet'
expandForIMS(obj, ...)

# S4 method for class 'componentsOpenMS'
expandForIMS(obj, ...)

# S4 method for class 'componentsRC'
expandForIMS(obj, ...)
```

## Arguments

- obj, object, x:

  The `component` object.

- fGroup:

  The name (thus a character) of the feature group that should be
  searched for.

- ...:

  For `delete`: passed to the function specified as `j`.

  For `plotChroms`: Further (optional) arguments passed to the
  `plotChroms` method for the
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  class. Note that the `groupBy`, `showPeakArea`, `showFGroupRect` and
  `topMost` arguments cannot be set as these are set by this method.

  For `plotSpectrum`: Further arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

  For `consensus`: `components` objects that should be used to generate
  the consensus.

  For [sets
  workflow](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
  methods: further arguments passed to the base `components` method.

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  components/feature groups by their index or name, respectively (for
  the order/names see `names()/groupNames()`).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all components/feature groups are
  selected.\
  \
  For `[[`: should be a scalar value. `j` is optional.\
  \
  For `delete`: The data to remove from. `i` are the components as
  numeric index, logical or character, `j` the feature groups as numeric
  index/logical (relative to component) or character. If either is
  `NULL` then data for all is removed. `j` may also be a function: it
  will be called for each component, with the component (a
  `data.table`), the component name and any other arguments passed as
  `...` to `delete`. The return value of this function specifies the
  feature groups to be removed (same format as `j`).

- drop:

  ignored.

- name:

  The component name (partially matched).

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object that was used to generate the components.

  For `expandForIMS`, this should be the `featureGroups` object that
  contains the post IMS feature groups.

- size:

  Should be a two sized vector with the minimum/maximum size of a
  component. Set to `NULL` to ignore.

- adducts:

  Remove any feature groups within components that do not match given
  adduct rules. If `adducts` is a logical then only results are kept
  when an adduct is assigned (`adducts=TRUE`) or not assigned
  (`adducts=FALSE`). Otherwise, if `adducts` contains one or more
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md)
  objects (or something that can be converted to it with
  [`as.adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-utils.md))
  then only results are kept that match the given adducts. Set to `NULL`
  to ignore this filter.

- isotopes:

  Only keep results that match a given isotope rule. If `isotopes` is a
  logical then only results are kept with (`isotopes=TRUE`) or without
  (`isotopes=FALSE`) isotope assignment. Otherwise `isotopes` should be
  a numeric vector with isotope identifiers to keep (*e.g.* `0` for
  monoisotopic results, `1` for `M+1` results etc.). Set to `NULL` to
  ignore this filter.

- rtIncrement, mzIncrement:

  Should be a two sized vector with the minimum/maximum retention or mz
  increment of a homologous series. Set to `NULL` to ignore.

- checkComponentsSession:

  If set then components and/or feature groups are removed that were
  selected for removal (see
  [check-GUI](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  and the
  [`checkComponents`](https://rickhelmus.github.io/patRoon/reference/check-GUI.md)
  function). The value of `checkComponentsSession` should either by a
  path to the session file or `TRUE`, in which case the default session
  file name is used. If `negate=TRUE` then all non-selected data is
  removed instead.

- negate:

  If `TRUE` then filters are applied in opposite manner.

- verbose:

  If set to `FALSE` then no text output is shown.

- index:

  The index of the component. Can be a numeric index or a character with
  its name.

- markFGroup:

  If specified (*i.e.* not `NULL`) this argument can be used to mark a
  feature group in the plotted spectrum. The value should be a character
  with the name of the feature group. Setting this to `NULL` will not
  mark any peak.

- xlim, ylim:

  Sets the plot size limits used by
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). Set to `NULL`
  for automatic plot sizing.

- EICParams:

  A named `list` with parameters used for extracted ion chromatogram
  (EIC) creation. See the [EIC
  parameters](https://rickhelmus.github.io/patRoon/reference/EIXParams.md)
  documentation for more details.

- sets:

  (**sets workflow**) A `character` with name(s) of the sets to keep (or
  remove if `negate=TRUE`).

- set:

  (**sets workflow**) The name of the set.

## Value

`delete` returns the object for which the specified data was removed.

`consensus` returns a `components` object that is produced by merging
multiple specified `components` objects.

## Details

`components` objects are obtained from
[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md).

## Methods (by generic)

- `componentTable(components)`: Accessor method for the `components`
  slot of a `components` class. Each component is stored as a
  `data.table`.

- `componentInfo(components)`: Accessor method for the `componentInfo`
  slot of a `components` class.

- `groupNames(components)`: returns a `character` vector with the names
  of the feature groups for which data is present in this object.

- `length(components)`: Obtain total number of components.

- `names(components)`: Obtain the names of all components.

- `show(components)`: Show summary information for this object.

- `x[i`: Subset on components/feature groups.

- `x[[i`: Extracts a component table, optionally filtered by a feature
  group.

- `$`: Extracts a component table by component name.

- `delete(components)`: Completely deletes specified (parts of)
  components.

- `as.data.table(components)`: Returns all component data in a table.

- `expandForIMS(components)`: Expands the components data for IMS
  feature groups. See the `IMS expansion` section below.

- `filter(components)`: Provides rule based filtering for components.

- `findFGroup(components)`: Returns the component id(s) to which a
  feature group belongs.

- `plotSpectrum(components)`: Plot a *pseudo* mass spectrum for a single
  component.

- `plotChroms(components)`: Plot an extracted ion chromatogram (EIC) for
  all feature groups within a single component.

- `consensus(components)`: Generates a consensus from multiple
  `components` objects. At this point results are simply combined and no
  attempt is made to merge similar components.

## Slots

- `components`:

  List of all components in this object. Use the `componentTable` method
  for access.

- `componentInfo`:

  A `data.table` containing general information for each component. Use
  the `componentInfo` method for access.

## Note

`filter` Applies only those filters for which a component has data
available. For instance, filtering by adduct will only filter any
results within a component if that component contains adduct
information.

For `plotChroms`: The `topMost` and `topMostByReplicate` EIC parameters
are ignored unless the components are from homologous series.

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`components`**

    - `componentsCamera`

    - `componentsFeatures`

      - `componentsCliqueMS`

      - `componentsOpenMS`

    - [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md)

      - [`componentsIntClust`](https://rickhelmus.github.io/patRoon/reference/componentsIntClust-class.md)

      - [`componentsSpecClust`](https://rickhelmus.github.io/patRoon/reference/componentsSpecClust-class.md)

    - `componentsSet`

      - [`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

    - `componentsUnset`

    - [`componentsNT`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

      - [`componentsNTUnset`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

    - `componentsRC`

    - [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)

## IMS expansion

In IMS workflows with post mobility assignment (see
[`assignMobilities`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md)),
specifically when the assignment occurs *after* generation of the
components, it may be desired to copy the results of the IMS precursors
to the IMS feature groups. For instance, for components from [intensity
clusters](https://rickhelmus.github.io/patRoon/reference/generateComponentsIntClust.md)
or
[TPs](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
one could assume that results for IMS feature groups will be largely the
same as their IMS precursors. The `expandForIMS` method function is used
to *expand* the original `components` object by adding in IMS feature
groups with data copied from their IMS precursors.

Currently, only components generated with
[`generateComponentsIntClust`](https://rickhelmus.github.io/patRoon/reference/generateComponentsIntClust.md),
[`generateComponentsSpecClust`](https://rickhelmus.github.io/patRoon/reference/generateComponentsSpecClust.md)
and
[`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md)
support this operation.

## Sets workflows

The `componentsSet` class is applicable for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).
This class is derived from `components` and therefore largely follows
the same user interface.

The following methods are specifically defined for sets workflows:

- `unset` Converts the object data for a specified set into a 'non-set'
  object (`componentsUnset`), which allows it to be used in 'regular'
  workflows. Only the components in the specified set are kept.

The following methods are changed or with new functionality:

- `filter` and the subset operator (`[`) Can be used to select
  components that are only present for selected sets.

## See also

[`generateComponents`](https://rickhelmus.github.io/patRoon/reference/generateComponents.md)
