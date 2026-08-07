# Compounds cluster class

Objects from this class are used to store hierarchical clustering data
of candidate structures within
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)
objects.

## Usage

``` r
# S4 method for class 'compoundsCluster'
clusters(obj)

# S4 method for class 'compoundsCluster'
cutClusters(obj)

# S4 method for class 'compoundsCluster'
clusterProperties(obj)

# S4 method for class 'compoundsCluster'
groupNames(obj)

# S4 method for class 'compoundsCluster'
length(x)

# S4 method for class 'compoundsCluster'
lengths(x, use.names = TRUE)

# S4 method for class 'compoundsCluster'
show(object)

# S4 method for class 'compoundsCluster,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'compoundsCluster'
treeCut(obj, k = NULL, h = NULL, groupName)

# S4 method for class 'compoundsCluster'
treeCutDynamic(obj, maxTreeHeight, deepSplit, minModuleSize, groupName)

# S4 method for class 'compoundsCluster,missing'
plot(
  x,
  ...,
  groupName,
  pal = "Paired",
  colourBranches = lengths(x)[groupName] < 50,
  showLegend = lengths(x)[groupName] < 20
)

# S4 method for class 'compoundsCluster'
getMCS(obj, groupName, cluster)

# S4 method for class 'compoundsCluster'
plotStructure(
  obj,
  groupName,
  cluster,
  width = 500,
  height = 500,
  withTitle = TRUE
)

# S4 method for class 'compoundsCluster'
plotSilhouettes(obj, kSeq, groupName, pch = 16, type = "b", ...)
```

## Arguments

- obj, x, object:

  A `compoundsCluster` object.

- use.names:

  A logical value specifying whether the returned vector should be named
  with the feature group names.

- i:

  For `[`: A numeric or character value which is used to select feature
  groups by their index or name, respectively (for the order/names see
  [`groupNames()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).
  Can also be logical to perform logical selection (similar to regular
  vectors). If missing all feature groups are selected.

- ...:

  Further arguments passed directly to the plotting function (`plot` or
  [`plot.dendrogram`](https://rdrr.io/r/stats/dendrogram.html)).

- drop, j:

  ignored.

- k, h:

  Desired number of clusters or tree height to be used for cutting the
  dendrogram, respecitively. One or the other must be specified.
  Analogous to [`cutree`](https://rdrr.io/r/stats/cutree.html).

- groupName:

  A character specifying the feature group name.

- maxTreeHeight, deepSplit, minModuleSize:

  Arguments used by `cutreeDynamicTree`.

- pal:

  Colour palette to be used from RColorBrewer.

- colourBranches:

  Whether branches from cut clusters (and their labels) should be
  coloured. Might be slow with large numbers of clusters, hence, the
  default is only `TRUE` when this is not the case.

- showLegend:

  If `TRUE` and `colourBranches` is also `TRUE` then a legend will be
  shown which outlines cluster numbers and their colours. By default
  `TRUE` for small amount of clusters to avoid overflowing the plot.

- cluster:

  A numeric value specifying the cluster.

- width, height:

  The dimensions (in pixels) of the raster image that should be plotted.

- withTitle:

  A logical value specifying whether a title should be added.

- kSeq:

  An integer vector containing the sequence that should be used for
  average silhouette width calculation.

- pch, type:

  Passed to [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`cutTree` and `cutTreeDynamic` return the modified `compoundsCluster`
object.

`getMCS` returns an [rcdk](https://CRAN.R-project.org/package=rcdk)
molecule object (`IAtomContainer`).

## Details

Objects from this type are returned by the `compounds` method for
[`makeHCluster`](https://rickhelmus.github.io/patRoon/reference/compounds-cluster.md).

## Methods (by generic)

- `clusters(compoundsCluster)`: Accessor method to the `clusters` slot.
  Returns a list that contains for each feature group an object as
  returned by [`hclust`](https://rdrr.io/r/stats/hclust.html).

- `cutClusters(compoundsCluster)`: Accessor method to the `cutClusters`
  slot. Returns a list that contains for each feature group a vector
  with cluster membership for each candidate (format as
  [`cutree`](https://rdrr.io/r/stats/cutree.html)).

- `clusterProperties(compoundsCluster)`: Returns a list with properties
  on how the clustering was performed.

- `groupNames(compoundsCluster)`: returns a `character` vector with the
  names of the feature groups for which data is present in this object.

- `length(compoundsCluster)`: Returns the total number of clusters.

- `lengths(compoundsCluster)`: Returns a `vector` with the number of
  clusters per feature group.

- `show(compoundsCluster)`: Show summary information for this object.

- `x[i`: Subset on feature groups.

- `treeCut(compoundsCluster)`: Manually (re-)cut a dendrogram that was
  generated for a feature group.

- `treeCutDynamic(compoundsCluster)`: Automatically (re-)cut a
  dendrogram that was generated for a feature group using the
  `cutreeDynamicTree` function from dynamicTreeCut.

- `plot(x = compoundsCluster, y = missing)`: Plot the dendrogram for
  clustered compounds of a feature group. Clusters are highlighted using
  [dendextend](https://CRAN.R-project.org/package=dendextend).

- `getMCS(compoundsCluster)`: Calculates the maximum common substructure
  (MCS) for all candidate structures within a specified cluster. This
  method uses the `get.mcs` function from
  [rcdk](https://CRAN.R-project.org/package=rcdk).

- `plotStructure(compoundsCluster)`: Plots the maximum common
  substructure (MCS) for all candidate structures within a specified
  cluster.

- `plotSilhouettes(compoundsCluster)`: Plots the average silhouette
  width when the clusters are cut by a sequence of k numbers. The k
  value with the highest value (marked in the plot) may be considered as
  the optimal number of clusters.

## Slots

- `clusters`:

  A `list` with [`hclust`](https://rdrr.io/r/stats/hclust.html) objects
  for each feature group.

- `dists`:

  A `list` with distance matrices for each feature group.

- `SMILES`:

  A `list` containing a vector with `SMILES` for all candidate
  structures per feature group.

- `cutClusters`:

  A `list` with assigned clusters for all candidates per feature group
  (same format as what [`cutree`](https://rdrr.io/r/stats/cutree.html)
  returns).

- `properties`:

  A list containing general properties and parameters used for
  clustering.
