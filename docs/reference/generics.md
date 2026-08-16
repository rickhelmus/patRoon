# Miscellaneous generics

Various (S4) generic functions providing a common interface for common
tasks such as plotting and filtering data. The actual functionality and
function arguments are often specific for the implemented methods, for
this reason, please refer to the linked method documentation for each
generic.

## Usage

``` r
adducts(obj, ...)

adducts(obj, ...) <- value

algorithm(obj)

analysisInfo(obj, df = FALSE)

analysisInfo(obj) <- value

analyses(obj)

annotatedPeakList(obj, ...)

annotations(obj, ...)

assignMobilities(obj, ...)

calculatePeakQualities(
  obj,
  weights = NULL,
  flatnessFactor = 0.05,
  featureQualities = NULL,
  featureGroupQualities = NULL,
  ...
)

clusterProperties(obj)

clusters(obj)

consensus(obj, ...)

convertToMFDB(TPs, out, ...)

convertToSuspects(obj, ...)

cutClusters(obj)

defaultExclNormScores(obj)

export(obj, type, out, ...)

featureTable(obj, ...)

filter(obj, ...)

fromIMS(obj)

getBPCs(obj, ...)

getFeatures(obj)

getFeatureQualityNames(obj, ...)

getMCS(obj, ...)

getTICs(obj, ...)

groupNames(obj)

hasIMS(obj)

plotBPCs(obj, ...)

plotChord(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...)

plotChroms(obj, ...)

plotChroms3D(obj, ...)

plotGraph(obj, ...)

plotInt(obj, ...)

plotScores(obj, ...)

plotSilhouettes(obj, kSeq, ...)

plotSpectrum(obj, ...)

plotStructure(obj, ...)

plotTICs(obj, ...)

plotVenn(obj, ...)

plotUpSet(obj, ...)

predictRespFactors(obj, ...)

predictTox(obj, ...)

delete(obj, ...)

plotVolcano(obj, ...)

replicates(obj)

setObjects(obj)

sets(obj)

treeCut(obj, k = NULL, h = NULL, ...)

treeCutDynamic(
  obj,
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1,
  ...
)

unset(obj, set)
```

## Arguments

- obj:

  The object the generic should be applied to.

- ...:

  Any further method specific arguments. See method documentation for
  details.

- value:

  The replacement value.

- df:

  If `TRUE` then a `data.frame` is returned, otherwise a `data.table` is
  returned.

- weights, flatnessFactor, featureQualities, featureGroupQualities:

  See method documentation.

- TPs:

  The
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)
  derived object.

- out:

  Output file.

- type:

  The export type.

- addSelfLinks:

  If `TRUE` then 'self-links' are added which represent non-shared data.

- addRetMzPlots:

  Set to `TRUE` to enable *m/z* *vs* retention time scatter plots.

- kSeq:

  An integer vector containing the sequence that should be used for
  average silhouette width calculation.

- k, h:

  Desired numbers of clusters. See
  [`cutree`](https://rdrr.io/r/stats/cutree.html).

- maxTreeHeight, deepSplit, minModuleSize:

  Arguments used by `cutreeDynamicTree`.

- set:

  The name of the set.

## Details

`adducts` returns assigned adducts of the object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

`adducts<-` sets adducts of the object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

`algorithm` returns the algorithm that was used to generate the object.

- Methods are defined for:
  [`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md);
  [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md).

`analysisInfo` returns the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
of an object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

`analysisInfo<-` modifies the [analysis
information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
of an object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsXCMS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresXCMS`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`analyses` returns a `character` vector with the analyses for which data
is present in this object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

`annotatedPeakList` returns an annotated MS peak list.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md).

`annotations` returns annotations.

- Methods are defined for:
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md).

`assignMobilities` assigns ion mobility and/or CCS values to workflow
data.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_comp.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_comp.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/assignMobilities_feat.md).

`calculatePeakQualities` calculates chromatographic peak qualities and
scores.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`clusterProperties` Obtain a list with properties of the generated
cluster(s).

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`clusters` Obtain clustering object(s).

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`consensus` combines and merges data from various algorithms to generate
a consensus.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`featureGroupsComparisonSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

`convertToMFDB` Exports the object to a local database that can be used
with `MetFrag`.

- Methods are defined for: .

`convertToSuspects` Converts an object to a suspect list.

- Methods are defined for:
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`cutClusters` Returns assigned cluster indices of a cut cluster.

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`defaultExclNormScores` Returns default scorings that are excluded from
normalization.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md).

`export` exports workflow data to a given format.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md).

`featureTable` returns feature information.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`filter` provides various functionality to do post-filtering of data.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md);
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/feature-filtering.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md);
  [`transformationProductsAnnComp`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnComp-class.md);
  [`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md);
  [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

`fromIMS` returns `TRUE` if the object was directly created from IMS
data.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`getBPCs` gets base peak chromatogram(s).

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`getFeatures` returns the object's
[`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)
object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md).

`getFeatureQualityNames` returns the object's feature quality names.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`getMCS` Calculates the maximum common substructure.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`getTICs` gets total ion chromatogram(s).

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`groupNames` returns a `character` vector with the names of the feature
groups for which data is present in this object.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

`hasIMS` returns `TRUE` if the object has ion mobility values

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`plotBPCs` plots base peak chromatogram(s).

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`plotChord` plots a Chord diagram to assess overlapping data.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md).

`plotChroms` plots extracted ion chromatogram(s).

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md).

`plotChroms3D` plots a three dimensional chromatogram.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md).

`plotGraph` Plots an interactive network graph.

- Methods are defined for:
  [`componentsNet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md);
  [`componentsNetSet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md);
  [`componentsNT`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md);
  [`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md);
  [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md);
  [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

`plotInt` plots the intensity of all contained features.

- Methods are defined for:
  [`componentsIntClust`](https://rickhelmus.github.io/patRoon/reference/componentsIntClust-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md).

`plotScores` plots candidate scorings.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md).

`plotSilhouettes` plots silhouette widths to evaluate the desired
cluster size.

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`plotSpectrum` plots a (annotated) spectrum.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

`plotStructure` plots a chemical structure.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`plotTICs` plots total ion chromatogram(s).

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`plotVenn` plots a Venn diagram to assess unique and overlapping data.

- Methods are defined for:
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

`plotUpSet` plots an UpSet diagram to assess unique and overlapping
data.

- Methods are defined for:
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md).

`predictRespFactors` Prediction of response factors.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md);
  [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/pred-quant.md).

`predictTox` Prediction of toxicity values.

- Methods are defined for:
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md);
  [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/pred-tox.md).

`delete` Deletes results.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsKPIC2`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsXCMS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsXCMS3`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresKPIC2`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresPiek`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresXCMS`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresXCMS3`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`plotVolcano` plots a volcano plot.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md).

`replicates` returns a `character` vector with the replicates for which
data is present in this object.

- Methods are defined for:
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md).

`setObjects` returns the *set objects* of this object. See the
documentation of
[`workflowStepSet`](https://rickhelmus.github.io/patRoon/reference/workflowStepSet-class.md).

- Methods are defined for:
  [`workflowStepSet`](https://rickhelmus.github.io/patRoon/reference/workflowStepSet-class.md).

`sets` returns the names of the sets inside this object. See the
documentation for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).

- Methods are defined for:
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`workflowStepSet`](https://rickhelmus.github.io/patRoon/reference/workflowStepSet-class.md).

`treeCut` Manually cut a cluster.

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`treeCutDynamic` Automatically cut a cluster.

- Methods are defined for:
  [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md).

`unset` Converts this object to a regular non-set object. See the
documentation for [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md).

- Methods are defined for:
  [`componentsNetSet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md);
  [`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md);
  [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compoundsConsensusSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulasConsensusSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md).

## Other generics

Below are methods that are defined for existing generics (*e.g.* defined
in `base`). Please see method specific documentation for more details.

`[` Subsets data within an object.

- Methods are defined for:
  [`components,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsSet,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compoundsCluster,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`compoundsSet,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureAnnotations,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md);
  [`featureGroupsScreening,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsScreeningSet,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsSet,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresSet,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulasSet,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet,ANY,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts,ANY,missing,missing`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`[[` Extract data from an object.

- Methods are defined for:
  [`components,ANY,ANY`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`featureAnnotations,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups,ANY,ANY`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md);
  [`features,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulas,ANY,ANY`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists,ANY,ANY`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts,ANY,missing`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`$` Extract data from an object.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`as.data.table` Converts an object to a table (`data.table`).

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/feature-table.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/feature-table.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/feature-table.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md);
  [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md).

`as.data.frame` Converts an object to a table (`data.frame`).

- Methods are defined for:
  [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md).

`length` Returns the length of an object.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`lengths` Returns the lengths of elements within this object.

- Methods are defined for:
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md).

`names` Return names for this object.

- Methods are defined for:
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsComparison`](https://rickhelmus.github.io/patRoon/reference/featureGroupsComparison-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md).

`plot` Generates a plot for an object.

- Methods are defined for:
  [`componentsClust,missing`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md);
  [`compoundsCluster,missing`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`featureGroups,missing`](https://rickhelmus.github.io/patRoon/reference/feature-plotting.md);
  [`featureGroupsComparison,missing`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md);
  [`optimizationResult,missing`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md).

`show` Prints information about this object.

- Methods are defined for:
  [`adduct`](https://rickhelmus.github.io/patRoon/reference/adduct-class.md);
  `C++Object`;
  [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsFeatures`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md);
  [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`compoundsCluster`](https://rickhelmus.github.io/patRoon/reference/compoundsCluster-class.md);
  [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md);
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md);
  [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md);
  [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md);
  [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md);
  [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md);
  [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md);
  [`optimizationResult`](https://rickhelmus.github.io/patRoon/reference/optimizationResult-class.md);
  [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md);
  [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md);
  [`workflowStepSet`](https://rickhelmus.github.io/patRoon/reference/workflowStepSet-class.md).
