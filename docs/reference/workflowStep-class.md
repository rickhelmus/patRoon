# (Virtual) Base class for all workflow objects.

All workflow objects (*e.g.*
[`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md),
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md),
etc) are derived from this class. Objects from this class are never
created directly.

## Usage

``` r
# S4 method for class 'workflowStep'
algorithm(obj)

# S4 method for class 'workflowStep'
as.data.table(x, keep.rownames = FALSE, ...)

# S4 method for class 'workflowStep'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S4 method for class 'workflowStep'
show(object)
```

## Arguments

- obj, x, object:

  An object (derived from) this class.

- keep.rownames:

  Ignored.

- ...:

  Method specific arguments. Please see the documentation of the derived
  classes.

- row.names, optional:

  Ignored.

## Methods (by generic)

- `algorithm(workflowStep)`: Returns the algorithm that was used to
  generate an object.

- `as.data.table(workflowStep)`: Summarizes the data in this object and
  returns this as a `data.table`.

- `as.data.frame(workflowStep)`: This method simply calls
  `as.data.table` and converts the result to a classic a `data.frame`.

- `show(workflowStep)`: Shows summary information for this object.

## Slots

- `algorithm`:

  The algorithm that was used to generate this object. Use the
  `algorithm` method for access.

## S4 class hierarchy

- **`workflowStep`**

  - [`transformationProducts`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)

    - [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsStructureConsensus`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsCTS`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsAnnComp`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnComp-class.md)

      - [`transformationProductsBT`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsLibrary`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

    - [`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)

      - [`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md)

      - [`transformationProductsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)

    - [`transformationProductsLogic`](https://rickhelmus.github.io/patRoon/reference/transformationProducts-class.md)

  - [`features`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresSet`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresUnset`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresFromFeatGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - [`featuresConsensus`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - [`featuresBruker`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresEnviPick`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresKPIC2`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresOpenMS`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresPiek`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresSAFD`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresSIRIUS`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresTable`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresBrukerTASQ`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerTASQ.md)

    - [`featuresXCMS`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

    - [`featuresXCMS3`](https://rickhelmus.github.io/patRoon/reference/features-class.md)

  - [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsSet`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

      - [`featureGroupsScreeningSet`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

    - [`featureGroupsUnset`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsScreening`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

      - [`featureGroupsSetScreeningUnset`](https://rickhelmus.github.io/patRoon/reference/featureGroupsScreening-class.md)

    - [`featureGroupsBruker`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsConsensus`](https://rickhelmus.github.io/patRoon/reference/featureGroups-compare.md)

    - [`featureGroupsEnviMass`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsGreedy`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsIMS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsKPIC2`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsOpenMS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsTable`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsBrukerTASQ`](https://rickhelmus.github.io/patRoon/reference/importFeatureGroupsBrukerTASQ.md)

    - [`featureGroupsXCMS`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

    - [`featureGroupsXCMS3`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)

  - [`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsCamera`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsFeatures`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

      - [`componentsCliqueMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

      - [`componentsOpenMS`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsClust`](https://rickhelmus.github.io/patRoon/reference/componentsClust-class.md)

      - [`componentsIntClust`](https://rickhelmus.github.io/patRoon/reference/componentsIntClust-class.md)

      - [`componentsSpecClust`](https://rickhelmus.github.io/patRoon/reference/componentsSpecClust-class.md)

    - [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

      - [`componentsNetSet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)

      - [`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

    - [`componentsUnset`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsNet`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)

      - [`componentsNetUnset`](https://rickhelmus.github.io/patRoon/reference/componentsNet-class.md)

    - [`componentsNT`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

      - [`componentsNTUnset`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

    - [`componentsRC`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsTPs`](https://rickhelmus.github.io/patRoon/reference/componentsTPs-class.md)

  - [`featureAnnotations`](https://rickhelmus.github.io/patRoon/reference/featureAnnotations-class.md)

    - [`formulas`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasConsensus`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasUnset`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

      - [`formulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/formulasSIRIUS-class.md)

    - [`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsConsensus`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsMF`](https://rickhelmus.github.io/patRoon/reference/compoundsMF-class.md)

      - [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsUnset`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

      - [`compoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/compoundsSIRIUS-class.md)

  - [`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)

    - [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)

    - [`MSPeakListsUnset`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)

  - [`MSLibrary`](https://rickhelmus.github.io/patRoon/reference/MSLibrary-class.md)
