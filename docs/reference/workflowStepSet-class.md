# (Virtual) base class for sets related workflow objects

This class is the base for many [sets
workflows](https://rickhelmus.github.io/patRoon/reference/sets-workflow.md)
related classes. This class is virtual, and therefore never created
directly.

## Usage

``` r
# S4 method for class 'workflowStepSet'
setObjects(obj)

# S4 method for class 'workflowStepSet'
sets(obj)

# S4 method for class 'workflowStepSet'
show(object)
```

## Arguments

- obj, object:

  An object that is derived from `workflowStepSet`.

## Details

The most important purpose of this class is to hold data that is
specific for a set. These *set objects* are typically objects with
classes from a regular non-sets workflow (*e.g.*
[`components`](https://rickhelmus.github.io/patRoon/reference/components-class.md),
[`compounds`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)),
and are used by the sets workflow object to *e.g.* form a consensus.
Since the set objects may contain additional data, such as algorithm
specific slots, it may in some cases be of interest to access them
directly with the `setObjects` method (described below).

## Methods (by generic)

- `setObjects(workflowStepSet)`: Accessor for the `setObjects` slot.

- `sets(workflowStepSet)`: Returns the names for each set in this
  object.

- `show(workflowStepSet)`: Shows summary information for this object.

## Slots

- `setObjects`:

  A `list` with the *set objects* (see the `Details` section). The
  `list` is named with the set names.

## S4 class hierarchy

- **`workflowStepSet`**

  - [`componentsSet`](https://rickhelmus.github.io/patRoon/reference/components-class.md)

    - [`componentsNTSet`](https://rickhelmus.github.io/patRoon/reference/componentsNT-class.md)

  - [`compoundsSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

    - [`compoundsConsensusSet`](https://rickhelmus.github.io/patRoon/reference/compounds-class.md)

  - [`formulasSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

    - [`formulasConsensusSet`](https://rickhelmus.github.io/patRoon/reference/formulas-class.md)

  - [`MSPeakListsSet`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
