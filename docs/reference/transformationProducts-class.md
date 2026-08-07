# Base transformation products (TP) class

Holds information for all TPs for a set of parents.

## Usage

``` r
parents(TPs)

products(TPs)

# S4 method for class 'transformationProducts'
parents(TPs)

# S4 method for class 'transformationProducts'
products(TPs)

# S4 method for class 'transformationProducts'
length(x)

# S4 method for class 'transformationProducts'
names(x)

# S4 method for class 'transformationProducts'
show(object)

# S4 method for class 'transformationProducts,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'transformationProducts,ANY,missing'
x[[i, j]]

# S4 method for class 'transformationProducts'
x$name

# S4 method for class 'transformationProducts'
as.data.table(x)

# S4 method for class 'transformationProducts'
convertToSuspects(obj, includeParents = FALSE)

# S4 method for class 'transformationProducts'
delete(obj, i = NULL, j = NULL, ...)

# S4 method for class 'transformationProducts'
filter(obj, properties = NULL, verbose = TRUE, negate = FALSE)
```

## Arguments

- TPs, x, obj, object:

  `transformationProducts` object to be accessed

- i, j:

  For `[`/`[[`: A numeric or character value which is used to select
  parents by their index or name, respectively (for the order/names see
  [`names()`](https://rickhelmus.github.io/patRoon/reference/generics.md)).\
  \
  For `[`: Can also be logical to perform logical selection (similar to
  regular vectors). If missing all parents are selected.\
  \
  For `[[`: should be a scalar value.\
  \
  For `delete`: The data to remove from. `i` are the parents as numeric
  index, logical or character, `j` the transformation products as
  numeric index (row) or name of the TP. If either is `NULL` then data
  for all is removed. `j` may also be a function: it will be called for
  each parent, with the TP info table (a `data.table`), the parent name
  and any other arguments passed as `...` to `delete`. The return value
  of this function specifies the TP indices (rows) (specified as an
  `integer` or `logical` vector) or names to be

- ...:

  For `delete`: passed to the function specified as `j`. Otherwise
  ignored.

- drop:

  ignored.

- name:

  The parent name (partially matched).

- includeParents:

  If `TRUE` then parents are also included in the returned suspect list.

- properties:

  A named `list` with properties to be filtered. Each item in the `list`
  should be named with the name of the property, and should be a vector
  with allowed values. To obtain the possible properties, run *e.g.*
  `names(TPs)[[1]]`. Example:
  `properties=list(likelihood=c("LIKELY","PROBABLE"))`. Set to `NULL` to
  ignore.

- verbose:

  If set to `FALSE` then no text output is shown.

- negate:

  If `TRUE` then filters are performed in opposite manner.

## Value

`delete` returns the object for which the specified data was removed.

`filter` returns a filtered `transformationProducts` object.

## Details

This class holds all generated data for transformation products for a
set of parents. The class is `virtual` and derived objects are created
by [TP
generators](https://rickhelmus.github.io/patRoon/reference/generateTPs.md).

## Methods (by generic)

- `parents(transformationProducts)`: Accessor method for the `parents`
  slot of a `transformationProducts` class.

- `products(transformationProducts)`: Accessor method for the `products`
  slot.

- `length(transformationProducts)`: Obtain total number of
  transformation products.

- `names(transformationProducts)`: Obtain the names of all parents in
  this object.

- `show(transformationProducts)`: Show summary information for this
  object.

- `x[i`: Subset on parents.

- `x[[i`: Extracts a table with TPs for a parent.

- `$`: Extracts a table with TPs for a parent.

- `as.data.table(transformationProducts)`: Returns all TP data in a
  table.

- `convertToSuspects(transformationProducts)`: Converts this object to a
  suspect list that can be used as input for
  [`screenSuspects`](https://rickhelmus.github.io/patRoon/reference/suspect-screening.md).

- `delete(transformationProducts)`: Completely deletes specified
  transformation product data.

- `filter(transformationProducts)`: Performs rule-based filtering.
  Useful to simplify and clean-up the data.

## Slots

- `parents`:

  A `data.table` with metadata for all parents that have TPs in this
  object. Use the `parents` method for access.

- `products`:

  A `list` with `data.table` entries with TP information for each
  parent. Use the `products` method for access.

## S4 class hierarchy

- [`workflowStep`](https://rickhelmus.github.io/patRoon/reference/workflowStep-class.md)

  - **`transformationProducts`**

    - [`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsStructureConsensus`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsCTS`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsAnnComp`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnComp-class.md)

      - [`transformationProductsBT`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

      - [`transformationProductsLibrary`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)

    - [`transformationProductsFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)

      - [`transformationProductsAnnForm`](https://rickhelmus.github.io/patRoon/reference/transformationProductsAnnForm-class.md)

      - [`transformationProductsLibraryFormula`](https://rickhelmus.github.io/patRoon/reference/transformationProductsFormula-class.md)

    - `transformationProductsLogic`

## See also

The derived
[`transformationProductsStructure`](https://rickhelmus.github.io/patRoon/reference/transformationProductsStructure-class.md)
class for more methods,
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md)
and
[`retDir`](https://rickhelmus.github.io/patRoon/reference/retDir.md).
