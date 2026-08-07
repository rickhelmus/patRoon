# Sets workflows

With sets workflows in patRoon a complete non-target (or suspect)
screening workflow is performed with sample analyses that were measured
with different MS methods (typically positive and negative ionization).

## Details

The analyses files that were measured with a different method are
grouped in *sets*. In the most typical case, there is a `"positive"` and
`"negative"` set, for the positively/negatively ionized data,
respectively. However, other distinctions than polarity are also
possible (although currently the chromatographic method should be the
same between sets). A sets workflow is typically initiated with the
[`makeSet`](https://rickhelmus.github.io/patRoon/reference/makeSet.md)
method. The handbook contains much more details about sets workflows.

## See also

[`makeSet`](https://rickhelmus.github.io/patRoon/reference/makeSet.md)
to initiate sets workflows,
[`workflowStepSet`](https://rickhelmus.github.io/patRoon/reference/workflowStepSet-class.md),
the `Sets workflows` sections in other documentation pages and the
patRoon handbook.
