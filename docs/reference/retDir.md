# Retention order direction

Calculation of the relative retention order between a parent and its
transformation product (TP).

## Details

The relative retention order between a parent and its TP (`retDir`) is
used throughout TP screening workflows for characterization and
prioritization purposes. These are `numeric` values that hint what the
the chromatographic retention order of a TP might be compared to its
parent: a value of `-1` means it will elute earlier, `1` it will elute
later and `0` that there is no significant difference or the direction
is unknown.

For TP data obtained with
[`generateTPs`](https://rickhelmus.github.io/patRoon/reference/generateTPs.md),
the missing `retDir` values are automatically calculated based on the
`log P` difference between the parent and TP. Here, a typical reversed
phase separation is assumed, *i.e.* compounds with (significantly) lower
log P values likely elute earlier. The `minLogPDiff` parameter of the
[TPStructParams](https://rickhelmus.github.io/patRoon/reference/getDefTPStructParams.md)
argument sets the minimum log P difference to be considered significant.

For TP feature candidates that were linked by
[`generateComponentsTPs`](https://rickhelmus.github.io/patRoon/reference/generateComponentsTPs.md),
the `retDir` values are calculated based on the retention time
difference between the parent and TP feature groups. The `minRTDiff`
argument sets the minimum difference to be considered significant.

## References

Helmus R, Bagdonaite I, de Voogt P, van Bommel MR, Schymanski EL, van
Wezel AP, ter Laak TL (2025). “Comprehensive Mass Spectrometry Workflows
to Systematically Elucidate Transformation Processes of Organic
Micropollutants: A Case Study on the Photodegradation of Four
Pharmaceuticals.” *Environmental Science & Technology*, **59**(7),
3723–3736. ISSN 1520-5851.
[doi:10.1021/acs.est.4c09121](https://doi.org/10.1021/acs.est.4c09121) .
<http://dx.doi.org/10.1021/acs.est.4c09121>.
