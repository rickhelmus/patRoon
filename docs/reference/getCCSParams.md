# Parameters for CCS calculation

Configuration and description of parameters used for CCS calculation.

## Usage

``` r
getCCSParams(method, ..., calibrant = NULL)
```

## Source

The calculation formulas was derived from (Haler et al. 2017) , (George
et al. 2024) and the implementation used by
[MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html)
(Tsugawa et al. 2020) . (`MobilityToCrossSection` method from the
`IonMobilityUtility` class).

## Arguments

- method, calibrant:

  Sets the CCS calculation method and the calibrant data (only required
  if `method="agilent"`). See details.

- ...:

  optional named arguments that override defaults.

## Details

The following parameters exist to configure the CCS calculation:

- `method` The CCS calculation method. Should be `"bruker"`,
  `"mason-schamp_k"`, `"mason-schamp_1/k"` or `"agilent"`. See details
  below.

- `defaultCharge` The default charge of the ions. This is used when no
  charge information is available.

- `temperature`,`massGas` The temperature (Kelvin) and exact mass of the
  drift gas. See calculation details below.

- `MasonSchampConstant` The Mason-Schamp constant. See calculation
  details below.

- `calibrant` If `method="agilent"`: the calibrant data to be used for
  CCS calculation. This should either be

  - A path to an Agilent `.d` file.

  - A path to an `OverrideImsCal.xml` file (found in
    `sample.d/AcqData`).

  - A named `list` with the elements `massGas`, `TFix` and `beta`.

The CCS calculation depends on the `method` parameter:

- `bruker`: uses the Bruker `TDF-SDK` for calculations. See
  [msdata](https://rickhelmus.github.io/patRoon/reference/msdata.md) for
  configuration options. Only applicable to TIMS data.

- `mason-schamp_k`: uses the Mason-Schamp equation: \$\$CCS = C \cdot
  \frac{charge}{\sqrt{u \cdot T}} \cdot \frac{1}{mobility}\$\$

  With

  - *C* the Mason-Schamp constant, can be changed by setting the
    `MasonSchampConstant` parameter. See (George et al. 2024) for
    details.

  - *u* the reduced mass of the drift gas and the ion: \$\$u =
    \frac{m\_{gas} \cdot m\_{ion}}{m\_{gas} + m\_{ion}}\$\$ The mass of
    the drift gas is defined by the `massGas` parameter.

  - *T* the temperature (Kelvin) as defined by the `temperature`
    parameter.

- `mason-schamp_1/k`: as `mason-schamp_k` but assuming an inversed
  mobility (\\\frac{1}{k}\\). This is meant for TIMS data. Compared to
  `method="bruker"`, this doesn't rely on the `TDF-SDK` but may produce
  results with very minor differences (George et al. 2024) .

- `agilent`: uses Agilent calibration data with the following equation:

  \$\$CCS = (mobility - t\_{fix}) \cdot \frac{charge}{\beta} \cdot
  \frac{1}{\sqrt{\frac{m\_{ion}}{m\_{ion} + m\_{gas}}}}\$\$

  With \\t\_{fix}\\ and \\\beta\\ the `TFix` and `beta` values from the
  calibration data. The `massGas` parameter sets the \\m\_{gas}\\ value.

The `getCCSParams` function generates such parameter list with defaults.

## References

George AC, Schmitz I, Rouviere F, Alves S, Colsch B, Heinisch S, Afonso
C, Fenaille F, Loutelier-Bourhis C (2024). “Interplatform comparison
between three ion mobility techniques for human plasma lipid collision
cross sections.” *Analytica Chimica Acta*, **1304**, 342535. ISSN
0003-2670.
[doi:10.1016/j.aca.2024.342535](https://doi.org/10.1016/j.aca.2024.342535)
. <http://dx.doi.org/10.1016/j.aca.2024.342535>.\
\
Haler JRN, Kune C, Massonnet P, Comby-Zerbino C, Jordens J, Honing M,
Mengerink Y, Far J, De Pauw E (2017). “Comprehensive Ion Mobility
Calibration: Poly(ethylene oxide) Polymer Calibrants and General
Strategies.” *Analytical Chemistry*, **89**(22), 12076–12086. ISSN
1520-6882.
[doi:10.1021/acs.analchem.7b02564](https://doi.org/10.1021/acs.analchem.7b02564)
. <http://dx.doi.org/10.1021/acs.analchem.7b02564>.\
\
Tsugawa H, Ikeda K, Takahashi M, Satoh A, Mori Y, Uchino H, Okahashi N,
Yamada Y, Tada I, Bonini P, Higashi Y, Okazaki Y, Zhou Z, Zhu Z, Koelmel
J, Cajka T, Fiehn O, Saito K, Arita M, Arita M (2020). “A lipidome atlas
in MS-DIAL 4.” *Nature Biotechnology*, **38**(10), 1159–1163. ISSN
1546-1696.
[doi:10.1038/s41587-020-0531-2](https://doi.org/10.1038/s41587-020-0531-2)
. <http://dx.doi.org/10.1038/s41587-020-0531-2>.
