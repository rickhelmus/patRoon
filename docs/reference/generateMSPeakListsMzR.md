# Generate peak lists with mzR (deprecated)

Uses the mzR package to read the MS data needed for MS peak lists. This
function is now deprecated, please use
[`generateMSPeakLists`](https://rickhelmus.github.io/patRoon/reference/generateMSPeakLists.md)
instead.

## Usage

``` r
generateMSPeakListsMzR(fGroups, ...)

# S4 method for class 'featureGroups'
generateMSPeakListsMzR(
  fGroups,
  maxMSRTWindow = 5,
  precursorMzWindow = 4,
  topMost = NULL,
  avgFeatParams = getDefAvgPListParams(clusterMzWindow = 0.005),
  avgFGroupParams = getDefAvgPListParams(clusterMzWindow = 0.005, withPrecursorMS =
    FALSE)
)

# S4 method for class 'featureGroupsSet'
generateMSPeakListsMzR(fGroups, ...)
```

## Arguments

- fGroups:

  The
  [`featureGroups`](https://rickhelmus.github.io/patRoon/reference/featureGroups-class.md)
  object for which MS peak lists should be generated.

- ...:

  (**sets workflow**) Further arguments passed to the non-sets workflow
  method.

- maxMSRTWindow:

  Maximum chromatographic peak window used for spectrum averaging (in
  seconds, +/- retention time). If `NULL` all spectra from a feature
  will be taken into account. Lower to decrease processing time.

- precursorMzWindow:

  The *m/z* window (in Da) to find MS/MS spectra of a precursor. This is
  typically used for Data-Dependent like MS/MS data and should
  correspond to the isolation *m/z* window (*i.e.* +/- the precursor
  *m/z*) that was used to collect the data. For Data-Independent MS/MS
  experiments, where precursor ions are not isolated prior to
  fragmentation (*e.g.* bbCID, MSe, all-ion, ...) the value should be
  `NULL`.

- topMost:

  Only extract MS peak lists from a maximum of `topMost` analyses with
  highest intensity. If `NULL` all analyses will be used.

- avgFeatParams:

  Parameters used for averaging MS peak lists of individual features.
  Analogous to `avgFGroupParams`.

- avgFGroupParams:

  A `list` with parameters used for averaging of peak lists for feature
  groups. See
  [`getDefAvgPListParams`](https://rickhelmus.github.io/patRoon/reference/getDefAvgPListParams.md)
  for more details.

## Value

An
[`MSPeakLists`](https://rickhelmus.github.io/patRoon/reference/MSPeakLists-class.md)
object.

## Details

The MS data files should be either in `.mzXML` or `.mzML` format.

The input MS data files need to be centroided. The
[`convertMSFiles`](https://rickhelmus.github.io/patRoon/reference/MSConversion.md)
function can be used to centroid data.

## References

Chambers, C. M, Maclean, Brendan, Burke, Robert, Amodei, Dario,
Ruderman, L. D, Neumann, Steffen, Gatto, Laurent, Fischer, Bernd, Pratt,
Brian, Egertson, Jarrett, Hoff, Katherine, Kessner, Darren, Tasman,
Natalie, Shulman, Nicholas, Frewen, Barbara, Baker, A. T, Brusniak,
Mi-Youn, Paulse, Christopher, Creasy, David, Flashner, Lisa, Kani, Kian,
Moulding, Chris, Seymour, L. S, Nuwaysir, M. L, Lefebvre, Brent,
Kuhlmann, Frank, Roark, Joe, Rainer, Paape, Detlev, Suckau, Hemenway,
Tina, Huhmer, Andreas, Langridge, James, Connolly, Brian, Chadick, Trey,
Holly, Krisztina, Eckels, Josh, Deutsch, W. E, Moritz, L. R, Katz, E. J,
Agus, B. D, MacCoss, Michael, Tabb, L. D, Mallick, Parag (2012). “A
cross-platform toolkit for mass spectrometry and proteomics.” *Nat
Biotech*, **30**(10), 918–920.
[doi:10.1038/nbt.2377](https://doi.org/10.1038/nbt.2377) .
<http://dx.doi.org/10.1038/nbt.2377>.\
\
Keller A, Eng J, Zhang N, Li X, Aebersold R (2005). “A uniform
proteomics MS/MS analysis platform utilizing open XML file formats.”
*Mol Syst Biol*.\
\
Kessner D, Chambers M, Burke R, Agus D, Mallick P (2008). “ProteoWizard:
open source software for rapid proteomics tools development.”
*Bioinformatics*, **24**(21), 2534–2536.
[doi:10.1093/bioinformatics/btn323](https://doi.org/10.1093/bioinformatics/btn323)
.\
\
Martens L, Chambers M, Sturm M, Kessner D, Levander F, Shofstahl J, Tang
WH, Rompp A, Neumann S, Pizarro AD, Montecchi-Palazzi L, Tasman N,
Coleman M, Reisinger F, Souda P, Hermjakob H, Binz P, Deutsch EW (2010).
“mzML - a Community Standard for Mass Spectrometry Data.” *Mol Cell
Proteomics*.
[doi:10.1074/mcp.R110.000133](https://doi.org/10.1074/mcp.R110.000133)
.\
\
Pedrioli PGA, Eng JK, Hubley R, Vogelzang M, Deutsch EW, Raught B, Pratt
B, Nilsson E, Angeletti RH, Apweiler R, Cheung K, Costello CE, Hermjakob
H, Huang S, Julian RK, Kapp E, McComb ME, Oliver SG, Omenn G, Paton NW,
Simpson R, Smith R, Taylor CF, Zhu W, Aebersold R (2004). “A common open
representation of mass spectrometry data and its application to
proteomics research.” *Nat Biotechnol*, **22**(11), 1459–1466.
[doi:10.1038/nbt1031](https://doi.org/10.1038/nbt1031) .
