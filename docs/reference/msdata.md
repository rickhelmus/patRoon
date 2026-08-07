# Interface for HRMS and IMS-HRMS raw data

Description, configuration and utilities for the raw (IMS-)HRMS data
interface of patRoon.

## Usage

``` r
availableBackends(
  anaInfo = NULL,
  needTypes = NULL,
  checkOption = TRUE,
  verbose = TRUE
)
```

## Arguments

- anaInfo:

  Optional. If not `NULL` then `anaInfo` should be a [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md)
  table, and only those backends that can read each of the analyses are
  returned.

- needTypes:

  Only applicable if `anaInfo` is set: should be `"centroid"`,
  `"profile"` and/or `"ims"` to filter file types.

- checkOption:

  Set to `TRUE` to only consider backends that are included in the
  `patRoon.MS.backends` option.

- verbose:

  Set to `TRUE` to print the status of each backend.

## Value

`availableBackends` returns (invisibly) a `character` vector with the
names of the available backends.

## Details

Version 3.0 of patRoon introduced an extensible and highly optimized
interface to read raw data from HRMS and IMS-HRMS instruments. This
interface supports chooseable 'backends' which perform the reading of
file data from various formats. Subsequent steps such as the formation
of extracted ion chromatograms, mobilograms and collection and averaging
of mass spectra are then performed in `patRoon`. The interface is
largely coded in C++ (using
[Rcpp](https://CRAN.R-project.org/package=Rcpp)), uses
[OpenMP](https://www.openmp.org/) parallelization and applies several
other optimization strategies to make it suitable to rapidly process
large amounts of raw data, *e.g.* as encountered in IMS-HRMS workflow.

The following backends for reading (IMS-)HRMS data are currently
available:

- `"opentims"`: uses the
  [OpenTIMS](https://github.com/michalsta/opentims) library to read
  Bruker TIMS data. This backends supports very fast reading of raw
  instrument `.d` data files directly, and therefore does not require
  any file conversions. This backend only supports 64 bit Windows and
  Linux systems. See the `Backend installation` section below
  installation instructions.

- `"mzr"`: uses the [mzR](https://CRAN.R-project.org/package=mzR)
  package to read `.mzML` and `.mzXML` files. This package was more or
  less the default in patRoon prior to 3.0, and due to its popularity
  and age is a stable and well tested option. The `mzr` backend
  currently reads the complete analysis file at once, which makes it
  more RAM intensive compared to other backends. The read data is cached
  to speed up any subsequent operations that require the file data. This
  backend currently does not support IMS data. Since mzR is a dependency
  of patRoon, no additional installation is necessary.

- `"mstoolkit"`: Uses the
  [MSToolKit](https://github.com/mhoopmann/mstoolkit) C++ library to
  read `.mzML` and `.mzXML` files, including IMS-HRMS data. The
  `MSToolKit` library has been developed for many years, and was
  recently updated with IMS-HRMS support. See the `Backend installation`
  section below installation instructions.

- `"streamcraft"`: Uses the
  [StreamCraft](https://github.com/odea-project/streamcraft) C++ library
  to read `.mzML` and `.mzXML` files, including IMS-HRMS data. The
  `StreamCraft` library is still young and somewhat experimental. The
  library is integrated within patRoon and therefore does not require
  any further installation.

The `availableBackends` function is used to query the available backends
on the system.

## Configuration

The following package options influence the behavior of raw data
interface:

- `patRoon.MS.backends`: A `character` vector with the names of the
  backends that may be choosen. The default is all backends. The first
  backend will be chosen that is available, is able to read at least one
  of the available analysis file types and formats (as configured by the
  [analysis
  information](https://rickhelmus.github.io/patRoon/reference/analysis-information.md))
  and supports IMS if needed.

- `patRoon.MS.preferIMS`: A `logical` value that indicates whether the
  IMS data should be preferred, even if the processing step does not
  require IMS data and non-IMS data is also available. Setting this to
  `TRUE` probably result in some additional computational overhead, but
  may avoid any inconsistencies between the IMS data and non-IMS data
  that may have been introduced during the conversion step of the
  latter. This option is only relevant for the `mstoolkit` and
  `streamcraft` backends (and if one of these backends is actually
  used).

- `patRoon.threads`: An `integer` value that indicates the number of
  threads to use for parallelization (multithreading). The default is
  determined from the number of physical cores of the system (obtained
  with the
  [parallel::detectCores](https://rdrr.io/r/parallel/detectCores.html)
  function).

- `patRoon.path.TDFSDK`: The file path to the Bruker `TDF-SDK` library
  file. See the `Backend installation` section below.

## Backend installation

The `opentims` backend requires the `win64/timsdata.dll` (Windows) or
`linux64/libtimsdata.so` (Linux) file from the `TDF-SDK` from
[Bruker](https://www.bruker.com/protected/en/services/software-downloads/mass-spectrometry/raw-data-access-libraries.html)
(requires login). The patRoonExt package makes these files automatically
available for patRoon. Otherwise the `patRoon.path.TDFSDK` option should
be manually set to the file path of the `timsdata.dll` or
`linux64/libtimsdata.so` file.

When patRoon is installed from source, *e.g.* on Linux/macOS systems or
when using
[`remotes::install_github`](https://remotes.r-lib.org/reference/install_github.html)
for installation, then the <https://github.com/rickhelmus/Rmstoolkitlib>
R package must be installed in advance.

The `availableBackends` function can be used to verify if the
dependencies for these backends are met.

## References

Eddelbuettel D (2013). *Seamless R and C++ Integration with Rcpp*.
Springer, New York.
[doi:10.1007/978-1-4614-6868-4](https://doi.org/10.1007/978-1-4614-6868-4)
. ISBN 978-1-4614-6867-7.\
\
Eddelbuettel D, Balamuta J (2018). “Extending R with C++: A Brief
Introduction to Rcpp.” *The American Statistician*, **72**(1), 28-36.
[doi:10.1080/00031305.2017.1375990](https://doi.org/10.1080/00031305.2017.1375990)
.\
\
Eddelbuettel D, François R (2011). “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software*, **40**(8), 1–18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08) .\
\
Eddelbuettel D, Francois R, Allaire J, Ushey K, Kou Q, Russell N, Ucar
I, Bates D, Chambers J (2026). *Rcpp: Seamless R and C++ Integration*. R
package version 1.1.2, <https://www.rcpp.org>.\
\
Dagum L, Menon R (1998). “OpenMP: an industry standard API for
shared-memory programming.” *IEEE Computational Science and
Engineering*, **5**(1), 46-55.
[doi:10.1109/99.660313](https://doi.org/10.1109/99.660313) .\
\
Łącki MK, Startek MP, Brehmer S, Distler U, Tenzer S (2021). “OpenTIMS,
TimsPy, and TimsR: Open and Easy Access to timsTOF Raw Data.” *Journal
of Proteome Research*, **20**(4), 2122–2129. ISSN 1535-3907.
[doi:10.1021/acs.jproteome.0c00962](https://doi.org/10.1021/acs.jproteome.0c00962)
. <http://dx.doi.org/10.1021/acs.jproteome.0c00962>.\
\
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
