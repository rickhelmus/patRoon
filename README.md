# patRoon
`patRoon` aims to provide a common interface to various (primarily
open-source) software solutions for mass spectrometry based non-target analysis.
The name is derived from a Dutch word that means _pattern_ and may also be an acronym for _hyPhenated
mAss-specTROmetry nOn-target aNalysis_.

`patRoon` is developed as an `R` package and leverages other `R` based packages and external software to provide the following workflow steps that are typically used within non-target (NT) analysis:

* **Data preparation**: conversion between and export to common open MS data formats (e.g. mzXML and mzML) ([OpenMS], [DataAnalysis]).
* **Extraction and grouping of features** ([XCMS], [OpenMS], [enviPick], [enviMass], [ProfileAnalysis]).
* **Data cleanup**: post filtering of data to improve its quality and aid prioritization.
* **Automatic extraction of MS and MS/MS** ([mzR], [DataAnalysis]).
* **Formula calculation**: automatic calculation of candidate formulae for detected features ([GenForm], [SIRIUS], [DataAnalysis]).
* **Compound identification**: automatic annotation of MS/MS spectra and retrieval of candidate structures ([MetFrag], [SIRIUS with CSI:FingerID][SIRIUS]).
* **Generation of components**: grouping of (chemically) related features such as isotopes, adducts and homologs ([RAMClustR], [CAMERA], [nontarget R package][nontarget]).
* **reporting** of all workflow steps in various formats: _CSV_, _PDF_ and _HTML_.

## Installation
You can install patRoon from github with:

``` r
install.packages("devtools")
devtools::install_github("rickhelmus/patRoon")
```

### R Dependencies
Besides 'regular' R package dependencies, which should be installed automatically with `install_github()`, the following optional R packages may also need to be installed:

``` r
install.packages("RDCOMClient") # required for DataAnalysis functionality
devtools::install_github("cbroeckl/RAMClustR", build_vignettes = TRUE, dependencies = TRUE)
devtools::install_github("c-ruttkies/MetFragR/metfRag") # only when using R interface (not by default)
```

The following [Bioconductor] packages dependencies may need to be installed manually:
```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("mzR", "xcms", "CAMERA"))
```

Note that in order to fulfill installation of the `rJava` package dependency you may need to setup a proper JDK (see for instance [the rJava website][rJava]).

### Other dependencies

Depending on which functionality is used, the following (optional) external dependencies should be installed:

* [OpenMS]
* [SIRIUS]
* [MetFrag CL][MetFrag-CL] (required if the command-line version of MetFrag is used, the default)
* [pngquant] (required for `reportMD` with `optimizePng` argument set to `TRUE`)

The location of the OpenMS, SIRIUS and pngquant should either be set within the _PATH_ environment variable (the OpenMS installer will do this automatically). Alternatively, one or more of the `R` options should be set below. Note that the location of the MetFrag CL _jar_ should always be set.

```r
options("patRoon.path.SIRIUS" = "~/C:/sirius-win64-3.5.1") # location where SIRIUS was extracted
options("patRoon.path.OpenMS" = "/usr/local/bin") # directory with the OpenMS binaries
options("patRoon.path.pngquant" = "~/pngquant") # directory containing pngquant binary
options("patRoon.path.metFragCL" = "~/MetFrag2.4.2-CL.jar") # full location to the jar file
```

## Getting started
``` r
library(patRoon)
newProject()
```

read the [tutorial](docs/articles/tutorial.html).

read the [reference](docs/reference/index.html).

## WIP
More documentation will folow..


[XCMS]: https://github.com/sneumann/xcms
[OpenMS]: http://openms.de/
[enviPick]: https://cran.r-project.org/web/packages/enviPick/index.html
[DataAnalysis]: https://www.bruker.com/
[enviMass]: http://www.looscomputing.ch/eng/enviMass/overview.htm
[ProfileAnalysis]: https://www.bruker.com/
[mzR]: https://github.com/sneumann/mzR/
[GenForm]: https://sourceforge.net/projects/genform
[SIRIUS]: https://bio.informatik.uni-jena.de/software/sirius/
[MetFrag]: http://c-ruttkies.github.io/MetFrag/
[RAMClustR]: https://github.com/sneumann/RAMClustR
[CAMERA]: http://msbi.ipb-halle.de/msbi/CAMERA/
[nontarget]: https://cran.r-project.org/web/packages/nontarget/index.html
[MetFrag-CL]: http://c-ruttkies.github.io/MetFrag/projects/metfragcl/
[pngquant]: https://pngquant.org/
[Bioconductor]: https://www.bioconductor.org
[rJava]: http://www.rforge.net/rJava/
