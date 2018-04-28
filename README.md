# patRoon

[![CircleCI](https://circleci.com/gh/rickhelmus/patRoon.svg?style=svg)](https://circleci.com/gh/rickhelmus/patRoon)
[![Build status](https://ci.appveyor.com/api/projects/status/52nnpq8kqpkjqc92/branch/master?svg=true)](https://ci.appveyor.com/project/rickhelmus/patroon/branch/master)
[![codecov](https://codecov.io/gh/rickhelmus/patRoon/branch/master/graph/badge.svg)](https://codecov.io/gh/rickhelmus/patRoon)
[![](https://images.microbadger.com/badges/image/patroonorg/patroon.svg)](https://microbadger.com/images/patroonorg/patroon)

`patRoon` aims to provide a common interface to various (primarily
open-source) software solutions for mass spectrometry based non-target analysis.
The name is derived from a Dutch word that means _pattern_ and may also be an acronym for _hyPhenated
mAss specTROmetry nOn-target aNalysis_.

Mass spectrometry based non-target analysis is used to screen large numbers of chemicals simultaneously. For this purpose, high resolution mass spectrometry instruments are used which are typically coupled (or _hyphenated_) with chromatography (_e.g._ LC or GC) to provide additional separation of analytes prior to mass detection. The size and complexity of resulting data makes manual processing impractical, however, many dedicated software tools were developed (and are generally still in active development) to facilitate a more automated approach. The aim of `patRoon` is to harmonize the many tools available to provide a consistent user interface without the need to know all the details of each individual software tool and remove the need for tedious conversion of data when multiple tools are used. Many of the available tools were developed using the [R language][R] or otherwise easily interface with `R`, hence, it was a straightforward decision to develop `patRoon` as an `R` package.

The workflow of non-target analysis is typically highly dependent on several factors such as the analytical instrumentation used and requirements of the study. For this reason, `patRoon` does not enforce a certain workflow. Instead, most workflow steps are optional, are highly configurable and algorithms can easily be mixed or even combined.

Below is an overview of the most common non-target workflow steps supported along with various software tools that were used to implement them:

* **Data preparation**: conversion between and export to common open MS data formats (e.g. mzXML and mzML) ([OpenMS], [DataAnalysis]).
* **Extraction and grouping of features** ([XCMS], [OpenMS], [enviPick], [enviMass], [ProfileAnalysis]).
* **Data cleanup**: post filtering of data to improve its quality and aid prioritization.
* **Automatic extraction of MS and MS/MS** data ([mzR], [DataAnalysis]).
* **Formula calculation**: automatic calculation of candidate formulae for detected features ([GenForm], [SIRIUS], [DataAnalysis]).
* **Compound identification**: automatic annotation of MS/MS spectra and retrieval of candidate structures ([MetFrag], [SIRIUS with CSI:FingerID][SIRIUS]).
* **Generation of components**: grouping of (chemically) related features such as isotopes, adducts and homologs ([RAMClustR], [CAMERA], [nontarget R package][nontarget]).
* **reporting** of all workflow steps in various formats: _CSV_, _PDF_ and _HTML_.

An example HTML report which gives an overview of what data can be produced can viewed [here][example].

Some implementation notes:

* Developed on both Windows and Linux
* `data.table` is used internally as a generally much more efficient alternative to `data.frame`.
* The [processx] R package is used to execute command line processes in parallel to reduce computation times.
* Results from workflow steps are cached within a [SQLite] database to avoid repeated computations.
* The [RDCOMClient] is used to provide several functionalities from DataAnalysis (Bruker).
* Depending on which software algorithms, the following MS data formats are supported: `mzXML`, `mzML` and `.d` (Bruker).
* The [Shiny] R package was used to implement several GUI tools.
* S4 classes and generics are used to implement a consistent interface to the many different supported software algorithms for each workflow step.

Currently most functionalities have been implemented. On the short term the main focus is to improve and extend documentation, improve robustness when dealing with user input and find & remove bugs by (implementing automated) testing.

## Installation

Prior to installation you have to make several software packages dependencies are installed.

### R Dependencies

Most `R` package dependencies are installed automatically. However, it may be necessary to manually install the following [Bioconductor] dependencies:

```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("mzR", "xcms", "CAMERA"))
```

In addition, the following optional `R` packages may also need to be installed:

``` r
install.packages("RDCOMClient") # required for DataAnalysis functionality

install.packages("devtools") # needed only if not already installed
devtools::install_github("cbroeckl/RAMClustR", build_vignettes = TRUE, dependencies = TRUE)
devtools::install_github("c-ruttkies/MetFragR/metfRag") # only when using the R interface (not by default)
```

Note that in order to fulfill installation of the `rJava` package dependency you may need to setup a proper JDK (see for instance [the rJava website][rJava]).

Finally, Windows users need to make sure that [Rtools] is installed.

### Other dependencies

Depending on which functionality is used, the following (optional) external dependencies should be installed:

* [OpenMS]
* [SIRIUS]
* [MetFrag CL][MetFrag-CL] (required if the command-line version of MetFrag is used, the default)
* [pngquant] (required for `reportMD` with `optimizePng` argument set to `TRUE`)

The location of the OpenMS, SIRIUS and pngquant binaries should either be set within the _PATH_ environment variable (the OpenMS installer will do this automatically). Alternatively, one or more of the `R` options should be set below.

```r
options("patRoon.path.SIRIUS" = "~/C:/sirius-win64-3.5.1") # location where SIRIUS was extracted
options("patRoon.path.OpenMS" = "/usr/local/bin") # directory with the OpenMS binaries
options("patRoon.path.pngquant" = "~/pngquant") # directory containing pngquant binary
options("patRoon.path.metFragCL" = "~/MetFrag2.4.2-CL.jar") # full location to the jar file
```

Please note that the location of the MetFrag CL jar file should always be set manually.


### patRoon installation

Finally, you can install patRoon from github with:

``` r
install.packages("devtools") # needed only if not already installed
devtools::install_github("rickhelmus/patRoon")
devtools::install_github("rickhelmus/patRoonData") # example data used by tutorial
```

In some cases the `patRoon` installation may fail with an error telling you that `SVN` is absent. In this case you may have better luck by changing to the [remotes] package:

``` r
install.packages("remotes")
remotes::install_github("rickhelmus/patRoon")
```


## Getting started

``` r
library(patRoon)
newProject()
```

The `newProject()` function will pop-up a dialog screen (requires [R Studio][RStudio]!) which will allow you to quickly select the analyses and common workflow options to subsequently generate a template `R` processing script.

However, for a better guide to get started it is recommended to read the [tutorial] (work in progress).

Finally, the [reference] outlines all the details of the `patRoon` package.

[R]: https://www.r-project.org/
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
[tutorial]: https://rickhelmus.github.io/patRoon/articles/tutorial.html
[reference]: https://rickhelmus.github.io/patRoon/reference/index.html
[remotes]: https://github.com/r-lib/remotes#readme
[Rtools]: https://cran.r-project.org/bin/windows/Rtools/
[RStudio]: https://www.rstudio.com/
[processx]: https://github.com/r-lib/processx
[SQLite]: https://www.sqlite.org/index.html
[RDCOMClient]: http://www.omegahat.net/RDCOMClient/
[Shiny]: https://shiny.rstudio.com/
[example]: https://rickhelmus.github.io/patRoon/examples/report.html


