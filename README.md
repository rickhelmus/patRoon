# patRoon <img src="https://raw.githubusercontent.com/rickhelmus/patRoon/master/logo.png" style="width: 25%;" align="right"></img>

[![CircleCI](https://circleci.com/gh/rickhelmus/patRoon.svg?style=svg)](https://circleci.com/gh/rickhelmus/patRoon)
[![btd](https://github.com/rickhelmus/patRoon/actions/workflows/btd.yml/badge.svg)](https://github.com/rickhelmus/patRoon/actions/workflows/btd.yml)
[![codecov](https://codecov.io/gh/rickhelmus/patRoon/branch/master/graph/badge.svg)](https://codecov.io/gh/rickhelmus/patRoon)
[![patRoon r-universe](https://rickhelmus.r-universe.dev/badges/patRoon)](https://rickhelmus.r-universe.dev/patRoon)
[![DOI:10.1186/s13321-020-00477-w](https://zenodo.org/badge/DOI/10.1186/s13321-020-00477-w.svg)](https://doi.org/10.1186/s13321-020-00477-w)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04029/status.svg)](https://doi.org/10.21105/joss.04029)
[![REUSE status](https://api.reuse.software/badge/github.com/rickhelmus/patRoon)](https://api.reuse.software/info/github.com/rickhelmus/patRoon)

`patRoon` aims to provide comprehensive mass spectrometry based non-target analysis (NTA) workflows for environmental
analysis. The name is derived from a Dutch word that means _pattern_ and may also be an acronym for _hyPhenated mAss
specTROmetry nOn-target aNalysis_.

## Project news

**October 2024** `patRoon 2.3.3` is now released. This is a small maintenance release which correct the SIRIUS login procedure and comes with several other bug fixes. Please see the [Project NEWS][NEWS] for more details.

**July 2024** `patRoon 2.3.2` is now released. This is a small maintenance release with several important bug fixes and new functions to inspect raw MS data (thanks to Ricardo Cunha). Please see the [Project NEWS][NEWS] for more details.

**December 2023** `patRoon 2.3.1` is now released. This is a small maintenance release with several bug fixes. The most important fix is a correction for the prediction of concentrations from SIRIUS fingerprints. Please see the [Project NEWS][NEWS] for more details and important notes on updating `patRoon`.

**November 2023** `patRoon 2.3` is now released. The most significant changes include improved installation strategies and integration of [MS2Tox] and [MS2Quant] to predict feature toxicities and concentrations. Please see the [Project NEWS][NEWS] for more details and important notes on updating `patRoon`.

**May 2023** `patRoon 2.2` is now released. The most significant change is the addition of a new reporting interface, which brings a much improved HTML interface, many optimizations and other important new functionality. Furthermore, `patRoon 2.2` introduces improved `SIRIUS 5` support, a new TP screening algorithm using formula libraries and many other improvements, of which many thanks to the great user feedback. Please see the [Project NEWS][NEWS] for details.

**March 2023** The Docker images moved to a new host. Please see the see [installation details in the handbook][handbook-inst] to obtain the latest images.

**May 2022** `patRoon 2.1` is now available. This new release integrates prediction of transformation products with
[CTS], adds several feature intensity normalization methods, adds new functionality and improvements for reporting TP
data and supports loading, processing and annotation with MS libraries such as MassBank. Please see the [Project
NEWS][NEWS] for details.

## Introduction

Mass spectrometry based non-target analysis is used to screen large numbers of chemicals simultaneously. For this
purpose, high resolution mass spectrometry instruments are used which are typically coupled (or _hyphenated_) with
chromatography (_e.g._ LC or GC). The size and complexity of resulting data makes manual processing impractical. Many
software tools were/are developed to facilitate a more automated approach. However, these tools are generally not
optimized for environmental workflows and/or only implement parts of the functionality required.

`patRoon` combines established software tools with novel functionality in order to provide comprehensive NTA workflows.
The different algorithms are provided through a consistent interface, which removes the need to know all the details of
each individual software tool and performing tedious data conversions during the workflow. The table below outlines the
major functionality of `patRoon`.

Functionality | Description | Algorithms
---------------------- | ------------------------------------------------------------------------ | -----------------------
Raw data pre-treatment | MS format conversion (e.g. vendor to `mzML`) and calibration.            | [ProteoWizard], [OpenMS], [DataAnalysis]
Feature extraction     | Finding features and grouping them across analyses.                      | [XCMS], [OpenMS], [enviPick], [DataAnalysis], [KPIC2], [SIRIUS], [SAFD]
Suspect screening      | Finding features with suspected presence by MS and chromatographic data. Estimation of identification confidence levels. | Native
MS data extraction     | Automatic extraction and averaging of feature MS(/MS) peak lists.        | Native, [mzR], [DataAnalysis]
Formula annotation     | Automatic calculation of formula candidates for features.                | [GenForm], [SIRIUS], [DataAnalysis]
Compound annotation    | Automatic (_in silico_) compound annotation of features.                 | [MetFrag], [SIRIUS], Native
Componentization & adduct annotation | Grouping of related features based on chemistry (e.g. isotopes, adducts and homologs), hierarchical clustering or MS/MS similarity into components. Using adduct and isotope annotations for prioritizing features and improving formula/compound annotations. | [RAMClustR], [CAMERA], [nontarget R package][nontarget], [OpenMS], [cliqueMS], Native
Combining algorithms   | Combine data from different algorithms (e.g. features, annotations) and generate a consensus. | Native
_Sets workflows_       | Simultaneous processing and combining +/- MS ionization data             | Native
Transformation product (TP) screening | Automatic screening of TPs using library/_in-silico_ data, MS similarities and classifications. Tools to improve compound TP annotation. | [BioTransformer], [PubChemLite][PubChemLiteTR], Native
Reporting              | Automatic reporting of all important workflow data. An example HTML report can be viewed [here][example]. | Native
Data clean-up & prioritization | Filters for blanks, replicates, intensity thresholds, neutral losses, annotation scores, identification levels and many more. | Native
Data curation          | Several graphical interactive tools and functions to inspect and remove unwanted data. | Native

The workflow of non-target analysis typically depends on the aims and requirements of the study and the instrumentation
and methodology used for sample analysis. For this reason, `patRoon` does not enforce a certain workflow. Instead, most
workflow steps are optional, fully configurable and algorithms can easily be mixed or even combined.

## Implementation details

* `patRoon` is implemented as an [R] package, which allows easy interfacing with the many other `R` based MS tools and other data processing functionality from `R`.
* Fully open-source (GPLv3).
* Developed on Windows, Linux and macOS
* S4 classes and generics are used to implement a consistent interface to all supported algorithms.
* Continuous integration is used to automatically perform unit tests, update the [Website][ghweb] and documentation, and maintaining installation resources such as the [patRoonDeps repository][patRoonDeps], [Docker image][DockerImg] and `patRoon` bundle (see [the handbook][handbook-inst] for more details).
* Supports all major instrument vendor input formats (through usage of [ProteoWizard] and [DataAnalysis]).
* Optimizations
    * `data.table` is used internally as a generally much more efficient alternative to `data.frame`.
    * The [processx] and [future] `R` packages are used for parallelization.
    * Results from workflow steps are cached within a [SQLite] database to avoid repeated computations.
    * Code for loading MS and EIC data, MS similarity calculations and others were implemented in `C++` to reduce computational times.
* The [RDCOMClient] `R` package is used to interface with Bruker DataAnalysis algorithms.
* The [Shiny] `R` package was used to implement several GUI tools.
* The reporting functionality relies on the excellent [R markdown][Rmd] and related packages such as [flexdashboard], [bslib] and [reactable].


## Installation

`patRoon` and its dependencies can be installed in various ways. Please see the [installation section in the handbook][handbook-inst] for more
information.


## Getting started

For a very quick start:

``` r
library(patRoon)
newProject()
```

The `newProject()` function will pop-up a dialog screen (requires [R Studio][RStudio]), which will allow you to quickly
select the analyses and common workflow options to subsequently generate a template `R` processing script.

However, for a better guide to get started it is recommended to read the [tutorial]. Afterwards the [handbook] is a
recommended read if you want to know more about advanced usage of `patRoon`. Finally, the [reference] outlines all the
details of the `patRoon` package.


## Citing

When you use `patRoon` please cite its publications:

Rick Helmus, Thomas L. ter Laak, Annemarie P. van Wezel, Pim de Voogt and Emma L. Schymanski. [patRoon: open source
software platform for environmental mass spectrometry based non-target
screening](https://doi.org/10.1186/s13321-020-00477-w). _Journal of Cheminformatics_ **13**, 1 (2021)

Rick Helmus, Bas van de Velde, Andrea M. Brunner, Thomas L. ter Laak, Annemarie P. van Wezel and Emma L. Schymanski.
[patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening]( https://doi.org/10.21105/joss.04029). _Journal of Open Source Software_, 7(71), 4029

`patRoon` builds on many open-source software tools and open data sources. Therefore, it is important to also cite their
work when using these algorithms via `patRoon`.

## Contributing

For bug reports, code contributions (pull requests), questions, suggestions and general feedback please use the [GitHub page](https://github.com/rickhelmus/patRoon).


[R]: https://www.r-project.org/
[NEWS]: https://github.com/rickhelmus/patRoon/blob/master/NEWS.md
[XCMS]: https://github.com/sneumann/xcms
[OpenMS]: http://openms.de/
[enviPick]: https://cran.r-project.org/web/packages/enviPick/index.html
[KPIC2]: https://github.com/hcji/KPIC2
[SAFD]: https://bitbucket.org/SSamanipour/safd.jl/src/master/
[DataAnalysis]: https://www.bruker.com/
[ProfileAnalysis]: https://www.bruker.com/
[mzR]: https://github.com/sneumann/mzR/
[GenForm]: https://sourceforge.net/projects/genform
[SIRIUS]: https://bio.informatik.uni-jena.de/software/sirius/
[MetFrag]: http://c-ruttkies.github.io/MetFrag/
[RAMClustR]: https://github.com/sneumann/RAMClustR
[CAMERA]: http://msbi.ipb-halle.de/msbi/CAMERA/
[nontarget]: https://cran.r-project.org/web/packages/nontarget/index.html
[cliqueMS]: https://github.com/osenan/cliqueMS
[BioTransformer]: https://bitbucket.org/djoumbou/biotransformer/src/master/
[PubChemLiteTR]: https://doi.org/10.5281/zenodo.5644560
[future]: https://github.com/HenrikBengtsson/future
[tutorial]: https://rickhelmus.github.io/patRoon/articles/tutorial.html
[handbook]: https://rickhelmus.github.io/patRoon/handbook_bd/index.html
[handbook-inst]: https://rickhelmus.github.io/patRoon/handbook_bd/installation.html
[reference]: https://rickhelmus.github.io/patRoon/reference/index.html
[RStudio]: https://www.rstudio.com/
[processx]: https://github.com/r-lib/processx
[SQLite]: https://www.sqlite.org/index.html
[RDCOMClient]: http://www.omegahat.net/RDCOMClient/
[Shiny]: https://shiny.rstudio.com/
[example]: https://rickhelmus.github.io/patRoon/examples/report.html
[ProteoWizard]: http://proteowizard.sourceforge.net/
[ghweb]: https://rickhelmus.github.io/patRoon/
[patRoonDeps]: https://github.com/rickhelmus/patRoonDeps
[miniCRAN]: https://cran.r-project.org/web/packages/miniCRAN/index.html
[DockerImg]: https://uva-hva.gitlab.host/R.Helmus/patroon/container_registry/2
[CTS]: https://qed.epa.gov/cts/
[Rmd]: https://rmarkdown.rstudio.com/
[flexdashboard]: https://pkgs.rstudio.com/flexdashboard/
[bslib]: https://rstudio.github.io/bslib/
[reactable]: https://glin.github.io/reactable/
[MS2Tox]: https://github.com/kruvelab/MS2Tox
[MS2Quant]: https://github.com/kruvelab/MS2Quant
