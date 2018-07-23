# patRoon 0.1.0.9000

## July 2018
* getMCS() method for compounds
* plotStructure() method for compounds will draw MCS when mutiple indices are specified


## June 2018
* Added removeRefAnalyses argument to filter() (featureGroups method) to easily remove e.g. analyses that are used as blanks after blank subtraction.
* Added filterBy() method which removes any feature groups from a featureGroups object of which no results are present in a specified object. Methods are defined for MSPeakLists, formulaConsenus, compounds and components. This method replaces some of the functionality of the filter() method for featureGroups (formConsensus and compounds arguments).
* Added mz and chromatographic peak width range options to filter() method for feature groups.
* Moved intensity clustering code (makeHCluster) to new component type (see componentsIntClust class documentation).


## May 2018
* Added compound clustering (see makeHCluster method for compounds). This is an useful tool to get an overview of all the candidate chemical structures after compound identification. The clustering will reduce data complexity. Furthermore, maximum common sucstructures (MCS) can be calculated and plotted for each cluster to get a quick impression of the different structures of candidates.
* Added function arguments checks using [checkmate]. This guards all exported functions and methods from wrong user input arguments. If any problems are found (e.g. a wrong data type or range was specified) then the user is informed about this and what kind of input is to be expected.
* Added workaround (removed latex dependencies added automatically by `kableExtra` package) that may cause memory leakage when `reportMD()` is called repeatedly.


## April 2018
* Added unit tests (using [testthat]) and fixed several small bugs that were revealed in the process.
* Continuous integration (CI) with running tests on [CircleCI] (Linux builds) and [AppVeyor] (Windows builds) and testing coverage on [Codecov]. Docker images with patRoon and all its dependencies are automatically pushed on [Docker Hub][DH].
* Many small bug fixes.




[checkmate]: https://github.com/mllg/checkmate
[testthat]: https://github.com/r-lib/testthat
[CircleCI]: https://circleci.com/gh/rickhelmus/patRoon
[AppVeyor]: https://ci.appveyor.com/project/rickhelmus/patroon/branch/master
[Codecov]: https://codecov.io/gh/rickhelmus/patRoon
[DH]: https://hub.docker.com/r/patroonorg/patroon/
