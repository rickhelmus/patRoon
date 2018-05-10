# patRoon 0.1.0.9000

## May 2018
* Added function arguments checks using [checkmate]. This change guards all exported functions and methods from wrong user input arguments. If any problems are found (e.g. a wrong data type or range was specified) then the user is informed about this and what kind of input is to be expected.
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
