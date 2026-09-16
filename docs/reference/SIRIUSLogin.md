# Log in to SIRIUS

Performs an account login to `SIRIUS`. This function starts the `SIRIUS`
API (if needed) and performs the login. It can be used to explicitly log
in before running
[`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)
or
[`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md).

## Usage

``` r
SIRIUSLogin(login = "check", alwaysLogin = FALSE, SIRIUSAPI = NULL)
```

## Arguments

- login, alwaysLogin:

  Specifies if and how account logging of SIRIUS should be handled:

  `login=FALSE`: no automatic login is performed and the active login
  status is not checked.

  `login="check"`: aborts if no active login is present.

  `login="interactive"`: interactively ask for login (using
  [getPass](https://CRAN.R-project.org/package=getPass)).

  `login=c(username="...", password="...")`: perform the login with the
  given details. For security reasons, please do not enter the details
  directly, but use e.g. environment variables or store/retrieve them
  with the [keyring](https://CRAN.R-project.org/package=keyring)
  package.

  if `alwaysLogin=TRUE` then a login is always performed, otherwise only
  if SIRIUS reports no active login.

  See the [SIRIUS
  website](https://v6.docs.sirius-ms.io/account-and-license/) and
  patRoon handbook for more information.

  **NOTE**: By loggin in you will accept the terms of the Service and
  Privacy Policy of the SIRIUS Webservice.

- SIRIUSAPI:

  An `rsirius_api` object for connecting to the `SIRIUS` API. If `NULL`,
  a new connection will be started automatically.
