# Obtain a SIRIUS job configuration

This function obtains a
[JobSubmission](https://rdrr.io/pkg/RSirius/man/JobSubmission.html)
configuration object for use with
[`generateFormulasSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateFormulasSIRIUS.md)
and
[`generateCompoundsSIRIUS`](https://rickhelmus.github.io/patRoon/reference/generateCompoundsSIRIUS.md).

## Usage

``` r
getSIRIUSConfig(
  config = NULL,
  import = NULL,
  login = "check",
  alwaysLogin = FALSE,
  SIRIUSAPI = NULL
)
```

## Arguments

- config:

  A configuration specification. Use `NULL` to obtain the server
  default, a `character` string to fetch a named configuration, `NA` to
  list available configuration names, or a `list` /
  [`RSirius::JobSubmission`](https://rdrr.io/pkg/RSirius/man/JobSubmission.html)
  object to use directly.

- import:

  Path to a JSON file to import the configuration from, *e.g.* exported
  from the `SIRIUS` GUI. Cannot be used together with `config`.

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

## Value

A
[`RSirius::JobSubmission`](https://rdrr.io/pkg/RSirius/man/JobSubmission.html)
object with the job configuration, or a `character` vector of
configuration names when `config = NA`.

## Details

The function can return the server default configuration, fetch a named
configuration, import a configuration from a JSON file, or use a
provided configuration object/list. If needed, the SIRIUS API is started
automatically and a login can be performed. There are *many*
configuration options available in `SIRIUS`. It is probably easiest to
explore and configure them through the `SIRIUS` GUI, save or export the
configuration profile and then load it via the `config` or `import`
arguments.
