# Default limits and tolerances

Get and configure default limits and tolerances used by patRoon.

## Usage

``` r
defaultLim(category, level)

getLimIMS()

genLimitsFile(out = "limits.yml", IMS = "bruker")
```

## Arguments

- category, level:

  The category and level of the limit to be returned. See the detail
  sections below. For mobility related limits, the `category="mobility"`
  (instead of instrument specific categories) should be used.

- out:

  The full output path for the new file.

- IMS:

  The IMS instrument type. This sets the `IMS` variable in the
  **general** section. Valid values are `"bruker"` and `"agilent"`.

## Details

Tolerances for retention times, *m/z* values and other numerical limits
and tolerances are widely used in patRoon to process data. Their
defaults used to be hardcoded directly as defaults to function
arguments. Since version 3.0 the defaults are centralized in a
`limits YAML` file. This simplifies their configuration and makes it
easier to switch the defaults between different HRMS instruments.

The limits configuration file is taken from one of the following
locations (in order):

- the path specified in the `patRoon.path.limits` option (if specified)

- the `limits.yml` file in the current working directory (if present)

- the default `limits.yml` file embedded in patRoon

`defaultLim` returns the limits for a specific category and tolerance
level.

`getLimIMS` returns the type of IMS instrument specified in the limits
configuration file. This is used to determine which mobility limits to
use, and is also used in some other functions to determine the default
behavior for IMS related processing.

`genLimitsFile` generates a new `limits.yml` configuration file in the
specified path. The file is created with the defaults embedded in
patRoon (see details below). Generating a custom limits is primarily
useful for non-Bruker IMS workflows.

## Note

Most of the defaults were derived with Bruker TOF HRMS instrumentation
in mind, but should be reasonable with minor adjustments for others.

## limits YAML file format

The limits are configured in a simple `YAML` file. A brief summary of
the format is given here.

The **general** section has the `IMS` field which specifies the type of
instrument used. Currently this is either `bruker` or `agilent`.

The next sections describe absolute and relative (suffixed by `_rel`)
limits and tolerances for retention, *m/z*, mobility and CCS data. These
are divided into several tolerance levels: `very_narrow`, `narrow`,
`medium` and `wide`. In general, the `very narrow` tolerances are used
to compare data that should be equivalent, but may be slightly different
due to *e.g.* small rounding errors. The `narrow` tolerances are
generally sufficient when only small deviations are expected, *e.g.*
comparing different *m/z* values from well calibrated HRMS instruments.
The `medium` tolerances are generally used when a slightly larger
tolerance is needed, *e.g.* when *m/z* values from raw (non-averaged)
spectra. The `wide` values are mainly used for plot limits, *e.g.* to
provide a reasonable zoom-out. The sections only define values for the
tolerance levels actually used in patRoon.

The `mobility_bruker` and `mobility_agilent` sections specify the
default mobility limits for Bruker and Agilent systems, respectively.
Which are used is set by the `IMS` variable in the **general** section.

To see which limits are used in which functions, please refer to the
`Usage` section of the respective functions, specifically how the
`defaultLim` function is used to assign function argument defaults.

**NOTE**: the choice between using the `narrow` and `medium` tolerance
is not always clear, and there is some inconsistency in its use
throughout patRoon (primarily due to legacy code or keeping defaults
from external algorithms).

## Default limits YAML file

    general:
        version: 1
        IMS: bruker
    retention:
        very_narrow: 2
        narrow: 6
        medium: 12
        wide: 30
    mz:
        very_narrow: 0.001
        narrow: 0.002
        medium: 0.005
        wide: 0.02
        narrow_rel: 5
        medium_rel: 10
    mobility_bruker:
        very_narrow: 0.01
        narrow: 0.02
        medium: 0.04
        wide: 0.2
        medium_rel: 0.05
    mobility_agilent:
        very_narrow: 0.1
        narrow: 0.2
        medium: 0.4
        wide: 2
        medium_rel: 0.05
    CCS:
        medium: 10
        medium_rel: 0.05
