# Miscellaneous DeeDeeExperiment methods

Miscellaneous methods for the [`DeeDeeExperiment`](DeeDeeExperiment.md)
class and its descendants that do not fit into any other documentation
category such as, for example, show and summary methods.

## Usage

``` r
# S4 method for class 'DeeDeeExperiment'
show(object)

# S4 method for class 'DeeDeeExperiment'
summary(object, FDR = 0.05, show_scenario_info = FALSE, ...)
```

## Arguments

- object:

  a [`DeeDeeExperiment`](DeeDeeExperiment.md) object

- FDR:

  Numeric, sets the significance threshold for subsetting differentially
  expressed genes based on adjusted p-values. Defaults to 0.05

- show_scenario_info:

  Logical; if TRUE, displays the associated scenario info for each DE
  contrast, if available. Defaults to FALSE

- ...:

  additional argument passed to the summary method. Currently supports:

  - `FDR`: Numeric, sets the significance threshold for subsetting
    differentially expressed genes based on adjusted p-values. Defaults
    to 0.05

  - `show_scenario_info`: Logical; if `TRUE`, displays the associated
    scenario info for each DE contrast, if available.

## Value

Returns NULL
