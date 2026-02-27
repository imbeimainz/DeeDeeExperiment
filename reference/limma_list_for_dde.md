# Convert a `MArrayLM` object with multiple contrasts into a list of DE results tables compatible with DeeDeeExperiment

This helper function extracts DE results for each contrast contained in
a `limma::MArrayLM` object and reformats them into a list of
standardized data frames suitable for integration in a
`DeeDeeExperiment` object. Each resulting data frame includes renamed
columns: `logFC` to `log2FoldChange` ,`P.Value` to `pvalue`, and
`adj.P.Val` to `padj`

## Usage

``` r
limma_list_for_dde(fit, number = nrow(fit), sort.by = "none")
```

## Arguments

- fit:

  A `MArrayLM` object, as produced by the `limma` workflow

- number:

  An integer specifying the maximum number of genes to extract per
  contrast

- sort.by:

  A character string specifying which statistic to rank the genes by. It
  must be one of the values accepted by the `sort.by` argument in the
  [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html)
  function. Defaults to "none".

## Value

A named list of DE results tables, one per contrast. Each table contains
standardized columns and the list is annotated with metadata indicating
its package origin (`limma`).

## Details

The function assumes that each column in `fit$coefficients` corresponds
to a contrast of interest. The names of the resulting list elements are
taken directly from `colnames(fit$coefficients)`. The names in the input
object must therefore accurately reflect the intended contrast names.

## See also

[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html)

## Examples

``` r
data("de_limma", package = "DeeDeeExperiment")
new_limma_list <- limma_list_for_dde(de_limma)
#> Error in topTable(fit, coef = x, number = number, sort.by = sort.by): could not find function "topTable"
```
