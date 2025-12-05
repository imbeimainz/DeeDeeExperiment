# Convert `muscat::pbDS()` results into a flat list of data frames compatible with `DeeDeeExperiment`

This helper function extracts and flattens the nested structure returned
by [`muscat::pbDS()`](https://rdrr.io/pkg/muscat/man/pbDS.html),
returning one table per contrast–cluster combination. Each resulting
data frame will have standardized column names (`log2FoldChange`,
`pvalue`, `padj`).

## Usage

``` r
muscat_list_for_dde(res, padj_col = c("p_adj.loc", "p_adj.glb"))
```

## Arguments

- res:

  A list, typically the output of
  [`muscat::pbDS()`](https://rdrr.io/pkg/muscat/man/pbDS.html) function,
  containing one or more contrasts

- padj_col:

  A character string specifying which adjusted p-value column to
  extract. It can be either "p_adj.loc" or "p_adj.glb".

## Value

A named list of data frames

## Details

The function is intended to simplify the integration of muscat results
into a `DeeDeeExperiment` object. It automatically renames columns
(`logFC` to `log2FoldChange`, `p_val` to `pvalue`, and the selected
adjusted p-value column to `padj`) and annotates the resulting list with
metadata about the originating package.

The function checks that each contrast entry contains a valid `table`
component as expected from `pbDS()` output. Invalid or empty contrasts
are skipped with a warning message. The names of the list elements in
`res` must match contrast names found in the `table` slot of each entry.

## Examples

``` r
data("muscat_pbDS_res", package = "DeeDeeExperiment")
new_muscat_list <- muscat_list_for_dde(list(`stim-ctrl` = muscat_res))
#> ℹ Returning 3 muscat contrast-cluster tables formatted for DeeDeeExperiment
```
