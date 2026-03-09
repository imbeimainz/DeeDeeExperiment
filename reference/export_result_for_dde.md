# Export DEA/FEA/ASSAY results from a `DeeDeeExperiment` to excel files

Extracts DEA and/or FEA results stored in a `DeeDeeExperiment` and
writes them to excel files. Each contrast/result is written to a
separate sheet. DEA and FEA are written to separate files.

## Usage

``` r
export_result_for_dde(
  x,
  res_type = c("assay", "dea", "fea", "all"),
  res_format = c("minimal", "original"),
  output_dir = getwd(),
  file_name = NULL,
  assay_type = c("counts"),
  force = FALSE
)
```

## Arguments

- x:

  A `DeeDeeExperiment` object containing DEA/FEA results to be extracted

- res_type:

  A character vector indicating the result to extract (i.e. `"dea"`,
  `"fea"`, `"assay"`) from the corresponding slot If set to `"all"`,
  DEAs, FEAs, as well as the specified assays are extracted into
  separate excel files

- res_format:

  A character string specifying the DEA/FEAs output format to be
  exported (i.e. `"minimal"` or `"original"`)

- output_dir:

  A character string specifying the directory where the excel file will
  be written

- file_name:

  Optional character string specifying a file name. If `NULL`, a default
  name is generated

- assay_type:

  Optional character vector specifying the assays to export. It defaults
  to `"counts"`

- force:

  Logical. If `TRUE`, an existing file with the same name will be
  overwritten

## Value

An (invisible) character vector of file paths to the exported excel
files. The vector may contain one path `res_type = "dea"` or
`res_type = "fea" `or three paths `res_type = "all"`

## Examples

``` r
data("de_named_list", package = "DeeDeeExperiment")
data("topGO_results_list", package = "DeeDeeExperiment")
dde <- DeeDeeExperiment(de_results = de_named_list, enrich_results = topGO_results_list)
#> ℹ creating a mock SCE from the rows of the DE result objects, if available
#> ℹ FEA "ifng_vs_naive" matched directly to DE contrast "ifng_vs_naive"
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> ℹ FEA "ifngsalmo_vs_naive" matched directly to DE contrast "ifngsalmo_vs_naive"
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> ℹ FEA "salmonella_vs_naive" matched directly to DE contrast "salmonella_vs_naive"
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> ℹ FEA "salmo_both" matched directly to DE contrast "salmo_both"
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...

export_result_for_dde(dde,
res_type = "dea",
res_format = "minimal", output_dir = tempdir(), force = TRUE)
#> Found 4 DEA results
#> Writing results to: /var/folders/5q/v_ms_h9x6mv05dzlf94g48d00000gn/T//Rtmph5BxFT/dde_minimal_DEA.xlsx

export_result_for_dde(dde,
res_type = c("dea", "fea"),
res_format = "original", output_dir = tempdir(), force = TRUE)
#> Error in match.arg(res_type): 'arg' must be of length 1

```
