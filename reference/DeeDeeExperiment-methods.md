# Methods for [DeeDeeExperiment](DeeDeeExperiment.md) objects

The [`DeeDeeExperiment()`](DeeDeeExperiment.md) class provides a family
of methods to get and set DE-related information and functional
enrichment results in [DeeDeeExperiment](DeeDeeExperiment.md) objects.

## Usage

``` r
# S4 method for class 'DeeDeeExperiment'
getDEAInfo(x)

# S4 method for class 'DeeDeeExperiment'
getDEAInfo(x) <- value

# S4 method for class 'DeeDeeExperiment'
getDEANames(x)

# S4 method for class 'DeeDeeExperiment'
renameDEA(x, old_name, new_name)

# S4 method for class 'DeeDeeExperiment'
addDEA(x, dea, force = FALSE)

# S4 method for class 'DeeDeeExperiment'
removeDEA(x, dea_name, remove_linked_fea = FALSE)

# S4 method for class 'DeeDeeExperiment'
getDEA(
  x,
  dea_name = NULL,
  format = c("minimal", "original"),
  extra_rd = NULL,
  type = c("DFrame", "data.frame"),
  verbose = FALSE
)

# S4 method for class 'DeeDeeExperiment'
getDEAList(x, format = c("minimal", "original"), verbose = FALSE)

# S4 method for class 'DeeDeeExperiment'
addScenarioInfo(x, dea_name, info = NULL, force = FALSE)

# S4 method for class 'DeeDeeExperiment'
getFEAInfo(x)

# S4 method for class 'DeeDeeExperiment'
getFEAInfo(x) <- value

# S4 method for class 'DeeDeeExperiment'
getFEANames(x)

# S4 method for class 'DeeDeeExperiment'
renameFEA(x, old_name, new_name)

# S4 method for class 'DeeDeeExperiment'
addFEA(
  x,
  fea,
  de_name = NA_character_,
  fe_name = NULL,
  fea_tool = "auto",
  force = FALSE,
  verbose = FALSE
)

# S4 method for class 'DeeDeeExperiment'
removeFEA(x, fea_name)

# S4 method for class 'DeeDeeExperiment'
getFEA(x, fea_name = NULL, format = c("minimal", "original"), verbose = FALSE)

# S4 method for class 'DeeDeeExperiment'
getFEAList(
  x,
  dea_name = NULL,
  format = c("minimal", "original"),
  verbose = FALSE
)

# S4 method for class 'DeeDeeExperiment'
linkDEAandFEA(x, dea_name, fea_name, force = FALSE)
```

## Arguments

- x:

  A [`DeeDeeExperiment()`](DeeDeeExperiment.md) object

- value:

  Replacement value for replacement methods.

- old_name:

  A character vector of existing DEA names to be renamed in a
  `DeeDeeExperiment` object

- new_name:

  A character vector with new names to assign to existing DEA names in a
  `DeeDeeExperiment` object. It must be the same length of `old_name`,
  and contains unique values that don't overlap with existing DEA names.

- dea:

  A named list of DE results, in any of the formats supported by the
  package (currently: results from DESeq2, edgeR, limma).

- force:

  A logical, indicating whether to overwrite results when introducing
  the same results name. It defaults to FALSE.

- dea_name:

  Character value, specifying the name of the DE analysis to get or
  remove, or match against (e.g., to fetch associated FEA results), or
  to which additional context and information can be attached

- remove_linked_fea:

  A logical, specifying whether to remove or not the linked FEA when a
  DEA results is removed

- format:

  A character string, specifying the DEA/FEAs output format. It takes
  either "minimal" to return only essential columns (e.g. log2FC,
  p-value, adjusted p-value for DEAs, or gs_id, gs_description,
  gs_pvalue, gs_genes... for FEAs), or "original" to return the full
  result object. It defaults to "minimal"

- extra_rd:

  A character vector of additional columns from rowData(x) to include.
  It defaults to c("gene_id", "SYMBOL").

- type:

  A character string referring to the type of object returned by
  `getDEA()`. It defaults to `DFrame`, but can also take the value of
  `data.frame`

- verbose:

  Logical, whether or not to display warnings. If TRUE,
  warnings/messages will be displayed. If FALSE, the function runs
  silently

- info:

  A character vector, containing contextual information about the
  specified DE analysis. It defaults to NULL

- fea:

  A named list of Functional Enrichment results. Each element can be
  either a data.frame (currently supports results from `topGO`,
  `enrichR`, `gProfiler`, `fgsea`, `gsea`, `DAVID`, and output of
  `GeneTonic` shakers), or an `enrichResult`/`gseaResult` objects
  (currently supports `clusterProfiler`)

- de_name:

  A character string to explicitly specify the name of the de result
  this fea should be linked to. If not provided, the function will
  attempt to match fea names to de results automatically.

- fe_name:

  A character string giving a name to the FE results.

- fea_tool:

  A character string indicating the FEA tool used. It can take any of
  the following values : "topGO", "clusterProfiler", "GeneTonic",
  "DAVID", "gsea", "fgsea", "enrichr", "gProfiler". When not specified,
  it defaults to "auto" and the tool is inferred automatically based on
  the input.

- fea_name:

  Character value, specifying the name of the functional enrichment
  result to add or remove

## Value

Return value varies depending on the individual methods, as described
below.

## Details

DEAs

- `getDEAInfo` and `getDEAInfo<-` are the methods to get and set the
  `dea` information as a whole. These methods return `DeeDeeExperiment`
  objects.

- `getDEANames` returns the names of the available DE contrasts in
  `DeeDeeExperiment` objects.

- `renameDEA` is the method to rename one or multiple DEAs stored in a
  `DeeDeeExperiment` object.

- `addDEA` and `removeDEA` are used to respectively add or remove
  DE-results items. These methods also return `DeeDeeExperiment`
  objects, with updated content in the `dea` slot.

- `getDEA` and `getDEAList` retrieve the DEA information, as well as
  some extra rowData information and provide this as a `DataFrame`
  object (for a specific analysis) or as a list, with one element for
  each reported analysis.

- `addScenarioInfo` is the method to add user defined contextual
  information for a specific DE analysis. It allows users to attach
  free-text notes to a specific DEA results that stored in a
  `DeeDeeExperiment` object. This information can include any other
  relevant information to help document that DEA scenario. This context
  is stored in the `dea` slot under the name `addScenarioInfo`, which is
  not a default element in `dea`.

FEAs

- `getFEAInfo` and `getFEAInfo<-` are the methods to get and set the
  `fea` information as a whole. These methods return `DeeDeeExperiment`
  objects.

- `getFEANames` returns the names of the available enrichment results in
  `DeeDeeExperiment` objects.

- `renameFEA` is the method to rename one or multiple FEAs stored in a
  `DeeDeeExperiment` object.

- `addFEA` and `removeFEA` are used to respectively add or remove
  functional enrichment results items. These methods also return
  `DeeDeeExperiment` objects, with updated content in the `fea` slot.

- `getFEA` is the method to retrieve FE results stored in a
  `DeeDeeExperiment` object for a specific contrast, as a standardized
  format similar to the output of `GeneTonic` shakers.

- `getFEAList` is the method that retrieves FEA results as a list. if
  the `dea_name` is indicated, the method will return only FEAs linked
  to that `dea_name`, otherwise it returns all FEAs in the `fea` slot.

- `linkDEAandFEA` is the method that allows the user to manually link a
  FEA result to a specific DEA result.

- `show` is the method to nicely print out the information of a
  `DeeDeeExperiment` object.

- `summary` is the method to print a summary of the available DE and FE
  results in a `DeeDeeExperiment` object.

## Examples

``` r
data("de_named_list", package = "DeeDeeExperiment")
data("topGO_results_list", package = "DeeDeeExperiment")
library("SummarizedExperiment")

rd_macrophage <- DataFrame(
  gene_id = rownames(de_named_list$ifng_vs_naive)
)
rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rd_macrophage
)

# creating a `DeeDeeExperiment`
dde <- DeeDeeExperiment(
  se_macrophage_noassays,
  de_results = de_named_list
)
dde
#> class: DeeDeeExperiment 
#> dim: 1000 0 
#> metadata(2): singlecontrast version
#> assays(0):
#> rownames(1000): ENSG00000164741 ENSG00000078808 ... ENSG00000273680
#>   ENSG00000073536
#> rowData names(13): gene_id ifng_vs_naive_log2FoldChange ...
#>   salmo_both_pvalue salmo_both_padj
#> colnames: NULL
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(4): ifng_vs_naive, ifngsalmo_vs_naive, salmonella_vs_naive, salmo_both 
#> fea(0): 

new_del <- list(
  ifng2 = de_named_list$ifng_vs_naive,
  ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
)

# add a new (set of) DE result(s)
dde_new <- addDEA(dde, new_del)
dde_new
#> class: DeeDeeExperiment 
#> dim: 1000 0 
#> metadata(2): singlecontrast version
#> assays(0):
#> rownames(1000): ENSG00000164741 ENSG00000078808 ... ENSG00000273680
#>   ENSG00000073536
#> rowData names(19): gene_id ifng_vs_naive_log2FoldChange ...
#>   ifngsalmo2_pvalue ifngsalmo2_padj
#> colnames: NULL
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(6): ifng_vs_naive, ifngsalmo_vs_naive, salmonella_vs_naive, salmo_both, ifng2, ifngsalmo2 
#> fea(0): 

# removing DEAs
dde_removed <- removeDEA(dde, "ifng_vs_naive")
dde_removed
#> class: DeeDeeExperiment 
#> dim: 1000 0 
#> metadata(2): singlecontrast version
#> assays(0):
#> rownames(1000): ENSG00000164741 ENSG00000078808 ... ENSG00000273680
#>   ENSG00000073536
#> rowData names(10): gene_id ifngsalmo_vs_naive_log2FoldChange ...
#>   salmo_both_pvalue salmo_both_padj
#> colnames: NULL
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(3): ifngsalmo_vs_naive, salmonella_vs_naive, salmo_both 
#> fea(0): 

# add a new (set of) FE result(s)
dde_new <- addFEA(dde, fea = topGO_results_list)
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...
#> Found 955 gene sets in `topGOtableResult` object.
#> Converting for usage within the DeeDeeExperiment framework...

# removing FEAs
dde_rem <- removeFEA(dde_new, "ifng_vs_naive")

# display available DEAs
getDEANames(dde)
#> [1] "ifng_vs_naive"       "ifngsalmo_vs_naive"  "salmonella_vs_naive"
#> [4] "salmo_both"         

# display available FEAs
getFEANames(dde)
#> NULL

# print a summary of the available DEAs and FEAs
summary(dde, FDR = 0.01)
#> Error in .Vector_summary(object, ...): unused argument (FDR = 0.01)

# rename DEA
dde_new <- renameDEA(dde_new,
  old_name = "salmonella_vs_naive",
  new_name = "Salmo_vs_Naive_renamed"
)
#> ✔ Renamed DEA entries: "salmonella_vs_naive" to "Salmo_vs_Naive_renamed"

# assign DEA to FEA

dde_new <- linkDEAandFEA(dde_new,
  dea_name = "ifngsalmo_vs_naive",
  fea_name = "ifngsalmo_vs_naive"
)
#> ✔ Assigning DEA: "ifngsalmo_vs_naive" to FEA "ifngsalmo_vs_naive"
```
