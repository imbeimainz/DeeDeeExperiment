# The DeeDeeExperiment class

The `DeeDeeExperiment` class is integrate and manage omics analysis
results. It inherits from the `SingleCellExperiment` class, and
additionally stores DE-related/functional enrichment information via
dedicated slots and `rowData`.

## Usage

``` r
DeeDeeExperiment(
  sce = SingleCellExperiment(),
  de_results = NULL,
  enrich_results = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object, that will be used as a scaffold to
  store the DE related information.

- de_results:

  A named list of DE results, in any of the formats supported by the
  package (currently: results from `DESeq2`, `edgeR`, `limma`,
  `muscat`).

- enrich_results:

  A named list of functional enrichment results. Each element can be
  either a data.frame (currently supports results from `topGO`,
  `enrichR`, `gProfiler`, `fgsea`, `gsea`, `DAVID`, and output of
  `GeneTonic` shakers), or an `enrichResult`/`gseaResult` objects
  (currently supports `clusterProfiler`)

## Value

A `DeeDeeExperiment` object.

## Details

The `sce` parameter can be optionally left unspecified. If this is the
case, the resulting `DeeDeeExperiment` object will contain as features
the ones specified by the provided components of the object supplied via
the `de_results` parameter.

The conversion of the components of the `de_results` list will be
handled via conversion functions to uniform the names and set of
information which will be stored in the returned `DeeDeeExperiment`
object. The names of the list will be used to define the `contrasts` for
the different DE analyses included, which will determine the way to
access the information stored in the `dea` slot of the
`DeeDeeExperiment` object.

The content of the `enrich_results` provided by the user will be
validated to ensure that it is properly formatted and correctly named.
The FE tool can be automatically detected, and based on that, the
appropriate shaking method is used to return a standardized format of
the FEA results. The names of the list will be used to attempt to
associate each enrichment result with a corresponding DE contrast stored
in the `DeeDeeExperiment` object, but it also can be defined by the
user.

Since a `DeeDeeExperiment` is also a `SummarizedExperiment` object, it
can be seamlessly provided downstream for visualization and in-depth
exploration to packages such as `iSEE` or similar.

## Slots

- `dea`:

  This slot is designed to hold the DE-related information. This is
  internally being created upon importing from the list of DE results
  objects, provided when instantiating the DeeDeeExperiment.

- `fea`:

  This slot is designed to hold Functional Enrichment related
  information.

## Author

Najla Abassi, Lea Schwarz, and Federico Marini

## Examples

``` r
data("de_named_list", package = "DeeDeeExperiment")

dde_onlyde <- DeeDeeExperiment(
  de_results = de_named_list
)
#> ℹ creating a mock SCE from the rows of the DE result objects, if available

# or, with a SE object as support - even without assay data available
library("SummarizedExperiment")

rd_macrophage <- DataFrame(
  gene_id = rownames(de_named_list$ifng_vs_naive)
)
rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rd_macrophage
)

dde <- DeeDeeExperiment(
  se_macrophage_noassays,
  de_results = de_named_list
)
```
