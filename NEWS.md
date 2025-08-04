# DeeDeeExperiment 0.99.2

* Addressed the points raised in the Bioc review

* Expended the supported DEAs (as simple data.frame object that contains at least
a column named `log2FoldChange`, `pvalue`, and `padj`)

# DeeDeeExperiment 0.99.1

* `DeeDeeExperiment` now extends the `SingleCellExperiment` class, accommodating
both bulk and single-cell analysis workflows

* Ready for Bioconductor review!

# DeeDeeExperiment 0.99.0

* Ready for Bioconductor submission!

# DeeDeeExperiment 0.6.0

* Full vignette now available for user guidance

# DeeDeeExperiment 0.5.0

* Reducing the dependency stack to make the class/package as slim as it can be,
with a focus on the full interoperability with SE (therefore now a full dependency)

# DeeDeeExperiment 0.4.0

* New methods to rename DEA and FEA elements in `dde` objects, get all FEAs for
a specific DEA, assign DEA to FEA

* Method dispatch is now only base on x being `DeeDeeExperiment`, other arguments
are no longer part of the method signature

* expended supported FEAs: now takes results from `topGO`, `clusterProfiler`,
`enrichR`, `gProfiler`, `fgsea`, `gsea`, `DAVID`, and output of `GeneTonic` shakers

# DeeDeeExperiment 0.3.0

## New features

* Methods to retrieve and manage the `fea` slot content are established

* supported FEAs: topGO (data.frame) & enrichResult objects

* A summary method is included to provide a quick overview of stored DEA and
FEA results

# DeeDeeExperiment 0.2.0

## New features

* `dea()` has been refactored for clarity:
  - Now it returns only the main results table, acting as `get_dea_df()`
  - A new function, `dea_info()`, has been added to return additional information
  previously returned by `dea()`

* Initial implementation of the `fea` slot

# DeeDeeExperiment 0.1.0

## New features

* The `DeeDeeExperiment` class is now being developed in its own standalone
package, following its separation from the original `DeeDee_legacy` repository:
https://github.com/imbeimainz/DeeDee_legacy
