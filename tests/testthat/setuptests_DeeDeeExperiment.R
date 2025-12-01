suppressPackageStartupMessages(
  library("SummarizedExperiment")
)

data("de_named_list", package = "DeeDeeExperiment")

rd_macrophage <- DataFrame(
  gene_id = rownames(de_named_list$ifng_vs_naive)
)
rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rd_macrophage
)

names(de_named_list)

data("de_limma", package = "DeeDeeExperiment")

data("DGEExact_IFNg_both", package = "DeeDeeExperiment")

data("topGO_results_list", package = "DeeDeeExperiment")

data("enrichr_res", package = "DeeDeeExperiment")

data("clusterPro_res", package = "DeeDeeExperiment")

data("gost_res", package = "DeeDeeExperiment")

data("fgseaRes", package = "DeeDeeExperiment")

data("gsea_res", package = "DeeDeeExperiment")

data("muscat_pbDS_res", package = "DeeDeeExperiment")
