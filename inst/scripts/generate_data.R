library("macrophage")
library("DESeq2")
library("SummarizedExperiment")
library("edgeR")
library("limma")

# load macrophage dataset
data(gse, "macrophage")

# DE with DESeq2 ----------------------------------------------------------------

dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]

# set seed for reproducibility
set.seed(42)
# sample randomly for 1k genes
selected_genes <- sample(rownames(dds_macrophage), 1000)
dds_macrophage <- dds_macrophage[selected_genes, ]

dds_macrophage <- DESeq(dds_macrophage)

IFNg_naive <- results(dds_macrophage,
                      contrast = c("condition", "IFNg", "naive"),
                      lfcThreshold = 1, alpha = 0.05
)

save(IFNg_naive, file = "data/dea_IFNg_naive.RData", compress = "xz")


IFNg_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "IFNg"),
                     lfcThreshold = 1, alpha = 0.05
)

save(IFNg_both, file = "data/dea_IFNg_both.RData", compress = "xz")


Salm_naive <- results(dds_macrophage,
                      contrast = c("condition", "SL1344", "naive"),
                      lfcThreshold = 1, alpha = 0.05
)

save(Salm_naive, file = "data/dea_Salm_naive.RData", compress = "xz")


Salm_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "SL1344"),
                     lfcThreshold = 1, alpha = 0.05
)

save(Salm_both, file = "data/dea_Salm_both.RData", compress = "xz")

de_named_list <- list(
  ifng_vs_naive = IFNg_naive,
  ifngsalmo_vs_naive = IFNg_both,
  salmonella_vs_naive = Salm_naive,
  salmo_both = Salm_both
)

save(de_named_list, file = "data/de_named_list.RData", compress = "xz")

# DE with limma ----------------------------------------------------------------

# extract counts and sample metadata
counts <- assays(gse)$counts
sample_info <- colData(gse)

# extract the conditions and cell line info
condition <- factor(sample_info$condition)
line <- factor(sample_info$line)

# create design for DE
design <- model.matrix(~ line + condition)

# create DGE list
dge <- DGEList(counts)

# normalize the counts
dge <- calcNormFactors(dge)

# filter by expression
keep <- rowSums(cpm(dge) >= 10) >= 6
dge <- dge[keep, , keep.lib.sizes=TRUE]

# sample
set.seed(42)
selected_genes <- sample(rownames(dge), 1000)

dge <- dge[selected_genes, , keep.lib.sizes=TRUE]

# transform counts into logCPM
v <- voom(dge, design)

# fitting linear models using weighted least squares for each gene
fit <- lmFit(v, design)

# setup comparisons
contrast_matrix <- makeContrasts(
  IFNgNaive     = conditionIFNg,
  IFNg_both = conditionIFNg_SL1344 - conditionIFNg,
  SalmNaive     = conditionSL1344,
  Salm_both = conditionIFNg_SL1344 - conditionSL1344,
  levels = colnames(design)
)

# apply contrast
fit2 <- contrasts.fit(fit, contrast_matrix)

# empirical Bayes smoothing of standard errors
fit2 <- eBayes(fit2)
de_limma <- fit2 # MArrayLM object

# save
save(de_limma, file = "data/de_limma.RData", compress = "xz")

# DE with edgeR ----------------------------------------------------------------

# extract counts and sample metadata
counts <- assays(gse)$counts
sample_info <- colData(gse)

# extract the conditions and cell line info
condition <- factor(sample_info$condition)
line <- factor(sample_info$line)

# create DGEList object
dge <- DGEList(counts=counts, group=condition)

# normalize the counts
dge <- calcNormFactors(dge)

# filter by expression
keep <- rowSums(cpm(dge) >= 10) >= 6
dge <- dge[keep, , keep.lib.sizes=TRUE]

set.seed(42)
selected_genes <- sample(rownames(dge), 1000)

dge <- dge[selected_genes, , keep.lib.sizes=TRUE]



# create design
design <- model.matrix(~ line + condition)

# estimate dispersion
dge <- estimateDisp(dge, design)

# perform likelihood ratio test
fit <- glmFit(dge, design)

# setup comparisons
contrast_matrix <- makeContrasts(
  IFNgNaive    = conditionIFNg,
  IFNg_both    = conditionIFNg_SL1344 - conditionIFNg,
  SalmNaive    = conditionSL1344,
  Salm_both    = conditionIFNg_SL1344 - conditionSL1344,
  levels = design
)

# DGELRT objects
dge_lrt_IFNg_naive <- glmLRT(fit, contrast = contrast_matrix[, "IFNgNaive"])
dge_lrt_IFNg_both  <- glmLRT(fit, contrast = contrast_matrix[, "IFNg_both"])
dge_lrt_Salm_naive <- glmLRT(fit, contrast = contrast_matrix[, "SalmNaive"])
dge_lrt_Salm_both  <- glmLRT(fit, contrast = contrast_matrix[, "Salm_both"])

# perform exact test
# exact test doesn't handle multi factor models, so we have to subset
# IFNg vs naive


keep_samples <- sample_info$condition %in% c("naive", "IFNg")
dge_sub <- dge[, keep_samples]
# droplevel
group <- droplevels(factor(sample_info[colnames(dge_sub), "condition"]))
dge_sub$samples$group <- group
# renormalizw
dge_sub <- calcNormFactors(dge_sub)
dge_sub <- estimateDisp(dge_sub,design = model.matrix(~group))
dge_exact_IFNg_naive <- exactTest(dge_sub, pair=c("naive", "IFNg")) # DGEExact object

# SL1344 vs naive
keep_samples <- sample_info$condition %in% c("naive", "SL1344")
dge_sub <- dge[, keep_samples]
group <- droplevels(factor(sample_info[colnames(dge_sub), "condition"]))
dge_sub$samples$group <- group
dge_sub <- calcNormFactors(dge_sub)
dge_sub <- estimateDisp(dge_sub,design = model.matrix(~group))
dge_exact_Salm_naive <- exactTest(dge_sub, pair=c("naive", "SL1344"))

# IFNg_SL1344 vs IFNg
keep_samples <- sample_info$condition %in% c("IFNg", "IFNg_SL1344")
dge_sub <- dge[, keep_samples]
group <- droplevels(factor(sample_info[colnames(dge_sub), "condition"]))
dge_sub$samples$group <- group
dge_sub <- calcNormFactors(dge_sub)
dge_exact_IFNg_both <- exactTest(dge_sub, pair=c("IFNg", "IFNg_SL1344"))

# IFNg_SL1344 vs SL1344
keep_samples <- sample_info$condition %in% c("SL1344", "IFNg_SL1344")
dge_sub <- dge[, keep_samples]
group <- droplevels(factor(sample_info[colnames(dge_sub), "condition"]))
dge_sub$samples$group <- group
dge_sub <- calcNormFactors(dge_sub)
dge_sub <- estimateDisp(dge_sub,design = model.matrix(~group))
dge_exact_Salm_both <- exactTest(dge_sub, pair=c("SL1344", "IFNg_SL1344"))

# save
save(dge_lrt_IFNg_naive, file = "data/dgeLRT_IFNg_naive.RData", compress = "xz")
save(dge_lrt_Salm_naive, file = "data/dgeLRT_Salm_naive.RData", compress = "xz")
save(dge_lrt_IFNg_both, file = "data/dgeLRT_IFNg_both.RData", compress = "xz")
save(dge_lrt_Salm_both, file = "data/dgeLRT_Salm_both.RData", compress = "xz")

save(dge_exact_IFNg_naive, file = "data/dgeExact_IFNg_naive.RData", compress = "xz")
save(dge_exact_Salm_naive, file = "data/dgeExact_Salm_naive.RData", compress = "xz")
save(dge_exact_IFNg_both, file = "data/dgeExact_IFNg_both.RData", compress = "xz")
save(dge_exact_Salm_both, file = "data/dgeExact_Salm_both.RData", compress = "xz")


# FE with topGO ----------------------------------------------------------------
library("org.Hs.eg.db")
library("topGO")
FDR = 0.05

de_named_list <- list(
  ifng_vs_naive = IFNg_naive,
  ifngsalmo_vs_naive = IFNg_both,
  salmonella_vs_naive = Salm_naive,
  salmo_both = Salm_both
)

topGO_results <- list()

for (name in names(de_named_list)) {
  de <- de_named_list[[name]]

  topGO_results[[name]] <-
    mosdef::run_topGO(de_container = dds_macrophage,
                      res_de = de,
                      FDR_threshold = FDR,
                      ontology = "BP",
                      add_gene_to_terms = TRUE,
                      mapping = "org.Hs.eg.db")
}

#save
save(topGO_results, file = "data/topGO_results_list.RData", compress = "xz")

# FE with clusterProfiler ----------------------------------------------------------------
library("clusterProfiler")

de_results_2 <- list(ifng_vs_naive = IFNg_naive,
                   salmonella_vs_naive = Salm_naive)

clusterPro_res <- list()

for (name in names(de_results_2)) {
  de <- de_results_2[[name]]

  de_genes <- mosdef::deresult_to_df(de,FDR = FDR)

  clusterPro_res[[name]] <-
    enrichGO(
      gene = rownames(de_genes),
      universe      = rownames(dds_macrophage),
      keyType       = "ENSEMBL",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2)
}

#save
save(clusterPro_res, file = "data/clusterPro_res.RData", compress = "xz")

# FE with enrichR ----------------------------------------------------------------
library("enrichR")
anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds_macrophage), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)

dbs <- c("GO_Biological_Process_2018",
         "KEGG_2019_Human",
         "Reactome_2016")

degenes <- rownames(mosdef::deresult_to_df(Salm_naive, FDR = FDR))
deg_symbols <- anno_df[degenes, "gene_name"]
# remove nas
deg_symbols <- deg_symbols[!is.na(deg_symbols)]


enrichr_res <- enrichr(deg_symbols, dbs)

save(enrichr_res, file = "data/enrichr_res.RData", compress = "xz")

# FE with gprofiler2 ----------------------------------------------------------------
library("gprofiler2")
degenes <- rownames(mosdef::deresult_to_df(Salm_naive, FDR = FDR))
deg_symbols <- anno_df[degenes, "gene_name"]
# remove nas
deg_symbols <- deg_symbols[!is.na(deg_symbols)]
gost_res <- gost(
  query = deg_symbols,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = TRUE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = "g_SCS",
  domain_scope = "annotated",
  numeric_ns = "",
  sources = "GO:BP",
  as_short_link = FALSE)

save(gost_res, file = "data/gost_res.RData", compress = "xz")

# fgseaRes object ---------------------------------------------------------
library("dplyr")
library("tibble")
library("fgsea")
IFNg_naive$SYMBOL <- anno_df[rownames(IFNg_naive), "gene_name"]
res2 <- IFNg_naive %>%
  as.data.frame() %>%
  dplyr::select(SYMBOL, stat)
de_ranks <- deframe(res2)
de_ranks <- de_ranks[!is.na(names(de_ranks))]
head(de_ranks, 20)
pathways_gmtfile <- gmtPathways("../msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt")
fgseaRes <- fgsea(pathways = pathways_gmtfile,
                  stats = de_ranks,
                  nperm=100000)
fgseaRes <- fgseaRes %>%
  arrange(desc(NES))

save(fgseaRes, file = "data/fgseaRes.RData", compress = "xz")



# gseaResult object ---------------------------------------------------------
sorted_genes <- sort(
    setNames(IFNg_naive$log2FoldChange,
             IFNg_naive$SYMBOL),
    decreasing = TRUE
  )
sorted_genes <- sorted_genes[!is.na(names(sorted_genes))]
gsea_res <- gseGO(
    geneList = sorted_genes,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = TRUE
  )

save(gsea_res, file = "data/gsea_res.RData", compress = "xz")
