# How to use DeeDeeExperiment with single-cell data

## Differential State (DS) analysis with `muscat` & `DeeDeeExperiment`

### `DeeDeeExperiment` on the `Kang` dataset

In this vignette, we illustrate how to integrate Differential Expression
Analysis (DEA) results generated with the
*[muscat](https://bioconductor.org/packages/3.24/muscat)* framework into
a
*[DeeDeeExperiment](https://bioconductor.org/packages/3.24/DeeDeeExperiment)*
object. As an example, we use a publicly available dataset from Kang, et
al. “Multiplexed droplet single-cell RNA-sequencing using natural
genetic variation”, published in Nature Biotechnology, December 2017
[Kang et al. (2017)](https://doi.org/10.1038/nbt.4042)

The data is made available via the
*[ExperimentHub](https://bioconductor.org/packages/3.24/ExperimentHub)*
Bioconductor package as a
*[SingleCellExperiment](https://bioconductor.org/packages/3.24/SingleCellExperiment)*
object containing scRNA-seq data from PBMCs obtained from 8 lupus
patients before and after IFNβ stimulation.

For demonstration purposes we adapt parts of the code from the original
*[muscat](https://bioconductor.org/packages/3.24/muscat)*
[vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html#differential-state-ds-analysis))

``` r

library("DeeDeeExperiment")
library("ExperimentHub")
library("scater")
library("muscat")
library("limma")
```

We begin by creating an `ExperimentHub` instance, which provides access
to curated datasets stored in the Bioconductor cloud. Using
[`query()`](https://rdrr.io/pkg/AnnotationHub/man/AnnotationHub-class.html),
we filter available records for entries matching the keyword “Kang”, and
then load the dataset of interest using its accession ID “EH2259”.

``` r

# retrieve the data
eh <- ExperimentHub()
query(eh, "Kang")
#> ExperimentHub with 1 record
#> # snapshotDate(): 2026-07-16
#> # names(): EH2259
#> # package(): muscData
#> # $dataprovider: GEO
#> # $species: Homo sapiens
#> # $rdataclass: SingleCellExperiment
#> # $rdatadateadded: 2019-04-15
#> # $title: Kang18_8vs8
#> # $description: Droplet-based scRNA-seq PBMC data from 8 Lupus patients befo...
#> # $taxonomyid: 9606
#> # $genome: NA
#> # $sourcetype: tar.gz
#> # $sourceurl: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
#> # $sourcesize: NA
#> # $tags: c("ExperimentHub", "ExperimentData", "ExpressionData",
#> #   "SingleCellData", "ImmunoOncologyData", "Homo_sapiens_Data") 
#> # retrieve record with 'object[["EH2259"]]'
sce <- eh[["EH2259"]]
```

Before running muscat, we perform standard single-cell preprocessing
steps: removing undetected genes, filtering low-quality cells using
*[scater](https://bioconductor.org/packages/3.24/scater)*, and removing
lowly expressed genes.

``` r

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
```

``` r

# calculate per-cell quality control (QC) metrics
qc <- perCellQCMetrics(sce)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)
#> [1] 18890 26820
```

``` r

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
#> [1]  7118 26820
```

``` r

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
```

We then prepare the data for `muscat`. The package expects a certain
format of the input SCE. Specifically, the following cell metadata
(colData) columns have to be provided:

- `sample_id` : unique sample identifiers
- `cluster_id` : subpopulation (cluster) assignments
- `group_id` : experimental group/condition

``` r

# data preparation
sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce,
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns
#> class: SingleCellExperiment 
#> dim: 7118 26820 
#> metadata(1): experiment_info
#> assays(2): counts logcounts
#> rownames(7118): NOC2L HES4 ... S100B PRMT2
#> rowData names(2): ENSEMBL SYMBOL
#> colnames(26820): AAACATACAATGCC-1 AAACATACATTTCC-1 ... TTTGCATGGTTTGG-1
#>   TTTGCATGTCTTAC-1
#> colData names(3): cluster_id sample_id group_id
#> reducedDimNames(1): TSNE
#> mainExpName: NULL
#> altExpNames(0):
```

We compute UMAP for visualization. The dataset already includes
precomputed TSNE coordinates, so we only run UMAP here.

``` r

# compute UMAP using 1st 20 PCs
sce <- runUMAP(sce, pca = 20)
```

Then we aggregate measurements for each sample (in each cluster) to
obtain pseudobulk data

``` r

# aggregate by cell type
pb <- aggregateData(
  sce,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)

assayNames(pb)
#> [1] "B cells"           "CD14+ Monocytes"   "CD4 T cells"      
#> [4] "CD8 T cells"       "Dendritic cells"   "FCGR3A+ Monocytes"
#> [7] "Megakaryocytes"    "NK cells"
```

And construct the contrast matrix

``` r

# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("stim-ctrl", levels = mm)
```

With the pseudobulk data assembled, we can now test for differential
state (DS) using `pbDS`

``` r

# run DS analysis
muscat_res <- pbDS(pb, design = mm, contrast = contrast)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%

names(muscat_res$table[["stim-ctrl"]])
#> [1] "B cells"           "CD14+ Monocytes"   "CD4 T cells"      
#> [4] "CD8 T cells"       "Dendritic cells"   "FCGR3A+ Monocytes"
#> [7] "Megakaryocytes"    "NK cells"
```

Now we integrate the output of muscat in
[`muscat_list_for_dde()`](../reference/muscat_list_for_dde.md) to
transform it into a format accepted by `DeeDeeExperiment`

``` r

# preparing the results as muscat list
muscat_list <- muscat_list_for_dde(res = list(`stim-ctrl` = muscat_res),
                                   padj_col = "p_adj.loc")
```

Finally the results can be directly added to a new `dde` object, or to
an existing one.

``` r

# create dde
dde <- DeeDeeExperiment(sce = sce,
                        de_results = muscat_list)
dde
#> class: DeeDeeExperiment 
#> dim: 7118 26820 
#> metadata(3): experiment_info singlecontrast version
#> assays(2): counts logcounts
#> rownames(7118): NOC2L HES4 ... S100B PRMT2
#> rowData names(26): ENSEMBL SYMBOL ... stim-ctrl_NK cells_pvalue
#>   stim-ctrl_NK cells_padj
#> colnames(26820): AAACATACAATGCC-1 AAACATACATTTCC-1 ... TTTGCATGGTTTGG-1
#>   TTTGCATGTCTTAC-1
#> colData names(3): cluster_id sample_id group_id
#> reducedDimNames(2): TSNE UMAP
#> mainExpName: NULL
#> altExpNames(0):
#> dea(8): stim-ctrl_B cells, stim-ctrl_CD14+ Monocytes, stim-ctrl_CD4 T cells, stim-ctrl_CD8 T cells, stim-ctrl_Dendritic cells, stim-ctrl_FCGR3A+ Monocytes, stim-ctrl_Megakaryocytes, stim-ctrl_NK cells 
#> fea(0):
```

As for the other `DeeDeeExperiment` objects, we can call some specific
methods to extract/retrieve/integrate some information.  
We can extract the names of the DEA included:

``` r

# check DEAs
getDEANames(dde)
#> [1] "stim-ctrl_B cells"           "stim-ctrl_CD14+ Monocytes"  
#> [3] "stim-ctrl_CD4 T cells"       "stim-ctrl_CD8 T cells"      
#> [5] "stim-ctrl_Dendritic cells"   "stim-ctrl_FCGR3A+ Monocytes"
#> [7] "stim-ctrl_Megakaryocytes"    "stim-ctrl_NK cells"
```

Also, we can directly retrieve the content itself of each DEA by typing

``` r

# retrieve results
getDEA(dde,
       dea_name = "stim-ctrl_NK cells")|> head()
#> DataFrame with 6 rows and 3 columns
#>          stim-ctrl_NK cells_log2FoldChange stim-ctrl_NK cells_pvalue
#>                                  <numeric>                 <numeric>
#> NOC2L                                   NA                        NA
#> HES4                                    NA                        NA
#> ISG15                             4.745959               5.97846e-20
#> TNFRSF18                         -0.555606               7.83911e-02
#> TNFRSF4                          -0.189651               5.60070e-01
#> SDF4                             -0.493880               2.42677e-02
#>          stim-ctrl_NK cells_padj
#>                        <numeric>
#> NOC2L                         NA
#> HES4                          NA
#> ISG15                2.58867e-17
#> TNFRSF18             1.93188e-01
#> TNFRSF4              7.22213e-01
#> SDF4                 8.27394e-02

getDEA(dde,
       dea_name = "stim-ctrl_CD14+ Monocytes",
       format = "original") |> head()
#>              gene      cluster_id log2FoldChange    logCPM          F
#> HES4         HES4 CD14+ Monocytes      6.5479774  7.994510 284.482163
#> ISG15       ISG15 CD14+ Monocytes      7.0760135 14.772034 235.655412
#> SDF4         SDF4 CD14+ Monocytes     -0.6895721  5.444335  13.433991
#> UBE2J2     UBE2J2 CD14+ Monocytes     -0.7693502  5.346941  21.212665
#> CPSF3L     CPSF3L CD14+ Monocytes     -0.6216370  4.334590   6.703139
#> AURKAIP1 AURKAIP1 CD14+ Monocytes     -0.5380074  7.386092  30.591918
#>                pvalue         padj    p_adj.glb  contrast
#> HES4     6.255203e-14 2.600377e-12 3.314853e-12 stim-ctrl
#> ISG15    5.981632e-13 1.571425e-11 2.315585e-11 stim-ctrl
#> SDF4     1.430674e-03 3.227334e-03 5.824784e-03 stim-ctrl
#> UBE2J2   1.507560e-04 4.218817e-04 8.010871e-04 stim-ctrl
#> CPSF3L   1.713350e-02 2.922635e-02 4.740595e-02 stim-ctrl
#> AURKAIP1 1.682349e-05 6.161011e-05 1.169786e-04 stim-ctrl
```

General info on the DEA performed can be shown with
[`getDEAInfo()`](../reference/DeeDeeExperiment-methods.md)

``` r

# retrieving the DEA information
dea_name <- "stim-ctrl_B cells"
getDEAInfo(dde)[[dea_name]][["package"]]
#> [1] "muscat"
getDEAInfo(dde)[[dea_name]][["package_version"]]
#> [1] "1.27.4"
```

To add some information on the scenario under investigation, we can use
[`addScenarioInfo()`](../reference/DeeDeeExperiment-methods.md). This
can be e.g. later processed as a contextually relevant bit if a Large
Language Model is used to interact with this object.

``` r

# adding info on the scenario under investigation
dde <- addScenarioInfo(dde,
                       dea_name = "stim-ctrl_Dendritic cells",
                       info = "This result contains the output of pseudobulk DE analysis performed on dendritic cells, comparing untreated samples to those stimulated with IFNβ")
```

As usual, the [`summary()`](https://rdrr.io/r/base/summary.html) method
can be called to obtain a quick overview on all performed analyses.

``` r

summary(dde, show_scenario_info = TRUE)
#> DE Results Summary:
#>                     DEA_name   Up Down  FDR
#>            stim-ctrl_B cells  346  269 0.05
#>    stim-ctrl_CD14+ Monocytes 1157 1202 0.05
#>        stim-ctrl_CD4 T cells  668  455 0.05
#>        stim-ctrl_CD8 T cells  182  104 0.05
#>    stim-ctrl_Dendritic cells  122   46 0.05
#>  stim-ctrl_FCGR3A+ Monocytes  502  399 0.05
#>     stim-ctrl_Megakaryocytes   28    3 0.05
#>           stim-ctrl_NK cells  269  174 0.05
#> 
#> No FEA results stored.
#> 
#> Scenario Info:
#>  - stim-ctrl_Dendritic cells :
#>  This result contains the output of pseudobulk DE analysis performed on
#>   dendritic cells, comparing untreated samples to those stimulated with IFNβ 
#>  
#> 
#> No scenario info for: stim-ctrl_B cells, stim-ctrl_CD14+ Monocytes, stim-ctrl_CD4 T cells, stim-ctrl_CD8 T cells, stim-ctrl_FCGR3A+ Monocytes, stim-ctrl_Megakaryocytes, stim-ctrl_NK cells
```

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.4
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: UTC
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] limma_3.69.2                muscat_1.27.4              
#>  [3] scater_1.41.2               ggplot2_4.0.3              
#>  [5] scuttle_1.23.1              ExperimentHub_3.3.1        
#>  [7] AnnotationHub_4.3.2         BiocFileCache_3.3.0        
#>  [9] dbplyr_2.6.0                DeeDeeExperiment_1.3.0     
#> [11] SingleCellExperiment_1.35.2 SummarizedExperiment_1.43.0
#> [13] Biobase_2.73.1              GenomicRanges_1.65.1       
#> [15] Seqinfo_1.3.0               IRanges_2.47.2             
#> [17] S4Vectors_0.51.5            BiocGenerics_0.59.10       
#> [19] generics_0.1.4              MatrixGenerics_1.25.0      
#> [21] matrixStats_1.5.0           BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23         splines_4.6.1            bitops_1.0-9            
#>   [4] filelock_1.0.3           tibble_3.3.1             lifecycle_1.0.5         
#>   [7] httr2_1.3.0              Rdpack_2.6.6             edgeR_4.11.4            
#>  [10] doParallel_1.0.17        lattice_0.22-9           MASS_7.3-66             
#>  [13] backports_1.5.1          magrittr_2.0.5           sass_0.4.10             
#>  [16] rmarkdown_2.31           jquerylib_0.1.4          yaml_2.3.12             
#>  [19] otel_0.2.0               DBI_1.3.0                minqa_1.2.8             
#>  [22] RColorBrewer_1.1-3       abind_1.4-8              EnvStats_3.1.0          
#>  [25] glmmTMB_1.1.14           purrr_1.2.2              rappdirs_0.3.4          
#>  [28] sandwich_3.1-2           circlize_0.4.18          pbkrtest_0.5.5          
#>  [31] ggrepel_0.9.8            irlba_2.3.7              RSpectra_0.16-2         
#>  [34] pkgdown_2.2.1.9000       codetools_0.2-20         DelayedArray_0.39.3     
#>  [37] tidyselect_1.2.1         shape_1.4.6.1            farver_2.1.2            
#>  [40] lme4_2.0-6               ScaledMatrix_1.21.0      viridis_0.6.5           
#>  [43] jsonlite_2.0.0           GetoptLong_1.1.1         BiocNeighbors_2.7.2     
#>  [46] iterators_1.0.14         systemfonts_1.3.2        foreach_1.5.2           
#>  [49] tools_4.6.1              progress_1.2.3           ragg_1.5.2              
#>  [52] Rcpp_1.1.2               blme_1.0-7               glue_1.8.1              
#>  [55] gridExtra_2.3.1          SparseArray_1.13.2       BiocBaseUtils_1.15.1    
#>  [58] xfun_0.60                mgcv_1.9-4               DESeq2_1.53.2           
#>  [61] dplyr_1.2.1              withr_3.0.3              numDeriv_2016.8-1.1     
#>  [64] BiocManager_1.30.27      fastmap_1.2.0            boot_1.3-32             
#>  [67] caTools_1.18.4           digest_0.6.39            rsvd_1.0.5              
#>  [70] R6_2.6.1                 textshaping_1.0.5        colorspace_2.1-3        
#>  [73] gtools_3.9.5             RSQLite_3.53.3           RhpcBLASctl_0.23-42     
#>  [76] tidyr_1.3.2              variancePartition_1.43.1 corpcor_1.6.10          
#>  [79] prettyunits_1.2.0        httr_1.4.8               htmlwidgets_1.6.4       
#>  [82] S4Arrays_1.13.0          uwot_0.2.4               pkgconfig_2.0.3         
#>  [85] gtable_0.3.6             blob_1.3.0               ComplexHeatmap_2.29.0   
#>  [88] S7_0.2.2                 XVector_0.53.0           remaCor_0.0.20          
#>  [91] htmltools_0.5.9          bookdown_0.47            TMB_1.9.23              
#>  [94] clue_0.3-68              scales_1.4.0             png_0.1-9               
#>  [97] fANCOVA_0.6-1            reformulas_0.4.4         knitr_1.51              
#> [100] reshape2_1.4.5           rjson_0.2.23             nlme_3.1-170            
#> [103] curl_7.1.0               nloptr_2.2.1             cachem_1.1.0            
#> [106] zoo_1.8-15               GlobalOptions_0.1.4      stringr_1.6.0           
#> [109] KernSmooth_2.23-26       BiocVersion_3.24.0       parallel_4.6.1          
#> [112] vipor_0.4.7              AnnotationDbi_1.75.2     desc_1.4.3              
#> [115] pillar_1.11.1            grid_4.6.1               vctrs_0.7.3             
#> [118] gplots_3.3.0             BiocSingular_1.29.0      beachmat_2.29.0         
#> [121] cluster_2.1.8.2          beeswarm_0.4.0           evaluate_1.0.5          
#> [124] mvtnorm_1.4-2            cli_3.6.6                locfit_1.5-9.12         
#> [127] compiler_4.6.1           rlang_1.3.0              crayon_1.5.3            
#> [130] plyr_1.8.9               fs_2.1.0                 ggbeeswarm_0.7.3        
#> [133] stringi_1.8.7            writexl_1.5.4            viridisLite_0.4.3       
#> [136] BiocParallel_1.47.0      lmerTest_3.2-1           Biostrings_2.81.5       
#> [139] aod_1.3.3                scrapper_1.7.3           Matrix_1.7-5            
#> [142] hms_1.1.4                bit64_4.8.2              KEGGREST_1.53.5         
#> [145] statmod_1.5.2            rbibutils_2.4.1          broom_1.0.13            
#> [148] memoise_2.0.1            bslib_0.11.0             bit_4.6.0
```

## References

Kang, Hyun Min, Meena Subramaniam, Sasha Targ, et al. 2017. “Multiplexed
Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation.”
*Nature Biotechnology* 36 (1): 89–94.
<https://doi.org/10.1038/nbt.4042>.
