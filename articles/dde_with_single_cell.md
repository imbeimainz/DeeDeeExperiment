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
#> # snapshotDate(): 2026-05-19
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
#> ISG15                             4.715080               6.94235e-20
#> TNFRSF18                         -0.583846               6.49706e-02
#> TNFRSF4                          -0.219549               5.00289e-01
#> SDF4                             -0.523241               1.80659e-02
#>          stim-ctrl_NK cells_padj
#>                        <numeric>
#> NOC2L                         NA
#> HES4                          NA
#> ISG15                3.00604e-17
#> TNFRSF18             1.63560e-01
#> TNFRSF4              6.67095e-01
#> SDF4                 6.37857e-02

getDEA(dde,
       dea_name = "stim-ctrl_CD14+ Monocytes",
       format = "original") |> head()
#>              gene      cluster_id log2FoldChange    logCPM          F
#> HES4         HES4 CD14+ Monocytes      6.5041181  7.975074 282.968142
#> ISG15       ISG15 CD14+ Monocytes      7.0258237 14.751515 229.037506
#> SDF4         SDF4 CD14+ Monocytes     -0.7256656  5.446124  15.236040
#> UBE2J2     UBE2J2 CD14+ Monocytes     -0.8083617  5.351035  22.886511
#> CPSF3L     CPSF3L CD14+ Monocytes     -0.6600343  4.337862   7.343274
#> AURKAIP1 AURKAIP1 CD14+ Monocytes     -0.5728129  7.389181  34.211958
#>                pvalue         padj    p_adj.glb  contrast
#> HES4     6.912599e-14 2.959294e-12 3.687094e-12 stim-ctrl
#> ISG15    8.146888e-13 2.139849e-11 3.052753e-11 stim-ctrl
#> SDF4     8.128918e-04 1.921981e-03 3.528924e-03 stim-ctrl
#> UBE2J2   9.898603e-05 2.898329e-04 5.545317e-04 stim-ctrl
#> CPSF3L   1.314382e-02 2.297739e-02 3.751614e-02 stim-ctrl
#> AURKAIP1 8.127731e-06 3.167114e-05 6.047946e-05 stim-ctrl
```

General info on the DEA performed can be shown with
[`getDEAInfo()`](../reference/DeeDeeExperiment-methods.md)

``` r

# retrieving the DEA information
dea_name <- "stim-ctrl_B cells"
getDEAInfo(dde)[[dea_name]][["package"]]
#> [1] "muscat"
getDEAInfo(dde)[[dea_name]][["package_version"]]
#> [1] "1.25.4"
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
#>            stim-ctrl_B cells  344  270 0.05
#>    stim-ctrl_CD14+ Monocytes 1118 1282 0.05
#>        stim-ctrl_CD4 T cells  640  511 0.05
#>        stim-ctrl_CD8 T cells  178  115 0.05
#>    stim-ctrl_Dendritic cells  113   70 0.05
#>  stim-ctrl_FCGR3A+ Monocytes  497  411 0.05
#>     stim-ctrl_Megakaryocytes   29    2 0.05
#>           stim-ctrl_NK cells  254  200 0.05
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
#> R version 4.6.0 (2026-04-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Sequoia 15.7.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] muscData_1.27.0             limma_3.69.0               
#>  [3] muscat_1.25.4               scater_1.41.1              
#>  [5] ggplot2_4.0.3               scuttle_1.21.6             
#>  [7] ExperimentHub_3.3.0         AnnotationHub_4.1.0        
#>  [9] BiocFileCache_3.3.0         dbplyr_2.5.2               
#> [11] DeeDeeExperiment_1.3.0      SingleCellExperiment_1.35.0
#> [13] SummarizedExperiment_1.43.0 Biobase_2.73.1             
#> [15] GenomicRanges_1.65.0        Seqinfo_1.3.0              
#> [17] IRanges_2.45.0              S4Vectors_0.51.1           
#> [19] BiocGenerics_0.59.0         generics_0.1.4             
#> [21] MatrixGenerics_1.25.0       matrixStats_1.5.0          
#> [23] BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23         splines_4.6.0            bitops_1.0-9            
#>   [4] filelock_1.0.3           tibble_3.3.1             lifecycle_1.0.5         
#>   [7] httr2_1.2.2              Rdpack_2.6.6             edgeR_4.9.9             
#>  [10] doParallel_1.0.17        lattice_0.22-9           MASS_7.3-65             
#>  [13] backports_1.5.1          magrittr_2.0.5           sass_0.4.10             
#>  [16] rmarkdown_2.31           jquerylib_0.1.4          yaml_2.3.12             
#>  [19] otel_0.2.0               DBI_1.3.0                minqa_1.2.8             
#>  [22] RColorBrewer_1.1-3       multcomp_1.4-30          abind_1.4-8             
#>  [25] EnvStats_3.1.0           glmmTMB_1.1.14           purrr_1.2.2             
#>  [28] TH.data_1.1-5            rappdirs_0.3.4           sandwich_3.1-1          
#>  [31] circlize_0.4.18          pbkrtest_0.5.5           ggrepel_0.9.8           
#>  [34] irlba_2.3.7              RSpectra_0.16-2          pkgdown_2.2.0           
#>  [37] codetools_0.2-20         DelayedArray_0.39.1      tidyselect_1.2.1        
#>  [40] shape_1.4.6.1            farver_2.1.2             lme4_2.0-1              
#>  [43] ScaledMatrix_1.21.0      viridis_0.6.5            jsonlite_2.0.0          
#>  [46] GetoptLong_1.1.1         BiocNeighbors_2.7.1      survival_3.8-6          
#>  [49] iterators_1.0.14         emmeans_2.0.3            systemfonts_1.3.2       
#>  [52] foreach_1.5.2            tools_4.6.0              progress_1.2.3          
#>  [55] ragg_1.5.2               Rcpp_1.1.1-1.1           blme_1.0-7              
#>  [58] glue_1.8.1               gridExtra_2.3            SparseArray_1.11.13     
#>  [61] xfun_0.57                mgcv_1.9-4               DESeq2_1.51.7           
#>  [64] dplyr_1.2.1              withr_3.0.2              numDeriv_2016.8-1.1     
#>  [67] BiocManager_1.30.27      fastmap_1.2.0            boot_1.3-32             
#>  [70] caTools_1.18.3           digest_0.6.39            rsvd_1.0.5              
#>  [73] R6_2.6.1                 estimability_1.5.1       textshaping_1.0.5       
#>  [76] colorspace_2.1-2         gtools_3.9.5             dichromat_2.0-0.1       
#>  [79] RSQLite_3.52.0           RhpcBLASctl_0.23-42      tidyr_1.3.2             
#>  [82] variancePartition_1.41.5 corpcor_1.6.10           prettyunits_1.2.0       
#>  [85] httr_1.4.8               htmlwidgets_1.6.4        S4Arrays_1.13.0         
#>  [88] uwot_0.2.4               pkgconfig_2.0.3          gtable_0.3.6            
#>  [91] blob_1.3.0               ComplexHeatmap_2.29.0    S7_0.2.2                
#>  [94] XVector_0.53.0           remaCor_0.0.20           htmltools_0.5.9         
#>  [97] bookdown_0.46            TMB_1.9.21               clue_0.3-68             
#> [100] scales_1.4.0             png_0.1-9                fANCOVA_0.6-1           
#> [103] reformulas_0.4.4         knitr_1.51               rstudioapi_0.18.0       
#> [106] reshape2_1.4.5           rjson_0.2.23             coda_0.19-4.1           
#> [109] nlme_3.1-169             curl_7.1.0               nloptr_2.2.1            
#> [112] cachem_1.1.0             zoo_1.8-15               GlobalOptions_0.1.4     
#> [115] stringr_1.6.0            KernSmooth_2.23-26       BiocVersion_3.24.0      
#> [118] parallel_4.6.0           vipor_0.4.7              AnnotationDbi_1.75.0    
#> [121] desc_1.4.3               pillar_1.11.1            grid_4.6.0              
#> [124] vctrs_0.7.3              gplots_3.3.0             BiocSingular_1.29.0     
#> [127] beachmat_2.29.0          xtable_1.8-8             cluster_2.1.8.2         
#> [130] beeswarm_0.4.0           evaluate_1.0.5           mvtnorm_1.3-7           
#> [133] cli_3.6.6                locfit_1.5-9.12          compiler_4.6.0          
#> [136] rlang_1.2.0              crayon_1.5.3             plyr_1.8.9              
#> [139] fs_2.1.0                 ggbeeswarm_0.7.3         stringi_1.8.7           
#> [142] writexl_1.5.4            viridisLite_0.4.3        BiocParallel_1.47.0     
#> [145] lmerTest_3.2-1           Biostrings_2.81.1        aod_1.3.3               
#> [148] Matrix_1.7-5             hms_1.1.4                bit64_4.8.0             
#> [151] KEGGREST_1.53.0          statmod_1.5.1            rbibutils_2.4.1         
#> [154] broom_1.0.12             memoise_2.0.1            bslib_0.10.0            
#> [157] bit_4.6.0
```

## References

Kang, Hyun Min, Meena Subramaniam, Sasha Targ, et al. 2017. “Multiplexed
Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation.”
*Nature Biotechnology* 36 (1): 89–94.
<https://doi.org/10.1038/nbt.4042>.
