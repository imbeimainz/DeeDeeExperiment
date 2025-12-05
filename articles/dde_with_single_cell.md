# How to use DeeDeeExperiment with single-cell data

## Differential State (DS) analysis with `muscat` & `DeeDeeExperiment`

### `DeeDeeExperiment` on the `Kang` dataset

In this vignette, we illustrate how to integrate Differential Expression
Analysis (DEA) results generated with the
*[muscat](https://bioconductor.org/packages/3.23/muscat)* framework into
a
*[DeeDeeExperiment](https://bioconductor.org/packages/3.23/DeeDeeExperiment)*
object. As an example, we use a publicly available dataset from Kang, et
al. “Multiplexed droplet single-cell RNA-sequencing using natural
genetic variation”, published in Nature Biotechnology, December 2017
[Kang et al. (2017)](https://doi.org/10.1038/nbt.4042)

The data is made available via the
*[ExperimentHub](https://bioconductor.org/packages/3.23/ExperimentHub)*
Bioconductor package as a
*[SingleCellExperiment](https://bioconductor.org/packages/3.23/SingleCellExperiment)*
object containing scRNA-seq data from PBMCs obtained from 8 lupus
patients before and after IFNβ stimulation.

For demonstration purposes we adapt parts of the code from the original
*[muscat](https://bioconductor.org/packages/3.23/muscat)*
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
#> ExperimentHub with 3 records
#> # snapshotDate(): 2025-12-03
#> # $dataprovider: NCI_GDC, GEO
#> # $species: Homo sapiens
#> # $rdataclass: character, SingleCellExperiment, BSseq
#> # additional mcols(): taxonomyid, genome, description,
#> #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#> #   rdatapath, sourceurl, sourcetype 
#> # retrieve records with, e.g., 'object[["EH1661"]]' 
#> 
#>            title                                               
#>   EH1661 | Whole Genome Bisulfit Sequencing Data for 47 samples
#>   EH1662 | Whole Genome Bisulfit Sequencing Data for 47 samples
#>   EH2259 | Kang18_8vs8
sce <- eh[["EH2259"]]
```

Before running muscat, we perform standard single-cell preprocessing
steps: removing undetected genes, filtering low-quality cells using
*[scater](https://bioconductor.org/packages/3.23/scater)*, and removing
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
#> ISG15                             4.714773               6.66734e-20
#> TNFRSF18                         -0.583992               6.47700e-02
#> TNFRSF4                          -0.220578               4.98028e-01
#> SDF4                             -0.524322               1.77759e-02
#>          stim-ctrl_NK cells_padj
#>                        <numeric>
#> NOC2L                         NA
#> HES4                          NA
#> ISG15                2.88696e-17
#> TNFRSF18             1.63055e-01
#> TNFRSF4              6.65061e-01
#> SDF4                 6.29608e-02

getDEA(dde,
       dea_name = "stim-ctrl_CD14+ Monocytes",
       format = "original") |> head()
#>              gene      cluster_id log2FoldChange    logCPM          F
#> HES4         HES4 CD14+ Monocytes      6.4921975  7.975074 305.569411
#> ISG15       ISG15 CD14+ Monocytes      7.0240366 14.751515 232.635530
#> SDF4         SDF4 CD14+ Monocytes     -0.7242092  5.446124  15.386517
#> UBE2J2     UBE2J2 CD14+ Monocytes     -0.8153555  5.351035  23.424243
#> CPSF3L     CPSF3L CD14+ Monocytes     -0.6688228  4.337862   7.582583
#> AURKAIP1 AURKAIP1 CD14+ Monocytes     -0.5741367  7.389181  34.562176
#>                pvalue         padj    p_adj.glb  contrast
#> HES4     2.883448e-14 1.398472e-12 1.686302e-12 stim-ctrl
#> ISG15    6.200271e-13 1.606550e-11 2.344791e-11 stim-ctrl
#> SDF4     7.676618e-04 1.818450e-03 3.340543e-03 stim-ctrl
#> UBE2J2   8.502381e-05 2.510891e-04 4.817526e-04 stim-ctrl
#> CPSF3L   1.187149e-02 2.083016e-02 3.426077e-02 stim-ctrl
#> AURKAIP1 7.378159e-06 2.880451e-05 5.506716e-05 stim-ctrl
```

General info on the DEA performed can be shown with
[`getDEAInfo()`](../reference/DeeDeeExperiment-methods.md)

``` r
# retrieving the DEA information
dea_name <- "stim-ctrl_B cells"
getDEAInfo(dde)[[dea_name]][["package"]]
#> [1] "muscat"
getDEAInfo(dde)[[dea_name]][["package_version"]]
#> [1] "1.25.0"
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
#>            stim-ctrl_B cells  345  271 0.05
#>    stim-ctrl_CD14+ Monocytes 1121 1289 0.05
#>        stim-ctrl_CD4 T cells  642  516 0.05
#>        stim-ctrl_CD8 T cells  178  115 0.05
#>    stim-ctrl_Dendritic cells  113   70 0.05
#>  stim-ctrl_FCGR3A+ Monocytes  498  413 0.05
#>     stim-ctrl_Megakaryocytes   29    2 0.05
#>           stim-ctrl_NK cells  254  201 0.05
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
#> R Under development (unstable) (2025-11-24 r89053)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.7.2
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
#>  [1] limma_3.67.0                muscat_1.25.0              
#>  [3] scater_1.39.0               ggplot2_4.0.1              
#>  [5] scuttle_1.21.0              ExperimentHub_3.1.0        
#>  [7] AnnotationHub_4.1.0         BiocFileCache_3.1.0        
#>  [9] dbplyr_2.5.1                DeeDeeExperiment_1.1.2     
#> [11] SingleCellExperiment_1.33.0 SummarizedExperiment_1.41.0
#> [13] Biobase_2.71.0              GenomicRanges_1.63.0       
#> [15] Seqinfo_1.1.0               IRanges_2.45.0             
#> [17] S4Vectors_0.49.0            BiocGenerics_0.57.0        
#> [19] generics_0.1.4              MatrixGenerics_1.23.0      
#> [21] matrixStats_1.5.0           BiocStyle_2.39.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.22         splines_4.6.0            bitops_1.0-9            
#>   [4] filelock_1.0.3           tibble_3.3.0             lpsymphony_1.39.0       
#>   [7] lifecycle_1.0.4          httr2_1.2.1              Rdpack_2.6.4            
#>  [10] edgeR_4.9.0              doParallel_1.0.17        globals_0.18.0          
#>  [13] lattice_0.22-7           MASS_7.3-65              backports_1.5.0         
#>  [16] magrittr_2.0.4           sass_0.4.10              rmarkdown_2.30          
#>  [19] jquerylib_0.1.4          yaml_2.3.11              sctransform_0.4.2       
#>  [22] DBI_1.2.3                minqa_1.2.8              RColorBrewer_1.1-3      
#>  [25] abind_1.4-8              EnvStats_3.1.0           glmmTMB_1.1.13          
#>  [28] purrr_1.2.0              rappdirs_0.3.3           sandwich_3.1-1          
#>  [31] circlize_0.4.16          pbkrtest_0.5.5           ggrepel_0.9.6           
#>  [34] irlba_2.3.5.1            listenv_0.10.0           RSpectra_0.16-2         
#>  [37] parallelly_1.45.1        pkgdown_2.2.0.9000       codetools_0.2-20        
#>  [40] DelayedArray_0.37.0      tidyselect_1.2.1         shape_1.4.6.1           
#>  [43] farver_2.1.2             lme4_1.1-38              ScaledMatrix_1.19.0     
#>  [46] viridis_0.6.5            jsonlite_2.0.0           GetoptLong_1.1.0        
#>  [49] BiocNeighbors_2.5.0      iterators_1.0.14         systemfonts_1.3.1       
#>  [52] foreach_1.5.2            tools_4.6.0              progress_1.2.3          
#>  [55] ragg_1.5.0               Rcpp_1.1.0               blme_1.0-6              
#>  [58] glue_1.8.0               gridExtra_2.3            SparseArray_1.11.8      
#>  [61] xfun_0.54                mgcv_1.9-4               DESeq2_1.51.6           
#>  [64] dplyr_1.1.4              withr_3.0.2              numDeriv_2016.8-1.1     
#>  [67] BiocManager_1.30.26      fastmap_1.2.0            boot_1.3-32             
#>  [70] caTools_1.18.3           digest_0.6.39            rsvd_1.0.5              
#>  [73] R6_2.6.1                 textshaping_1.0.4        colorspace_2.1-2        
#>  [76] gtools_3.9.5             RSQLite_2.4.5            RhpcBLASctl_0.23-42     
#>  [79] tidyr_1.3.1              variancePartition_1.41.0 corpcor_1.6.10          
#>  [82] data.table_1.17.8        prettyunits_1.2.0        httr_1.4.7              
#>  [85] htmlwidgets_1.6.4        S4Arrays_1.11.1          uwot_0.2.4              
#>  [88] pkgconfig_2.0.3          gtable_0.3.6             blob_1.2.4              
#>  [91] ComplexHeatmap_2.27.0    S7_0.2.1                 XVector_0.51.0          
#>  [94] remaCor_0.0.20           htmltools_0.5.9          bookdown_0.46           
#>  [97] TMB_1.9.18               clue_0.3-66              scales_1.4.0            
#> [100] png_0.1-8                fANCOVA_0.6-1            reformulas_0.4.2        
#> [103] knitr_1.50               reshape2_1.4.5           rjson_0.2.23            
#> [106] nlme_3.1-168             curl_7.0.0               nloptr_2.2.1            
#> [109] cachem_1.1.0             zoo_1.8-14               GlobalOptions_0.1.3     
#> [112] stringr_1.6.0            KernSmooth_2.23-26       BiocVersion_3.23.1      
#> [115] parallel_4.6.0           vipor_0.4.7              AnnotationDbi_1.73.0    
#> [118] desc_1.4.3               pillar_1.11.1            grid_4.6.0              
#> [121] vctrs_0.6.5              gplots_3.3.0             slam_0.1-55             
#> [124] BiocSingular_1.27.1      IHW_1.39.0               beachmat_2.27.0         
#> [127] cluster_2.1.8.1          beeswarm_0.4.0           evaluate_1.0.5          
#> [130] mvtnorm_1.3-3            cli_3.6.5                locfit_1.5-9.12         
#> [133] compiler_4.6.0           rlang_1.1.6              crayon_1.5.3            
#> [136] future.apply_1.20.0      fdrtool_1.2.18           plyr_1.8.9              
#> [139] fs_1.6.6                 ggbeeswarm_0.7.3         stringi_1.8.7           
#> [142] viridisLite_0.4.2        BiocParallel_1.45.0      lmerTest_3.1-3          
#> [145] Biostrings_2.79.2        aod_1.3.3                Matrix_1.7-4            
#> [148] hms_1.1.4                bit64_4.6.0-1            future_1.68.0           
#> [151] KEGGREST_1.51.1          statmod_1.5.1            rbibutils_2.4           
#> [154] broom_1.0.11             memoise_2.0.1            bslib_0.9.0             
#> [157] bit_4.6.0
```

## References

Kang, Hyun Min, Meena Subramaniam, Sasha Targ, Michelle Nguyen, Lenka
Maliskova, Elizabeth McCarthy, Eunice Wan, et al. 2017. “Multiplexed
Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation.”
*Nature Biotechnology* 36 (1): 89–94.
<https://doi.org/10.1038/nbt.4042>.
