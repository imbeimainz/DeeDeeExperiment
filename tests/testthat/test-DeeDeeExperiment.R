test_that("creating", {
  ## create dde objects --------------------------------------------------------
  expect_error(DeeDeeExperiment())

  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  print(dde)

  expect_s4_class(dde, "DeeDeeExperiment")

  dde_only_de <- DeeDeeExperiment(
    de_results = de_named_list
  )
  expect_s4_class(dde_only_de, "DeeDeeExperiment")

  expect_type(dea_names(dde_only_de), "character")

  expect_equal(
    dea_names(dde_only_de),
    c("ifng_vs_naive", "ifngsalmo_vs_naive", "salmonella_vs_naive", "salmo_both")
  )

  dde_nodd <- DeeDeeExperiment(
    se = se_macrophage_noassays,
  )
  expect_s4_class(dde_nodd, "DeeDeeExperiment")


  expect_s4_class(
    dea(dde, "ifng_vs_naive"), "DataFrame"
  )

  expect_error(
    dea(dde, "wrong_name", verbose = TRUE)
  )

  dde_gone_wrong <- dde
  rowData(dde_gone_wrong)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(
    dea(dde_gone_wrong, "ifng_vs_naive", verbose = TRUE)
  )

  expect_error(
    DeeDeeExperiment(
      rowData(se_macrophage_noassays),
      de_results = de_named_list
    )
  )

  expect_error(
    DeeDeeExperiment(
      assay(se_macrophage),
      de_results = de_named_list
    )
  )

  salmo_both <- de_named_list$salmo_both

  dde_one <- DeeDeeExperiment(se = se_macrophage_noassays,
                              de_results = salmo_both)

  expect_false(is.list(salmo_both))

  expect_s4_class(salmo_both, "DESeqResults")

  expect_type(dea_info(dde_one), "list")

  expect_length(dea_info(dde_one), 1)

  expect_true("salmo_both" == names(dea_info(dde_one)))

  de_results_mismatch <- list(
    contrast = de_named_list$ifng_vs_naive
  )
  rownames(de_results_mismatch$contrast) <- paste0(
    "gene", 101:(100 + nrow(de_results_mismatch$contrast)))

  expect_warning(
    DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_results_mismatch
  ))

  broken_de_res <- salmo_both
  rownames(broken_de_res) <- NULL

  expect_error(DeeDeeExperiment(de_results = list(salmo_both = salmo_both,
                                                  broken_salmo_both = broken_de_res)))

  expect_true(is(.check_de_results(list(salmo_both = salmo_both)), "list"))
  expect_error(.check_de_results(
    list(salmo_both = salmo_both,
         simple_list = list())
    ),
    "All elements in the list must be of type"
  )

  dea1 <- de_limma

  de_res_list <- list(de_deseq = salmo_both,
                      dge_lrt = dea1)

  dde_list <- DeeDeeExperiment(de_results = de_res_list)

  expect_warning(get_dea_list(dde_list, verbose =TRUE), regexp = NULL)

  expect_warning(dea(dde_list, dea_name = "dge_lrt", verbose =TRUE) , regexp = NULL)

  expect_warning(DeeDeeExperiment(se = se_macrophage_noassays,
                                   de_results = dea1))

  dea2 <- dge_exact_IFNg_both
  expect_warning(DeeDeeExperiment(se = se_macrophage_noassays,
                                   de_results = dea2))

  fea1 <- topGO_results

  expect_warning(DeeDeeExperiment(enrich_results = fea1))

  contrast1 <- topGO_results$ifng_vs_naive

  expect_warning(DeeDeeExperiment(de_results = de_named_list,
                                  enrich_results = contrast1))



  expect_message({DeeDeeExperiment(
    de_results = list(ifng_vs_naive = de_named_list$ifng_vs_naive,
                      salmonella_vs_naive = de_named_list$salmonella_vs_naive),
    enrich_results = list(
      topGO_ifng_vs_naive = topGO_results$ifng_vs_naive,
      salmonella_vs_naive = topGO_results$salmonella_vs_naive))}
    )


  broken_limma <- de_limma
  broken_limma$coefficients <- broken_limma$coefficients[, "Salm_both", drop = FALSE]
  broken_limma$t           <- broken_limma$t[, "Salm_both", drop = FALSE]
  broken_limma$p.value     <- broken_limma$p.value[, "Salm_both", drop = FALSE]
  broken_limma$lods        <- broken_limma$lods[, "Salm_both", drop = FALSE]

  expect_error(DeeDeeExperiment(se = se_macrophage_noassays,
                               de_results = broken_limma))

  expect_error(DeeDeeExperiment(enrich_results = "1st enrich res"))

  expect_error(DeeDeeExperiment(de_results = list(de_limma)))

  dde5 <- DeeDeeExperiment(se = se_macrophage_noassays,
                           enrich_results = list(enrichr_salmo_vs_naive = enrichr_res$Reactome_2016))

  expect_s3_class(fea(dde5,"enrichr_salmo_vs_naive",verbose =TRUE), "data.frame")

  expect_length(fea_info(dde5), 1)


  expect_error(DeeDeeExperiment(se = se_macrophage_noassays,
                     de_results = de_named_list,
                     enrich_results = list(clusterPro_res$ifng_vs_naive)))

  expect_error(DeeDeeExperiment(se = se_macrophage_noassays,
                                de_results = de_named_list,
                                enrich_results = gost_res))

  new_dde <- DeeDeeExperiment(se = se_macrophage_noassays,
                              de_results = de_named_list,
                              enrich_results = list(
                                clusterPro_res = clusterPro_res$salmonella_vs_naive,
                                gPro_res = gost_res$result,
                                fgsea = fgseaRes,
                                gsea = gsea_res))

  expect_s4_class(new_dde, "DeeDeeExperiment")

  expect_equal(fea_info(new_dde)$fgsea$fe_tool, "fgsea")

  expect_equal(fea_info(new_dde)$gPro_res$fe_tool, "gProfiler")

  expect_equal(fea_info(new_dde)$clusterPro_res$fe_tool, "clusterProfiler")

  failing_fgsea <- fgseaRes
  expect_error({
    .DeeDeefy_fgseaResult(matrix())
  }, regexp = "should be a data.frame")

  expect_error({
    fail1 <- failing_fgsea[, -1]
    .DeeDeefy_fgseaResult(fail1)
  }, "I could not find some of the usual column names")

  expect_error({
    fail2 <- failing_fgsea
    fail2$leadingEdge <- as.character(fail2$leadingEdge)
    .DeeDeefy_fgseaResult(fail2)
  }, "Expecting 'leadingEdge' column to be a list")

  expect_error({
    fail_enrichr <- enrichr_res$Reactome_2016
    fail_enrichr$Genes <- NULL
    .DeeDeefy_enrichr(fail_enrichr)
  }, "I could not find some of the usual column names")


  expect_message(DeeDeeExperiment(se = se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = topGO_results), )

  expect_message(DeeDeeExperiment(se = se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = list(
                                    topGO_Salm_naive = topGO_results$salmonella_vs_naive,
                                    topGO_IFNg_naive = topGO_results$ifng_vs_naive)))



})
