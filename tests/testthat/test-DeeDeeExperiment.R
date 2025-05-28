test_that("creating", {
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


}

)


test_that("adding and removing", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  new_del <- list(
    ifng2 = de_named_list$ifng_vs_naive,
    ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
  )
  # add a new (set of) DE result(s)
  dde_new <- add_dea(dde, new_del)
  expect_s4_class(dde_new, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde)), 4)
  expect_equal(length(dea_info(dde_new)), 6)

  expect_error({add_dea(x = dde, dea = list(de_named_list$ifng_vs_naive))})

  expect_error({dde_add <- add_dea(x = "sthg else", dea = new_del)})

  expect_error({dde_add <- add_dea(x = dde, dea = list(
    ifng2 = de_named_list$ifng_vs_naive,
    ifng2 = de_named_list$ifngsalmo_vs_naive
  ))})

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_s4_class(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_removed)), 3)

  expect_warning(dde_removed <- remove_dea(dde, "lol"))

  expect_error(remove_dea(dde))



  dde_edgeR <- add_dea(dde, dea = list(DGEExact_IFNg_both = dge_exact_IFNg_both))
  expect_s4_class(dde_edgeR, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_edgeR)), 5)

  dde_limma <- add_dea(dde, dea = list(de_limma = de_limma))
  expect_s4_class(dde_limma, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_limma)), 5)

  expect_no_error(add_dea(dde_limma,
                          dea = list(de_limma = de_limma),
                          force = TRUE))

  genes <- rownames(de_limma)[1:20]

  de_custom <- data.frame(p_val = rep(0.5,20),
                          adj_pvalue = rep(0.5,20),
                          gene = genes)
  expect_error(add_dea(dde, dea = list(de_custom = de_custom)))


  topGO_Salm_naive <- topGO_results$salmonella_vs_naive
  topGO_IFNg_naive <- topGO_results$ifng_vs_naive

  dde2 <- add_fea(dde, fea = list(topGO_Salm_naive = topGO_Salm_naive,
                            topGO_IFNg_naive = topGO_IFNg_naive),
                  fea_tool = "topGO")

  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 2)

  dde2 <- remove_fea(dde2, "topGO_IFNg_naive")
  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 1)


  expect_error(remove_fea(dde2,"IFNgVSnaive"))

  expect_s3_class(fea(dde2, "topGO_Salm_naive"), "data.frame")

  dde3 <- DeeDeeExperiment(se = se_macrophage_noassays)
  expect_warning(add_fea(dde3, fea = list(topGO_Salm_naive = topGO_Salm_naive,
                                      topGO_IFNg_naive = topGO_IFNg_naive), verbose = TRUE))

  expect_error(add_fea(dde3, fea = list(topGO_Salm_naive = topGO_Salm_naive,
                                              topGO_IFNg_naive)))

  dde3 <- add_fea(dde3, fea = list(topGO_Salm_naive = topGO_Salm_naive,
                                       topGO_IFNg_naive = topGO_IFNg_naive), verbose = TRUE)
  expect_error({
    add_fea(dde3, fea = list(topGO_Salm_naive = topGO_Salm_naive),
            force = FALSE)
  })

  expect_error(add_fea(dde3, fea = list(FE1 = topGO_Salm_naive,
                                        FE1 = topGO_IFNg_naive)))

  expect_message(DeeDeeExperiment(se = se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = topGO_results), )

  expect_message(DeeDeeExperiment(se = se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = list(topGO_Salm_naive = topGO_Salm_naive,
                                                        topGO_IFNg_naive = topGO_IFNg_naive)))

  fea_name <- character(0)
  expect_error({
    remove_fea(dde3, fea_name)})

  dde3 <- add_fea(dde3, fea = list(gPro_salmonella_vs_naive = gost_res$result))

  expect_equal(fea_info(dde3)$gPro_salmonella_vs_naive$fe_tool, "gProfiler")

  expect_message(dde3 <- add_fea(dde3,
                                 fea =
                                   list(salmonella_vs_naive = enrichr_res$KEGG_2019_Human), verbose = TRUE))


  expect_error(fea_rename(dde3, old_name = "contrast1",
                          new_name = "topGO_IFNg_naive"))

  expect_error(fea_rename(dde, old_name = "topGO_Salm_naive",
                          new_name = c("salmonella_vs_naive","new_column")))

  expect_warning(dea(dde, verbose = TRUE))

  expect_error(dea(dde3, dea_name = "contrast1"))

  expect_error(dea(dde, dea_name = c("salmonella_vs_naive","salmo_both")))

  expect_error(link_dea_and_fea(dde,
                                 dea_name = "contrast1",
                                 fea_name = "topGO_IFNg_naive"))

  expect_error(link_dea_and_fea(dde,
                                 dea_name = "salmonella_vs_naive",
                                 fea_name = "topGO"))

  expect_error(get_fea_list(dde3,dea_name = c("salmonella_vs_naive", "salmo_both")))
  expect_error(get_fea_list(dde))

  expect_error(add_dea(dde, "ifng_vs_naive"))

  expect_error(dea(dde, format = "simple"))
  expect_error(fea(dde, format = "simple"))

  expect_error(dea(dde3))
  expect_error(fea(dde))

  dde_overlap <- DeeDeeExperiment(se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = topGO_results)

  expect_error(add_dea(dde_overlap,
                       dea = list(ifng_vs_naive = de_named_list$ifng_vs_naive)))

  dde_de_empty <- DeeDeeExperiment(se_macrophage_noassays)

  expect_error(dea_rename(dde_de_empty,
                          old_name = "ifng_vs_naive", new_name = "IFNgvsNaive"))

  expect_error(dea(dde_overlap,
      dea_name = "ifng_vs_naive",
      extra_rd = "ifng_vs_naive_pvalue"))

  extract_dea <- dea(dde_overlap, dea_name = "ifng_vs_naive", format = "original")

  expect_s4_class(extract_dea, "DESeqResults")

  expect_error(get_dea_list(dde_overlap, format = "simple"))

  expect_error(fea_rename(dde_overlap,
                          old_name = c("ifng_vs_naive", "ifngsalmo_vs_naive"),
                          new_name = "naive"))

  expect_error(fea_rename(dde_overlap,
                          old_name = c("ifng_vs_naive", "ifngsalmo_vs_naive"),
                          new_name = c("naive", "naive")))

  expect_error(fea_rename(dde_overlap,
                          old_name = c("ifng_vs_naive"),
                          new_name = c("ifngsalmo_vs_naive")))


  expect_error(add_fea(dde_de_empty, fea = list(topGO_results$salmo_both)))

  expect_warning(fea(dde_overlap, verbose = TRUE))

  expect_error(fea(dde_overlap, fea_name = c("salmonella_vs_naive", "ifng_vs_naive")))

  expect_error(fea(dde_overlap, fea_name = "sthg else"))

  expect_s3_class(fea(dde_overlap, fea_name = "salmonella_vs_naive",
                      format = "original"), "data.frame")

  expect_error(link_dea_and_fea(dde_de_empty,
                                 "not there",
                                 "not there"))

  expect_error(link_dea_and_fea(dde_overlap,
                                 "new_name",
                                 "salmonella_vs_naive"))

  expect_error(link_dea_and_fea(dde_overlap,
                                 "ifng_vs_naive",
                                 "salmonella_vs_naive",
                                 force = FALSE))

  dde_overlap_add <- add_fea(dde_overlap, fea = list(INFg_vs_Naive = topGO_results$ifng_vs_naive))

  expect_warning(link_dea_and_fea(dde_overlap_add,
                                "ifng_vs_naive",
                                "salmonella_vs_naive",
                                force = TRUE))

  expect_error(link_dea_and_fea(dde_overlap,
                                dea_name = 2,
                                fea_name = "salmonella_vs_naive"))

  expect_error(link_dea_and_fea(dde_overlap,
                                dea_name = "salmonella_vs_naive",
                                fea_name = 2))

  expect_error(add_scenario_info(dde_overlap,
                    dea_name = "i dont exist"))

  expect_error(add_scenario_info(dde_overlap,
                                 dea_name = 2))

  dde_with_info <- add_scenario_info(dde_overlap,
                                     dea_name = "ifng_vs_naive",
                                     info = "here goes some txt about the contrast ifng_vs_naive")

  expect_s4_class(dde_with_info, "DeeDeeExperiment")

  expect_error(add_scenario_info(dde_with_info,
                                 dea_name = "ifng_vs_naive",
                                 info = "sthg else about ifng_vs_naive"))

  expect_error(add_scenario_info(dde_overlap,
                                 dea_name = c("i dont exist", "ifng_vs_naive")))


  expect_error(get_fea_list(dde_overlap,
               dea_name = "ifng_vs_naive",
               format = "simple"))

  expect_warning(get_fea_list(dde_overlap_add, dea_name = "INFg_vs_Naive"))



  new_remove_dea <- remove_dea(dde_overlap_add,
                               dea_name = "ifng_vs_naive",
                               remove_linked_fea = TRUE)

  expect_equal(length(fea_info(new_remove_dea)), 4)

  expect_error(remove_dea(dde, dea_name = NULL))

  expect_error(add_fea(dde_overlap_add,
          fea = list(gPro_res = gost_res)))


  dde_overlap_add <- add_fea(dde_overlap_add,
                             fea = list(gPro_res = gost_res$result))

  expect_equal(length(fea_info(dde_overlap_add)), 6)


  expect_error(dea(dde_overlap_add, extra_rd = NA))

  expect_warning(dea(dde_overlap_add, extra_rd = c("guiga","other"),
                     verbose = TRUE))

  expect_error(add_scenario_info(dde_overlap_add,
                                 dea_name = "salmo_both",
                                 info = NA ))

  expect_error(add_scenario_info(dde_overlap_add,
                                 dea_name = "salmo_both",
                                 info = data.frame(Info = "here is some context") ))

  expect_error(add_fea(dde_overlap_add,
                       fea = topGO_results$ifng_vs_naive,
                       de_name  = NA)
               )

  expect_error(add_fea(dde_overlap_add,
                       fea = topGO_results$ifng_vs_naive,
                       fea_tool = NA)
  )

  expect_no_error(add_fea(dde_with_info,
                       fea = topGO_results$ifng_vs_naive,
                       de_name = "ifng_vs_naive",
                       verbose = TRUE))

  expect_warning(add_fea(dde_with_info,
                         fea = topGO_results$ifng_vs_naive,
                         de_name = "IFNg_vs_naive"))

  expect_message(add_fea(dde_with_info,
                         fea = list(topGO_ifng_vs_naive = topGO_results$ifng_vs_naive),
                         verbose = TRUE))

  expect_message(add_fea(dde_with_info,
                         fea = list(ifng_vs_naive = topGO_results$ifng_vs_naive),
                         verbose = TRUE,
                         force = TRUE))

  expect_error(add_fea(dde_with_info,
                       fea = topGO_results$ifng_vs_naive,
                       fea_tool = c("topGO","topGO")
                       ))

  expect_message(add_fea(dde_with_info,
          fea = list(FEA1= clusterPro_res$ifng_vs_naive,
                    FEA2= gost_res$result,
                    FEA3=gsea_res,
                    FEA4=fgseaRes),
          fea_tool = c("clusterProfiler",
                       "gProfiler",
                       "gsea",
                       "fgsea"
                       )))



  expect_warning(fea(dde_overlap_add, verbos = TRUE, format = "original"))


  expect_length(get_fea_list(dde_overlap_add, format = "minimal", verbose = TRUE), 6)

  expect_length(get_fea_list(dde_overlap_add, dea_name = "ifng_vs_naive", format = "original", verbose = TRUE), 1)


})


test_that("renaming", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )


  dde_rename <- dea_rename(dde,old_name = "salmonella_vs_naive" ,
                           new_name = "SalmvsNaive")

  expect_s4_class(dde_rename, "DeeDeeExperiment")

  expect_equal(names(dea_info(dde_rename)), c("ifng_vs_naive", "ifngsalmo_vs_naive",
                                              "SalmvsNaive","salmo_both"))

  expect_error(dea_rename(dde, old_name = "salmonella_vs_naive",
                          new_name = "ifng_vs_naive"))

  expect_error(dea_rename(dde, old_name = NULL,
                          new_name = "ifng_vs_naive"))

  expect_error(dea_rename(dde, old_name = "salmonella_vs_naive",
                          new_name = NULL))

  expect_error(dea_rename(dde, old_name = "contrast1",
                          new_name = "ifng_vs_naive"))

  expect_error(dea_rename(dde, old_name = 1,
                          new_name = "1"))

  dde3 <- add_fea(dde_rename,fea = topGO_results)

  expect_error(dea_rename(dde3, old_name = 1,
                          new_name = "1"))

  expect_error(fea_rename(dde3, old_name = NULL,
                          new_name = "ifng_vs_naive"))

  expect_error(fea_rename(dde3, old_name = "salmonella_vs_naive",
                          new_name = NULL))

  expect_error(dea_rename(dde, old_name = "salmonella_vs_naive",
                          new_name = c("ifng_vs_naive","new_column")))

  expect_error(dea_rename(dde, old_name = c("salmonella_vs_naive","salmo_both"),
                          new_name = c("salmonel_vs_naive", "salmonel_vs_naive")))


  dde3_rename <- dea_rename(dde3,old_name = "ifngsalmo_vs_naive" ,
                           new_name = "IFNg_SalmvsNaive")

  expect_s4_class(dde3_rename, "DeeDeeExperiment")

  expect_equal(names(dea_info(dde3_rename)), c("ifng_vs_naive", "IFNg_SalmvsNaive",
                                              "SalmvsNaive","salmo_both"))

  expect_equal(fea_info(dde3_rename)[["ifngsalmo_vs_naive"]][["de_name"]],
               "IFNg_SalmvsNaive")



  fea_rename <- fea_rename(dde3,old_name = "salmonella_vs_naive" ,
                           new_name = "topGO_SalmonellavsNaive")
  expect_s4_class(fea_rename, "DeeDeeExperiment")

  expect_equal(fea_names(fea_rename), c("ifng_vs_naive",
                                        "ifngsalmo_vs_naive",
                                        "topGO_SalmonellavsNaive",
                                        "salmo_both"))


})


test_that("validity and so", {
  dde2 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  expect_true(validObject(dde2))

  rowData(dde2)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(validObject(dde2))

  dde3 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  dde3@dea <- list("foo", "bar")
  expect_error(validObject(dde3))

  dea_not_list <- dde3

  expect_error({
    dea_info(dea_not_list) <- data.frame()
  })

  expect_error({
    fea_info(dea_not_list) <- data.frame()
  })

  dde4 <- DeeDeeExperiment(
    se_macrophage_noassays,
    enrich_results = topGO_results
  )
  dde4@fea <- list("foo", "bar")
  expect_error(validObject(dde4))

  expect_error({
    add_fea(dde4, fea = "meow", fea_tool = "fujitsu")
  })


  expect_error({
    fea_info(dde4) <- list(list(de_name = "c1",
                           fe_name = "c1",
                           shaken_results = NULL,
                           original_object = "meow",
                           fe_tool = "topGO"))
  }
  )

  expect_error({
    fea_info(dde4) <- list(list(de_name = "c1",
                                fe_name = "c1",
                                shaken_results = NULL,
                                original_object = NULL,
                                fe_tool = "chico"))
  }
  )


})


test_that("misc", {

  dde_no_fea <- DeeDeeExperiment(se = se_macrophage_noassays,
                                 de_results = de_named_list)
  expect_no_error(summary(dde_no_fea))

  dde_no_dea <- DeeDeeExperiment(se = se_macrophage_noassays,
                                 enrich_results =  topGO_results)
  expect_no_error(summary(dde_no_dea))

  dde_empty <- DeeDeeExperiment(se = se_macrophage_noassays)
  expect_no_error(summary(dde_empty))



})

