# GSVA scores for genes associated with energy metabolism and TRP degradation

  insert_head()

# container -----

  meta_gene <- list()

# globals: reactions of interest ------

  insert_msg('Globals: reactions of interest')

  meta_gene$reacts[c('faox', 'oxphos', 'tca', 'trp')] <-
    meta_hg[c('faox', 'oxphos', 'tca', 'trp')] %>%
    map(~.x$reacts)

# searching for genes ------

  insert_msg('Seraching for genes')

  ## entrez IDs

  meta_gene$entrez_id <- meta_gene$reacts %>%
    map(react_to_gene, geneSBML = meta_pool$model) %>%
    map(~.x$entrez_id) %>%
    map(unlist) %>%
    unique %>%
    set_names(names(meta_gene$reacts))

  ## gene symbols

  meta_gene$gene_symbols <- dge$dict %>%
    map(function(cohort) meta_gene$entrez_id %>%
          map(~filter(cohort, entrez_id %in% .x)) %>%
          map(~.$gene_symbol))

# analysis tables with the gene expression values only -------

  insert_msg('Analysis tables')

  meta_gene$analysis_tbl <-
    map2(dge$analysis_tbl,
         map(dge$dict, ~.x$gene_symbol),
         ~.x[.y])

# Calculating GSVA scores for each gene signature ------

  insert_msg('Calculating GSVA scores')

  plan('multisession')

  meta_gene$score_tbl <-
    future_map2(meta_gene$gene_symbols,
                meta_gene$analysis_tbl,
                ~gseaTools::calculate(.x, data = .y))

  plan('sequential')

# Appending with the CXCL9 strata information and infiltration estimates -----

  insert_msg('Appanding with the CXCL9 strata and immune infiltration data')

  meta_gene$infil_tbl <- dge$analysis_tbl %>%
    map(select, cxcl9_group, CXCL9, all_of(globals$imm_estimates)) %>%
    map(mutate,
        cd8_treg_ratio = `T cell CD8+`/`T cell regulatory (Tregs)`,
        cd8_treg_ratio = ifelse(is.na(cd8_treg_ratio) | is.infinite(cd8_treg_ratio),
                                0, cd8_treg_ratio),
        cd8_m2_ratio = `T cell CD8+`/`Macrophage M2`,
        cd8_m2_ratio = ifelse(is.na(cd8_m2_ratio) | is.infinite(cd8_m2_ratio),
                              0, cd8_m2_ratio))

  meta_gene$score_tbl <-
    map2(meta_gene$score_tbl,
         meta_gene$infil_tbl,
         cbind) %>%
    map(as_tibble)

# END ------

  meta_gene$analysis_tbl <- NULL
  meta_gene$infil_tbl <- NULL

  meta_gene <- compact(meta_gene)
