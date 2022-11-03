# Checks expression of signature genes connected to activation or suppression
# of T cells in the tumor

  insert_head()

# container list ------

  dge_sig <- list()

# signature database -----

  insert_msg('Signature databases')

  dge_sig$db <- load_dbsig('./input data/signatures/msigdb.v7.5.1.symbols.gmt')

  dge_sig$sigs[c('cytotoxicity',
                 'exhaustion',
                 'treg',
                 'cd8')] <-
    c('CYTOTOXIC|CYTOTOXICITY',
      'EXHAUST*',
      'TREG|(REGULATORY_T_CELL*)|(REGULATORY_T_LYMPHOCYTE*)|(REGULATORY_TCELL*)',
      '(CD8_TCELL*)|(CD8_T_CELL*)|(CD8_T_LYMPHOCYTE*)') %>%
    map(~filter(dge_sig$db, stri_detect(sign_name, regex = .x)))

  dge_sig$sigs[c('mitosis',
                 'cell_cycle')] <- c('MITOSIS|MITOTIC',
                                     'CELLCYCLE|(CELL_CYCLE)') %>%
    map(~filter(dge_sig$db, stri_detect(sign_name, regex = .x)))

  dge_sig$db <- NULL

  dge_sig <- compact(dge_sig)

# calculation of the signatures ------

  insert_msg('GSVA estimates of the signatures')

  plan('multisession')

  dge_sig[c('gsva_cytotoxicity',
            'gsva_exhaustion',
            'gsva_treg',
            'gsva_cd8',
            'gsva_mitosis',
            'gsva_cellcycle')] <- dge_sig$sigs %>%
    map(function(sig_db) list(x = dge$analysis_tbl,
                              y = dge$dict) %>%
          future_pmap(function(x, y) calculate(sig_db,
                                               data = x[y$gene_symbol]) %>%
                        mutate(patient_id = x$patient_id,
                               cxcl9_group = x$cxcl9_group),
                      .options = furrr_options(seed = TRUE)))

  plan('sequential')

# Testing for the signature differences between the CXCL9 expression strata ------

  insert_msg('Testing for the signature differences between the strata')

  plan('multisession')

  dge_sig$test_results <- dge_sig[c('gsva_cytotoxicity',
                                    'gsva_exhaustion',
                                    'gsva_treg',
                                    'gsva_cd8',
                                    'gsva_mitosis',
                                    'gsva_cellcycle')] %>%
    map(function(sig_est) sig_est %>%
          future_map(~test_two_groups(.x,
                                      split_fct = 'cxcl9_group',
                                      variables = names(.x)[!names(.x) %in% c('patient_id', 'cxcl9_group')],
                                      type = 't',
                                      parallel = FALSE),
                     .options = furrr_options(seed = TRUE))) %>%
    unlist(recursive = FALSE)

  plan('sequential')

# formatting the results -----

  insert_msg('Formatting the results')

  dge_sig$test_results <- dge_sig$test_results %>%
    map(mutate,
        regulation = ifelse(significant == 'no',
                            'ns',
                            ifelse(estimate > 0, 'upregulated', 'downregulated')),
        regulation = factor(regulation, c('upregulated', 'downregulated', 'ns')))

# plot titles -----

  insert_msg('Generating plot titles')

  dge_sig$plot_titles <- names(dge_sig$test_results) %>%
    tibble(study = stri_split_fixed(., pattern = '.', simplify = TRUE)[, 2],
           family = stri_split_fixed(., pattern = '.', simplify = TRUE)[, 1]) %>%
    mutate(study = globals$cohorts[study],
           family = stri_replace(family, fixed = 'gsva_', replacement = ''),
           family = car::recode(family,
                                "'cytotoxicity' = 'Cytotoxicity';
                                'exhaustion' = 'Exhaustion';
                                'treg' = 'TReg';
                                'cd8' = 'CD8+ T cells';
                                'mitosis' = 'Mitosis';
                                'cellcycle' = 'Cell cycle'"),
           plot_title = paste(family, study, sep = ', ')) %>%
    .$plot_title

# Numbers of significantly regulated signatures ------

  insert_msg('Numbers of significantly regulated signatures')

  dge_sig$sign_n <- dge_sig$test_results %>%
    map(filter, regulation != 'ns') %>%
    map(mutate,
        regulation = factor(regulation, c('upregulated', 'downregulated'))) %>%
    map(count, regulation, .drop = FALSE) %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# plotting the estimates -------

  insert_msg('Plotting the estimates')

  dge_sig$est_plots <-
    list(data = dge_sig$test_results,
         plot_title = dge_sig$plot_titles,
         plot_subtitle = dge_sig$sign_n,
         n_sign = map_dbl(dge_sig$test_results, nrow)) %>%
    pmap(function(data, plot_title, n_sign, plot_subtitle) ggplot(data,
                                                                  aes(x = estimate,
                                                                      y = reorder(response, estimate),
                                                                      fill = regulation)) +
           geom_bar(stat = 'identity') +
           geom_vline(xintercept = 0,
                      linetype = 'solid',
                      size = 0.75) +
           scale_fill_manual(values = c(upregulated = 'firebrick',
                                        downregulated = 'steelblue',
                                        ns = 'gray60'),
                             name = 'CXCL9 high vs low') +
           scale_x_continuous(limits = c(-0.8, 0.8)) +
           globals$common_theme +
           theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank(),
                 panel.grid.major = element_blank()) +
           labs(title = plot_title,
                subtitle = plot_subtitle,
                x = 'Regulation, CXCL9 high vs low tumors',
                y = paste('Signature, total n =', n_sign)))

# END ----

  insert_tail()
