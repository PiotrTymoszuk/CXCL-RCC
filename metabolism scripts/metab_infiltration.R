# Correlation of the metabolic gene GSVA scores with infiltration estimates
# by QuanTISeq

  insert_head()

# globals ----

  meta_infil <- list()

# globals -----

  insert_msg('Globals')

  ## metabolic variables

  meta_infil$metab_vars <-
    names(meta_gene$gene_symbols$tcga)

  meta_infil$metab_labs <- c('FAOx',
                             'OxPhos',
                             'TCA',
                             'TRP') %>%
    set_names(meta_infil$metab_vars)

  ## variable pairs

  meta_infil$pairs <- meta_infil$metab_vars %>%
    map(function(metab) globals$imm_estimates %>%
          map(~c(metab, .x))) %>%
    unlist(recursive = FALSE)

  ## parallel backend

  plan('multisession')

# Serial correlation analysis, Spearman ------

  insert_msg('Serial correlation analysis, Spearman')

  meta_infil$test_results <- meta_gene$score_tbl %>%
    future_map(function(cohort) meta_infil$pairs %>%
                 map_dfr(~correlate_variables(cohort,
                                              variables = .x,
                                              what = 'correlation',
                                              type = 'spearman',
                                              ci = TRUE,
                                              pub_styled = FALSE)) %>%
                 mutate(p_adjusted = p.adjust(p_value, 'BH'),
                        est_lab = paste0(signif(estimate, 2),
                                         ' [', signif(lower_ci, 2),
                                         ' - ', signif(upper_ci, 2), ']'),
                        significance = ifelse(p_adjusted >= 0.05,
                                              paste0('ns (',
                                                     signif(p_adjusted, 2),
                                                     ')'),
                                              paste('p =', signif(p_adjusted, 2)))),
               .options = furrr_options(seed = TRUE))

  meta_infil$test_results <- meta_infil$test_results %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 'bold', 'plain'),
        plot_lab = signif(estimate, 2),
        plot_cap = paste(est_lab, significance, sep = ', '),
        plot_cap = paste(plot_cap, n, sep = ', n = '))

# Bubble plots with the correlation coefficients ------

  insert_msg('Bubble plots with the correlation coefficients')

  meta_infil$bubble <-
    map2(meta_infil$test_results,
         globals$cohorts,
         ~ggplot(.x,
                 aes(x = variable1,
                     y = variable2,
                     size = abs(estimate),
                     fill = estimate)) +
           geom_point(shape = 21) +
           geom_text(aes(label = plot_lab,
                         fontface = significant),
                     size = 2.75,
                     hjust = 0.5,
                     vjust = -1.5) +
           scale_y_discrete(labels = c(globals$imm_labels,
                                       c(cd8_treg_ratio = 'CD8+/Treg',
                                         cd8_m2_ratio = 'CD8+/M2 TAM'))) +
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'white',
                                high = 'firebrick',
                                midpoint = 0,
                                limits = c(-1, 1),
                                name = expression(rho)) +
           scale_size_area(limits = c(0, 1),
                           max_size = 6,
                           name = expression('abs('*rho*')')) +
           scale_x_discrete(limits = meta_infil$metab_vars,
                            labels = meta_infil$metab_labs) +
           globals$common_theme +
           theme(axis.title = element_blank()) +
           labs(title = .y,
                subtitle = paste('Spearman correlation, n =',
                                 .x$n[[1]])))

# Plotting the correlations of TRegs and FAOX -----

  insert_msg('Plotting correlations for T Regs and FAOX')

  ## plot captions

  meta_infil$treg_faox$captions <- meta_infil$test_results %>%
    map(filter,
        variable1 == 'faox',
        variable2 == 'T cell regulatory (Tregs)') %>%
    map(~.x$plot_cap)

  ## plots

  meta_infil$treg_faox$plots <-
    list(data = meta_gene$score_tbl,
         plot_title = paste('FAOx and Tregs',
                            globals$cohorts,
                            sep = ', '),
         plot_subtitle = meta_infil$treg_faox$captions) %>%
    pmap(plot_correlation,
         variables = c('faox', 'T cell regulatory (Tregs)'),
         show_trend = TRUE,
         point_color = 'orangered2',
         x_lab = 'FAOx geneSign, GSVA score',
         y_lab = 'Tregs, QuanTIseq estimate',
         point_hjitter = 0,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(plot.tag = element_blank()) +
          scale_y_continuous(trans = 'log2',
                             labels = function(x) signif(x, 2)))

# Plotting the correlations between Tregs and OxPhos ------

  insert_msg('Plotting correlations between the T regs and OxPhos')

  ## plot captions

  meta_infil$treg_oxphos$captions <- meta_infil$test_results %>%
    map(filter,
        variable1 == 'oxphos',
        variable2 == 'T cell regulatory (Tregs)') %>%
    map(~.x$plot_cap)

  ## plots

  meta_infil$treg_oxphos$plots <-
    list(data = meta_gene$score_tbl,
         plot_title = paste('Oxphos and Tregs',
                            globals$cohorts,
                            sep = ', '),
         plot_subtitle = meta_infil$treg_faox$captions) %>%
    pmap(plot_correlation,
         variables = c('oxphos', 'T cell regulatory (Tregs)'),
         show_trend = TRUE,
         point_color = 'gray60',
         x_lab = 'OxPhos geneSign, GSVA score',
         y_lab = 'Tregs, QuanTIseq estimate',
         point_hjitter = 0,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(plot.tag = element_blank()) +
          scale_y_continuous(trans = 'log2',
                             labels = function(x) signif(x, 2)))

# Plotting the correlations between Tregs and TCA ------

  insert_msg('Plotting correlations between the T regs and TCA')

  ## plot captions

  meta_infil$treg_tca$captions <- meta_infil$test_results %>%
    map(filter,
        variable1 == 'tca',
        variable2 == 'T cell regulatory (Tregs)') %>%
    map(~.x$plot_cap)

  ## plots

  meta_infil$treg_tca$plots <-
    list(data = meta_gene$score_tbl,
         plot_title = paste('TCA and Tregs',
                            globals$cohorts,
                            sep = ', '),
         plot_subtitle = meta_infil$treg_faox$captions) %>%
    pmap(plot_correlation,
         variables = c('tca', 'T cell regulatory (Tregs)'),
         show_trend = TRUE,
         point_color = 'darkolivegreen3',
         x_lab = 'TCA geneSign, GSVA score',
         y_lab = 'Tregs, QuanTIseq estimate',
         point_hjitter = 0,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(plot.tag = element_blank()) +
          scale_y_continuous(trans = 'log2',
                             labels = function(x) signif(x, 2)))

# END -----

  plan('sequential')

  insert_tail()
