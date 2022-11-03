# Inter-correlation of the CXCL9, 10 and 11 expression in the tumor tissue.
# Pearson correlation.

  insert_head()

# container list ----

  correl <- list()

# globals ----

  insert_msg('Globals setup')

  correl$analysis_tbl <- list(tcga = tcga$expression %>%
                                filter(tissue_type == 'Tumor'),
                              emtab1980 = emtab1980$expression,
                              gse73731 = gse73731$expression,
                              gse167093 = gse167093$expression %>%
                                filter(tissue_type == 'Tumor'),
                              reca = reca$expression,
                              cm10 = cm10$expression,
                              cm25ev = cm25ev$expression,
                              cm25ni = cm25ni$expression) %>%
    map(~.x[globals$cxcl_genes])

  correl$pairs <- combn(globals$cxcl_genes, m = 2, simplify = FALSE)

# testing -----

  insert_msg('Correlation')

  correl$test_results <- set_names(names(correl$analysis_tbl),
                                   names(correl$analysis_tbl)) %>%
    map(function(cohort) correl$pairs %>%
          map_dfr(~correlate_variables(correl$analysis_tbl[[cohort]],
                                       variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       pub_styled = TRUE,
                                       simplify_p = FALSE)) %>%
          mutate(correl_id = paste(variable1, variable2, sep = '_'),
                 p_adjusted = p.adjust(p_value, 'BH'),
                 significance = paste('p =', signif(p_adjusted, 2)),
                 plot_cap = paste(eff_size, significance, sep = ', '),
                 plot_cap = paste0(plot_cap, ', n = ', n),
                 plot_cap = stri_replace(plot_cap, fixed = 'rho', replacement = '\u03C1')))

# plotting -----

  insert_msg('Plotting')

  correl$plots <- list(cohort = set_names(names(correl$analysis_tbl),
                                          names(correl$analysis_tbl)),
                       suffix = globals$cohorts) %>%
    pmap(function(cohort, suffix) list(variables = correl$pairs,
                                      plot_title = map_chr(correl$pairs,
                                                           paste,
                                                           collapse = ': ') %>%
                                        paste(suffix, sep = ', '),
                                      plot_subtitle = correl$test_results[[cohort]]$plot_cap,
                                      x_lab = map(correl$pairs,
                                                  ~paste0(.x[1], ', log2 expression')),
                                      y_lab = map(correl$pairs,
                                                  ~paste0(.x[2], ', log2 expression'))) %>%
          pmap(plot_correlation,
               correl$analysis_tbl[[cohort]],
               type = 'correlation',
               cust_theme = globals$common_theme,
               point_color = 'indianred3') %>%
          map(~.x +
                theme(plot.tag = element_blank(),
                      plot.title = element_text(face = 'bold.italic'))) %>%
          set_names(map(correl$pairs, paste, collapse = '_')))

# Summary bubble plots ------

  insert_msg('Summary bubble plots')

  ## correlation coefficients

  correl$summ_plots <- set_names(names(correl$analysis_tbl),
                                 names(correl$analysis_tbl)) %>%
    map(function(cohort) correl$pairs %>%
          map_dfr(~correlate_variables(correl$analysis_tbl[[cohort]],
                                       variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       pub_styled = FALSE,
                                       simplify_p = FALSE))) %>%
    map(select, variable1, variable2, estimate) %>%
    map2(., map(correl$test_results,
                ~.x[c('variable1', 'variable2', 'p_adjusted', 'n')]),
         left_join,
         by = c('variable1', 'variable2')) %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 'bold', 'plain'),
        plot_lab = signif(estimate, 2))

  ## plotting

  correl$summ_plots <- map2(correl$summ_plots,
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
                              scale_fill_gradient2(low = 'steelblue',
                                                   mid = 'white',
                                                   high = 'firebrick',
                                                   midpoint = 0,
                                                   limits = c(-1, 1),
                                                   name = expression(rho)) +
                              scale_size_area(max_size = 6,
                                              limits = c(0, 1),
                                              name = expression('abs(' * rho * ')')) +
                              scale_x_discrete(limits = c('CXCL9',
                                                          'CXCL10',
                                                          'CXCL11',
                                                          'CXCR3',
                                                          'CCL11')) +
                              scale_y_discrete(limits = c('CXCL9',
                                                          'CXCL10',
                                                          'CXCL11',
                                                          'CXCR3',
                                                          'CCL11')) +
                              globals$common_theme +
                              theme(axis.title = element_blank(),
                                    axis.text = element_text(face = 'italic')) +
                              labs(title = .y,
                                   subtitle = paste('Spearman correlation, n =',
                                                    .x$n[1])))

# END ----

  insert_tail()
