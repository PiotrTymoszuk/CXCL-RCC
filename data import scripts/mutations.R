# Compares the differences in total mutation burden (TMB) and frequency of
# the most common mutations between the high and low expressors of CXCL10
# and CCL11. Tools: Mann-Whitney U test and Chi-squared test

  insert_head()

# container list -----

  mut_strata <- list()

# globals: analysis table -----

  insert_msg('Globals setup')

  mut_strata$analysis_tbl <- left_join(somut$mut_summary,
                                       dge$anlysis_tbl[c('patient_id',
                                                         'cxcl10_group',
                                                         'ccl11_group')],
                                       by = 'patient_id') %>%
    filter(complete.cases(.))

# TMB comparison ----

  insert_msg('Comaring the TMB between the expression strata')

  ## testing

  mut_strata$tmb_results <- c(cxcl10 = 'cxcl10_group',
                              ccl11 = 'ccl11_group') %>%
    map(~compare_variables(mut_strata$analysis_tbl,
                           variables = 'mut_count',
                           split_factor = .x,
                           what = 'eff_size',
                           types = 'wilcoxon_r',
                           ci = FALSE,
                           pub_styled = TRUE,
                           simplify_p = FALSE)) %>%
    map2_dfr(., names(.),
             ~mutate(.x,
                     split_factor = .y,
                     plot_cap = paste(eff_size, significance, sep = ', ')))

  ## plots

  mut_strata$tmb_plots <- list(split_factor = c(cxcl10 = 'cxcl10_group',
                                                ccl11 = 'ccl11_group'),
                               x_lab = c('CXCL10 strata',
                                         'CCL11 strata'),
                               plot_subtitle = mut_strata$tmb_results$plot_cap,
                               plot_title = c('CXCL10: TMB',
                                              'CCL11: TMB')) %>%
    pmap(plot_variable,
         mut_strata$analysis_tbl,
         variable = 'mut_count',
         type = 'violin',
         cust_theme = globals$common_theme,
         y_lab = '# mutations') %>%
    map(~.x +
          scale_y_continuous(trans = 'log10') +
          scale_fill_manual(values = c(high = 'firebrick',
                                       low = 'steelblue')) +
          labs(tag = stri_replace(.x$labels$tag,
                                  fixed = '\n',
                                  replacement = ', ') %>%
                 paste0('\n', .)))

# Comparing the mutation frequencies between the strata -----

  insert_msg('Comparing the mutation frequencies')


