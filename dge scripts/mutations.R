# Compares the differences in total mutation burden (TMB) and frequency of
# the most common mutations between the high and low expressors of CXCL10
# and CCL11. Tools: Mann-Whitney U test and Chi-squared test

  insert_head()

# container list -----

  mut_strata <- list()

# globals: analysis table -----

  insert_msg('Globals setup')

  mut_strata$analysis_tbl <- left_join(somut$mut_summary,
                                       dge$analysis_tbl$tcga[c('patient_id',
                                                               'cxcl9_group',
                                                               'ccl11_group')],
                                       by = 'patient_id') %>%
    filter(complete.cases(.)) %>%
    select(patient_id,
           ends_with('group'),
           mut_count,
           somut$top_mutations$variable)

  ## n numbers

  mut_strata$n_numbers <- c(cxcl9 = 'cxcl9_group',
                            ccl11 = 'ccl11_group') %>%
    map(~count(mut_strata$analysis_tbl, .data[[.x]]))

  mut_strata$n_tags <- mut_strata$n_numbers %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ') %>%
    map(~paste0('\n', .x))

# TMB comparison ----

  insert_msg('Comaring the TMB between the expression strata')

  ## testing

  mut_strata$tmb_results <- c(cxcl9 = 'cxcl9_group',
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

  mut_strata$tmb_plots <- list(split_factor = c(cxcl9 = 'cxcl9_group',
                                                ccl11 = 'ccl11_group'),
                               x_lab = c('CXCL9 strata',
                                         'CCL11 strata'),
                               plot_subtitle = mut_strata$tmb_results$plot_cap,
                               plot_title = c('CXCL9: TMB',
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

  mut_strata$freq_test <- c(cxcl9 = 'cxcl9_group',
                            ccl11 = 'ccl11_group') %>%
    map(~compare_variables(mut_strata$analysis_tbl,
                           variables = somut$top_mutations$variable,
                           split_factor = .x,
                           what = 'test',
                           types = 'chisq_test',
                           ci = FALSE,
                           pub_styled = FALSE,
                           adj_method = 'BH',
                           simplify_p = FALSE,
                           .parallel = TRUE,
                           .paropts = furrr_options(seed = TRUE,
                                                    globals = c('mut_strata')))) %>%
    map(mutate,
        variable = stri_replace(variable,
                                fixed = '_mut',
                                replacement = ''),
        plot_lab = ifelse(p_value < 0.05, variable, NA),
        plot_cap = paste0('V = ', signif(estimate, 2),
                         '\nraw p = ', signif(p_value, 2),
                         ', pFDR = ', signif(p_adjusted, 2)))

# Summary plots -------

  insert_msg('Summary plots of the testing results')

  mut_strata$freq_test_plot <- list(x = mut_strata$freq_test,
                                    plot_title = c('CXCL9: mutations',
                                                   'CCL11: mutations'),
                                    plot_subtitle = list(expression('CXCL9 high vs low tumors, '*chi^2*' test'),
                                                         expression('CCL11 high vs low tumors, '*chi^2*' test'))) %>%
    pmap(plot,
         show_labels = 'none',
         fontface = 'italic',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          geom_hline(yintercept = -log10(0.05),
                     linetype = 'dashed') +
          geom_text_repel(aes(label = plot_lab),
                          size = 2.5,
                          fontface = 'italic',
                          max.overlaps = 25))

  mut_strata$freq_test_plot <- map2(mut_strata$freq_test_plot,
                                    mut_strata$n_tags,
                                    ~.x +
                                      labs(tag = .y,
                                           x = 'Effect size, Cramer V',
                                           y = expression('-log'[10]*' pFDR')))

# Stack bar plots with the mutation frequencies in the expression strata -----

  insert_msg('Plots for particular alleles')

  plan('multisession')

  mut_strata$freq_plots_cxcl9 <- list(variable = somut$top_mutations$variable,
                                       plot_title = somut$top_mutations$gene,
                                       plot_subtitle = mut_strata$freq_test$cxcl9$plot_cap) %>%
    pmap(plot_variable,
         mut_strata$analysis_tb,
         split_factor = 'cxcl9_group',
         type = 'stack',
         scale = 'percent',
         cust_theme = globals$common_theme,
         x_lab = 'CXCL9 strata',
         y_lab = '% of expression strata')

  mut_strata$freq_plots_ccl11 <- list(variable = somut$top_mutations$variable,
                                      plot_title = somut$top_mutations$gene,
                                      plot_subtitle = mut_strata$freq_test$ccl11$plot_cap) %>%
    pmap(plot_variable,
         mut_strata$analysis_tb,
         split_factor = 'ccl11_group',
         type = 'stack',
         scale = 'percent',
         cust_theme = globals$common_theme,
         x_lab = 'CCL11 strata',
         y_lab = '% of expression strata')

  plan('sequential')

  ## adjustment

  mut_strata[c('freq_plots_cxcl9',
               'freq_plots_ccl11')] <- mut_strata[c('freq_plots_cxcl9',
                                                    'freq_plots_ccl11')] %>%
    map(~map(.x,
             ~.x +
               scale_fill_manual(values = c('cornsilk2', 'darkolivegreen3'),
                                 labels = c('WT', 'mutant')) +
               labs(tag = stri_replace(.x$labels$tag,
                                       fixed = '\n',
                                       replacement = ', ') %>%
                      paste0('\n', .)) +
               theme(plot.title = element_text(face = 'bold.italic')))) %>%
    map(set_names, somut$top_mutations$gene)

# For information: frequencies of the most common mutated protein coded alleles ----

  insert_msg('Mutated allele frequency')

  mut_strata$cohort_freq_plot <- somut$mut_freq %>%
    top_n(20, percent) %>%
    ggplot(aes(x = percent,
               y = reorder(gene, percent))) +
    geom_bar(stat = 'identity',
             color = 'black',
             fill = 'darkolivegreen3') +
    geom_text(aes(label = signif(percent, 3)),
              size = 2.75,
              hjust = -0.4,
              vjust = 0.5) +
    scale_x_continuous(limits = c(0, 55),
                       breaks = seq(0, 50, by = 10)) +
    globals$common_theme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = 'italic')) +
    labs(title = 'Top 20 mutated genes',
         subtitle = paste('KIRC cohort, n =', nrow(mut_strata$analysis_tbl)),
         x = '% of complete samples')

# END -----

  insert_tail()
