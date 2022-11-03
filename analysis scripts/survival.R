# Survival in the low and high expressors of CXCL9/10/11 and CCL11
# Stratification by Mentel-Henszel test. TCGA and E-MTAB 1980 cohorts

  insert_head()

# container list ----

  surv_cut <- list()

# globals -----

  insert_msg('Globals')

  ### analysis tables

  surv_cut$analysis_tbl <-
    list(tcga_os = tcga$expression[c('patient_id',
                                     'os_days',
                                     'death',
                                     'tissue_type',
                                     globals$cxcl_genes)],
         tcga_tfs = tcga$expression[c('patient_id',
                                      'os_days',
                                      'tumor_death',
                                      'tissue_type',
                                      globals$cxcl_genes)],
         tcga_rfs = tcga$expression[c('patient_id',
                                      'rfs_days',
                                      'relapse',
                                      'tissue_type',
                                      globals$cxcl_genes)],
         emtab1980_os = emtab1980$expression[c('patient_id',
                                               'os_days',
                                               'death',
                                               globals$cxcl_genes)],
         reca_os = reca$expression[c('patient_id',
                                     'os_days',
                                     'death',
                                     globals$cxcl_genes)],
         cm10_os = cm10$expression[c('patient_id',
                                     'os_days',
                                     'death',
                                     globals$cxcl_genes)],
         cm10_rfs = cm10$expression[c('patient_id',
                                      'rfs_days',
                                      'relapse',
                                      globals$cxcl_genes)],
         cm25ev_os = cm25ev$expression[c('patient_id',
                                         'os_days',
                                         'death',
                                         globals$cxcl_genes)],
         cm25ev_rfs = cm25ev$expression[c('patient_id',
                                          'rfs_days',
                                          'relapse',
                                          globals$cxcl_genes)],
         cm25ni_os = cm25ni$expression[c('patient_id',
                                         'os_days',
                                         'death',
                                         globals$cxcl_genes)],
         cm25ni_rfs = cm25ni$expression[c('patient_id',
                                          'rfs_days',
                                          'relapse',
                                          globals$cxcl_genes)]) %>%
    map(~filter(.x, complete.cases(.x)))

  ## x titles and plot titles

  surv_cut$x_labs <- names(surv_cut$analysis_tbl) %>%
    map_chr(~ifelse(stri_detect(.x, fixed = 'rfs'),
                    'RFS, days',
                    ifelse(stri_detect(.x, fixed = 'tfs'),
                           'TFS, days',
                           'OS, days')))

  surv_cut$plot_titles <- names(surv_cut$analysis_tbl) %>%
    map(~list(stri_split_fixed(.x,
                               pattern = '_',
                               simplify = TRUE)[, 1],
              stri_split_fixed(.x,
                               pattern = '_',
                               simplify = TRUE)[, 2])) %>%
    map_chr(~paste(globals$cohorts[.x[[1]]],
                   toupper(.x[[2]]),
                   sep = ', '))

  ## time and event vars and minimal n numbers (20% of the cohort size)

  surv_cut$time_vars <- names(surv_cut$analysis_tbl) %>%
    map_chr(~ifelse(stri_detect(.x, fixed = 'rfs'), 'rfs_days', 'os_days')) %>%
    set_names(names(surv_cut$analysis_tbl))

  surv_cut$event_vars <- names(surv_cut$analysis_tbl) %>%
    map_chr(~ifelse(stri_detect(.x, fixed = 'rfs'),
                    'relapse',
                    ifelse(stri_detect(.x, fixed = 'tfs'),
                           'tumor_death',
                           'death'))) %>%
    set_names(names(surv_cut$analysis_tbl))

  surv_cut$min_n <- surv_cut$analysis_tbl %>%
    map_dbl(~nrow(.x) * 0.2)

  ## parallel backend

  plan('multisession')

# stratification ------

  insert_msg('Stratification')

  surv_cut[c('cxcl9_cutoff',
             'cxcl10_cutoff',
             'cxcl11_cutoff',
             'ccl11_cutoff',
             'cxcr3_cutoff')] <- globals$cxcl_genes %>%
    map(function(gene) list(data = surv_cut$analysis_tbl,
                            time = surv_cut$time_vars,
                            event = surv_cut$event_vars,
                            min_n = surv_cut$min_n) %>%
          future_pmap(find_cutoff,
                      variable = gene,
                      rho = 0,
                      .parallel = FALSE))

# Optimal cutoff stats, p value correction within the cohort ------

  insert_msg('Cutoff summary, p correction')

  surv_cut$test_results <-
    surv_cut[c('cxcl9_cutoff',
               'cxcl10_cutoff',
               'cxcl11_cutoff',
               'ccl11_cutoff',
               'cxcr3_cutoff')] %>%
    map(~map2_dfr(.x, names(.x),
                  ~mutate(summary(.x)[1, ], model = .y)))

  ## stripping the information to be shown in the plot captions
  ## gene-wise correction form multiple testing

  surv_cut$test_results <- surv_cut$test_results %>%
    map(mutate,
        p_adjusted = p.adjust(p_value, 'BH'),
        plot_cap = paste0('pFDR = ', signif(p_adjusted, 2),
                          ', p raw = ', signif(p_value, 2)),
        plot_cap = ifelse(p_adjusted < 0.05,
                          plot_cap,
                          paste('ns,', plot_cap)),
        plot_cap = paste0('cutoff = ',
                          signif(cutoff, 2),
                          ', ',
                          plot_cap),
        study = stri_split_fixed(model, pattern = '_', simplify = TRUE)[, 1],
        response = stri_split_fixed(model, pattern = '_', simplify = TRUE)[, 2],
        x_ax = paste(globals$cohorts[study],
                     toupper(response),
                     sep = ', '))

# significance plots -------

  insert_msg('Summary bar plots with p values for each gene')

  surv_cut$summ_plots <-
    list(data = surv_cut$test_results,
         plot_title = list(expression(italic('CXCL9') * ', survival'),
                           expression(italic('CXCL10') * ', survival'),
                           expression(italic('CXCL11') * ', survival'),
                           expression(italic('CCL11') * ', survival'),
                           expression(italic('CXCR3') * ', survival'))) %>%
    pmap(plot_signifcant,
         p_variable = 'p_adjusted',
         label_variable = 'x_ax',
         top_significant = 20,
         plot_subtitle = 'Optimal cutoff, Mentel-Henszel test',
         x_lab = expression('-log'[10] * ' pFDR'),
         cust_theme = globals$common_theme)

# Kaplan-Meier plots -----

  insert_msg('Kaplan-Meier plots')

  surv_cut[c('cxcl9_plots',
             'cxcl10_plots',
             'cxcl11_plots',
             'ccl11_plots',
             'cxcr3_plots')] <- surv_cut[c('cxcl9_cutoff',
                                           'cxcl10_cutoff',
                                           'cxcl11_cutoff',
                                           'ccl11_cutoff',
                                           'cxcr3_cutoff')] %>%
    map(function(cut_obj) list(x = cut_obj,
                               title = surv_cut$plot_titles,
                               xlab = surv_cut$x_labs) %>%
          pmap(plot,
               type = 'km') %>%
          map(~.x$plot +
                globals$common_theme))

  ## appending with the raw and p values

  surv_cut[c('cxcl9_plots',
             'cxcl10_plots',
             'cxcl11_plots',
             'ccl11_plots',
             'cxcr3_plots')] <- list(x = surv_cut[c('cxcl9_plots',
                                                    'cxcl10_plots',
                                                    'cxcl11_plots',
                                                    'ccl11_plots',
                                                    'cxcr3_plots')],
                                     y = surv_cut$test_results) %>%
    pmap(function(x, y) map2(x,
                             y$plot_cap,
                             ~.x + labs(subtitle = .y)))


# END -----

  plan('sequential')

  insert_tail()
