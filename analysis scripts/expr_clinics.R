# Expression of CXCL9, 10 and 11 and clinical features.

  insert_head()

# container list -----

  exp_clinic <- list()

# globals -----

  insert_msg('Globals setup')

  ## analysis tables

  exp_clinic$analysis_tbl <- list(tcga = tcga$expression %>%
                                    filter(tissue_type == 'Tumor'),
                                  emtab1980 = emtab1980$expression,
                                  gse73731 = gse73731$expression,
                                  gse167093 = gse167093$expression %>%
                                    filter(tissue_type == 'Tumor'),
                                  reca = reca$expression,
                                  cm10 = cm10$expression,
                                  cm25ev = cm25ev$expression,
                                  cm25ni = cm25ni$expression) %>%
    map(~select(.x,
                any_of(globals$clin_vars[!globals$clin_vars %in% c('neoadjuvant')]),
                all_of(globals$cxcl_genes)))


  ## manual re-coding of the variables
  ## stratification of age, exclusion of samples with undetermined
  ## pm, pn, pt stage and grade
  ## classification of the tumors with and without shrinkage
  ## in the ChecMate cohorts

  exp_clinic$analysis_tbl[c('tcga',
                            'emtab1980',
                            'gse167093',
                            'reca',
                            'cm10',
                            'cm25ev',
                            'cm25ni')] <-
    exp_clinic$analysis_tbl[c('tcga',
                              'emtab1980',
                              'gse167093',
                              'reca',
                              'cm10',
                              'cm25ev',
                              'cm25ni')] %>%
    map(mutate, age = cut(age, c(0, 50, 60, 70, 100)))

  exp_clinic$analysis_tbl[c('tcga',
                            'emtab1980',
                            'reca')] <-
    exp_clinic$analysis_tbl[c('tcga',
                              'emtab1980',
                              'reca')] %>%
    map(mutate,
        pn_stage = ifelse(pn_stage == 'NX',
                          NA, as.character(pn_stage)),
        pm_stage = ifelse(pm_stage == 'MX',
                          NA, as.character(pm_stage)),
        pn_stage = factor(pn_stage),
        pm_stage = factor(pm_stage))

  exp_clinic$analysis_tbl[c('tcga',
                            'emtab1980',
                            'gse73731',
                            'gse167093',
                            'reca')] <-
    exp_clinic$analysis_tbl[c('tcga',
                              'emtab1980',
                              'gse73731',
                              'gse167093',
                              'reca')] %>%
    map(mutate,
        tumor_grade = ifelse(tumor_grade %in% c('GX', 'undetermined'),
                             NA, as.character(tumor_grade)),
        tumor_grade = factor(tumor_grade))

  exp_clinic$analysis_tbl[c('cm10', 'cm25ev', 'cm25ni')] <-
    exp_clinic$analysis_tbl[c('cm10', 'cm25ev', 'cm25ni')] %>%
    map(mutate,
        tumor_delta_vol = ifelse(tumor_delta_vol < 0, 'yes', 'no'),
        tumor_delta_vol = factor(tumor_delta_vol, c('no', 'yes')),
        biresponse = car::recode(as.character(response_type),
                                 "'CR' = 'CR/PR';
                                 'PR' = 'CR/PR';
                                 'SD' = 'SD/PD';
                                 'PD' = 'SD/PD';
                                 'CRPR' = 'CR/PR'"),
        biresponse = factor(biresponse, c('SD/PD', 'CR/PR')))

  ## manual recoding of the tumor stage information: simplifying

  exp_clinic$analysis_tbl[c('tcga',
                            'emtab1980',
                            'reca',
                            'gse73731',
                            'gse167093')] <-
    exp_clinic$analysis_tbl[c('tcga',
                              'emtab1980',
                              'reca',
                              'gse73731',
                              'gse167093')] %>%
    map(mutate,
        pt_stage = stri_extract(pt_stage, regex = '^T\\d{1}'),
        pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')))

  ## manual recoding of the grade: common levels

  exp_clinic$analysis_tbl[c('emtab1980',
                            'reca',
                            'gse73731',
                            'gse167093')] <-
    exp_clinic$analysis_tbl[c('emtab1980',
                              'reca',
                              'gse73731',
                              'gse167093')] %>%
    map(mutate,
        tumor_grade = ifelse(!is.na(tumor_grade),
                             paste0('G', as.character(tumor_grade)),
                             NA),
        tumor_grade = factor(tumor_grade))

  ## independent clinical variables
  ## there's only one arm in the CheckMate 010 cohort,
  ## an rhabdoid differentiation
  ## the variable is hence removed

  exp_clinic$variables <- exp_clinic$analysis_tbl %>%
    map(names) %>%
    map(~.x[.x %in% c(globals$clin_vars, 'biresponse')])

  exp_clinic$variables$cm10 <-
    exp_clinic$variables$cm10[!exp_clinic$variables$cm10 %in% c('rhab_diff')]

  ## effect size statistics, dependent on the variable type
  ## and number of levels

  exp_clinic$types <-
    map2(exp_clinic$analysis_tbl,
         exp_clinic$variables,
         ~map(.x[.y], levels) %>%
           map_dbl(length)) %>%
    map(~ifelse(.x > 2, 'kruskal_etasq', 'wilcoxon_r'))

  ## parallel backend

  plan('multisession')

# Testing ----

  insert_msg('Testing')

  exp_clinic$test_results <- names(exp_clinic$variables) %>%
    future_map(function(cohort) list(split_factor = exp_clinic$variables[[cohort]],
                                     types = exp_clinic$types[[cohort]]) %>%
                 pmap(function(split_factor,
                               types) compare_variables(exp_clinic$analysis_tbl[[cohort]] %>%
                                                          filter(!is.na(.data[[split_factor]])),
                                                        variables = globals$cxcl_genes,
                                                        split_factor = split_factor,
                                                        what = 'eff_size',
                                                        types = types,
                                                        pub_styled = TRUE,
                                                        ci = FALSE,
                                                        simplify_p = FALSE,
                                                        adj_method = 'BH',
                                                        exact = FALSE)) %>%
                 set_names(exp_clinic$variables[[cohort]]) %>%
                 map2_dfr(., names(.), ~mutate(.x, split_factor = .y)),
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', ')) %>%
    set_names(names(exp_clinic$variables))

# Plotting ----

  insert_msg('Plotting: violin plots')

  exp_clinic$plots <- list(cohort = set_names(names(exp_clinic$variables),
                                              names(exp_clinic$variables)),
                           suffix = globals$cohorts) %>%
    future_pmap(function(cohort, suffix) set_names(globals$cxcl_genes,
                                            globals$cxcl_genes) %>%
           map(function(gene) list(split_factor = exp_clinic$variables[[cohort]],
                                   x_lab = globals$var_labs[exp_clinic$variables[[cohort]]],
                                   plot_subtitle = exp_clinic$test_results[[cohort]] %>%
                                     filter(variable == gene) %>%
                                     .$plot_cap) %>%
                 pmap(function(split_factor,
                               x_lab,
                               plot_subtitle) plot_variable(exp_clinic$analysis_tbl[[cohort]] %>%
                                                              filter(!is.na(.data[[split_factor]])),
                                                            variable = gene,
                                                            split_factor = split_factor,
                                                            type = 'violin',
                                                            plot_title = paste(gene, suffix, sep = ', '),
                                                            plot_subtitle = plot_subtitle,
                                                            x_lab = x_lab,
                                                            y_lab = expression('log'[2]*' expression'),
                                                            cust_theme = globals$common_theme)) %>%
                 map(~.x +
                       theme(plot.title = element_text(face = 'bold.italic')) +
                       scale_fill_brewer(palette = 1) +
                       labs(tag = .x$labels$tag %>%
                              stri_replace_all(fixed = '\n',
                                               replacement = ', ') %>%
                              paste0('\n', .))) %>%
                 set_names(exp_clinic$variables[[cohort]])),
           .options = furrr_options(seed = TRUE))

# Summary plots: effect sizes and significance for tumor, node and meta stages ------

  insert_msg('Summary plots for tumor stages')

  ## etest objects

  exp_clinic$etest_stages <- names(exp_clinic$variables) %>%
    future_map(function(cohort) list(split_factor = exp_clinic$variables[[cohort]],
                                     types = exp_clinic$types[[cohort]]) %>%
                 pmap(function(split_factor,
                               types) compare_variables(exp_clinic$analysis_tbl[[cohort]] %>%
                                                          filter(!is.na(.data[[split_factor]])),
                                                        variables = globals$cxcl_genes,
                                                        split_factor = split_factor,
                                                        what = 'eff_size',
                                                        types = types,
                                                        pub_styled = FALSE,
                                                        ci = FALSE,
                                                        simplify_p = FALSE,
                                                        adj_method = 'BH',
                                                        exact = FALSE)) %>%
                 set_names(exp_clinic$variables[[cohort]]) %>%
                 map2_dfr(., names(.), ~mutate(.x, split_factor = .y)),
               .options = furrr_options(seed = TRUE)) %>%
    set_names(names(exp_clinic$variables))

  ## filtering and merging

  exp_clinic$etest_stages <- exp_clinic$etest_stages %>%
    map(filter,
        split_factor %in% c('tumor_grade', 'pt_stage', 'pm_stage', 'pn_stage')) %>%
    map2_dfr(., globals$cohorts[names(.)],
             ~mutate(.x, variable = paste(variable, .y, sep = ', ')))

  ## plotting

  exp_clinic$stages_summ_plots <- exp_clinic$etest_stages %>%
    dlply(.(split_factor),
          function(x) structure(x, class = c('etest', 'data.frame'))) %>%
    list(x = .,
         plot_title = globals$var_labs[names(.)]) %>%
    pmap(plot,
         plot_subtitle = 'significance vs effect size',
         cust_theme = globals$common_theme,
         show_labels = 'signif')

  exp_clinic$stages_summ_plots <-
    list(x = exp_clinic$stages_summ_plots,
         y = list('Effect size, Wilcoxon r',
                  'Effect size, Wilcoxon r',
                  expression('Effect size, ' * eta^2),
                  expression('Effect size, ' * eta^2)),
         z = c('Mann-Whitney U test, FDR correction',
               'Mann-Whitney U test, FDR correction',
               'Kruskal-Wallis test, FDR correction',
               'Kruskal-Wallis test, FDR correction')) %>%
    pmap(function(x, y, z) x +
           geom_hline(yintercept = -log10(0.05),
                      linetype = 'dashed') +
           labs(x = y,
                subtitle = z,
                y = expression('-log'[10] * ' pFDR')) +
           theme(plot.tag = element_blank()))

# END -----

  plan('sequential')

  insert_tail()
