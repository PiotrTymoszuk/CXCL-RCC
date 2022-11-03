# Differences in overall- and relapse-free survival between the clusters

  insert_head()

# container ------

  clust_surv <- list()

# globals -------

  insert_msg('Globals')

  ## analysis tables with he survival and lustering information

  clust_surv$analysis_tbl <-
    list(tcga_os = tcga$expression[c('patient_id',
                                     'os_days',
                                     'death')],
         tcga_tfs = tcga$expression[c('patient_id',
                                      'os_days',
                                      'tumor_death')],
         tcga_rfs = tcga$expression[c('patient_id',
                                      'rfs_days',
                                      'relapse')],
         emtab1980_os = emtab1980$expression[c('patient_id',
                                               'os_days',
                                               'death')],
         reca_os = reca$expression[c('patient_id',
                                     'os_days',
                                     'death')],
         cm10_os = cm10$expression[c('patient_id',
                                     'os_days',
                                     'death')],
         cm10_rfs = cm10$expression[c('patient_id',
                                      'rfs_days',
                                      'relapse')],
         cm25ev_os = cm25ev$expression[c('patient_id',
                                         'os_days',
                                         'death')],
         cm25ev_rfs = cm25ev$expression[c('patient_id',
                                          'rfs_days',
                                          'relapse')],
         cm25ni_os = cm25ni$expression[c('patient_id',
                                         'os_days',
                                         'death')],
         cm25ni_rfs = cm25ni$expression[c('patient_id',
                                          'rfs_days',
                                          'relapse')]) %>%
    map(~filter(.x, complete.cases(.x)))

  clust_surv$analysis_tbl <-
    map2(clust_surv$analysis_tbl,
         map(infil_clust$clust_obj[c(rep('tcga', 3),
                                     'emtab1980',
                                     'reca',
                                     rep('cm10', 2),
                                     rep('cm25ev', 2),
                                     rep('cm25ni', 2))],
             ~set_names(.x$clust_assignment[c('observation', 'clust_id')],
                        c('patient_id', 'clust_id'))),
         right_join, by = 'patient_id')

  ## model formulas

  clust_surv$formulas <- names(clust_surv$analysis_tbl) %>%
    map(function(x) switch(stri_extract(x, regex = 'os|tfs|rfs'),
                           os = Surv(os_days, death) ~ clust_id,
                           tfs = Surv(os_days, tumor_death) ~ clust_id,
                           rfs = Surv(rfs_days, relapse) ~ clust_id)) %>%
    set_names(names(clust_surv$analysis_tbl))

  ## x titles and plot titles

  clust_surv$x_labs <- names(clust_surv$analysis_tbl) %>%
    map_chr(~ifelse(stri_detect(.x, fixed = 'rfs'),
                    'RFS, days',
                    ifelse(stri_detect(.x, fixed = 'tfs'),
                           'TFS, days',
                           'OS, days')))

  clust_surv$plot_titles <- names(clust_surv$analysis_tbl) %>%
    map(~list(stri_split_fixed(.x,
                                   pattern = '_',
                                   simplify = TRUE)[, 1],
                  stri_split_fixed(.x,
                                   pattern = '_',
                                   simplify = TRUE)[, 2])) %>%
    map_chr(~paste(globals$cohorts[.x[[1]]],
                   toupper(.x[[2]]),
                   sep = ', '))

  ## n numbers per cluster

  clust_surv$n_numbers <- clust_surv$analysis_tbl %>%
    map(count, clust_id)

  ## plot tags

  clust_surv$n_tags <- clust_surv$n_numbers %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ') %>%
    map(~paste0('\n', .x))

# Construction of the survdiff and survfit objects, log-rank test, BH correction ------

  insert_msg('Testing for survival differences')

  clust_surv$survdiff_obj <-
    list(formula = clust_surv$formulas,
         data = clust_surv$analysis_tbl) %>%
    pmap(survdiff)

  clust_surv$test_results <- clust_surv$survdiff_obj %>%
    map2_dfr(., names(.), ~tibble(model = .y,
                                  chisq = .x$chisq,
                                  df = length(.x$n) - 1)) %>%
    mutate(p_value = pchisq(chisq, df, lower.tail = FALSE),
           p_adjusted = p.adjust(p_value, 'BH'),
           plot_cap = ifelse(p_adjusted < 0.05,
                             paste0('pFDR = ', signif(p_adjusted, 2),
                                    ', p raw = ', signif(p_value, 2)),
                             paste0('ns, pFDR = ', signif(p_adjusted, 2),
                                    ', p raw = ', signif(p_value, 2))))

  clust_surv$survfit_obj <-
    list(formula = clust_surv$formulas,
         data = clust_surv$analysis_tbl) %>%
    pmap(surv_fit)

# Kaplan-Meier plots -------

  insert_msg('Kaplan-Meier plots')

  clust_surv$plots <-
    list(fit = clust_surv$survfit_obj,
         title = clust_surv$plot_titles,
         xlab = clust_surv$x_labs) %>%
    pmap(ggsurvplot,
         palette = unname(globals$clust_colors),
         pval = NULL,
         legend.title = 'Infiltration cluster',
         legend.labs = names(globals$clust_colors))

  clust_surv$plots <- list(x = clust_surv$plots,
                           y = clust_surv$test_results$plot_cap,
                           z = clust_surv$n_tags) %>%
    pmap(function(x, y, z) x$plot +
           labs(subtitle = y,
                tag = z) +
           globals$common_theme)

# END -------

  insert_tail()
