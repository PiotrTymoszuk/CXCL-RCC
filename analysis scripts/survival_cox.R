# Univariable Cox modeling of the CXCL9 expression strata


  insert_head()

# container -----

  surv_cox <- list()

# globals ------

  insert_msg('Globals setup')

  ## the analysis tables are taken over from the survcut script

  surv_cox$analysis_tbl <- surv_cut$cxcl9_cutoff %>%
    map(~.x$data)

  ## model formulas

  surv_cox$formulas <- names(surv_cox$analysis_tbl) %>%
    map(function(x) switch(stri_extract(x, regex = 'os|tfs|rfs'),
                           os = Surv(os_days, death) ~ CXCL9_strata,
                           tfs = Surv(os_days, tumor_death) ~ CXCL9_strata,
                           rfs = Surv(rfs_days, relapse) ~ CXCL9_strata)) %>%
    set_names(names(surv_cox$analysis_tbl))

  ## x titles and plot titles

  surv_cox$x_labs <- names(surv_cox$analysis_tbl) %>%
    map_chr(~ifelse(stri_detect(.x, fixed = 'rfs'),
                    'RFS, days',
                    ifelse(stri_detect(.x, fixed = 'tfs'),
                           'TFS, days',
                           'OS, days')))

  surv_cox$plot_titles <- names(surv_cox$analysis_tbl) %>%
    map(~list(stri_split_fixed(.x,
                               pattern = '_',
                               simplify = TRUE)[, 1],
              stri_split_fixed(.x,
                               pattern = '_',
                               simplify = TRUE)[, 2])) %>%
    map_chr(~paste(globals$cohorts[.x[[1]]],
                   toupper(.x[[2]]),
                   sep = ', '))

# model construction -------

  insert_msg('Model construction')

  ## sorry: working with AST, because the coxph model needs to be bound to data

  surv_cox$models <- map2(surv_cox$formulas,
                          surv_cox$analysis_tbl,
                          ~rlang::call2('coxph', formula = .x, data = .y)) %>%
    map(eval) %>%
    map2(., surv_cox$analysis_tbl, as_coxex)

  ## fit stats

  surv_cox$fit_stats <- surv_cox$models %>%
    map(summary, type = 'fit') %>%
    map2_dfr(., names(.), ~mutate(.x, model = .y)) %>%
    mutate(response = stri_extract(model, regex = 'os|tfs|rfs'),
           study = stri_split_fixed(model,
                                    pattern = '_',
                                    simplify = TRUE)[, 1],
           x_ax = paste0(globals$cohorts[study],
                         '\ntotal n = ',
                         n_complete,
                         ', event n = ',
                         n_events))

  ## model assumptions

  surv_cox$assumptions <- surv_cox$models %>%
    map(summary, type = 'assumptions')

  ## model inference

  surv_cox$inference <- surv_cox$models %>%
    map(summary, type = 'inference') %>%
    map2_dfr(., names(.), ~mutate(.x, model = .y)) %>%
    mutate(response = stri_extract(model, regex = 'os|tfs|rfs'),
           study = stri_split_fixed(model,
                                    pattern = '_',
                                    simplify = TRUE)[, 1],
           p_adjusted = p.adjust(p_value, 'BH'),
           estimate = exp(estimate),
           lower_ci = exp(lower_ci),
           upper_ci = exp(upper_ci),
           x_ax = paste0(globals$cohorts[study],
                         '\ntotal n = ',
                         n_complete,
                         ', high n = ',
                         n))

# obtaining metaanalysis HR, inverted variance method ------

  insert_msg('Meta HR')

  surv_cox$meta_obj <- surv_cox$inference %>%
    mutate(estimate = log(estimate)) %>%
    dlply(.(response),
          function(x) metagen(TE = x$estimate,
                              seTE = abs(x$se),
                              sm = ''))

  ## extracting the stats

  surv_cox$meta_inference <- surv_cox$meta_obj %>%
    map(~tibble(parameter = 'CXCL9_stratahigh',
                variable = 'CXCL9_strata',
                level = 'high',
                n = 0,
                n_complete = 0,
                stat_name = 'z',
                stat = .x$zval.common,
                estimate = exp(.x$TE.common),
                se = .x$seTE.common,
                lower_ci = exp(.x$lower.common),
                upper_ci = exp(.x$upper.common),
                p_value = .x$pval.common,
                p_adjusted = .x$pval.common)) %>%
    map2_dfr(., names(.), ~mutate(.x, model = .y)) %>%
    mutate(response = stri_extract(model, regex = 'os|tfs|rfs'),
           study = 'Meta',
           x_ax = 'Meta, inv. variance')

  surv_cox$inference <- rbind(surv_cox$inference,
                              surv_cox$meta_inference) %>%
    mutate(type = ifelse(study == 'Meta', 'Meta-analysis', 'Study'))

# Plotting the C indexes for OS, TFS and RFS ------

  insert_msg('Plotting the C indexes')

  surv_cox$c_forest_plots <- surv_cox$fit_stats %>%
    mutate(p_value = ifelse(lower_ci > 0.5, 0.02, 1)) %>%
    dlply(.(response)) %>%
    list(data = .,
         plot_title = paste(toupper(names(.)), 'CXCL9 strata', sep = ', ')) %>%
    pmap(plot_top,
         regulation_variable = 'c_index',
         label_variable = 'x_ax',
         p_variable = 'p_value',
         lower_ci_variable = 'lower_ci',
         upper_ci_variable = 'upper_ci',
         regulation_level = 0.5,
         cust_theme = globals$common_theme,
         plot_subtitle = 'Concordance index, univariable Cox models, optimal cutoff',
         x_lab = 'C-index \u00B1 95% CI',
         show_txt = TRUE,
         show_ci_txt = TRUE) %>%
    map(~.x +
          scale_x_continuous(limits = c(0.45, 1)) +
          geom_vline(xintercept = 0.5,
                     linetype = 'dashed'))

# Plotting R-squares -----

  insert_msg('Rsq plots')

  surv_cox$rsq_plots <- surv_cox$fit_stats %>%
    mutate(p_value = 0.02) %>%
    dlply(.(response)) %>%
    list(data = .,
         plot_title = paste(toupper(names(.)), 'CXCL9 strata', sep = ', ')) %>%
    pmap(plot_top,
         regulation_variable = 'raw_rsq',
         label_variable = 'x_ax',
         p_variable = 'p_value',
         regulation_level = 0,
         cust_theme = globals$common_theme,
         plot_subtitle = expression('Raw R'^2 * 'univariable Cox models, optimal cutoff'),
         x_lab = expression('raw R'^2),
         show_txt = TRUE,
         show_ci_txt = TRUE)

# HR forest plots -------

  insert_msg('Forest plots with the HRs')

  surv_cox$hr_forest_plots <- surv_cox$inference %>%
    mutate(significance = ifelse(p_adjusted >= 0.05,
                                'ns',
                                ifelse(estimate > 1,
                                       'unfavorable',
                                       'favorable')),
           significance  = factor(significance,
                                  c('favorable', 'unfavorable', 'ns'))) %>%
    dlply(.(response)) %>%
    list(data = .,
         plot_title = paste(toupper(names(.)), 'CXCL9 strata', sep = ', ')) %>%
    pmap(function(data, plot_title) ggplot(data,
                                           aes(x = estimate,
                                               y = reorder(x_ax, estimate),
                                               color = significance,
                                               shape = type)) +
           geom_vline(xintercept = 1,
                      linetype = 'dashed') +
           geom_errorbarh(aes(xmin = lower_ci,
                              xmax = upper_ci),
                          height = 0) +
           geom_point(aes(size = type)) +
           geom_text(aes(label = paste0(signif(estimate, 2),
                                        ' [',
                                        signif(lower_ci, 2),
                                        ' - ',
                                        signif(upper_ci, 2),
                                        ']')),
                     size = 2.75,
                     hjust = 0.5,
                     vjust = -1,
                     show.legend = FALSE) +
           scale_shape_manual(values = c('Study' = 16,
                                         'Meta-analysis' = 18),
                              name = '') +
           scale_size_manual(values = c('Study' = 2,
                                         'Meta-analysis' = 4),
                              name = '') +
           scale_color_manual(values = c(favorable = 'steelblue',
                                         unfavorable = 'firebrick',
                                         ns = 'gray60'),
                              name = 'Risk') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = plot_title,
                subtitle = 'Hazard ratio, high vs low expressors, univariable Cox models, optimal cutoff',
                x = 'HR \u00B1 95% CI'))

# END ------

  insert_tail()
