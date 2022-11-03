# Compares clinical features of the CXCL10 and CCL11 expression strata

  insert_head()

# container list -----

  clin_strata <- list()

# globals ------

  insert_msg('Globals setup')

  ## analysis tables

  clin_strata$analysis_tbl <- dge$analysis_tbl %>%
    map(~select(.x,
                cxcl9_group,
                ccl11_group,
                any_of(globals$clin_vars)))

  clin_strata$analysis_tbl[c('tcga', 'emtab1980', 'reca')] <-
    clin_strata$analysis_tbl[c('tcga', 'emtab1980', 'reca')] %>%
    map(mutate,
        pn_stage = ifelse(pn_stage == 'NX', NA, as.character(pn_stage)),
        pn_stage = factor(pn_stage),
        pm_stage = ifelse(pm_stage == 'MX', NA, as.character(pm_stage)),
        pm_stage = factor(pm_stage))

  ## variables
  ## removing rhabdoid differentiation for Checkmate010

  clin_strata$variables <- dge$analysis_tbl %>%
    map(names) %>%
    map(~.x[.x %in% globals$clin_vars])

  clin_strata$variables$cm10 <-
    clin_strata$variables$cm10[clin_strata$variables$cm10 != 'rhab_diff']

  ## effect size types

  clin_strata$eff_sizes <- map2(clin_strata$variables,
                                clin_strata$analysis_tbl,
                                ~.y[.x]) %>%
    map(~map_lgl(.x, is.numeric)) %>%
    map(~ifelse(.x, 'cohen_d', 'cramer_v'))

  ## plot types

  clin_strata$plot_types <- map2(clin_strata$variables,
                                 clin_strata$analysis_tbl,
                                 ~.y[.x]) %>%
    map(~map_lgl(.x, is.numeric)) %>%
    map(~ifelse(.x, 'violin', 'stack'))

  ## parallel backend

  plan('multisession')

# serial testing ------

  insert_msg('Serial testing')

  clin_strata$test_results <- set_names(names(clin_strata$variables),
                                        names(clin_strata$variables)) %>%
    future_map(function(cohort) c(cxcl9 = 'cxcl9_group',
                                  ccl11 = 'ccl11_group') %>%
                 map(~compare_variables(clin_strata$analysis_tbl[[cohort]],
                                        variables = clin_strata$variables[[cohort]],
                                        split_factor = .x,
                                        what = 'eff_size',
                                        types = clin_strata$eff_sizes[[cohort]],
                                        ci = FALSE,
                                        exact = FALSE,
                                        adj_method = 'BH',
                                        pub_styled = TRUE,
                                        simplify_p = FALSE)) %>%
                 map(mutate, plot_cap = paste(eff_size, significance, sep = ', ')),
               .options = furrr_options(seed = TRUE))


# serial plotting -------

  insert_msg('Serial plotting')

  clin_strata$plots <- set_names(names(clin_strata$variables),
                                 names(clin_strata$variables)) %>%
    future_map(function(cohort) list(chemo = c(cxcl9 = 'cxcl9_group',
                                               ccl11 = 'ccl11_group'),
                                     chemo_lab = c(cxcl9 = 'CXCL9 strata',
                                                   ccl11 = 'CCL11 strata'),
                                     stats = clin_strata$test_results[[cohort]]) %>%
                 pmap(function(chemo, chemo_lab, stats) list(variable = clin_strata$variables[[cohort]],
                                                             type = clin_strata$plot_types[[cohort]],
                                                             plot_title = paste(globals$var_labs[clin_strata$variables[[cohort]]],
                                                                                globals$cohorts[[cohort]],
                                                                                sep = ', '),
                                                             plot_subtitle = stats$plot_cap) %>%
                        pmap(plot_variable,
                             clin_strata$analysis_tbl[[cohort]],
                             split_factor = chemo,
                             scale = 'percent',
                             x_lab = chemo_lab,
                             cust_theme = globals$common_theme) %>%
                        map(~.x +
                              scale_fill_brewer(palette = 1) +
                              labs(tag = .x$labels$tag %>%
                                     stri_replace_all(fixed = '\n', replacement = ', ') %>%
                                     paste0('\n', .))) %>%
                        set_names(clin_strata$variables[[cohort]])),
               .options = furrr_options(seed = TRUE))

# END -----

  plan('sequential')

  insert_tail()
