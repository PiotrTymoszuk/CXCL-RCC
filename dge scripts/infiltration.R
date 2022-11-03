# Differences in infiltration between the CXCL9 expression strata

  insert_head()

# containers -----

  cxcl9_infil <- list()

# globals -----

  insert_msg('Globals')

  cxcl9_infil$analysis_tbl <- dge$analysis_tbl %>%
    map(select, cxcl9_group, all_of(globals$imm_estimates)) %>%
    map(mutate,
        cd8_treg_ratio = `T cell CD8+`/`T cell regulatory (Tregs)`,
        cd8_treg_ratio = ifelse(is.na(cd8_treg_ratio) | is.infinite(cd8_treg_ratio),
                                0, cd8_treg_ratio),
        cd8_m2_ratio = `T cell CD8+`/`Macrophage M2`,
        cd8_m2_ratio = ifelse(is.na(cd8_m2_ratio) | is.infinite(cd8_m2_ratio),
                                0, cd8_m2_ratio))

  cxcl9_infil$variables <- c(globals$imm_estimates,
                             'cd8_treg_ratio',
                             'cd8_m2_ratio')

  cxcl9_infil$titles <- c(globals$imm_estimates,
                          'CD8+/Treg ratio',
                          'CD8+/M2 TAM ratio') %>%
    set_names(cxcl9_infil$variables)

  cxcl9_infil$var_labs <- globals$imm_labels

  ## parallel backend

  plan('multisession')

# serial testing ------

  insert_msg('Serial testing')

  cxcl9_infil$test_results <- cxcl9_infil$analysis_tbl %>%
    map(mutate,
        cxcl9_group = factor(cxcl9_group, c('high', 'low'))) %>% ## proper eff size sign
    future_map(~compare_variables(.x,
                                  variables = cxcl9_infil$variables,
                                  split_factor = 'cxcl9_group',
                                  what = 'eff_size',
                                  types = 'wilcoxon_r',
                                  ci = FALSE,
                                  exact = FALSE,
                                  pub_styled = TRUE,
                                  adj_method = 'BH'),
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

  ## panel plot labels

  cxcl9_infil$panel_labs <- cxcl9_infil$test_results %>%
    map(mutate,
        ribb_lab = cxcl9_infil$var_labs[variable],
        ribb_lab = paste(ribb_lab, significance, sep = ', '))

  cxcl9_infil$panel_labs <- cxcl9_infil$panel_labs %>%
    map(~set_names(.x$ribb_lab,
                   .x$variable))

  ## cohort plot labs

  cxcl9_infil$study_labs <- cxcl9_infil$test_results %>%
    map2_dfr(., names(.), ~mutate(.x, study = .y)) %>%
    mutate(study_lab = paste0(globals$cohorts[study],
                              '\n',
                              eff_size,
                              ', ',
                              significance)) %>%
    dlply(.(variable), function(x) set_names(x$study_lab, x$study))

# serial plotting --------

  insert_msg('Serial plotting')

  cxcl9_infil$plots <- list(cohort = cxcl9_infil$analysis_tbl,
                            stats = cxcl9_infil$test_results) %>%
    pmap(function(cohort, stats) list(variable = cxcl9_infil$variables,
                                      plot_title = cxcl9_infil$titles,
                                      plot_subtitle = stats$plot_cap) %>%
           pmap(plot_variable,
                cohort,
                split_factor = 'cxcl9_group',
                type = 'violin',
                x_lab = 'CXCL9 strata',
                y_lab = 'QuanTIseq estimate',
                point_hjitter = 0,
                cust_theme = globals$common_theme) %>%
           map(~.x +
                 scale_fill_manual(values = c('steelblue', 'firebrick')) +
                 labs(tag = .x$labels$tag %>%
                        stri_replace_all(fixed = '\n', replacement = ', ') %>%
                        paste0('\n', .))) %>%
           set_names(cxcl9_infil$variables))

# Ribbon plots -------

  insert_msg('Ribbon plots')

  ## min/max normalized data

  cxcl9_infil$panel_data <- cxcl9_infil$analysis_tbl %>%
    map(~cbind(.x[, 'cxcl9_group'],
               min_max(.x[, -1]))) %>%
    map(as_tibble)

  ## plots

  cxcl9_infil$panel_plots <- list(data = cxcl9_infil$panel_data,
                                   plot_title = globals$cohorts) %>%
    pmap(draw_violin_panel,
         variables = c('Macrophage M1',
                       'Macrophage M2',
                       'NK cell',
                       'T cell CD4+ (non-regulatory)',
                       'T cell regulatory (Tregs)',
                       'T cell CD8+'),
         split_factor = 'cxcl9_group',
         distr_geom  = 'violin',
         plot_subtitle = 'min/max normalized, Mann-Whitney test',
         x_lab = 'min/max QuanTIseq estimates',
         cust_theme = globals$common_theme,
         scale = 'width',
         dodge_w = 0.9,
         point_hjitter = 0,
         point_wjitter = 0.05,
         point_alpha = 0.3) %>%
    map2(., cxcl9_infil$panel_labs,
         ~.x +
           scale_y_discrete(labels = .y) +
           scale_fill_manual(values = c(low = 'steelblue',
                                        high = 'firebrick'),
                             name = 'CXCL9 strata'))

# plotting normalized CD8, TReg and CD8/Treg ratios for the cohorts ------

  insert_msg('CD8, Treg and CD8/Treg in the cohorts')

  ## data

  cxcl9_infil$cohort_data <- cxcl9_infil$analysis_tbl %>%
    map(~cbind(.x[, 'cxcl9_group'],
               min_max(.x[, -1]))) %>%
    map(as_tibble) %>%
    map2_dfr(., names(.), ~mutate(.x, study = .y)) %>%
    select(study,
           cxcl9_group,
           all_of(cxcl9_infil$variables)) %>%
    pivot_longer(cols = cxcl9_infil$variables,
                 names_to = 'imm_parameter',
                 values_to = 'imm_estimate') %>%
    dlply(.(imm_parameter), as_tibble)

  ## plots

  cxcl9_infil$cohort_panels <-
    list(data = cxcl9_infil$cohort_data,
         plot_title = cxcl9_infil$titles[names(cxcl9_infil$cohort_data)]) %>%
    pmap(function(data, plot_title) ggplot(data,
                                           aes(x = imm_estimate,
                                               y = study,
                                               fill = cxcl9_group)) +
           geom_boxplot(alpha = 0.25,
                        show.legend = FALSE,
                        position = position_dodge(width = 0.9)) +
           geom_point(shape = 21,
                      alpha = 0.3,
                      position = position_jitterdodge(jitter.width = 0.05,
                                                      dodge.width = 0.9)) +
           scale_fill_manual(values = c(low = 'steelblue',
                                        high = 'firebrick'),
                             name = 'CXCL9 strata') +
           labs(title = plot_title,
                subtitle = 'min/max normalized, Mann-Whitney test',
                x = 'min/max QuanTIseq estimates') +
           globals$common_theme +
           theme(axis.title.y = element_blank())) %>%
    map2(., cxcl9_infil$study_labs,
         ~.x +
           scale_y_discrete(limits = c('cm10', 'cm25ni', 'cm25ev', 'tcga',
                                       'reca', 'emtab1980'),
                            labels = .y))

# END ----

  plan('sequential')

  insert_tail()
