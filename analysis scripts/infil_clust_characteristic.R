# Levels of infiltration and expression of CXCL9/10/11 and CCL11
# in the infiltration clusters

  insert_head()

# container list ----

  clust_chara <- list()

# globals: variables and analysis tables ------

  insert_msg('Globals setup')

  ## analysis tables, appending with the cluster information

  clust_chara$analysis_tbl <- list(tcga = tcga$expression %>%
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
                patient_id,
                any_of(globals$clin_vars[globals$clin_vars != 'neoadjuvant']),
                all_of(globals$imm_estimates),
                all_of(globals$cxcl_genes))) %>%
    map2(.,
         map(infil_clust$clust_obj,
             ~set_names(.x$clust_assignment[c('observation', 'clust_id')],
                        c('patient_id', 'clust_id'))),
         right_join, by = 'patient_id')

  ## clinical variables
  ## Checkmate 010: there's only one arm and no rhabdoid differentiation
  ## hence, the variables are removed

  clust_chara$variables <- clust_chara$analysis_tbl %>%
    map(names) %>%
    map(~.x[.x %in% globals$clin_vars])

  clust_chara$variables$cm10 <-
    clust_chara$variables$cm10[!clust_chara$variables$cm10 %in% c('arm', 'rhab_diff')]

  ## immune infiltration variables, chemokine genes

  clust_chara$variables <- clust_chara$variables %>%
    map(~c(.x, globals$imm_estimates, globals$cxcl_genes))

  ## variable labels

  clust_chara$var_labs <- c(globals$var_labs,
                            set_names(globals$imm_estimates,
                                      globals$imm_estimates),
                            set_names(globals$cxcl_genes,
                                      globals$cxcl_genes))

  ## test type

  clust_chara$test_type <- map2(clust_chara$analysis_tbl,
                                clust_chara$variables,
                                ~.x[.y]) %>%
    map(~map_lgl(.x, is.numeric)) %>%
    map(~ifelse(.x, 'kruskal_test', 'chisq_test'))

  ## plot type

  clust_chara$plot_type <- map2(clust_chara$analysis_tbl,
                                clust_chara$variables,
                                ~.x[.y]) %>%
    map(~map_lgl(.x, is.numeric)) %>%
    map(~ifelse(.x, 'violin', 'stack'))

  ## parallel backend

  plan('multisession')

# Descriptive stats -------

  insert_msg('Descriptive statistic')

  clust_chara$desc_stats <- future_map2(clust_chara$analysis_tbl,
                                        clust_chara$variables,
                                        ~explore(.x,
                                                 split_factor = 'clust_id',
                                                 variables = .y,
                                                 what = 'table',
                                                 pub_styled = TRUE),
                                        .options = furrr_options(seed = TRUE)) %>%
    map(reduce, left_join, by = 'variable') %>%
    map(set_names, c('variable', 'clust_1', 'clust_2', 'clust_3'))

# Testing -------

  insert_msg('Testing for differences between the clusters')

  clust_chara$test_results <- list(x = clust_chara$analysis_tbl,
                                   y = clust_chara$variables,
                                   z = clust_chara$test_type) %>%
    future_pmap(function(x, y, z) compare_variables(x,
                                                    variables = y,
                                                    split_factor = 'clust_id',
                                                    what = 'test',
                                                    types = z,
                                                    ci = FALSE,
                                                    adj_method = 'BH',
                                                    pub_styled = TRUE,
                                                    simplify_p = FALSE),
                .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_lab = paste(eff_size, significance, sep = ', '))

# plotting -------

  insert_msg('Plotting')

  clust_chara$plots <- list(cohort = set_names(names(clust_chara$variables),
                                               names(clust_chara$variables)),
                            suffix = globals$cohorts) %>%
    future_pmap(function(cohort, suffix) list(variable = clust_chara$variables[[cohort]],
                                              type = clust_chara$plot_type[[cohort]],
                                              plot_title = paste(clust_chara$var_labs[clust_chara$variables[[cohort]]],
                                                                 suffix, sep = ', '),
                                              plot_subtitle = clust_chara$test_results[[cohort]]$plot_lab) %>%
                  pmap(plot_variable,
                       clust_chara$analysis_tbl[[cohort]],
                       split_factor = 'clust_id',
                       scale = 'percent',
                       cust_theme = globals$common_theme,
                       x_lab = 'infiltration cluster',
                       point_hjitter = 0) %>%
                  map(~.x +
                        labs(tag = .x$labels$tag %>%
                               stri_replace_all(fixed = '\n', replacement = ', ') %>%
                               paste0('\n', .))) %>%
                  set_names(clust_chara$variables[[cohort]]),
                .options = furrr_options(seed = TRUE))

  ## manual adjustment: immune estimate plots, gene expression plots

  for(i in names(clust_chara$plots)) {

    clust_chara$plots[[i]][globals$cxcl_genes] <-
      clust_chara$plots[[i]][globals$cxcl_genes] %>%
      map(~.x +
            theme(plot.title = element_text(face = 'bold.italic')) +
            labs(y = expression('log'[2]*' expression')) +
            scale_fill_manual(values = globals$clust_colors))

    clust_chara$plots[[i]][globals$imm_estimates] <-
      clust_chara$plots[[i]][globals$imm_estimates] %>%
      map(~.x +
            labs(y = 'quanTIseq estimate') +
            scale_fill_manual(values = globals$clust_colors))

  }

# END -----

  plan('sequential')

  rm(i)

  insert_tail()
