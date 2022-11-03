# Clinical cohort characteristic.

  insert_head()

# container list -----

  cohort <- list()

# globals ------

  insert_msg('Globals setup')

  ## variables

  cohort$variables <- c(globals$clin_vars,
                        c('relapse', 'death', 'os_days'))

  ## analysis tables

  cohort$analysis_tbl <- list(tcga = tcga$expression %>%
                                filter(tissue_type == 'Tumor'),
                              emtab1980 = emtab1980$expression,
                              gse73731 = gse73731$expression,
                              gse167093 = gse167093$expression %>%
                                filter(tissue_type == 'Tumor'),
                              reca = reca$expression,
                              cm10 = cm10$expression,
                              cm25ev = cm25ev$expression,
                              cm25ni = cm25ni$expression) %>%
    map(~select(.x, any_of(cohort$variables)))

  ## data set-specific variables

  cohort$variables <- cohort$analysis_tbl %>%
    map(names) %>%
    map(~.x[.x %in% cohort$variables])

  ## manual re-coding: grade and tumor stage

  cohort$analysis_tbl[c('tcga',
                        'emtab1980',
                        'reca',
                        'cm10',
                        'cm25ev',
                        'cm25ni')] <-
    cohort$analysis_tbl[c('tcga',
                          'emtab1980',
                          'reca',
                          'cm10',
                          'cm25ev',
                          'cm25ni')] %>%
    map(mutate,
        death = ifelse(death == 1, 'yes', 'no'),
        death = factor(death),
        relapse = ifelse(relapse == 1, 'yes', 'no'),
        relapse = factor(relapse))

  cohort$analysis_tbl[c('emtab1980',
                        'reca',
                        'gse73731',
                        'gse167093')] <-
    cohort$analysis_tbl[c('emtab1980',
                          'reca',
                          'gse73731',
                          'gse167093')] %>%
    map(mutate,
        tumor_grade = ifelse(!is.na(tumor_grade),
                             paste0('G', as.character(tumor_grade)),
                             NA),
        tumor_grade = factor(tumor_grade))

  cohort$analysis_tbl[c('tcga',
                        'emtab1980',
                        'reca',
                        'gse73731',
                        'gse167093')] <-
    cohort$analysis_tbl[c('tcga',
                          'emtab1980',
                          'reca',
                          'gse73731',
                          'gse167093')] %>%
    map(mutate,
        pt_stage = stri_extract(pt_stage, regex = '^T\\d{1}'),
        pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')))

# analysis -----

  insert_msg('Creating the feature table')

  cohort$ft_table <- list(data = cohort$analysis_tbl,
                          variables = cohort$variables) %>%
    pmap(explore,
         what = 'table',
         pub_styled = TRUE) %>%
    map(mutate,
        variable = globals$var_labs[variable],
        statistic = stri_replace(statistic,
                                 regex = 'no.*\\nyes:\\s{1}',
                                 replacement = ''),
        statistic = stri_replace_all(statistic,
                                     fixed = '% (',
                                     replacement = '% (n = '))

  cohort$ft_table <- cohort$ft_table %>%
    reduce(full_join, by = 'variable') %>%
    set_names(c('Variable',
                globals$cohorts))

# Differences in stage and grade between TCGA, E-MTAB 1980 and EURECA ------

  insert_msg('Stage distribution for TCGA, E-MTAB 1980 and EURECA')

  cohort$comparison_plots <-
    list(variable = c('tumor_grade',
                      'pt_stage'),
         plot_title = globals$var_labs[c('tumor_grade',
                                         'pt_stage')]) %>%
    pmap(plot_variable,
         cohort$analysis_tbl$tcga,
         cohort$analysis_tbl$emtab1980,
         cohort$analysis_tbl$reca,
         cohort$analysis_tbl$gse73731,
         cohort$analysis_tbl$gse167093,
         data_names = globals$cohorts[c('tcga', 'emtab1980', 'reca',
                                        'gse73731', 'gse167093')],
         type = 'stack',
         scale = 'percent',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          labs(tag = .x$labels$tag %>%
                 stri_replace_all(fixed = '\n',
                                  replacement = ', ') %>%
                 paste0('\n', .)) +
          scale_fill_viridis_d()) %>%
    set_names(c('tumor_grade',
                'pt_stage'))

# END -----

  insert_tail()
