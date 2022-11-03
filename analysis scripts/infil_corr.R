# Correlation of the CXCL9/10/11 and CXCL11 with immune cell infiltration.
# Spearman correlation.

  insert_head()

# container list ----

  infil_cor <- list()

# globals -----

  insert_msg('Globals setup')

  ## analysis table, including the CD8:Treg ratio

  infil_cor$analysis_tbl <- list(tcga = tcga$expression %>%
                                   filter(tissue_type == 'Tumor'),
                                 emtab1980 = emtab1980$expression,
                                 gse73731 = gse73731$expression,
                                 gse167093 = gse167093$expression %>%
                                   filter(tissue_type == 'Tumor'),
                                 reca = reca$expression,
                                 cm10 = cm10$expression,
                                 cm25ev = cm25ev$expression,
                                 cm25ni = cm25ni$expression) %>%
    map(~.x[c(globals$cxcl_genes,
              globals$imm_estimates)]) %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate,
        cd8_treg_ratio = `T cell CD8+`/`T cell regulatory (Tregs)`,
        cd8_treg_ratio = ifelse(is.nan(cd8_treg_ratio) | is.infinite(cd8_treg_ratio),
                                0,
                                cd8_treg_ratio),
        cd8_m2_ratio = `T cell CD8+`/`Macrophage M2`)

  ## feature pairs

  infil_cor$pairs <- globals$cxcl_genes %>%
    map(function(cxc) map(c(globals$imm_estimates,
                            'cd8_treg_ratio',
                            'cd8_m2_ratio'),
                          ~c(cxc, .x)) %>%
          set_names(c(globals$imm_estimates,
                      'cd8_treg_ratio',
                      'cd8_m2_ratio'))) %>%
    set_names(globals$cxcl_genes)

# correlations ----

  insert_msg('Serial correlation')

  infil_cor$test_results <- set_names(names(infil_cor$analysis_tbl),
                                      names(infil_cor$analysis_tbl)) %>%
    map(function(cohort) infil_cor$pairs %>%
          map(function(cxc) map_dfr(cxc,
                                    ~correlate_variables(infil_cor$analysis_tbl[[cohort]],
                                                         variables = .x,
                                                         what = 'correlation',
                                                         type = 'spearman',
                                                         pub_styled = FALSE,
                                                         simplify_p = FALSE))) %>%
          map_dfr(mutate, p_adjusted = p.adjust(p_value, 'BH')))

  infil_cor$test_results <- infil_cor$test_results %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 'bold', 'plain'),
        plot_lab = signif(estimate, 2))

# Summary bubble plots -----

  insert_msg('Summary plots')

  infil_cor$bubble <-
    map2(infil_cor$test_results,
         globals$cohorts,
         ~ggplot(.x,
                 aes(x = variable1,
                     y = variable2,
                     size = abs(estimate),
                     fill = estimate)) +
           geom_point(shape = 21) +
           geom_text(aes(label = plot_lab,
                         fontface = significant),
                     size = 2.75,
                     hjust = 0.5,
                     vjust = -1.5) +
           scale_y_discrete(labels = c(globals$imm_labels,
                                       c(cd8_treg_ratio = 'CD8+/Treg',
                                         cd8_m2_ratio = 'CD8+/M2 TAM'))) +
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'white',
                                high = 'firebrick',
                                midpoint = 0,
                                limits = c(-1, 1),
                                name = expression(rho)) +
           scale_size_area(limits = c(0, 1),
                           max_size = 6,
                           name = expression('abs('*rho*')')) +
           scale_x_discrete(limits = c('CXCL9',
                                       'CXCL10',
                                       'CXCL11',
                                       'CXCR3',
                                       'CCL11')) +
           globals$common_theme +
           theme(axis.title = element_blank(),
                 axis.text.x = element_text(face = 'italic')) +
           labs(title = .y,
                subtitle = paste('Spearman correlation, n =',
                                 .x$n[[1]])))

# Correlation between the CXCL9/10/11/CCL11 and the IF count of CD8 T cells -----

  insert_msg('Correlation with CD8+ IF count')

  ## data

  infil_cor$cd8_data <- list(cm10 = cm10$expression,
                             cm25ev = cm25ev$expression,
                             cm25ni = cm25ni$expression) %>%
    map(select, TM_CD8_Density, all_of(globals$cxcl_genes)) %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate, TM_CD8_Density = log2(TM_CD8_Density))

  ## Spearman correlations

  infil_cor$cd8_test_results <- infil_cor$cd8_data %>%
    map(function(cohort) globals$cxcl_genes %>%
          map_dfr(~correlate_variables(cohort,
                                   variables = c(.x, 'TM_CD8_Density'),
                                   what = 'correlation',
                                   type = 'spearman',
                                   pub_styled = TRUE,
                                   simplify_p = FALSE)) %>%
          mutate(p_adjusted = p.adjust(p_value, 'BH'),
                 significance = ifelse(p_adjusted < 0.05,
                                       paste('p =',
                                             signif(p_adjusted, 2)),
                                       paste0('ns (p = ',
                                              signif(p_adjusted, 2), ')')),
                 eff_size = stri_replace(eff_size, fixed = 'rho',
                                         replacement = '\u03C1'),
                 plot_cap = paste(eff_size, significance, sep = ', '),
                 plot_cap = paste0(plot_cap, ', n = ', n)))

  ## correlation plots

  infil_cor$cd8_plots <- list(x = infil_cor$cd8_data,
                              y = globals$cohorts[names(infil_cor$cd8_data)],
                              z = infil_cor$cd8_test_results ) %>%
    pmap(function(x, y, z) list(variables = map(globals$cxcl_genes,
                                                ~c(.x, 'TM_CD8_Density')),
                                plot_subtitle = z$plot_cap,
                                plot_title = paste(globals$cxcl_genes, y, sep = ', ')) %>%
           pmap(plot_correlation,
                data = x,
                x_lab = expression('log'[2] * ' expression'),
                y_lab = expression('log'[2] * ' CD8'^{'+'} * ' density, ID'),
                cust_theme = globals$common_theme,
                show_trend = FALSE) %>%
           map(~.x +
                 geom_smooth(method = 'gam')) %>%
           set_names(globals$cxcl_genes))

# Forest plots for the key immune populations: CD8, Treg, TAMS and the ratios -----

  insert_msg('Forest plots with the correlations')

  ## data

  infil_cor$forest_data <- infil_cor$test_results %>%
    map2_dfr(., names(.), ~mutate(.x, study = .y)) %>%
    filter(variable2 %in% c('T cell CD8+', 'T cell regulatory (Tregs)',
                            'Macrophage M1', 'Macrophage M2',
                            'cd8_treg_ratio', 'cd8_m2_ratio'),
           variable1 == 'CXCL9') %>%
    mutate(variable2 = car::recode(variable2,
                                   "'cd8_treg_ratio' = 'CD8+/Treg';
                                   'cd8_m2_ratio' = 'CD8+/M2 TAM'")) %>%
    dlply(.(variable2))

  ## forest plots

  infil_cor$forest_plots <- list(data = infil_cor$forest_data,
                                 plot_title = names(infil_cor$forest_data)) %>%
    pmap(plot_top,
         regulation_variable = 'estimate',
         label_variable = 'study',
         p_variable = 'p_adjusted',
         lower_ci_variable = 'lower_ci',
         upper_ci_variable = 'upper_ci',
         plot_subtitle = 'Spearman correlation, CXCL9 expression',
         x_lab = expression(rho * ' \u00B1 95% CI'),
         cust_theme = globals$common_theme,
         show_txt = TRUE) %>%
    map(~.x +
          scale_color_manual(values = c(upregulated = 'firebrick',
                                        downregulated = 'steelblue',
                                        ns = 'gray60'),
                             labels = c(upregulated = 'positive',
                                        downregulated = 'negative',
                                        ns = 'ns'),
                             name = 'Correlation') +
          guides(fill = 'none') +
          scale_y_discrete(labels = globals$cohorts))

# Bubble plots for CXCL9 and CXCR3 only -------

  insert_msg('Bubble-correlograms for CXCL9 and CXCR3 only')

  ## plot data

  infil_cor[c('CXCl9_data',
              'CXCR3_data')] <- c('CXCL9', 'CXCR3') %>%
    map(function(gene) map(infil_cor$test_results,
                           ~filter(.x, variable1 == gene)) %>%
          map2_dfr(., names(.),
                   ~mutate(.x, cohort = .y)) %>%
          mutate(signif_alpha = ifelse(significant == 'plain', 'ns', 'signif')))

  ## labels

  infil_cor$CXCl9_n_numbers <-
    map2_chr(globals$cohorts,
             infil_cor$analysis_tbl[names(globals$cohorts)],
             ~paste(.x, nrow(.y), sep = '\nn = '))

  ## bubble plots

  infil_cor[c('CXCl9_bubble',
              'CXCR3_bubble')] <-
    map2(infil_cor[c('CXCl9_data',
                     'CXCR3_data')],
         list(expression(bold(italic('CXCL9') * ' expression and immune infiltration')),
              expression(bold(italic('CXCR3') * ' expression and immune infiltration'))),
         ~ggplot(.x,
                 aes(x = cohort,
                     y = variable2,
                     size = abs(estimate),
                     fill = estimate)) +
           geom_point(shape = 21) +
           geom_text(aes(label = plot_lab,
                         fontface = significant,
                         alpha = signif_alpha),
                     size = 2.5,
                     hjust = 0.5,
                     vjust = -1.5) +
           scale_x_discrete(limits = c('tcga',
                                       'cm10',
                                       'cm25ev',
                                       'cm25ni',
                                       'gse73731',
                                       'gse167093',
                                       'emtab1980',
                                       'reca'),
                            labels =   infil_cor$CXCl9_n_numbers) +
           scale_y_discrete(limits = rev(globals$imm_estimates),
                            labels = globals$imm_labels) +
           scale_alpha_manual(values = c(ns = 0.5,
                                         signif = 1)) +
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'white',
                                high = 'firebrick',
                                midpoint = 0,
                                limits = c(-1, 1)) +
           guides(fill = 'none',
                  size = 'none',
                  alpha = 'none') +
           globals$common_theme +
           theme(axis.title = element_blank()) +
           labs(title = .y,
                subtitle = 'Spearman correlation with QuanTIseq infiltration estimates'))

# Point plots for the CD8 and CXCL9/CXCR3 correlations -------

  insert_msg('Correlation point plots')

  ## plot captions

  infil_cor$point_captions <- infil_cor[c('CXCl9_data',
                                          'CXCR3_data')] %>%
    map(filter, variable2 == 'T cell CD8+') %>%
    map(mutate,
        plot_cap = paste0('\u03C1 = ',
                          signif(estimate, 2),
                          ' [', signif(lower_ci, 2),
                          ' - ', signif(upper_ci, 2), ']'),
        plot_cap = paste0(plot_cap, ', p = ', signif(p_adjusted, 2)),
        plot_cap = paste0(plot_cap, ', n = ', n)) %>%
    set_names(c('CXCl9', 'CXCR3'))

  ## point plots

  infil_cor$point_plots$CXCl9 <-
    list(data = infil_cor$analysis_tbl[infil_cor$point_captions$CXCl9$cohort],
         plot_title = globals$cohorts[infil_cor$point_captions$CXCl9$cohort],
         plot_subtitle = infil_cor$point_captions$CXCl9$plot_cap) %>%
    pmap(plot_correlation,
         variables = c('CXCL9', 'T cell CD8+'),
         type = 'correlation',
         point_hjitter = 0,
         x_lab = expression(italic('CXCL9') * ', log'[2] * ' expression'),
         y_lab = expression('CD8'^'+' * ' T cells, QuanTIseq'),
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_y_continuous(trans = 'log2',
                             labels = function(x) signif(x, 2)) +
          theme(plot.tag = element_blank()))

  infil_cor$point_plots$CXCR3 <-
    list(data = infil_cor$analysis_tbl[infil_cor$point_captions$CXCR3$cohort],
         plot_title = globals$cohorts[infil_cor$point_captions$CXCR3$cohort],
         plot_subtitle = infil_cor$point_captions$CXCR3$plot_cap) %>%
    pmap(plot_correlation,
         variables = c('CXCR3', 'T cell CD8+'),
         type = 'correlation',
         point_hjitter = 0,
         x_lab = expression(italic('CXCR3') * ', log'[2] * ' expression'),
         y_lab = expression('CD8'^'+' * ' T cells, QuanTIseq'),
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_y_continuous(trans = 'log2',
                             labels = function(x) signif(x, 2)) +
          theme(plot.tag = element_blank()))

# END ------

  insert_tail()
