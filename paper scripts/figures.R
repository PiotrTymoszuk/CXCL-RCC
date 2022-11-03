# Main paper figures

  insert_head()

# container ------

  figs <- list()

# Figure 1: CXCL9 in normal kidney and RCC -------

  insert_msg('Figure 1: CXCL9 in the normal kidney and RCC')

  figs$norm_tumor <- norm_tumor$plots %>%
    map(~.x$CXCL9) %>%
    map2(., globals$cohorts[names(norm_tumor$plots)],
         ~.x +
           labs(title = .y,
                y = expression(italic('CXCL9') * ', log'[2] * ' expression')) +
           theme(legend.position = 'none',
                 plot.title = element_text(face = 'bold'))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv') %>%
    as_figure(label = 'figure_1_norm_tumor_CXCL9',
              ref_name = 'norm_tumor',
              caption = 'Expression of CXCL9 in the normal kindney and RCC.',
              w = 120,
              h = 90)

# Figure 2: correlation of CXCL9 with infiltration --------

  insert_msg('Figure 2: correlation of CXCL9 with infiltration')

  figs$infil <- infil_cor$CXCl9_bubble %>%
    as_figure(label = 'figure_2_infiltration_CXCL9',
              ref_name = 'infil',
              caption = 'Correlation of CXCL9 expression in RCC with immune cell infiltration.',
              w = 180,
              h = 150)

# Figure 4: GO enrichment, common top terms -------

  insert_msg('Figure 4: GO enrichment')

  figs$go <- enrich$go_cxcl9_cmm_plots[c('tcga',
                                         'cm25ev',
                                         'cm25ni',
                                         'gse73731',
                                         'gse167093',
                                         'emtab1980',
                                         'reca')] %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.text.y = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv') %>%
    as_figure(label = 'figure_4_go_cxcl9',
              ref_name = 'go',
              caption = 'Biological process GO term enrichment within the genes differentially regulated between CXCL9 high and CXCL9 low tumors.',
              w = 180,
              h = 220)

# Figure 5: signaling, common regulated pathways -----

  insert_msg('Figure 5: common regulated pathways')

  figs$spia <- plot_grid(pathways$bubble_plot +
                           guides(size = 'none') +
                           theme(legend.position = 'bottom')) %>%
    as_figure(label = 'figure_5_signaling_spia',
              ref_name = 'spia',
              caption = 'Predicted modulation of signaling pathways in CXCL9 high tumors.',
              w = 180,
              h = 190)

# Figure 6: regulation of FAOx -------

  insert_msg('Figure 6: regulation of FAOx')

  figs$faox <- meta_reg$faox$plots[c('tcga',
                                     'cm10',
                                     'cm25ev',
                                     'cm25ni',
                                     'gse73731',
                                     'gse167093',
                                     'emtab1980',
                                     'reca')] %>%
    map(~.x + theme(legend.position = 'none')) %>%
    c(list(get_legend(meta_reg$faox$plots[[1]]))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_6_faox',
              ref_name = 'faox',
              caption = 'Predicted regulation of fatty acid oxidation reactions in CXCL9 high tumors.',
              w = 180,
              h = 190)

# Figure 7: regulation of OxPhos ------

  insert_msg('Figure 7: Regulation of OxPhos')

  figs$oxphos <- meta_reg$oxphos$plots[c('tcga',
                                         'cm10',
                                         'cm25ev',
                                         'cm25ni',
                                         'gse73731',
                                         'gse167093',
                                         'emtab1980',
                                         'reca')] %>%
    map(~.x + theme(legend.position = 'none')) %>%
    c(list(get_legend(meta_reg$oxphos$plots[[1]]))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_7_oxphos',
              ref_name = 'oxphos',
              caption = 'Predicted regulation of oxidative phosphorylation reactions in CXCL9 high tumors.',
              w = 180,
              h = 190)

# Figure 8: survival, OS -------

  insert_msg('Figure 8: survival')

  figs$os <- surv_cut$cxcl9_plots[c('tcga_os',
                                    'cm10_os',
                                    'cm25ev_os',
                                    'cm25ni_os',
                                    'emtab1980_os',
                                    'reca_os')] %>%
    map(~.x +
          labs(subtitle = .x$labels$subtitle %>%
                 stri_replace(regex = ',\\s{1}p\\s{1}raw.*$',
                              replacement = '') %>%
                 stri_replace(fixed = 'FDR', replacement = '') %>%
                 paste(.x$labels$tag, ., sep = '\n')) +
          theme(legend.position = 'none',
                plot.tag = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(surv_cut$cxcl9_plots$tcga_os +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_8_cxcl_os',
              ref_name = 'os',
              caption = 'Overall survival of patients bearing CXCL9 high and CXCL low RCC.',
              w = 180,
              h = 150)

# Saving the figures ------

  insert_msg('Saving the figures')

  figs %>%
    walk(pickle,
         path = './paper/figures',
         format = 'pdf',
         device = cairo_pdf)

# END -----

  insert_tail()
