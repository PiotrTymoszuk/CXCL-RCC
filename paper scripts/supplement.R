# Supplementary figures

  insert_head()

# container -----

  suppl <- list()

# Figure S1: CXCL10, CXCL11 and CXCR3 in the normal kidney and RCC -----

  insert_msg('Figure S1: CXCL10/11 and CXCR3, normal - tumor')

  suppl$norm_tumor <-
    list(cxcl10 = map(norm_tumor$plots, ~.x$CXCL10),
         cxcl11 = map(norm_tumor$plots, ~.x$CXCL11),
         cxcl11 = map(norm_tumor$plots, ~.x$CXCR3)) %>%
    map(~map2(.x, globals$cohorts[c('tcga', 'gse167093')],
              ~.x +
                labs(title = .y) +
                theme(plot.title = element_text(face = 'bold'),
                      plot.title.position = 'plot',
                      legend.position = 'none'))) %>%
    map2(., list(expression(italic('CXCL10') * ', log'[2] * ' expression'),
                 expression(italic('CXCL11') * ', log'[2] * ' expression'),
                 expression(italic('CXCR3') * ', log'[2] * ' expression')),
         function(plot, y_lab) map(plot, ~.x + labs(y = y_lab))) %>%
    unlist(recursive = FALSE) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              labels = c('A', '', 'B', '', 'C'),
              label_size = 10) %>%
    as_figure(label = 'figure_s1_norm_tumor',
              ref_name = 'norm_tumor',
              caption = 'Expression of CXCL10, CXCL11 and CXCR3 in the normal kidney and RCC.',
              w = 180,
              h = 160)

# Figure S2: co-regulation of CXCL9/10/11 and CXCR3 -----

  insert_msg('Co-regulation of CXCL9/10/11')

  suppl$corr_exp <- correl$summ_plots[c('tcga',
                                        'cm10',
                                        'cm25ev',
                                        'cm25ni',
                                        'gse73731',
                                        'gse167093',
                                        'emtab1980',
                                        'reca')] %>%
    map(~.x +
          labs(subtitle = stri_replace(.x$labels$subtitle,
                                       fixed = 'Spearman correlation, ',
                                       replacement = '')) +
          scale_x_discrete(limits = c('CXCL9', 'CXCL10', 'CXCL11', 'CXCR3')) +
          scale_y_discrete(limits = c('CXCL9', 'CXCL10', 'CXCL11', 'CXCR3')) +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(correl$summ_plots[[1]] +
                           guides(size = 'none')),
              ncol = 2,
              rel_widths = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_s2_coexpression',
              ref_name = 'corr_exp',
              caption = 'Co-regulation of CXCL9, CXCL10, CXCL11 and CXCR3 expression in RCC.',
              w = 180,
              h = 220)

# Figure S3: correlation of CXCL9 and CD8+ T cells ------

  insert_msg('Figure S3: correlation of CXCL9 and T cells')

  suppl$cd8 <- infil_cor$point_plots$CXCl9[c('tcga',
                                             'cm10',
                                             'cm25ev',
                                             'cm25ni',
                                             'gse73731',
                                             'gse167093',
                                             'emtab1980',
                                             'reca')] %>%
    map(~.x + theme(plot.title.position = 'plot')) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv') %>%
    as_figure(label = 'figure_s3_cxcl9_cd8',
              ref_name = 'cd8',
              caption = 'Correlation of CXCL9 expression in tumor tissue with predicted infiltration by CD8+ T cells.',
              w = 180,
              h = 190)

# Figure S4: correlation of CXCR3 with QuanTIseq infiltration data ------

  insert_msg('Figure S4: QuanTI seq and CXCR3')

  suppl$infil_cxcr3 <- infil_cor$CXCR3_bubble %>%
    as_figure(label = 'figure_s4_infiltration_CXCR3',
              ref_name = 'infil',
              caption = 'Correlation of CXCR3 expression in RCC with immune cell infiltration.',
              w = 180,
              h = 150)

# Figure S5: correlation of CXCR3 and CD8+ T cells ------

  insert_msg('Figure S5: correlation of CXCR3 and T cells')

  suppl$cd8_cxcr3 <- infil_cor$point_plots$CXCR3[c('tcga',
                                             'cm10',
                                             'cm25ev',
                                             'cm25ni',
                                             'gse73731',
                                             'gse167093',
                                             'emtab1980',
                                             'reca')] %>%
    map(~.x + theme(plot.title.position = 'plot')) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv') %>%
    as_figure(label = 'figure_s5_cxcr3_cd8',
              ref_name = 'cd8_cxcr3',
              caption = 'Correlation of CXCR3 expression in tumor tissue with predicted infiltration by CD8+ T cells.',
              w = 180,
              h = 190)

# Figures S6 - S9: signatures in CXCL9 expression strata -------

  insert_msg('Figure S6 - S9: signatures in the CXCL9 high RCC')

  suppl[c('cd8_sign',
          'cytotox_sign',
          'exhaust_sign',
          'treg_sign')] <-
    list(dge_sig$est_plots[c('gsva_cd8.tcga',
                             'gsva_cd8.cm10',
                             'gsva_cd8.cm25ev',
                             'gsva_cd8.cm25ni',
                             'gsva_cd8.gse73731',
                             'gsva_cd8.gse167093',
                             'gsva_cd8.emtab1980',
                             'gsva_cd8.reca')],
         dge_sig$est_plots[c('gsva_cytotoxicity.tcga',
                             'gsva_cytotoxicity.cm10',
                             'gsva_cytotoxicity.cm25ev',
                             'gsva_cytotoxicity.cm25ni',
                             'gsva_cytotoxicity.gse73731',
                             'gsva_cytotoxicity.gse167093',
                             'gsva_cytotoxicity.emtab1980',
                             'gsva_cytotoxicity.emtab1980')],
         dge_sig$est_plots[c('gsva_exhaustion.tcga',
                             'gsva_exhaustion.cm10',
                             'gsva_exhaustion.cm25ev',
                             'gsva_exhaustion.cm25ni',
                             'gsva_exhaustion.gse73731',
                             'gsva_exhaustion.gse167093',
                             'gsva_exhaustion.emtab1980',
                             'gsva_exhaustion.reca')],
         dge_sig$est_plots[c('gsva_treg.tcga',
                             'gsva_treg.cm10',
                             'gsva_treg.cm25ev',
                             'gsva_treg.cm25ni',
                             'gsva_treg.gse73731',
                             'gsva_treg.gse167093',
                             'gsva_treg.emtab1980',
                             'gsva_treg.reca')]) %>%
    map(~map(.x,
             ~.x +
               theme(legend.position = 'none',
                     plot.title.position = 'plot') +
               labs(x = expression(Delta * ' GSVA score, ' * italic('CXCL9') * ' high vs low')))) %>%
    map(~c(.x,
           list(get_legend(dge_sig$est_plots[[1]] +
                             labs(fill = 'Regulation status'))))) %>%
    map(~plot_grid(plotlist = .x,
                   ncol = 3,
                   align = 'hv',
                   axis = 'tblr'))

  suppl[c('cd8_sign',
          'cytotox_sign',
          'exhaust_sign',
          'treg_sign')] <- suppl[c('cd8_sign',
                                   'cytotox_sign',
                                   'exhaust_sign',
                                   'treg_sign')] %>%
    list(x = .,
         label = c('figure_s6_cd8_signatures',
                   'figure_s7_cytotox_signatures',
                   'figure_s8_exhaust_signatures',
                   'figure_s9_treg_signatures'),
         ref_name = c('cd8_sign',
                      'cytotox_sign',
                      'exhaust_sign',
                      'treg_sign'),
         caption = c('Gene set variation analysis for CD8+ T cell signatures in CXCL9 high and CXCL9 low RCC.',
                     'Gene set variation analysis for cytotoxicity signatures in CXCL9 high and CXCL9 low RCC.',
                     'Gene set variation analysis for T cell exhaustion signatures in CXCL9 high and CXCL9 low RCC.',
                     'Gene set variation analysis for regulatory T cell signatures in CXCL9 high and CXCL9 low RCC.')) %>%
    pmap(as_figure,
         w = 180,
         h = 200)

# Figures S10 - S112: pathview diagrams -------

  insert_msg('Figures S10 - S12: pathview diagrams')

  suppl[c('pathview_cytotox',
          'pathview_chemo',
          'pathview_cyto')] <-
    c('./report/pathview/hsa04650.pathview.png',
      './report/pathview/hsa04062.pathview.png',
      './report/pathview/hsa04060.pathview.png') %>%
    map(~plot_grid(ggdraw() +
                     draw_image(.x)))

  suppl[c('pathview_cytotox',
          'pathview_chemo',
          'pathview_cyto')] <- suppl[c('pathview_cytotox',
                                       'pathview_chemo',
                                       'pathview_cyto')] %>%
    list(x = .,
         label = c('figure_s10_pathview_cytotoxicity',
                   'figure_s11_pathview_chemokine',
                   'figure_s12_pathview_cytokine'),
         ref_name = c('pathview_cytotox',
                      'pathview_chemo',
                      'pathview_cyto'),
         caption = c('Pathview diagram of modulation of the KEGG NK-cell mediate cytotoxicity pathway for CXCL9 high RCC.',
                     'Pathview diagram of modulation of the KEGG chemokine signaling pathway for CXCL9 high RCC.',
                     'Pathview diagram of modulation of the KEGG cytokine-cytokine receptor interaction pathway for CXCL9 high RCC.'),
         h = c(1025/1295,
               927/1246,
               1157/1795) * 180) %>%
    pmap(as_figure,
         w = 180)

# Figure S13: metabolic pathways, pooled estimates --------

  insert_msg('Figure S13: pooled estmates, metabolic reations')

  suppl$meta_pool <- plot_grid(meta_pool$cmm_react_forest +
                                 theme(axis.text.y = element_text(size = 5))) %>%
    as_figure(label = 'figure_s13_common_metabolic_reactions',
              ref_name = 'meta_pool',
              caption = 'Prediction of metabolic reactions modulation based on pooled estimates of differential gene expression in CXCL9 high tumors.',
              w = 180,
              h = 230)

# Figure S14: predicted regulation of TCA -----

  insert_msg('Figure S14: predicted regulation of TCA')

  suppl$tca <- plot_grid(meta_reg$tca$plot +
                           theme(legend.position = 'bottom')) %>%
    as_figure(label = 'figure_s14_tca',
              ref_name = 'tca',
              caption = 'Predicted regulation of TCA reactions in CXCL9 high tumors.',
              w = 180,
              h = 200)

# Figure S15: CXCL9 and tumor grade -----

  insert_msg('Figure S15: CXCL9 and tumor grade')

  suppl$grade <- exp_clinic$plots[c('tcga',
                                    'gse73731',
                                    'gse167093',
                                    'emtab1980',
                                    'reca')] %>%
    map(~.x$CXCL9$tumor_grade) %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title = element_text(face = 'bold')) +
          labs(title = .x$labels$title %>%
                 stri_replace(fixed = 'CXCL9, ',
                              replacement = ''),
               y = expression(italic('CXCL9') * ', log'[2] * ' expression'))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_s15_tumor_grade',
              ref_name = 'grade',
              caption = 'Tumor grade and CXCL9 expression.',
              w = 180,
              h = 220)

# Figure S16: CXCL9 and tumor stage -------

  insert_msg('Figure S16: CXCL9 and tumor stage')

  suppl$stage <- exp_clinic$plots[c('tcga',
                                    'gse73731',
                                    'gse167093',
                                    'emtab1980',
                                    'reca')] %>%
    map(~.x$CXCL9$pt_stage) %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title = element_text(face = 'bold')) +
          labs(title = .x$labels$title %>%
                 stri_replace(fixed = 'CXCL9, ',
                              replacement = ''),
               y = expression(italic('CXCL9') * ', log'[2] * ' expression'))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_s16_tumor_stage',
              ref_name = 'stage',
              caption = 'Tumor stage and CXCL9 expression.',
              w = 180,
              h = 220)

# Figure S17: CXCL9 and metastasis status --------

  insert_msg('Figure S17: CXCL9 and distant metastases')

  suppl$metastases <- exp_clinic$plots[c('tcga',
                                         'emtab1980',
                                         'reca')] %>%
    map(~.x$CXCL9$pm_stage) %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title = element_text(face = 'bold')) +
          labs(title = .x$labels$title %>%
                 stri_replace(fixed = 'CXCL9, ',
                              replacement = ''),
               y = expression(italic('CXCL9') * ', log'[2] * ' expression'))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_s17_metastases',
              ref_name = 'metastases',
              caption  = 'Distant metastases at RCC diagnosis and CXCL9 expression.',
              w = 180,
              h = 90)

# Figure S18: MSKCC and CXCL9 expression -------

  insert_msg('Figure S18: MSKCC and CXCL9')

  suppl$mskcc <- exp_clinic$plots[c('cm10',
                                    'cm25ev',
                                    'cm25ni')] %>%
    map(~.x$CXCL9$MSKCC_risk) %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title = element_text(face = 'bold')) +
          labs(title = .x$labels$title %>%
                 stri_replace(fixed = 'CXCL9, ',
                              replacement = ''),
               y = expression(italic('CXCL9') * ', log'[2] * ' expression'))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_s18_mskcc',
              ref_name = 'mskcc',
              caption = 'MSKCC risk classification and CXCL9 expression.',
              w = 180,
              h = 160)

# Figure S19: CXCL9 and relapse-free survival ------

  insert_msg('Figure S19: CXCL9 and relapse-free survival')

  suppl$rfs <- surv_cut$cxcl9_plots[c('tcga_rfs',
                                      'cm10_rfs',
                                      'cm25ev_rfs',
                                      'cm25ni_rfs')] %>%
    map(~.x +
          labs(subtitle = .x$labels$subtitle %>%
                 stri_replace(regex = ',\\s{1}p\\s{1}raw.*$',
                              replacement = '') %>%
                 stri_replace(fixed = 'FDR', replacement = '') %>%
                 paste(.x$labels$tag, ., sep = '\n')) +
          theme(legend.position = 'none',
                plot.tag = element_blank())) %>%
    c(list(get_legend(surv_cut$cxcl9_plots$tcga_rfs))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_s19_rfs',
              ref_name = 'rfs',
              caption = 'Relapse-free survival of patients bearing CXCL9 high and CXCL9 low RCC.',
              w = 180,
              h = 150)

# Figure S20: CXCL9 expression and response ------

  insert_msg('Figure S20: CXCL9 and therapy response')

  suppl$biresp <-
    exp_clinic$plots[c('cm10',
                       'cm25ev',
                       'cm25ni')] %>%
    map(~.x$CXCL9$biresponse) %>%
    map2(., globals$cohorts[c('cm10',
                              'cm25ev',
                              'cm25ni')],
         ~.x +
           theme(axis.title.x = element_blank(),
                 plot.title = element_text(face = 'bold'),
                 legend.position = 'none') +
           labs(title = paste('Overall response', .y, sep = ', '),
                y = expression(italic('CXCL9') * ', log'[2] * ' expression'))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv') %>%
    as_figure(label = 'figure_s20_clinical response',
              ref_name = 'biresp',
              caption = 'Expression of CXCL9 in patients with and without clinical reponse to everolimus and nivolumab therapy.',
              w = 180,
              h = 90)

# Saving the figures -----

  insert_msg('Saving the figures')

  suppl %>%
    walk(pickle,
         path = './paper/supplementary figures',
         format = 'pdf',
         device = cairo_pdf)

# END ------

  insert_tail()
