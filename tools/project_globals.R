# This script contains project globals

# libraries ----

  library(plyr)
  library(tidyverse)
  library(stringi)

# data container ------

  globals <- list()

# graphics -----

  globals$common_text <- element_text(size = 8,
                                      face = 'plain',
                                      color = 'black')

  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = 'mm')

  globals$common_theme <- theme_classic() + theme(axis.text = globals$common_text,
                                                  axis.title = globals$common_text,
                                                  plot.title = element_text(size = 8,
                                                                            face = 'bold'),
                                                  plot.subtitle = globals$common_text,
                                                  plot.tag = element_text(size = 8,
                                                                          face = 'plain',
                                                                          color = 'black',
                                                                          hjust = 0,
                                                                          vjust = 1),
                                                  plot.tag.position = 'bottom',
                                                  legend.text = globals$common_text,
                                                  legend.title = globals$common_text,
                                                  strip.text = globals$common_text,
                                                  strip.background = element_rect(fill = 'gray95',
                                                                                  color = 'gray80'),
                                                  plot.margin = globals$common_margin,
                                                  panel.grid.major = element_line(color = 'gray90'))

# genes of interest ------

  globals$cxcl_genes <- c('CXCL9', 'CXCL10', 'CXCL11', 'CCL11', 'CXCR3')

# clinical variables -----

  globals$clin_vars <- c('sex', 'age', 'race', 'neoadjuvant',
                        'laterality', 'tumor_grade', 'p_stage',
                        'pt_stage', 'pm_stage', 'pn_stage',
                        'MSKCC_risk', 'IMDC_risk',
                        'sarc_diff', 'rhab_diff',
                        'tumor_delta_vol', 'response_type',
                        'benefit')

  globals$surv_vars <- c('relapse', 'rfs_days',
                         'death', 'os_days',
                         'tumor_death')

  globals$var_labs <- c(globals$var_labs,
                        set_names(c('Sex', 'Age, years', 'Race',
                                    'Neoadjuvant treatment', 'Laterality',
                                    'Tumor grade', 'Pathological stage',
                                    'Tumor stage, pt', 'Metastasis stage, pm',
                                    'Node stage, pn',
                                    'MSKCC risk group',
                                    'IMDC risk group',
                                    'Sarcomatoid differentiation',
                                    'Rhabdoid differentiation',
                                    'Tumor shrinkage',
                                    'Therapy response',
                                    'Clinical benefit'),
                                  globals$clin_vars),
                        set_names(c('Relapse', 'RFS, days',
                                    'Death', 'OS, days', 'Tumor-related death'),
                                  globals$surv_vars),
                        set_names(c('Overall response',
                                    'biresponse')))

# infiltration estimates, cell types and labels -----

  globals$imm_estimates <- c('T cell CD4+ (non-regulatory)',
                             'T cell CD8+',
                             'T cell regulatory (Tregs)',
                             'NK cell',
                             'B cell',
                             'Macrophage M1',
                             'Macrophage M2',
                             'Monocyte',
                             'Neutrophil',
                             'Myeloid dendritic cell',
                             'uncharacterized cell')

  globals$imm_labels <-
    c('T cell CD4+ (non-regulatory)' = 'CD4+ T',
      'T cell CD8+' = 'CD8+ T',
      'T cell regulatory (Tregs)' = 'Treg',
      'NK cell' = 'NK',
      'B cell' = 'B',
      'Macrophage M1' = 'M1 TAM',
      'Macrophage M2' = 'M2 TAM',
      'Monocyte' = 'Mono',
      'Neutrophil' = 'Neutro',
      'Myeloid dendritic cell' = 'mDC',
      'uncharacterized cell' = 'rest',
      'cd8_treg_ratio' = 'CD8+ T : Treg',
      'cd8_m2_ratio' = 'CD8+ T : M2 TAM')

# study cohorts ------

  globals$cohorts <-c(tcga = 'TCGA',
                      emtab1980 = 'E-MTAB 1980',
                      gse73731 = 'GSE73731',
                      gse167093 = 'GSE167093',
                      reca = 'RECA-EU',
                      cm10 = 'CM 010',
                      cm25ev = 'CM 025 EVER',
                      cm25ni = 'CM 025 NIVO')

# infiltration cluster colors -----

  globals$clust_colors <- c('#1' = 'darkolivegreen4',
                            '#2' = 'coral3',
                            '#3' = 'steelblue')

# END -----
