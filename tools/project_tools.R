# This script provides tools for project data analyses

# tools ----

  library(plyr)
  library(tidyverse)


# annotation functions -----

  translate_gene <- function(id,
                             id_type = 'gene_symbol',
                             output_type = 'entrez_id',
                             dictionary = tcga$annotation) {

    ## gets gene identifier of interest

    naming_vec <- dictionary[[output_type]] %>%
      set_names(dictionary[[id_type]])

    naming_vec[id]

  }

  mm_inch <- function(x) x * 0.0393700787

# Customized Forest plot -----

  cust_forest <- function(data,
                          plot_title = 'Differences between TCGA/Checkmate and E-MTAB 1980/EURECA',
                          plot_subtitle = NULL) {

    ggplot(data,
           aes(x = 2^estimate,
               y = reorder(response, group_diff),
               color = study)) +
      geom_errorbarh(aes(xmin = 2^lower_ci,
                         xmax = 2^upper_ci),
                     height = 0,
                     position = position_dodge(width = 0.7)) +
      geom_point(shape = 16,
                 size = 2,
                 position = position_dodge(width = 0.7)) +
      scale_color_manual(values = c(tcga = 'coral3',
                                    cm25ev = 'orangered2',
                                    cm25ni = 'orangered3',
                                    cm10 = 'firebrick4',
                                    emtab1980 = 'steelblue2',
                                    reca = 'steelblue3',
                                    gse73731 = 'steelblue4',
                                    gse167093 = 'darkslateblue'),
                         labels = globals$cohorts,
                         name = '') +
      globals$common_theme +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(face = 'italic')) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = 'Fold regulation, CXC9 high vs low tumors')

  }

