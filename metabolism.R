# Differences in metabolism based on differential gene expression
# between the CXCL9 low and CXCL9 high kidney tumors

# Tools -------

  library(plyr)
  library(tidyverse)
  library(exda)
  library(microViz)
  library(soucer)
  library(ggrepel)
  library(stringi)
  library(furrr)
  library(BiGGR)
  library(biggrExtra)
  library(meta)
  library(readxl)
  library(gseaTools)

  explore <- exda::explore

  source_all('./tools/project_tools.R')

  insert_head()

# scripts ------

  insert_msg('Metabolism scripts')

  c('./metabolism scripts/metab_testing.R', ## investigating regulation of metabolic reactions
    './metabolism scripts/metab_plots.R', ## reaction regulation plots
    './metabolism scripts/metab_meta.R', ## pooled reaction regulation estimates
    './metabolism scripts/metab_hg.R', ## hypergraphs with pooled estimates for selected pathways
    './metabolism scripts/metab_genes.R', ## GSVA scores for genes associated with metabolic pathways
    './metabolism scripts/metab_cxcl9.R', ## differences in GSVA scores for metabolism genes in the CXCL9 expression strata
    './metabolism scripts/metab_infiltration.R', ## correlation of GSVA scores for metabolism genes with infiltration
    './metabolism scripts/metab_reg_plots.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
