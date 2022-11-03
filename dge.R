# Differential gene expression.

# tools -----

  library(plyr)
  library(tidyverse)
  library(exda)
  library(microViz)
  library(soucer)
  library(limma)
  library(SPIA)
  library(ggrepel)
  library(stringi)
  library(furrr)
  library(gseaTools)
  library(clustTools)
  library(meta)
  library(pathview)

  source_all('./tools/project_tools.R')

  insert_head()

# analysis scripts -----

  ## working with a cache with the DGE results

  if(file.exists('./input data/dge.RData')) {

    load('./input data/dge.RData')

  } else {

    source_all('./dge scripts/dge_analysis.R',
               message = TRUE, crash = TRUE)

  }

  c('./dge scripts/enrichment.R', ## GO and KEGG term enrichment
    './dge scripts/spia.R', ## signaling pathway modulation,
    './dge scripts/clinical.R', ## clinical characteristic of the expression strata
    './dge scripts/mutations.R', ## mutation frequency changes between the expression strata, TCGA only
    './dge scripts/signatures.R', ## signatures of T cell activation and suppression
    './dge scripts/selected_genes.R', ## genes differentially regulated between the cohorts by CXCL9
    './dge scripts/infiltration.R', ## differences in infiltration between the CXCL9 expression strata
    './dge scripts/pathview.R') %>% ## pathwview diagrams for selected signaling pathways
    source_all(message = TRUE, crash = TRUE)

# END ----

  insert_tail()
