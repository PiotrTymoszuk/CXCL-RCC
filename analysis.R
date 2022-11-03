# Sourcing of the analysis scripts.

# tools ------

  library(plyr)
  library(tidyverse)
  library(soucer)
  library(exda)
  library(kmOptimizer)
  library(survival)
  library(survminer)
  library(somKernels)
  library(clustTools)
  library(furrr)
  library(microViz)
  library(stringi)
  library(coxExtensions)
  library(survival)
  library(survminer)
  library(microViz)
  library(meta)

  insert_head()

  source_all('./tools/project_tools.R')

# analysis scripts ------

  c('./analysis scripts/cohort_features.R', ## characteristic of the study cohorts
    './analysis scripts/norm_tumor.R', ## normal vs tumor expression
    './analysis scripts/expr_clinics.R', ## expression and clinical features
    './analysis scripts/coexpression.R', ## co-expression of CXCL9/10/11 and CCL11 in the tumor tissue
    './analysis scripts/infil_corr.R',  ## correlation of QuanTIseq infiltration estimates with expression
    './analysis scripts/infil_pca.R',  ## PCA, UMAP and clustering tendency of the infiltration data
    './analysis scripts/infil_clust.R',  ## semi-supervised clustering
    './analysis scripts/infil_clust_characteristic.R', ## characteristic of the immune infiltration clusters
    './analysis scripts/clust_survival.R', ## survival in the infiltration clusters
    './analysis scripts/survival.R', ## survival differences between high and low expressors
    './analysis scripts/survival_cox.R') %>%  ## survival differences between high and low expressors, Cox modeling
    source_all(crash = TRUE, message = TRUE)

# END ------

  insert_tail()
