# This script imports TCGA RNA Seq and clinical data.

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(stringi)
  library(soucer)
  library(xena)
  library(readxl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(furrr)
  library(GEOquery)
  library(immunedeconv)

  reduce <- purrr::reduce
  select <- dplyr::select

  insert_head()

# data containers -----

  tcga_expr <- list()
  tcga_clinics <- list()
  gene_sign <- list()

  timer_est <- list() ## infiltration estimates
  somut <- list() ## mutations

  tcga <- list() ## for the final data container

  ## containers for validation data sets

  emtab1980 <- list()
  gse73731 <- list()
  gse167093 <- list()
  cm10 <- list() ## Checkmate 10
  cm25ev <- list() ## checkmate 25 evrolimus
  cm25ni <- list() ## checkmate 25 nivolumab
  reca <- list() ## EU-RECA

# globals setup -----

  insert_msg('Globals setup')

  source_all('./tools/project_globals.R',
             message = TRUE, crash = TRUE)

# downloading TCGA clinical and expression data if not present on the disc ------

  insert_msg('Reading the TCGA clinical and expression data')

  if(length(list.files('./input data/TCGA', pattern = 'KIRC__gene*')) == 0) {

    source_all('./data import scripts/data_download_TCGA.R',
               crash = TRUE, message = TRUE)

    stop('Copy the expression and clinical files into ./input data/TCGA and run the pipeline again')

  }

# clearing the clinical and expression data sets ------

  insert_msg('Clearning the clinical and expression data')

  c('./data import scripts/data_clearing_tcga_clinics.R',
    './data import scripts/data_clearing_tcga_expression.R',
    './data import scripts/emtab1980.R',
    './data import scripts/gse73731.R',
    './data import scripts/gse167093.R',
    './data import scripts/checkmate.R',
    './data import scripts/reca.R') %>%
    source_all(source, message = TRUE, crash = TRUE)

# merging the expression data with the essential set of clinical information -----

  insert_msg('Appending the expression data with the essential clinical information')

  tcga[c('clinical',
         'drugs',
         'radiation')] <- tcga_clinics[c('clinical',
                                         'drugs',
                                         'radiation')]

  tcga[c('expression',
         'annotation')] <- tcga_expr[c('exprs',
                                       'annotation')]

  tcga$expression <- right_join(tcga$clinical,
                                tcga$expression,
                                by = 'patient_id')

# Duplicate handling: samples with duplicated analyte ID are removed ------

  insert_msg('Duplicate handling, removing duplicated analyte id')

  ## duplicated analyte id

  tcga$expression <- tcga$expression %>%
    dlply(.(tissue_type)) %>%
    map_dfr(filter, !duplicated(analyte_id)) %>%
    as_tibble()

# XENA mutations ------

  insert_msg('XENA mutations')

  source_all('./data import scripts/xena.R',
             crash = TRUE, message = TRUE)

# infiltration estimates ------

  insert_msg('Infiltration estimates')

  ## working primarily with cache

  if(!file.exists('./input data/infiltration/quantiseq.RDa')) {

    source_all('./data import scripts/infiltration.R',
               message = TRUE, crash = TRUE)

  } else {

    load('./input data/infiltration/quantiseq.RDa')

  }

  deconv <- deconv %>%
    map(column_to_rownames, 'cell_type') %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'patient_id') %>%
    map(as_tibble)

  tcga$expression <- left_join(tcga$expression,
                               deconv$tcga,
                               by = 'patient_id')

  emtab1980$expression <- left_join(emtab1980$expression,
                                    deconv$emtab1980,
                                    by = 'patient_id')

  gse73731$expression <- left_join(gse73731$expression,
                                   deconv$gse73731,
                                   by = 'patient_id')

  gse167093$expression <- left_join(gse167093$expression,
                                    deconv$gse167093,
                                    by = 'patient_id')

  reca$expression <- left_join(reca$expression,
                               deconv$reca,
                               by = 'patient_id')

  cm10$expression <- left_join(cm10$expression,
                               deconv$cm10,
                               by = 'patient_id')

  cm25ev$expression <- left_join(cm25ev$expression,
                                 deconv$cm25ev,
                                 by = 'patient_id')

  cm25ni$expression <- left_join(cm25ni$expression,
                                 deconv$cm25ni,
                                 by = 'patient_id')

  for(i in globals$imm_estimates) {

    tcga$expression <- tcga$expression %>%
      mutate(!!i := ifelse(tissue_type == 'Normal', NA, .data[[i]]))

    gse167093$expression <- gse167093$expression %>%
      mutate(!!i := ifelse(tissue_type == 'Normal', NA, .data[[i]]))

  }

# saving the workspace -------

  insert_msg('Saving the workspace')

  save(tcga, emtab1980, gse73731, gse167093, somut, reca, cm10, cm25ev, cm25ni,
       file = './input data/processed_data.RDa')

# END ----

  rm(tcga_clinics,
     tcga_expr,
     gene_sign,
     timer_est,
     deconv,
     i)

  insert_tail()
