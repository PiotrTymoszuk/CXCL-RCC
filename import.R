# This script imports TCGA RNA Seq and clinical data, works with cache.

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(soucer)

  insert_head()

  source_all('./tools/project_globals.R',
             message = TRUE, crash = TRUE)

# loading the processed data ------

  insert_msg('Loading the processed data')

  if(file.exists('./input data/processed_data.RDa')) {

    load('./input data/processed_data.RDa')

  } else {

    source_all('./data import scripts/import_all.R')

  }

# END ----

  insert_tail()
