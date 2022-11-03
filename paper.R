# Paper scripts

# tools -------

  library(soucer)
  library(plyr)
  library(tidyverse)
  library(writexl)
  library(knitr)
  library(rmarkdown)
  library(bookdown)
  library(figur)
  library(cowplot)
  library(flextable)
  library(rmdformats)
  library(exda)
  library(stringi)
  library(biggrExtra)

  insert_head()

  explora <- exda::explore

  source_all('./tools/project_tools.R')

# scripts -------

  insert_msg('Paper scripts')

  c('./paper scripts/tables.R',
    './paper scripts/figures.R',
    './paper scripts/supplement.R',
    './paper scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
