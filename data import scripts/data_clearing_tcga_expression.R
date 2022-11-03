# This script reads and wrangles the expression data
# transcript counts will be log2(x + 1) transformed
# the gene data set is constrained to the transcripts with known symbols

  insert_head()

# reading the text table -----

  insert_msg('Reading the text table')

  tcga_expr$file_path <- list.files('./input data/TCGA',
                                    pattern = 'KIRC__gene*',
                                    full.names = TRUE) %>%
    sort

  tcga_expr$file_path <-  tcga_expr$file_path[length( tcga_expr$file_path)]

  tcga_expr$exprs <- read_tsv(tcga_expr$file_path)

# Clearing and transforming ----

  insert_msg('Clearing the data and transforming')

  tcga_expr$exprs <- tcga_expr$exprs %>%
    mutate(gene_id = stri_replace(gene_id,
                                  fixed = '|',
                                  replacement = '_')) %>%
    filter(!stri_detect(gene_id,
                        fixed = '?'))

  ## generating the annotation table

  tcga_expr$annotation <- tibble(gene_symbol = stri_split_fixed(tcga_expr$exprs$gene_id,
                                                                pattern = '_',
                                                                simplify = TRUE)[, 1],
                                 entrez_id = stri_split_fixed(tcga_expr$exprs$gene_id,
                                                              pattern = '_',
                                                              simplify = TRUE)[, 2])

  ## clearing the issue with the duplicated symbol entry for SLC35E2

  tcga_expr$annotation[tcga_expr$annotation$entrez_id == '9906', 'gene_symbol'] <- 'SLC35E2A'
  tcga_expr$annotation[tcga_expr$annotation$entrez_id == '728661', 'gene_symbol'] <- 'SLC35E2B'

  ## log2 transformation

  tcga_expr$exprs <- tcga_expr$exprs %>%
    select( - gene_id) %>%
    map_dfc(~log2(.x + 1))

  ## transposing

  tcga_expr$exprs <- tcga_expr$exprs %>%
    as.data.frame

  rownames(tcga_expr$exprs) <- tcga_expr$annotation$gene_symbol

  tcga_expr$exprs <- tcga_expr$exprs %>%
    t %>%
    as.data.frame %>%
    rownames_to_column('analyte_id') %>%
    as_tibble

  ## extracting the patient_id, sample_id and tissue type

  tcga_expr$exprs <- tcga_expr$exprs %>%
    mutate(patient_id = stri_extract(analyte_id,
                                     regex = 'TCGA-\\w{2}-\\w{4}'),
           sample_id = stri_extract(analyte_id,
                                    regex = 'TCGA-\\w{2}-\\w{4}-\\d{2}'),
           tissue_type = stri_extract(analyte_id, regex = '-\\d{2}\\w{1}-'),
           tissue_type = stri_extract(tissue_type, regex = '\\d{2}'),
           tissue_type = as.numeric(tissue_type),
           tissue_type = ifelse(tissue_type == 1,
                                'Tumor',
                                ifelse(tissue_type > 10,
                                       'Normal', NA)) %>%
             factor)

# END ----

  tcga_expr <- tcga_expr[c('exprs',
                           'annotation')]

  insert_tail()
