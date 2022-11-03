# Imports and clears the EURECA study data

  insert_head()

# reading the clinical, sample, specimen and expression data -------

  insert_msg('Reading the data sets')

  reca[c('clinical',
         'specimen',
         'sample')] <- c('./input data/RECA/donor.RECA-EU.tsv',
                         './input data/RECA/specimen.RECA-EU.tsv',
                         './input data/RECA/sample.RECA-EU.tsv') %>%
    map(read_tsv)

  reca$expression <- read_tsv('./input data/RECA/exp_seq.RECA-EU.tsv')

# clearing the clinical, sample and specimen data ------

  insert_msg('Clearing the clinical data')

  ## sample and specimen, selecting primary tumors

  reca$specimen <- reca$specimen[c('icgc_donor_id',
                                   'icgc_specimen_id',
                                   'specimen_type',
                                   'specimen_type_other',
                                   'tumour_grade',
                                   'level_of_cellularity')] %>%
    filter(specimen_type == 'Primary tumour - solid tissue')

  reca$sample <- reca$sample %>%
    filter(icgc_specimen_id %in% reca$specimen$icgc_specimen_id) %>%
    select('icgc_donor_id',
           'icgc_sample_id',
           'icgc_specimen_id',
           'analyzed_sample_interval')

  reca$sample <- left_join(reca$sample,
                           reca$specimen,
                           by = c('icgc_donor_id',
                                  'icgc_specimen_id'))

  ## clinical data

  reca$clinical <- reca$clinical %>%
    left_join(reca$sample, by = 'icgc_donor_id') %>%
    transmute(patient_id = icgc_donor_id,
              sample_id = icgc_specimen_id,
              sex = car::recode(donor_sex, "'male' = 'Male'; 'female' = 'Female'"),
              sex = factor(sex, c('Female', 'Male')),
              death = ifelse(donor_vital_status == 'alive', 0, 1),
              age = donor_age_at_diagnosis,
              relapse = ifelse(is.na(donor_relapse_interval), 0, 1),
              rfs_days = ifelse(relapse == 1,
                                donor_relapse_interval,
                                donor_interval_of_last_followup),
              pt_stage = stri_replace(donor_tumour_stage_at_diagnosis, regex = 'N.*$', replacement = ''),
              pt_stage = ifelse(is.na(pt_stage) | stri_detect(pt_stage, fixed = 'T'),
                                pt_stage, paste0('T', pt_stage)),
              pt_stage = factor(pt_stage),
              pn_stage = stri_extract(donor_tumour_stage_at_diagnosis, regex = 'N.{1}'),
              pn_stage = factor(pn_stage),
              pm_stage = stri_extract(donor_tumour_stage_at_diagnosis, regex = 'M.{1}'),
              pm_stage = factor(pm_stage),
              tumor_grade = factor(tumour_grade, c('1', '2', '3', '4')),
              cellularity = factor(level_of_cellularity),
              os_days = donor_survival_time)

# annotation table ------

  insert_msg('Annotation table')

  reca$annotation <- unique(reca$expression$gene_id)

  reca$annotation  <- AnnotationDbi::select(org.Hs.eg.db,
                                            keys = set_names(reca$annotation,
                                                             reca$annotation),
                                            column = c('ENTREZID', 'SYMBOL'),
                                            keytype = 'ENSEMBL') %>%
    set_names(c('ensembl_id', 'entrez_id', 'gene_symbol')) %>%
    filter(!is.na(entrez_id),
           !is.na(gene_symbol),
           !duplicated(entrez_id)) %>%
    as_tibble

# clearing the expression data ---------

  insert_msg('Clearing the expression data, log2 conversion')

  reca$expression <- reca$expression[c('icgc_specimen_id', 'gene_id', 'normalized_read_count')] %>%
    set_names(c('sample_id', 'gene_id', 'count')) %>%
    filter(gene_id %in% reca$annotation$ensembl_id,
           sample_id %in% reca$sample$icgc_specimen_id) %>%
    mutate(gene_symbol = set_names(reca$annotation$gene_symbol,
                                   reca$annotation$ensembl_id)[gene_id])

  reca$expression <- reca$expression %>%
    dlply(.(sample_id)) %>%
    map2(., names(.),
         ~matrix(log2(.x$count + 1),
                 ncol = 1,
                 dimnames = list(.x$gene_symbol, .y))) %>%
    map(t) %>%
    do.call('rbind', .) %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

# merging with the clinical data -------

  insert_msg('merging with the clinical data')

  reca$expression <- right_join(reca$clinical,
                                reca$expression,
                                by = 'sample_id')

# END -----

  reca <- reca[c('clinical',
                 'expression',
                 'annotation')]

  insert_tail()
