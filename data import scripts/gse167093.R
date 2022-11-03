# Data from GSE167093

  insert_head()

  options(future.globals.maxSize = 1900*1024^2)

# reading the data -------

  insert_msg('Reading the data')

  gse167093$expression_pairs <- read_excel('./input data/GSE167093/GSE167093_Pairs_normalization_matrix.xlsx')

  gse167093$expression <- read_excel('./input data/GSE167093/GSE167093_Tumor_only_norm_matrix.xlsx')

  gse167093$clinical <- read_tsv('./input data/GSE167093/clinical.csv')

  gse167093$annotation <- read_tsv('./input data/GSE167093/GPL10558-50081.txt',
                                   skip = 30)

# clearing the clinical information ------

  insert_msg('Clinical data')

  gse167093$clinical <- gse167093$clinical %>%
    transmute(patient_id = ifelse(stri_detect(sample, fixed = '2020_'),
                                  sample, `GEO accession`),
              description = title,
              tissue_type = stri_extract(type, regex = 'normal|tumor'),
              tissue_type = car::recode(tissue_type,
                                        "'normal' = 'Normal'; 'tumor' = 'Tumor'"),
              sex = stri_extract(sex, regex = 'male|female'),
              sex = car::recode(sex, "'female' = 'Female'; 'male' = 'Male'"),
              sex = factor(sex),
              age = age,
              histology = histo,
              pt_stage = factor(paste0('T', stage)),
              tumor_grade = factor(grade))

# Annotation data -----

  insert_msg('Annotation data')

  gse167093$annotation <- gse167093$annotation %>%
    transmute(probe_ID = ID,
              gene_symbol = ILMN_Gene,
              entrez_id = Entrez_Gene_ID) %>%
    filter(complete.cases(.)) %>%
    arrange(gene_symbol)

# Expression data -------

  insert_msg('Expression data')

  gse167093$expression <- gse167093$expression %>%
    column_to_rownames('...1')

  ## unique genes, skipping the evidently non-gene symbols

  gse167093$unique_genes <- gse167093$annotation[26:nrow(gse167093$annotation), ] %>%
    filter(probe_ID %in% rownames(gse167093$expression)) %>%
    dlply(.(gene_symbol), function(x) x$probe_ID)

  ## integration by geometric mean

  plan('multisession')

  gse167093$expression <- rownames(gse167093$expression) %>%
    future_map(~unlist(gse167093$expression[.x, ])) %>%
    set_names(rownames(gse167093$expression))

  gse167093$expression  <- gse167093$unique_genes %>%
    future_map(~reduce(compact(gse167093$expression[.x]), `*`)^(1/length(.x)))

  gse167093$expression  <- gse167093$expression %>%
    do.call('cbind', .) %>%
    as.data.frame %>%
    rownames_to_column('patient_id') %>%
    as_tibble

# Normal- Tumor pair expression ------

  insert_msg('Normal - tumor pairs')

  gse167093$expression_pairs <- gse167093$expression_pairs %>%
    column_to_rownames('...1')

  ## probes of interest

  gse167093$unique_genes <- gse167093$annotation[26:nrow(gse167093$annotation), ] %>%
    filter(probe_ID %in% rownames(gse167093$expression_pairs)) %>%
    dlply(.(gene_symbol), function(x) x$probe_ID)

  gse167093$unique_genes <- gse167093$unique_genes[globals$cxcl_genes]

  ## integration by geometric mean

  gse167093$expression_pairs <-
    gse167093$expression_pairs[reduce(gse167093$unique_genes, c), ]

  gse167093$expression_pairs <- rownames(gse167093$expression_pairs) %>%
    map(~unlist(gse167093$expression_pairs[.x, ])) %>%
    set_names(rownames(gse167093$expression_pairs))

  gse167093$expression_pairs  <- gse167093$unique_genes %>%
    future_map(~reduce(compact(gse167093$expression_pairs[.x]), `*`)^(1/length(.x)))

  gse167093$expression_pairs  <- gse167093$expression_pairs %>%
    do.call('cbind', .) %>%
    as.data.frame %>%
    rownames_to_column('patient_id') %>%
    as_tibble

# merging with clinical information -------

  insert_msg('Joining with the clinical information')

  gse167093$expression  <- left_join(gse167093$expression,
                                     gse167093$clinical,
                                     by = 'patient_id')

  gse167093$expression_pairs <- left_join(gse167093$expression_pairs,
                                          gse167093$clinical,
                                          by = 'patient_id')

# END -----

  plan('sequential')

  gse167093 <- gse167093[c('clinical',
                           'expression',
                           'annotation',
                           'expression_pairs')]


  insert_tail()
