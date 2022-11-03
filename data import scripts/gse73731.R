# Import and wrangling of the GSE73731 data set

  insert_head()

  options(future.globals.maxSize = 1500*1024^2)

# Fetching the GEO data ------

  insert_msg('Fetching the GEO data')

  gse73731$geo_object <- getGEO(GEO = 'GSE73731',
                                destdir = './input data/GSE73731',
                                GSEMatrix = TRUE,
                                AnnotGPL = TRUE,
                                getGPL = TRUE)

# clinical data -----

  insert_msg('Cleaning the clinical data')

  gse73731$clinical <- pData(gse73731$geo_object[[1]]) %>%
    rownames_to_column('patient_id') %>%
    transmute(patient_id = patient_id,
              sex = car::recode(`gender:ch1`, "'F' = 'Female'; 'M' = 'Male'") %>%
                factor(c('Female', 'Male')),
              tumor_grade = factor(ifelse(`grade:ch1` == 'NA', NA, `grade:ch1`)),
              pt_stage = ifelse(`Stage:ch1` == 'NA', NA, `Stage:ch1`),
              pt_stage = ifelse(!is.na(pt_stage),
                                paste0('T', pt_stage),
                                pt_stage),
              pt_stage = factor(pt_stage)) %>%
    as_tibble

# Annotation data ------

  insert_msg('Annotation data')

  gse73731$annotation <- fData(gse73731$geo_object[[1]]) %>%
    transmute(probe_ID = ID,
              entrez_id = `Gene ID`,
              gene_symbol = `Gene symbol`,
              description = `Gene title`) %>%
    mutate(entrez_id = stri_split_fixed(entrez_id,
                                        pattern = '///',
                                        simplify = TRUE)[, 1],
           gene_symbol = stri_split_fixed(gene_symbol,
                                          pattern = '///',
                                          simplify = TRUE)[, 1]) %>%
    as_tibble

  ## unique gene symbol list

  gse73731$unique_genes <- gse73731$annotation %>%
    filter(gene_symbol != '') %>%
    dlply(.(gene_symbol), function(x) x$probe_ID)

# Expression data ------

  insert_msg('Expression data')

  gse73731$expression <- gse73731$geo_object[[1]] %>%
    exprs

  ## integration by geometric mean

  plan('multisession')

  gse73731$expression <- rownames(gse73731$expression) %>%
    map(~gse73731$expression[.x, ]) %>%
    set_names(rownames(gse73731$expression))

  gse73731$expression  <- gse73731$unique_genes %>%
    future_map(~reduce(compact(gse73731$expression[.x]), `*`)^(1/length(.x)))

  gse73731$expression  <- gse73731$expression %>%
    do.call('cbind', .) %>%
    as.data.frame %>%
    rownames_to_column('patient_id') %>%
    as_tibble

# Joining with the clinical information ------

  insert_msg('Joining the expression data with clinical information')

  gse73731$expression <- left_join(gse73731$expression,
                                   gse73731$clinical,
                                   by = 'patient_id')

# END -----

  plan('sequential')

  gse73731 <- gse73731[c('clinical',
                         'expression',
                         'annotation')]

  insert_tail()
