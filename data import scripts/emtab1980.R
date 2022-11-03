# Reading and clearing the E-MTAB 1980 data set

  insert_head()

# Reading the data ------

  insert_msg('Reading the data')

  emtab1980$clinical <- read_excel('./input data/EMTAB1980/clinical.xlsx',
                                   skip = 1)

  emtab1980$expression <- read_tsv('./input data/EMTAB1980/expression.txt')

  emtab1980$annotation <- read_tsv('./input data/EMTAB1980/annotation.txt',
                                   skip = 19) %>%
    filter(!`Reporter Name` %in% c('GE_BrightCorner',
                                   'DarkCorner'),
           `Reporter Group[role]` == 'Experimental')

# clearing the clinical information ------

  insert_msg('Cleaning the clinical information')

  emtab1980$clinical <- emtab1980$clinical %>%
    transmute(patient_id = `sample ID`,
              sex = car::recode(Sex, "'F' = 'Female'; 'M' = 'Male'") %>%
                factor(c('Female', 'Male')),
              age = Age,
              tumor_grade = factor(`Fuhrman grade`, c('1', '2', '3', '4')),
              pt_stage = stri_replace(`Stage at diagnosis`,
                                      regex = 'N.*',
                                      replacement = '') %>%
                stri_replace(fixed = 'pT', replacement = 'T') %>%
                factor,
              pn_stage = stri_extract(`Stage at diagnosis`,
                                      regex = 'N\\w{1}') %>%
                factor,
              pm_stage = stri_extract(`Stage at diagnosis`,
                                      regex = 'M\\w{1}') %>%
                factor,
              relapse = ifelse(`metastases during the observation period` == '-',
                               0, 1),
              relapse_type = ifelse(`metastases during the observation period` == '-',
                                    NA, `metastases during the observation period`) %>%
                factor,
              os_months = `observation period (month)`,
              death = ifelse(outcome == 'alive', 0, 1),
              os_days = 30.437 * os_months)

# clearing the annotation data set -----

  insert_msg('Cleaning the annotation')

  emtab1980$annotation <- emtab1980$annotation %>%
    transmute(probe_ID = `Reporter Name`) %>%
    full_join(emtab1980$expression[c('REF', 'SystematicName')] %>%
                set_names(c('probe_ID', 'db_entry')),
              by = 'probe_ID')

  ## fetching the gene numbers and symbols

  emtab1980$annot_ensembl <- AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = set_names(emtab1980$annotation$db_entry,
                                                                    emtab1980$annotation$probe_ID),
                                                   column = c('ENTREZID', 'SYMBOL'),
                                                   keytype = 'ENSEMBLTRANS') %>%
    filter(complete.cases(.)) %>%
    filter(!duplicated(ENSEMBLTRANS)) %>%
    set_names(c('db_entry', 'entrez_id', 'gene_symbol')) %>%
    as_tibble

  emtab1980$annot_refseq <- AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = set_names(emtab1980$annotation$db_entry,
                                                                   emtab1980$annotation$probe_ID),
                                                  column = c('ENTREZID', 'SYMBOL'),
                                                  keytype = 'REFSEQ') %>%
    filter(complete.cases(.)) %>%
    filter(!duplicated(REFSEQ)) %>%
    set_names(c('db_entry', 'entrez_id', 'gene_symbol')) %>%
    as_tibble

  ## common annotation data frame

  emtab1980$annotation <- left_join(emtab1980$annotation,
                                    rbind(emtab1980$annot_ensembl,
                                          emtab1980$annot_refseq),
                                    by = 'db_entry') %>%
    filter(!duplicated(probe_ID),
           complete.cases(.))

  ## unique gene symbols

  emtab1980$unique_genes <- emtab1980$annotation %>%
    dlply(.(gene_symbol), function(x) x$probe_ID)

# Clearing the expression data set -------

  insert_msg('Clearing the expression data set')

  ## integration by mean for multiplicated probes

  plan('multisession')

  emtab1980$expression <- emtab1980$expression %>%
    select(REF, starts_with('ccRCC')) %>%
    dlply(.(REF), function(x) x[, -1]) %>%
    future_map(colMeans)

  ## integration by mean for multiple probes per gene

  emtab1980$expression <- emtab1980$unique_genes %>%
    future_map(~reduce(compact(emtab1980$expression[.x]), `*`)^(1/length(.x)))

  emtab1980$expression <- emtab1980$expression %>%
    do.call('cbind', .) %>%
    as.data.frame %>%
    rownames_to_column('patient_id') %>%
    as_tibble

# Joining the expression data set with the clinical information ------

  insert_msg('Appending the expression data with the clinical information')

  emtab1980$expression <- left_join(emtab1980$expression,
                                    emtab1980$clinical,
                                    by = 'patient_id')

# END -----

  plan('sequential')

  emtab1980 <- emtab1980[c('clinical',
                           'expression',
                           'annotation')]

  insert_tail()
