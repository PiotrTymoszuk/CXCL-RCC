# Imports mutations (Mutect 2 algorithm) stored in a datase derived from XENA
# (https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Kidney%20Clear%20Cell%20Carcinoma%20(KIRC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# calculates mutation count and mutation frequencies
# mutations affecting UTRs, exons and introns are analyzed

  insert_head()

# loading the database ------

  insert_msg('Loading the database')

  somut$db <- load_xena('./input data/xena/TCGA-KIRC.mutect2_snv.tsv') %>%
    filter(effect != 'synonymous_variant')

# counting the mutations -----

  insert_msg('Total mutation burden')

  somut$tmb <- count_mutations(somut$db) %>%
    mutate(patient_id = stri_replace(sample_id,
                                     regex = '-01(A|B)$',
                                     replacement = ''))

# mutation table -------

  insert_msg('Mutation table')

  somut$mut_tbl <- tab_mutations(somut$db)

  rownames(somut$mut_tbl) <- stri_replace(rownames(somut$mut_tbl),
                                          regex = '-01(A|B)$',
                                          replacement = '') ## conversion to patient_ID

# Mutation frequencies, identifying top mutations present in at least 2% samples ----

  insert_msg('Mutation frequency and top most frequent mutations')

  somut$mut_freq <- freq_mutations(somut$db)

  somut$top_mutations <- somut$mut_freq %>%
    filter(percent > 2) %>%
    mutate(variable = paste0(gene, '_mut'))

# A convenience table with the most common mutations and TMB ----

  insert_msg('Common mutation and TMB table')

  somut$mut_summary <- somut$mut_tbl[, somut$top_mutations$gene] %>%
    as.matrix %>%
    as.data.frame %>%
    rownames_to_column('patient_id') %>%
    as_tibble %>%
    set_names(c('patient_id', somut$top_mutations$variable)) %>%
    left_join(somut$tmb, by = 'patient_id')

# END ----

  insert_tail()
