# Paper and supplement tables

  insert_head()

# container ------

  tbls <- list()
  suppl_tbls <- list()

# Table 1: cohort characteristic ------

  insert_msg('Table 1: cohort characteristic')

  tbls$cohort <- cohort$ft_table %>%
    filter(!Variable %in% c('Race',
                            'Laterality',
                            'Pathological stage',
                            'Sarcomatoid differentiation',
                            'Rhabdoid differentiation',
                            'Tumor shrinkage',
                            'Clinical benefit',
                            'IMDC risk group',
                            'Neoadjuvant treatment')) %>%
    mutate(Variable = ifelse(Variable == 'OS, days',
                             'Observation time, days',
                             Variable)) %>%
    map_dfc(stri_replace_all,
            fixed = 'Complete',
            replacement = 'complete') %>%
    map_dfc(stri_replace_all,
            fixed = 'Range',
            replacement = 'range') %>%
    map_dfc(stri_replace_all,
            regex = 'Mean.*\\nMedian\\s{1}=\\s{1}',
            replacement = '') %>%
    as_mdtable(label = 'table_1_cohort_features',
               ref_name = 'cohort',
               caption = 'Characteristic of the study cohorts. Numeric variables are presented as medians with interquartile ranges and ranges. Categorical variables are presented as percentage and total number within the complete observation set.')

# Supplementary Table S1: tissue donors for flow cytometry ------

  insert_msg('Tissue donors fro flow cytometry')

  suppl_tbls$fc <- tibble(Variable = c('Age, years',
                                       'Sex',
                                       'Tumor grade',
                                       'Tumor stage',
                                       'Lymph node stage',
                                       'Previous RCC treatment',
                                       'Tumor entity',
                                       'Peripheral leukocytes, cells/nL'),
                          `Patient 1` = c('61',
                                          'male',
                                          '2',
                                          'T1',
                                          'N0',
                                          'no',
                                          'ccRCC',
                                          '6.38'),
                          `Patient 2` = c('55',
                                          'male',
                                          '1',
                                          'T3',
                                          'N0',
                                          'no',
                                          'ccRCC',
                                          '9.24'),
                          `Patient 3` = c('81',
                                          'male',
                                          '2',
                                          'T1',
                                          'N0',
                                          'no',
                                          'ccRCC',
                                          '6.34'),
                          `Patient 4` = c('77',
                                          'male',
                                          '1',
                                          'T1',
                                          'N0',
                                          'no',
                                          'ccRCC',
                                          '4.24'))

  suppl_tbls$fc <- mdtable(suppl_tbls$fc,
                           label = 'table_s1_flow_cytometry',
                           ref_name = 'fc',
                           caption = 'Descriptive characteristics of tissue donors for flow cytometry analysis.')

# Supplementary Table S2: numbers of differentially regulated genes -----

  insert_msg('Table S2: numbers of differentially regulated genes')

  suppl_tbls$n_dge <- dge$n_genes %>%
    filter(expression_strata == 'CXCL9') %>%
    select(cohort,
           downregulated,
           upregulated) %>%
    set_names(c('Cohort',
                'N downregulated',
                'N upregulated')) %>%
    as_mdtable(label = 'table_s2_number_dge',
               ref_name = 'n_dge',
               caption = 'Numbers of genes differentially regulated in in CXCL9 high versus CXCL9 low tumor samples.')

# Supplementary Table S3: listing of differentially regulated genes -------

  insert_msg('Table S3: differentially regulated gene list')

  suppl_tbls$dge_genes <- dge$test_results %>%
    map(~.x$cxcl9) %>%
    map(filter, regulation != 'ns') %>%
    map(mutate, entrez_id = as.character(entrez_id)) %>%
    map2_dfr(., names(.),
             ~mutate(.x, study = globals$cohorts[.y])) %>%
    select(study, regulation,
           gene_symbol, entrez_id,
           estimate, lower_ci, upper_ci, p_adjusted) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>%
    set_names(c('Cohort',
                'Regulation',
                'Gene symbol',
                'Entrez ID',
                'Log2 fold-regulation estimate',
                'lower CI',
                'upper CI',
                'pFDR')) %>%
    as_mdtable(label = 'table_s3_dge',
               ref_name = 'dge_gene',
               caption = 'Genes differentially regulated in CXCL9 high versus CXCL9 low tumor samples. The table is available as a supplementary Excel file.')

# Supplementary Table S4: listing of significantly enriched GOs -------

  insert_msg('Table S4: Significantly enriched GOs')

  suppl_tbls$go <- enrich$go_cxcl9 %>%
    map(filter, p_adjusted < 0.05) %>%
    map(arrange, p_adjusted) %>%
    map2_dfr(., names(.),
             ~mutate(.x, study = .y)) %>%
    mutate(study = stri_split_fixed(study, pattern = '.', simplify = TRUE)[, 1],
           study = globals$cohorts[study]) %>%
    select(study,
           ID,
           Term,
           p_adjusted) %>%
    set_names(c('Cohort',
                'GO ID',
                'GO name',
                'pFDR')) %>%
    as_mdtable(label = 'table_s4_go',
               ref_name = 'go',
               caption = 'Biological process gene ontology (GO) terms significantly enriched in the set of genes differentially regulated in CXCL9 high versus CXCL9 low tumors. The table is available as a supplementary Excel file.')

# Supplementary Table S5: signaling pathways -------

  insert_msg('Table S5: signaling pathways')

  suppl_tbls$spia <- pathways$spia[c('tcga.cxcl9',
                                     'cm25ev.cxcl9',
                                     'cm25ni.cxcl9',
                                     'gse73731.cxcl9',
                                     'gse167093.cxcl9',
                                     'emtab1980.cxcl9',
                                     'reca.cxcl9')] %>%
    map(filter, pGFdr < 0.05) %>%
    map(arrange, -tA) %>%
    map2_dfr(., names(.),
             ~mutate(.x, study = .y)) %>%
    mutate(study = stri_split_fixed(study, pattern = '.', simplify = TRUE)[, 1],
           study = globals$cohorts[study]) %>%
    select(study,
           ID,
           Name,
           tA,
           pGFdr) %>%
    set_names(c('Cohort',
                'KEGG pathway ID',
                'KEGG pathway name',
                'Magnitude of regulation, tA',
                'pFDR')) %>%
    as_mdtable(label = 'table_s5_spia',
               ref_name = 'spia',
               caption = 'KEGG (Kyoto Encyclopedia of Genes and Genomes) signaling pathways predicted to be significantly regulated in CXCL9 high versus CXCL9 low tumors by the SPIA tool. The table is available as a supplementary Excel file.')

# Supplementary Table S6: metabolic pathways ------

  insert_msg('Table S6: predicted regulation of metabolic reactions')

  suppl_tbls$reacts <- meta$models %>%
    map(components, 'regulation') %>%
    map(filter, p_value < 0.05) %>%
    map(mutate,
        react_name = annotate_bigg(react_id,
                                   annotation_db = biggrExtra::Recon2D),
        react_name = unlist(react_name)) %>%
    map2_dfr(., names(.),
             ~mutate(.x, study = .y)) %>%
    mutate(study = globals$cohorts[study]) %>%
    select(study,
           react_id,
           react_name,
           fold_reg,
           lower_ci,
           upper_ci,
           p_value) %>%
    set_names(c('Cohort',
                'BiGG reaction ID',
                'BiGG reaction name',
                'fold-regulation',
                'lower CI',
                'upper CI',
                'pFDR')) %>%
    as_mdtable(label = 'table_s6_reactions',
               ref_name = 'reacts',
               caption = 'Metabolic reactions predicted to be significantly modulated in CXCL9 high versus CXCL9 low tumors. The table is available as a supplementary Excel file.')


# Saving the Supplementary tables -----

  insert_msg('Saving the supplementary tables')

  suppl_tbls %>%
    set_names(paste0('Table S', 1:length(suppl_tbls))) %>%
    write_xlsx('./paper/supplementary_tables.xlsx')

# END ------

  insert_tail()
