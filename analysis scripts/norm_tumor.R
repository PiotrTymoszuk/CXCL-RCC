# Normal and tumor tissue expression of CXCL9/10/11.
# Paired T test, TCGA and GSE167093 cohorts

  insert_head()

# container list ------

  norm_tumor <- list()

# globals -----

  insert_msg('Globals setup')

  ## TCGA: IDs of the normal tissue donors

  norm_tumor$pat_ids <- tcga$expression %>%
    filter(tissue_type == 'Normal') %>%
    .$patient_id

  norm_tumor$analysis_tbl$tcga <- tcga$expression %>%
    filter(patient_id %in% norm_tumor$pat_ids) %>%
    dlply(.(tissue_type), as_tibble) %>%
    map(select, patient_id, tissue_type, all_of(globals$cxcl_genes))

  ## GSE167093, paired samples only

  norm_tumor$analysis_tbl$gse167093 <- gse167093$expression_pairs %>%
    select(description, tissue_type, all_of(globals$cxcl_genes)) %>%
    mutate(patient_id = stri_extract(description, regex = 'patient\\s{1}\\d{1,3}')) %>%
    dlply(.(patient_id)) %>%
    map_dfr(function(x) if(nrow(x) == 2) x else NULL) %>%
    dlply(.(tissue_type), as_tibble)

# testing ------

  insert_msg('Testing')

  norm_tumor$test_results <- norm_tumor$analysis_tbl %>%
    map(~compare_variables(.x$Tumor,
                           .x$Normal,
                           variables = globals$cxcl_genes,
                           what = 'eff_size',
                           types = 'paired_cohen_d',
                           ci = FALSE,
                           pub_styled = TRUE,
                           adj_method = 'BH',
                           simplify_p = FALSE)) %>%
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))

# plotting -----

  insert_msg('Paired plots')

  norm_tumor$plots <- list(cohort = norm_tumor$analysis_tbl,
                           stats = norm_tumor$test_results,
                           suffix = c('TCGA', 'GSE167093')) %>%
    pmap(function(cohort, stats, suffix)  list(variable = globals$cxcl_genes,
                                               plot_title = paste(globals$cxcl_genes,
                                                                  suffix,
                                                                  sep = ', '),
                                               plot_subtitle = stats$plot_cap) %>%
           pmap(plot_variable,
                cohort$Normal,
                cohort$Tumor,
                data_names = c('Normal', 'Tumor'),
                cust_theme = globals$common_theme,
                type = 'paired',
                y_lab = expression('log'[2]*' expression')) %>%
           map(~.x +
                 scale_fill_manual(values = c('steelblue', 'coral3')) +
                 labs(tag = paste('\nn =', nrow(norm_tumor$analysis_tbl[[1]]))) +
                 theme(plot.title = element_text(face = 'bold.italic'))) %>%
           set_names(globals$cxcl_genes)) %>%
    map2(., map(norm_tumor$analysis_tbl, ~nrow(.x[[1]])),
         function(cohort, n) cohort %>%
           map(~.x + labs(tag = paste('\nn =', n))))

# END -----

  insert_tail()
