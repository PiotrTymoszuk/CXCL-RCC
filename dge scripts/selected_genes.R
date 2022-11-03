# Checks for genes regulated in cohorts, where CXCL9 was associated
# with favorable prognosis (E-MTAB 1980 and EU RECA) but not in the cohorts
# where CXCL9 was an unfavorable marker (TCGA, CM10, CM25)
# CheckMate 010 is left ut, hence only few genes were regulated there

  insert_head()

# Containers ----

  gene_comp <- list()

# Globals: lists of differentially regulated genes -----

  insert_msg('globals')

  gene_comp$studies <- c('tcga',
                         'cm25ev',
                         'cm25ni',
                         'cm10',
                         'reca',
                         'emtab1980')

  gene_comp$all_genes <- dge$dict[gene_comp$studies ] %>%
    map(~.x$gene_symbol) %>%
    reduce(intersect)

  gene_comp$signif_genes <- dge$signif_genes[gene_comp$studies] %>%
    map(~.x$cxcl9) %>%
    transpose

  gene_comp$signif_genes <- gene_comp$signif_genes %>%
    map(~map(.x, ~.x[.x %in% gene_comp$all_genes])) %>%
    map(reduce, union)

# Expression estimates of the significant genes ----

  insert_msg('Regulation estimates')

  gene_comp$estimates <- dge$test_results[gene_comp$studies] %>%
    map(~.x$cxcl9) %>%
    map(filter, response %in% reduce(gene_comp$signif_genes, union)) %>%
    map(mutate,
        entrez_id = as.character(entrez_id)) %>%
    map2_dfr(., names(.),
             ~mutate(.x,
                     study = .y,
                     study = factor(study,
                                    gene_comp$studies))) %>%
    mutate(study_group = ifelse(study %in% c('emtab1980',
                                             'reca',
                                             'gse73731',
                                             'gse167093'),
                                'favorable_cohort',
                                'unfavorable_cohort'),
           study_group = factor(study_group, c('favorable_cohort',
                                               'unfavorable_cohort')))

# Checking for significant differences in regulation between the cohort groups ------

  insert_msg('Significance')

  gene_comp$test_results <-
    test_two_groups(data = gene_comp$estimates %>%
                      pivot_wider(names_from = response,
                                  values_from = estimate),
                    variables = reduce(gene_comp$signif_genes, union),
                    split_fct = 'study_group',
                    type = 't',
                    adj_method = 'BH',
                    .parallel = FALSE) %>%
    mutate(group_diff = estimate,
           group_p = p_value)

  gene_comp$estimates <- left_join(gene_comp$estimates,
                                   gene_comp$test_results[c('response',
                                                            'group_diff',
                                                            'group_p')],
                                   by = 'response')

# Top genes with the larges differences between the study cohorts ------

  insert_msg('Significant and top genes')

  ## significant genes

  gene_comp$signif_diff <- gene_comp$estimates %>%
    filter(group_p < 0.05) %>%
    filter(!duplicated(response)) %>%
    dlply(.(sign(group_diff))) %>%
    map(~set_names(.x$entrez_id, .x$response)) %>%
    set_names(c('down', 'up'))

  ## log2 fold regulation vectors

  gene_comp$reg_vector <-  gene_comp$estimates %>%
    filter(group_p < 0.05) %>%
    filter(!duplicated(response))

  gene_comp$reg_vector <- set_names(gene_comp$reg_vector$group_diff,
                                    gene_comp$reg_vector$entrez_id)

  gene_comp$all_vector <- dge$dict[gene_comp$studies] %>%
    map(~.x$entrez_id) %>%
    reduce(intersect)

  ## top genes

  gene_comp$top_genes <- gene_comp$estimates %>%
    filter(group_p < 0.05) %>%
    filter(!duplicated(response)) %>%
    dlply(.(sign(group_diff)), top_n, 20, abs(group_diff)) %>%
    map(~.x$response) %>%
    set_names(c('down', 'up'))

# Volcano plot with the regulation estimates -------

  insert_msg('Volcano plot')

  gene_comp$volcano_plot <-
    plot_volcano(data = gene_comp$test_results,
                 regulation_variable = 'group_diff',
                 p_variable = 'group_p',
                 signif_level = 0.05,
                 regulation_level = 0,
                 x_lab = expression('log'[2] * ' fold regulation, unfavorable vs favorable marker cohorts'),
                 y_lab = expression('-log'[10] * ' p raw'),
                 top_significant = 20,
                 label_variable = 'response',
                 label_type = 'label',
                 txt_face = 3,
                 plot_title = 'Differences in CXCL9 strata transcriptome',
                 plot_subtitle = 'Unfavorable vs favorable marker cohorts, T test',
                 cust_theme = globals$common_theme) +
    geom_vline(xintercept = 0,
               linetype = 'dashed')

  gene_comp$volcano_plot <-
    gene_comp$volcano_plot +
    labs(tag = gene_comp$volcano_plot$labels$tag %>%
           paste0('\n', .))

# Forest plot with the regulation estimates ----
  # for the genes with the highest differences between the cohort groups

  insert_msg('Forest plot with the top regulated genes')

  gene_comp$top_forest <- gene_comp$top_genes %>%
    map(~filter(gene_comp$estimates, response %in% .x)) %>%
    map2(., c('Upregulated in favorable marker cohorts, top 20 genes',
              'Downregulated in favorable marker cohorts, top 20 genes'),
         cust_forest)

# GO enrichment analysis ------

  insert_msg('GO enrichment analysis')

  plan('multisession')

  gene_comp$go_enrichment <- gene_comp$signif_diff %>%
    future_map(~goana(de = unname(.x)),
               .options = furrr_options(seed = TRUE))

  plan('sequential')

  gene_comp$go_enrichment <- gene_comp$go_enrichment %>%
    map(rownames_to_column, 'ID') %>%
    map(filter, Ont == 'BP') %>%
    map(mutate, p_adjusted = p.adjust(P.DE)) %>%
    map(filter, P.DE < 0.05) %>%
    map(as_tibble)

# Plotting the significant GO enrichment, upregulated features only -----

  insert_msg('Significant GOs bar plot, upregulated features only')

  gene_comp$go_plot <-
    plot_signifcant(data = gene_comp$go_enrichment$up,
                    p_variable = 'p_adjusted',
                    label_variable = 'Term',
                    top_significant = 20,
                    plot_title = 'BP GO enrichment, top 20 terms',
                    plot_subtitle = 'Genes unfavorable vs favorable marker cohorts',
                    x_lab = expression('-log'[10] * ' pFDR'),
                    cust_theme = globals$common_theme)

# SPIA analysis of the signaling pathways ------

  insert_msg('SPIA analysis')

  gene_comp$spia_results <- spia(de = gene_comp$reg_vector,
                                 all = gene_comp$all_vector,
                                 verbose = FALSE) %>%
    as_tibble

# Significant pathways and genes of interest -----

  insert_msg('Significant pathways and their member genes')

  gene_comp$signif_pathways <- gene_comp$spia_results %>%
    filter(pG < 0.05) %>%
    .$ID %>%
    paste0('path:hsa', .)

  gene_comp$signif_path_genes <- enrich$gene_pathway %>%
    filter(PathwayID %in% gene_comp$signif_pathways,
           GeneID %in% names(gene_comp$reg_vector)) %>%
    set_names(c('entrez_id', 'kegg_id')) %>%
    filter(!duplicated(entrez_id)) %>%
    mutate(gene_symbol = translate_gene(entrez_id,
                                       id_type = 'entrez_id',
                                       output_type = 'gene_symbol'),
           ID = stri_replace(kegg_id, fixed = 'path:hsa', replacement = ''),
           path_name = set_names(gene_comp$spia_results$Name,
                                 gene_comp$spia_results$ID)[ID])

# A volcano plot with the SPIA results -------

  insert_msg('SPIA results volcano plot')

  gene_comp$spia_volcano <-
    plot_volcano(data = gene_comp$spia_results,
                 regulation_variable = 'tA',
                 p_variable = 'pG',
                 signif_level = 0.05,
                 regulation_level = 0,
                 x_lab = 'Pathway modulation (tA), unfavorable vs favorable marker cohorts',
                 y_lab = expression('-log'[10] * ' pG raw'),
                 top_significant = 10,
                 label_variable = 'Name',
                 label_type = 'label',
                 plot_title = 'Signaling pathway modulation',
                 plot_subtitle = 'CXCL9 transcriptome, unfavorable vs favorable marker cohorts, SPIA',
                 cust_theme = globals$common_theme) +
    geom_vline(xintercept = 0,
               linetype = 'dashed')

  gene_comp$spia_volcano <- gene_comp$spia_volcano  +
    labs(tag = gene_comp$spia_volcano$labels$tag %>%
           paste0('\n', .))

# Plots of tA estimates for the top pathways ------

  insert_msg('Plots of TA for top pathways')

  gene_comp$top_spia_forest <- gene_comp$spia_results %>%
    filter(tA != 0) %>%
    ddply(.(sign(tA)), top_n, 20, abs(tA)) %>%
    as_tibble %>%
    plot_top(regulation_variable = 'tA',
             label_variable = 'Name',
             p_variable = 'pG',
             signif_level = 0.05,
             regulation_level = 0,
             top_regulated = 40,
             plot_title = 'Signaling pathway modulation, top 20 pathways',
             plot_subtitle = 'CXCL9 transcriptome, unfavorable vs favorable marker cohorts, SPIA',
             x_lab = 'Pathway modulation (tA), unfavorable vs favorable marker cohorts',
             cust_theme = globals$common_theme)


# Forest plots of regulation estimates of the genes associated with signif. pathways ----

  insert_msg('Signif pathways genes')

  gene_comp$path_gene_plots <- gene_comp$estimates %>%
    filter(entrez_id %in% gene_comp$signif_path_genes$entrez_id) %>%
    left_join(gene_comp$signif_path_genes, by = c('entrez_id', 'gene_symbol')) %>%
    dlply(.(path_name), as_tibble) %>%
    list(data = .,
         plot_title = names(.)) %>%
    pmap(cust_forest,
         plot_subtitle = 'unfavorable vs favorable marker cohorts')

# END ----

  insert_tail()
