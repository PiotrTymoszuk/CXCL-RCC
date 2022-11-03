# Pathview plots for particular signaling pathways, pooled expression estimates

  insert_head()

# container -----

  dge_path <- list()

# globals ------

  insert_msg('Globals setup')

  ## KEGG IDs of the pathways of interest

  dge_path$kegg_ids <- pathways$cmm_pathways

# pooled gene expression regulation estimates -------

  insert_msg('Gene expression estimates for each cohort')

  dge_path$cmm_genes <- dge$test_results %>%
    map(~.x$cxcl9$entrez_id) %>%
    reduce(intersect)

  dge_path$analysis_tbl <- dge$test_results %>%
    map(~.x$cxcl9) %>%
    map(filter, entrez_id %in% dge_path$cmm_genes) %>%
    map(mutate,
        error = abs(estimate/stat),
        error = log(2) * error * 2^estimate,
        estimate = 2^estimate) %>%
    map(~.x[c('entrez_id', 'estimate', 'error')])

  plan('multisession')

  dge_path$gene_expression <- dge_path$analysis_tbl %>%
    future_map(dlply, 'entrez_id', function(x) x$estimate) %>%
    transpose %>%
    map(unlist) %>%
    map(unname)

  dge_path$gene_error <- dge_path$analysis_tbl %>%
    future_map(dlply, 'entrez_id', function(x) x$error) %>%
    transpose %>%
    map(unlist) %>%
    map(unname)

# calculation of the pooled estimates ------

  insert_msg('Computing pooled estimates')

  dge_path$cmm_expression <-
    future_map2(dge_path$gene_expression,
                dge_path$gene_error,
                safely(metagen)) %>%
    map(~.x$result) %>%
    compact %>%
    map2(., names(.),
         ~tibble(entrez_id = .y,
                 estimate = .x$TE.common,
                 error = .x$seTE.common,
                 lower_ci = .x$lower.common,
                 upper_ci = .x$upper.common))

  plan('sequential')

  dge_path$cmm_expression <-
    do.call('rbind', dge_path$cmm_expression)


# pooled regulation estimates in a linear scale -----

  insert_msg('Pooled regulation estimates, linear scale')

  dge_path$reg_vec <- set_names(dge_path$cmm_expression$estimate,
                                dge_path$cmm_expression$entrez_id)

# Generating the pathview images ------

  insert_msg('pathview images')

  enter_directory('./report/pathview')

  dge_path$path_images <- dge_path$kegg_ids %>%
    map(~pathview(gene.data = dge_path$reg_vec,
                  pathway.id = .x,
                  low = list(gene = 'steelblue',
                             cpd = 'steelblue'),
                  mid = list(gene = 'white',
                             cpd = 'white'),
                  high = list(gene = 'firebrick',
                              cpd = 'firebrick'),
                  limit = list(gene = 5, cpd = 5)))

  go_proj_directory()

# END ------

  insert_tail()
