# Infiltration estimates by quanTIseq, PCA

  insert_head()

# containr list ------

  infil_pca <- list()

# globals: analysis table -------

  insert_msg('Globals setup')

  infil_pca$analysis_tbl <- list(tcga = tcga$expression %>%
                                   filter(tissue_type == 'Tumor'),
                                 emtab1980 = emtab1980$expression,
                                 gse73731 = gse73731$expression,
                                 gse167093 = gse167093$expression %>%
                                   filter(tissue_type == 'Tumor'),
                                 reca = reca$expression,
                                 cm10 = cm10$expression,
                                 cm25ev = cm25ev$expression,
                                 cm25ni = cm25ni$expression) %>%
    map(~.x[c('patient_id', globals$imm_estimates)]) %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(column_to_rownames, 'patient_id') %>%
    map(center_data, type = 'median')

# PCA: 4 dimensions ------

  insert_msg('4d PCA')

  ## objects

  infil_pca$pca_obj <- infil_pca$analysis_tbl %>%
    map(reduce_data,
        kdim = 4,
        red_fun = 'pca')

  ## plots

  infil_pca$pca_scree_plots <- infil_pca$pca_obj %>%
    map(plot,
        type = 'scree',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts,
         ~.x +
           expand_limits(y = 0) +
           labs(title = .y,
                subtitle = 'PCA variance'))

  infil_pca$pca_score_plots <- infil_pca$pca_obj %>%
    map(plot,
        type = 'scores',
        cust_theme = globals$common_theme,
        point_color = 'coral3') %>%
    map2(., globals$cohorts,
         ~.x +
           labs(title = .y,
                subtitle = 'PCA scores') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = 0,
                      linetype = 'dashed'))

  infil_pca$pca_loadings_plots <- infil_pca$pca_obj %>%
    map(plot,
        type = 'loadings',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts,
         ~.x +
           labs(title = .y,
                subtitle = 'PCA loadings') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = 0,
                      linetype = 'dashed'))

# UMAP ------

  insert_msg('UMAP')

  ## objects

  infil_pca$umap_obj <- infil_pca$analysis_tbl %>%
    map(reduce_data,
        distance_method = 'cosine',
        kdim = 2,
        red_fun = 'umap')

  ## score plots

  infil_pca$umap_score_plots <- infil_pca$umap_obj %>%
    map(plot,
        type = 'scores',
        cust_theme = globals$common_theme,
        point_color = 'coral3') %>%
    map2(., globals$cohorts,
         ~.x +
           labs(title = .y,
                subtitle = 'UMAP scores, cosine distance') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = 0,
                      linetype = 'dashed'))

# assessment of clustering tendencies -------

  insert_msg('Clustering tendencies')

  plan('multisession')

  infil_pca$clust_tend <- list(data = infil_pca$analysis_tbl,
                               n = map(infil_pca$analysis_tbl,
                                       ~nrow(.x)/2)) %>%
    future_pmap(get_clust_tendency,
                .options = furrr_options(seed = TRUE))

  plan('sequential')

# END -----

  insert_tail()
