# Infiltration cluster development in the TCGA training cohort.
# SOM and HCl clustering.

  insert_head()

# container list -----

  infil_clust_dev <- list()

# globals ------

  insert_msg('Analysis tables')

  infil_clust_dev$analysis_tbl <- list(tcga = tcga$expression %>%
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
    map(min_max)

  infil_clust_dev$test_dist <- c('euclidean', 'manhattan', 'cosine')

# algorithms -------

  insert_msg('Algorithms to test')

  plan('multisession')

  infil_clust_dev$algos[c('HCL_Euclidean',
                          'HCL_Manhattan',
                          'HCL_Cosine')] <- list(distance_som = infil_clust_dev$test_dist,
                                                 k = c(3, 3, 3)) %>%
    future_pmap(combi_cluster,
                data = infil_clust_dev$analysis_tbl$tcga,
                xdim = 10,
                ydim = 10,
                topo = 'hexagonal',
                neighbourhood.fct = 'gaussian',
                rlen = 2000,
                toroidal = FALSE,
                node_clust_fun = hcluster,
                distance_nodes = 'euclidean',
                .options = furrr_options(seed = TRUE,
                                         packages = c('clustTools',
                                                      'somKernels')))

  infil_clust_dev$algos[c('KMEANS_Euclidean',
                          'KMEANS_Manhattan',
                          'KMEANS_Cosine')] <- list(distance_som = infil_clust_dev$test_dist,
                                                    k = c(3, 3, 3)) %>%
    future_pmap(combi_cluster,
                data = infil_clust_dev$analysis_tbl$tcga,
                xdim = 10,
                ydim = 10,
                topo = 'hexagonal',
                neighbourhood.fct = 'gaussian',
                rlen = 2000,
                toroidal = FALSE,
                node_clust_fun = kcluster,
                distance_nodes = 'euclidean',
                clust_fun = 'kmeans',
                .options = furrr_options(seed = TRUE,
                                         packages = c('clustTools',
                                                      'somKernels')))

  plan('sequential')

# clustering variances ------

  insert_msg('Clustering variances')

  infil_clust_dev$variance <- infil_clust_dev$algos %>%
    map(var) %>%
    map2_dfr(., names(.),
             ~tibble(method = .y,
                     frac_var = .x$frac_var))

# cross-validation --------

  insert_msg('Clustering cross-validation errors')

  infil_clust_dev$cv <- infil_clust_dev$algos %>%
    map(cv,
        nfolds = 10,
        kNN = 5,
        simple_vote = TRUE,
        resolve_ties = TRUE,
        .parallel = TRUE)

  infil_clust_dev$cv <- infil_clust_dev$cv %>%
    map2_dfr(., names(.),
             ~mutate(.x$summary, method = .y))

# Plotting the results ------

  insert_msg('Plotting the testing results')

  infil_clust_dev$test_results <- left_join(infil_clust_dev$variance,
                                            infil_clust_dev$cv,
                                            by = 'method') %>%
    mutate(cv_accuracy = 1 - mean_error)

  infil_clust_dev$test_plot <- infil_clust_dev$test_results %>%
    select(method, frac_var, cv_accuracy) %>%
    gather(key = 'stat',
           value = 'value',
           frac_var,
           cv_accuracy) %>%
    mutate(method = stri_replace(method, fixed = '_', replacement = ' ')) %>%
    ggplot(aes(x = value,
               y = reorder(method, value),
               fill = stat)) +
    geom_bar(stat = 'identity',
             color = 'black',
             position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c('frac_var' = 'steelblue',
                                 'cv_accuracy' = 'coral3'),
                      labels = c('frac_var' = 'expl. variance',
                                 'cv_accuracy' = 'CV accuracy'),
                      name = '') +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Candidate clustering algorithms',
         subtitle = paste0('min/max normalization, k = 3 clusters, n = ',
                           nrow(infil_clust_dev$analysis_tbl$tcga),
                           ' observations'),
         x = 'statistic value')

# saving an image with the testing results and trained cluster structures -----

  insert_msg('Saving the image')

  save(infil_clust_dev, file = './input data/clust_training.RDa')

# END -----

  insert_tail()
