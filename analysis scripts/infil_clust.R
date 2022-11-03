# Semi-supervised SOM/HCL clustering.
# Cluster assignment in the test cohorts (E-MTAB 1980, GSE73731 and GSE167093)
# is predicted by k-NN label propagation based on the clustering structure
# developed in the TCGA collective

  insert_head()

# container list ------

  infil_clust <- list()

# loading the cache with trained cluster structures ------

  insert_msg('Loading the trained clustering structures')

  if(file.exists('./input data/clust_training.RDa')) {

    load('./input data/clust_training.RDa')

  } else {

    source_all('./analysis scripts/infil_clust_dev.R',
               message = TRUE, crash = TRUE)

  }

# Globals ------

  insert_msg('Globals setup')

  infil_clust$analysis_tbl <- infil_clust_dev$analysis_tbl

  infil_clust$clust_obj$tcga <- infil_clust_dev$algos$KMEANS_Euclidean

  infil_clust$clust_obj$tcga$clust_assignment <-
    infil_clust$clust_obj$tcga$clust_assignment %>%
    mutate(clust_id = paste0('#', clust_id),
           clust_id = factor(clust_id))

  ## scale limits

  infil_clust$scale_limits <- rev(c('T cell CD8+',
                                    'Macrophage M1',
                                    'Macrophage M2',
                                    'T cell regulatory (Tregs)',
                                    'Neutrophil',
                                    'Myeloid dendritic cell',
                                    'Monocyte',
                                    'T cell CD4+ (non-regulatory)',
                                    'NK cell',
                                    'B cell',
                                    'uncharacterized cell'))

# characteristic of the training clustering object ------

  insert_msg('Training clustering object')

  infil_clust$train_plots$som_training <- plot(infil_clust$clust_obj$tcga,
                                               type = 'training',
                                               cust_theme = globals$common_theme)$observation

  infil_clust$train_plots[c('wss', 'silhouette')] <- plot(infil_clust$clust_obj$tcga,
                                                          type = 'diagnostic',
                                                          cust_theme = globals$common_theme)$node

  infil_clust$train_plots$heat_map <- plot(infil_clust$clust_obj$tcga,
                                           type = 'heat_map',
                                           cust_theme = globals$common_theme)$node +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.text.x = element_blank())

# Cluster assignment prediction ------

  insert_msg('Cluster prediction')

  plan('multisession')

  infil_clust$clust_obj[c('emtab1980',
                          'gse73731',
                          'gse167093',
                          'reca',
                          'cm10',
                          'cm25ev',
                          'cm25ni')] <-
    list(newdata = infil_clust$analysis_tbl[c('emtab1980',
                                              'gse73731',
                                              'gse167093',
                                              'reca',
                                              'cm10',
                                              'cm25ev',
                                              'cm25ni')],
         kNN = c(9, 9, 11, 5, 5, 5, 5)) %>%
    pmap(predict,
         object = infil_clust$clust_obj$tcga,
         type = 'propagation',
         simple_vote = TRUE,
         resolve_ties = TRUE)

  plan('sequential')

# Clustering variances ------

  insert_msg('Clustering variances')

  infil_clust$variances <- infil_clust$clust_obj %>%
    map(var) %>%
    map2_dfr(., names(.),
             ~tibble(cohort = globals$cohorts[.y],
                     frac_var = .x$frac_var)) %>%
    mutate(n = map_dbl(infil_clust$analysis_tbl, nrow),
           axis_lab = paste(cohort, n, sep = '\nn = '))

  ## plotting

  infil_clust$variance_plot <- infil_clust$variances %>%
    ggplot(aes(x = frac_var,
               y = reorder(cohort, frac_var))) +
    geom_bar(color = 'black',
             fill = 'steelblue',
             stat = 'identity') +
    geom_text(aes(label = signif(frac_var, 2)),
              size = 2.75,
              hjust = -0.3) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Explained variance',
         subtitle = paste0('training: TCGA, test: ',
                           paste(globals$cohorts[-1],
                                 collapse = ', ')),
         x = 'explained clusetring variance')

# Infiltration heat maps -------

  insert_msg('Infiltration heat maps')

  infil_clust$feature_heat_maps <- list(x_object = infil_clust$clust_obj,
                                        plot_title = globals$cohorts,
                                        plot_subtitle = c('training',
                                                          rep('test', 7))) %>%
    pmap(plot_clust_hm,
         cust_theme = globals$common_theme,
         fill_lab = 'min/max level',
         x_lab = 'tumor sample') %>%
    map(~.x +
          scale_y_discrete(labels = globals$imm_labels,
                           limits = infil_clust$scale_limits) +
          labs(tag = paste0('\n', .x$labels$tag)))

# Tables for plotting ribbon plots, testing for differences in clust. features ------

  insert_msg('Tables for plotting ribbon panels')

  ## tables

  infil_clust$panel_tbl <-
    map2(infil_clust$analysis_tbl,
         infil_clust$clust_obj,
         ~left_join(rownames_to_column(.x, 'patient_id'),
                    set_names(.y$clust_assignment[c('observation', 'clust_id')],
                              c('patient_id', 'clust_id')),
                    by = 'patient_id')) %>%
    map(as_tibble)

  ## testing for differences between the clusters

  infil_clust$test_results <- infil_clust$panel_tbl %>%
    map(~compare_variables(.x,
                           variables = globals$imm_estimates,
                           split_factor = 'clust_id',
                           what = 'test',
                           types = 'kruskal_test',
                           ci = FALSE,
                           pub_styled = TRUE,
                           adj_method = 'BH')) %>%
    map(mutate,
        plot_lab = globals$imm_labels[variable],
        plot_lab = paste(plot_lab, significance, sep = '\n'))

  ## plot labels

  infil_clust$panel_labs <- infil_clust$test_results %>%
    map(~set_names(.x$plot_lab, .x$variable))

# Infiltration radial plots -------

  insert_msg('Infiltration radial plots')

  infil_clust$panel_plots <-
    list(data = infil_clust$panel_tbl,
         plot_title = globals$cohorts,
         plot_subtitle = 'min/max normalized mean \u00B1 2\u00D7SEM') %>%
    pmap(draw_stat_panel,
         variables = globals$imm_estimates,
         split_factor = 'clust_id',
         stat = 'mean',
         err_stat = '2se',
         form = 'line',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(axis.title = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank()) +
          scale_x_continuous(limits = c(0, 1),
                             breaks = seq(0, 1, by = 0.2)) +
          scale_fill_manual(values = globals$clust_colors,
                            name = 'Infiltration cluster') +
          scale_color_manual(values = globals$clust_colors,
                             name = 'Infiltration cluster')) %>%
    map2(., infil_clust$panel_labs,
         ~.x +
           scale_y_discrete(limits = infil_clust$scale_limits,
                            labels = .y) +
           coord_polar(theta = 'y',
                       clip = 'off'))

  for(i in seq(0, 1, by = 0.2)) {

    infil_clust$panel_plots <- infil_clust$panel_plots %>%
      map(~.x +
            annotate('text',
                     label = i,
                     y = 0.5,
                     x = i,
                     color ='gray60',
                     size = 2.75))

  }

  infil_clust$panel_plots <- map2(infil_clust$panel_plots ,
                                  infil_clust$feature_heat_maps,
                                  ~.x + labs(tag = .y$labels$tag))

# END -----

  rm(i)

  insert_tail()

