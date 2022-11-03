# Pathway perturbation analysis with SPIA

  insert_head()

# container list ------

  pathways <- list()

# globals -------

  insert_msg('Globals setup')

  ## regulation vectors

  pathways$reg_vectors <- dge$reg_vectors %>%
    compact %>%
    map(~map(compact(.x), reduce, c)) %>%
    unlist(recursive = FALSE) %>%
    map(function(x) if(length(x) < 50) NULL else x) %>%
    compact

  ## all ID vectors

  pathways$all_vectors <- dge$test_results %>%
    map(~map(.x, filter, !is.na(entrez_id))) %>%
    map(~map(.x, ~unique(.x$entrez_id))) %>%
    unlist(recursive = FALSE) %>%
    map(as.character)

  pathways$all_vectors <- pathways$all_vectors[names(pathways$reg_vectors)]

  ## plot specifications

  pathways$plot_specs <- tibble(plot_name = names(pathways$reg_vectors)) %>%
    mutate(cohort = stri_replace(plot_name, regex = '\\..*', replacement = ''),
           cohort = globals$cohorts[cohort],
           plot_volcano_title = paste('Signaling modulation,', cohort),
           plot_forest_title = paste('Top 20 pathways,', cohort),
           plot_subtitle = ifelse(stri_detect(plot_name, fixed = 'cxcl9'),
                                  'CXCL9 high versus low tumors',
                                  'CCL11 high versus low tumors'))

  ## parallel backend

  plan('multisession')

# Serial analysis -----

  insert_msg('Serial SPIA')

  pathways$spia <- list(de = pathways$reg_vectors,
                        all = pathways$all_vectors) %>%
    future_pmap(spia,
                organism = 'hsa',
                verbose = FALSE,
                .options = furrr_options(seed = TRUE))

  ## some formatting

  pathways$spia <- pathways$spia %>%
    map(as_tibble) %>%
    map(arrange, pG)

# Volcano plots -----

  insert_msg('Volcano plots')

  pathways$volcano_plot <- list(data = pathways$spia ,
                                plot_title = pathways$plot_specs$plot_volcano_title,
                                plot_subtitle = pathways$plot_specs$plot_subtitle) %>%
    pmap(plot_volcano,
         regulation_variable = 'tA',
         p_variable = 'pGFdr',
         label_variable = 'Name',
         signif_level = 0.05,
         regulation_level = 0,
         x_lab = 'Pathway regulation estimate, tA',
         y_lab = expression('-log'[10]*' pFDR'),
         cust_theme = globals$common_theme,
         txt_size = 2.3,
         txt_face = 1,
         top_significant = 20) %>%
    map(~.x +
          geom_vline(xintercept = 0, linetype = 'dashed') +
          theme(legend.position = 'bottom') +
          scale_fill_manual(values = c(upregulated = 'coral3',
                                       downregulated = 'steelblue',
                                       ns = 'gray60'),
                            labels = c(upregulated = 'activated',
                                       downregulated = 'inhibited',
                                       ns = 'ns')))

  ## tag adjustment

  pathways$volcano_plot  <- pathways$volcano_plot %>%
    map(~.x +
          labs(tag = stri_replace(pathways$volcano_plot$labels$tag,
                                  fixed = 'upregulated',
                                  replacement = 'activated'))) %>%
    map(~.x +
          labs(tag = stri_replace(pathways$volcano_plot$labels$tag,
                                  fixed = 'downregulated',
                                  replacement = 'inhibited')))

# Top regulated plots -----

  insert_msg('Top pathway plots')

  pathways$top_pathway_plot <- list(data = pathways$spia,
                                    plot_title = pathways$plot_specs$plot_forest_title,
                                    plot_subtitle = pathways$plot_specs$plot_subtitle) %>%
    pmap(plot_top,
         regulation_variable = 'tA',
         label_variable = 'Name',
         p_variable = 'pGFdr',
         signif_level = 0.05,
         regulation_level = 0,
         top_regulated = 20,
         x_lab = 'Pathway regulation estimate, tA',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_color_manual(values = c(upregulated = 'coral3',
                                        downregulated = 'steelblue',
                                        ns = 'gray60'),
                             labels = c(upregulated = 'activated',
                                        downregulated = 'inhibited',
                                        ns = 'ns'),
                             name = '') +
          guides(fill = 'none'))

# common regulated pathways ------

  insert_msg('Common regulated pathways')

  pathways$cmm_pathways <- pathways$spia[c('tcga.cxcl9',
                                           'cm25ev.cxcl9',
                                           'cm25ni.cxcl9',
                                           'gse73731.cxcl9',
                                           'gse167093.cxcl9',
                                           'emtab1980.cxcl9',
                                           'reca.cxcl9')] %>%
    map(filter, pGFdr < 0.05, tA != 0) %>%
    map(~.x$ID) %>%
    reduce(intersect)

# Bubble plot with the pathway regulation estimates ------

  insert_msg('Common pathway bubble plot')

  ## plotting data

  pathways$bubble_tbl <- pathways$spia[c('tcga.cxcl9',
                                         'cm25ev.cxcl9',
                                         'cm25ni.cxcl9',
                                         'gse73731.cxcl9',
                                         'gse167093.cxcl9',
                                         'emtab1980.cxcl9',
                                         'reca.cxcl9')] %>%
    map(filter, ID %in% pathways$cmm_pathways) %>%
    map2_dfr(., names(.),
             ~mutate(.x, cohort = .y)) %>%
    mutate(cohort = stri_split_fixed(cohort,
                                     pattern = '.',
                                     simplify = TRUE)[, 1],
           Name = stri_replace(Name,
                               regex = '\\s{1}\\(.*\\)',
                               replacement = ''),
           Name = stri_replace(Name,
                               fixed = ' for ',
                               replacement = '\nfor '))

  ## X axis labels

  pathways$bubble_labs <- globals$cohorts %>%
    stri_replace(fixed = ' NIVO', replacement = '\nNIVO') %>%
    stri_replace(fixed = ' EVER', replacement = '\nEVER') %>%
    set_names(names(globals$cohorts))

  ## bubble plot

  pathways$bubble_plot <- pathways$bubble_tbl %>%
    ggplot(aes(x = cohort,
               y = reorder(Name, tA),
               fill = Status,
               size = abs(tA))) +
    geom_point(shape = 21) +
    geom_text(aes(label = signif(tA, 2)),
              size = 2.5,
              hjust = 0.5,
              vjust = -1.1) +
    scale_x_discrete(limits = c('tcga',
                                'cm25ev',
                                'cm25ni',
                                'gse73731',
                                'gse167093',
                                'emtab1980',
                                'reca'),
                     labels = pathways$bubble_labs) +
    scale_fill_manual(values = c(Activated = 'firebrick',
                                 Inhibited = 'steelblue'),
                      labels = function(x) tolower(x),
                      name = 'Pathway status') +
    scale_size(range = c(0.5, 5)) +
    globals$common_theme +
    theme(axis.title = element_blank()) +
    labs(title = 'Modulation of signaling pathways, CXCL9 high vs CXCL9 low tumors',
         subtitle = 'SPIA, common significant pathways')

# END -----

  plan('sequential')

  insert_tail()
