# Identifies differentially regulated genes (T test, BH adjustment)
# The patients are stratified by the CXCL9 and CCL11 cutoffs for OS
# (TCGA, E-MTAB 1980, RECA and CheckMate)
# or, for the GSE73731 and GSE167031, by median split
# Significance cutoff: at least 1.5-fold regulation, pFDR < 0.05 (BH)

  insert_head()

# container list -----

  dge <- list()

# globals ------

  insert_msg('Globals setup')

  ## analysis tables

  dge$analysis_tbl <- list(tcga = tcga$expression %>%
                             filter(tissue_type == 'Tumor'),
                           emtab1980 = emtab1980$expression,
                           gse73731 = gse73731$expression,
                           gse167093 = gse167093$expression %>%
                             filter(tissue_type == 'Tumor'),
                           reca = reca$expression,
                           cm10 = cm10$expression,
                           cm25ev = cm25ev$expression,
                           cm25ni = cm25ni$expression)

  ## sample stratification by the defined cutoffs for

  dge$analysis_tbl[c('tcga',
                     'emtab1980',
                     'reca',
                     'cm10',
                     'cm25ev',
                     'cm25ni')] <-
    list(x = dge$analysis_tbl[c('tcga',
                                'emtab1980',
                                'reca',
                                'cm10',
                                'cm25ev',
                                'cm25ni')],
         y = surv_cut$cxcl9_cutoff[c('tcga_os',
                                     'emtab1980_os',
                                     'reca_os',
                                     'cm10_os',
                                     'cm25ev_os',
                                     'cm25ni_os')],
         z = surv_cut$ccl11_cutoff[c('tcga_os',
                                     'emtab1980_os',
                                     'reca_os',
                                     'cm10_os',
                                     'cm25ev_os',
                                     'cm25ni_os')]) %>%
    pmap(function(x, y, z) mutate(x,
                                  cxcl9_group = cut(CXCL9,
                                                    c(-Inf, y$cutoff, Inf),
                                                    c('low', 'high')),
                                  ccl11_group = cut(CCL11,
                                                    c(-Inf, z$cutoff, Inf),
                                                    c('low', 'high'))))

  ## sample stratification by median split

  dge$analysis_tbl[c('gse73731', 'gse167093')] <-
    dge$analysis_tbl[c('gse73731', 'gse167093')] %>%
    map(~mutate(.x,
                cxcl9_group = cut(CXCL9,
                                  c(-Inf, median(.x$CXCL9), Inf),
                                  c('low', 'high')),
                ccl11_group = cut(CCL11,
                                  c(-Inf, median(.x$CCL11), Inf),
                                  c('low', 'high'))))

  ## n numbers

  dge$n_numbers <- dge$analysis_tbl %>%
    map(~list(cxcl9 = count(.x, cxcl9_group),
              ccl11 = count(.x, ccl11_group)))

  ## plot tags

  dge$n_tags <- dge$n_numbers %>%
    map(function(cohort) cohort %>%
          map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
          map(paste, collapse = ', ') %>%
          map(~paste0('\n', .x)))

  ## variable lists

  dge$dict <- list(tcga = tcga$annotation,
                   emtab1980 = emtab1980$annotation,
                   gse73731 = gse73731$annotation,
                   gse167093 = gse167093$annotation,
                   reca = reca$annotation,
                   cm10 = cm10$annotation,
                   cm25ev = cm25ev$annotation,
                   cm25ni = cm25ni$annotation) %>%
    map2(., dge$analysis_tbl,
         ~filter(.x, gene_symbol %in% names(.y))) %>%
    map(~filter(.x, !duplicated(gene_symbol)))

# serial analysis ------

  insert_msg('Serial analysis')

  dge$test_results <- list(cohort = dge$analysis_tbl,
                           annotation = dge$dict) %>%
    pmap(function(cohort, annotation) c(cxcl9 = 'cxcl9_group',
                                        ccl11 = 'ccl11_group') %>%
           map(~test_two_groups(data = cohort,
                                split_fct = .x,
                                variables = annotation$gene_symbol,
                                type = 't',
                                adj_method = 'BH',
                                .parallel = TRUE)))

# Formatting the results -------

  insert_msg('Formatting the results')

  dge$test_results <- list(cohort = dge$test_results,
                           dict =  list(tcga$annotation,
                                        emtab1980$annotation,
                                        gse73731$annotation,
                                        gse167093$annotation,
                                        reca$annotation,
                                        cm10$annotation,
                                        cm25ev$annotation,
                                        cm25ni$annotation)) %>%
    pmap(function(cohort, dict) cohort %>%
           map(mutate,
               regulation = ifelse(significant == 'no',
                                   'ns',
                                   ifelse(estimate >= log2(1.5),
                                          'upregulated',
                                          ifelse(estimate <= -log2(1.5),
                                                 'downregulated',
                                                 'ns'))),
               gene_symbol = response,
               entrez_id = translate_gene(response,
                                          dictionary = dict)))

# Vectors of significant genes and regulation -----

  insert_msg('Signif gene vectors')

  ## symbols of significantly regulated genes

  dge$signif_genes <- dge$test_results %>%
    map(~map(.x, filter, regulation != 'ns')) %>%
    map(~map(.x, ~dlply(.x, .(regulation), function(x) x$gene_symbol)))

  ## Entrez IDs of significantly regulated genes, used later
  ## in GO enrichment analysis

  dge$signif_entrez <- dge$test_results %>%
    map(~map(.x, filter, regulation != 'ns', !is.na(entrez_id))) %>%
    map(~map(.x, ~dlply(.x, .(regulation), function(x) x$entrez_id)))

  ## log2-fold regulation vectors named with Entrez IDs used later
  ## in SPIA analysis

  dge$reg_vectors <- dge$test_results %>%
    map(~map(.x, filter, regulation != 'ns', !is.na(entrez_id))) %>%
    map(~map(.x,
             ~dlply(.x,
                    .(regulation),
                    function(x) set_names(x$estimate,
                                          x$entrez_id))))

# Numbers of differentially regulated features -----

  insert_msg('Numbers of diferentially regulated genes')

  dge$n_genes <- c(cxcl9 = 'cxcl9',
                   ccl11 = 'ccl11') %>%
    map(function(gene) dge$signif_entrez %>%
          map(~.x[[gene]]) %>%
          map2_dfr(., names(.),
                   ~tibble(cohort = globals$cohorts[.y],
                           downregulated = length(.x$downregulated),
                           upregulated = length(.x$upregulated)))) %>%
    map2_dfr(., c('CXCL9', 'CCL11'),
             ~mutate(.x, expression_strata = .y))


# Volcano plots -----

  insert_msg('Volcano plots')

  dge$volcano_plots <- set_names(names(dge$test_results),
                                 names(dge$test_results)) %>%
    map(function(cohort) list(data = dge$test_results[[cohort]],
                              plot_title = c('CXCL9: differential gene expression',
                                             'CCL11: differential gene expression') %>%
                                paste(globals$cohorts[[cohort]], sep = ', '),
                              plot_subtitle = list(expression('CXCL9'^{hi}*' vs CXCL9'^{lo}*' tumors'),
                                                   expression('CCL11'^{hi}*' vs CCL11'^{lo}*' tumors'))) %>%
          pmap(plot_volcano,
               regulation_variable = 'estimate',
               p_variable = 'p_adjusted',
               signif_level = 0.05,
               regulation_level = log2(1.5),
               x_lab = expression('log'[2]*' fold regulation'),
               y_lab = expression('-log'[10]*' pFDR'),
               top_significant = 20,
               label_variable = 'gene_symbol',
               txt_size = 2.3,
               txt_face = 3,
               cust_theme = globals$common_theme,
               fill_scale = c(upregulated = 'coral3',
                              downregulated = 'steelblue',
                              ns = 'gray60')) %>%
          map(~.x +
                theme(legend.position = 'bottom')) %>%
          map2(., dge$n_tags[[cohort]],
               ~.x +
                 labs(tag = paste0(.x$labels$tag, .y))))

# Top 20 regulated genes ----

  insert_msg('Top regulated genes')

  dge$top_gene_plots <- set_names(names(dge$test_results),
                                  names(dge$test_results)) %>%
    map(function(cohort) list(data = dge$test_results[[cohort]],
                              plot_title = c('CXCL9: top 20 genes',
                                             'CCL11: top 20 genes') %>%
                                paste(globals$cohorts[[cohort]], sep = ', '),
                              plot_subtitle = list(expression('CXCL9'^{hi}*' vs CXCL9'^{lo}*' tumors'),
                                                   expression('CCL11'^{hi}*' vs CCL11'^{lo}*' tumors')),
                              plot_tag = dge$n_tags[[cohort]]) %>%
          pmap(plot_top,
               regulation_variable = 'estimate',
               label_variable = 'gene_symbol',
               p_variable = 'p_adjusted',
               signif_level = 0.05,
               regulation_level = 0,
               lower_ci_variable = 'lower_ci',
               upper_ci_variable = 'upper_ci',
               top_regulated = 20,
               x_lab = expression('log'[2]*' fold regulation'),
               cust_theme = globals$common_theme,
               fill_scale = c(upregulated = 'coral3',
                              downregulated = 'steelblue',
                              ns = 'gray60')) %>%
          map(~.x + theme(axis.text.y = element_text(face = 'italic'))))

# saving the results ------

  insert_msg('Saving the results')

  save(dge, file = './input data/dge.RData')

# END -----

  insert_tail()
