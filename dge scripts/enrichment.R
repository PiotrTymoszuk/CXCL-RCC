# GO and KEGG term enrichment analysis.

  insert_head()

# container list -----

  enrich <- list()

# globals ------

  insert_msg('Globals setup')

  ## gene identifier list,
  ## eliminating entries with less than 50 genes

  enrich$entrez_list <- dge$signif_entrez %>%
    map(unlist, recursive = FALSE) %>%
    unlist(recursive = FALSE) %>%
    map(function(x) if(length(x) < 50) NULL else x) %>%
    compact

  ## plot subtitles and colors

  enrich$plot_specs <- tibble(plot_name = names(enrich$entrez_list))

  enrich$plot_specs$fill_scale <- enrich$plot_specs$plot_name %>%
    map(function(x) if(stri_detect(x, fixed = 'upregulated')) c(significant = 'coral3', ns = 'gray60') else c(significant = 'coral3', ns = 'gray60'))

  enrich$plot_specs <- enrich$plot_specs %>%
    mutate(regulation = ifelse(stri_detect(plot_name, fixed = 'downregulated'),
                               'Downregulated in ',
                               'Upregulated in '),
           gene = ifelse(stri_detect(plot_name, fixed = 'cxcl9'),
                         'CXCL9 high tumors',
                         'CCL11 high tumors'),
           plot_subtitle = paste0(regulation, gene),
           cohort = stri_replace(plot_name, regex = '\\..*', replacement = ''),
           cohort = globals$cohorts[cohort],
           plot_go_title = paste('Top 20 GOs,', cohort),
           plot_kegg_title = paste('Top 20 KEGG pathways, ', cohort))

  ## KEGG database fetch

  enrich$gene_pathway <- read_tsv('./input data/KEGGA input/pathway_links.csv',
                                  col_names = FALSE)  %>%
    set_names(c('GeneID', 'PathwayID')) %>%
    mutate(GeneID = stri_replace(GeneID, fixed = 'hsa:', replacement = ''))

  enrich$pathway_names <- read_tsv('./input data/KEGGA input/pathway_names.csv',
                                   col_names = FALSE) %>%
    set_names(c('PathwayID', 'Description'))

  ## parallel backend

  plan('multisession')

# serial analysis -----

  insert_msg('Serial analysis')

  enrich$go <- enrich$entrez_list %>%
    future_map(~goana(de = unname(.x)),
               .options = furrr_options(seed = TRUE))

  #enrich$kegg <-  enrich$entrez_list %>%
   # future_map(~kegga(de = unname(.x),
    #                  gene.pathway = enrich$gene_pathway,
    #                  pathway.names = enrich$pathway_names),
     #          .options = furrr_options(seed = TRUE))

# Clearing the results ------

  insert_msg('Clearing the results')

  enrich$go <- enrich$go %>%
    map(rownames_to_column, 'ID') %>%
    map(filter, Ont == 'BP') %>%
    map(mutate, p_adjusted = p.adjust(P.DE)) %>%
    map(filter, P.DE < 0.05) %>%
    map(as_tibble)

  enrich$kegg <- enrich$kegg %>%
    map(rownames_to_column, 'ID') %>%
    map(mutate, p_adjusted = p.adjust(P.DE)) %>%
    map(filter, P.DE < 0.05) %>%
    map(as_tibble)

# Plotting the top 20 most significantly enriched GOs and pathways ----

  insert_msg('Enrichment plots')

  enrich$go_plots <- list(data = map(enrich$go, top_n, n = 20, -P.DE),
                          plot_title = enrich$plot_specs$plot_go_title,
                          plot_subtitle = enrich$plot_specs$plot_subtitle,
                          fill_scale = enrich$plot_specs$fill_scale) %>%
    pmap(plot_signifcant,
         p_variable = 'p_adjusted',
         label_variable = 'Term',
         signif_level = 0.05,
         top_significant = 20,
         cust_theme = globals$common_theme,
         x_lab = expression('-log'[10]*' pFDR'))

 # enrich$kegg_plots <- list(data = map(enrich$kegg, top_n, n = 20, -P.DE) %>%
  #                            map(mutate,
   #                               Pathway = stri_replace(Pathway,
    #                                                     fixed = ' - Homo sapiens (human)',
    #                                                     replacement = '')),
     #                       plot_title = enrich$plot_specs$plot_kegg_title,
     #                       plot_subtitle = enrich$plot_specs$plot_subtitle,
      #                      fill_scale = enrich$plot_specs$fill_scale) %>%
  #  pmap(plot_signifcant,
   #      p_variable = 'p_adjusted',
   #      label_variable = 'Pathway',
   #      signif_level = 0.05,
    #     top_significant = 20,
    #     cust_theme = globals$common_theme,
    #     x_lab = expression('-log'[10]*' pFDR'))

# serial analysis, GO CXCL9, no discrimination for up- and downregulated ------

  insert_msg('GO enrichment, CXCL9-regulated transcriptome')

  enrich$de_cxcl9 <- dge$signif_entrez %>%
    map(~.x$cxcl9) %>%
    map(reduce, c)

  enrich$go_cxcl9 <- enrich$de_cxcl9 %>%
    future_map(~goana(de = unname(.x)),
               .options = furrr_options(seed = TRUE))

  enrich$go_cxcl9 <- enrich$go_cxcl9 %>%
    map(rownames_to_column, 'ID') %>%
    map(mutate, p_adjusted = p.adjust(P.DE)) %>%
    map(filter, Ont == 'BP') %>%
    map(as_tibble)

# Common significant GOs/top 20, CXCL9 high tumors -----

  insert_msg('Common significant GOs, CXCL9 high tumors')

  enrich$cmm_go_cxcl9 <- enrich$go_cxcl9[names(enrich$go_cxcl9) != 'cm10'] %>%
    map(filter, p_adjusted < 0.05) %>%
    map(top_n, n = 20, -p_adjusted) %>%
    map(~.x$ID) %>%
    reduce(intersect)

# Plotting the common enriched GO's p values -------

  insert_msg('Common isgnificant GOs: plotting')

  enrich$go_cxcl9_cmm_tbl <- enrich$go_cxcl9 %>%
    map(filter, ID %in% enrich$cmm_go_cxcl9) %>%
    map(arrange, p_adjusted) %>%
    map(mutate,
        term_lab = stri_replace(Term, fixed = ' of ', replacement = '\nof '))

  ## bar plots

  enrich$go_cxcl9_cmm_plots <-
    list(data = enrich$go_cxcl9_cmm_tbl,
         plot_title = globals$cohorts[names(enrich$go_cxcl9_cmm_tbl)]) %>%
    pmap(plot_signifcant,
         p_variable = 'p_adjusted',
         label_variable = 'term_lab',
         x_lab = expression('-log'[10] * ' pFDR'),
         fill_scale = c(significant = 'coral3', ns = 'gray70'),
         cust_theme = globals$common_theme)

# END -----

  plan('sequential')

  insert_tail()
