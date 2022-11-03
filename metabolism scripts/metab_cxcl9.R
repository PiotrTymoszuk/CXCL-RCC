# comparing the gene signatures of FAOx, OxPhos, TCA and TRP catabolism

  insert_head()

# container ------

  meta_cxcl9 <- list()

# some globals ------

  insert_msg('Globals setup')

  meta_cxcl9$variables <-
    names(meta_gene$gene_symbols$tcga)

  meta_cxcl9$labels <- c('FAOx geneSign',
                         'OxPhos geneSign',
                         'TCA geneSign',
                         'TRP geneSign')

  plan('multisession')

# descriptive stats ------

  insert_msg('Descriptive stats')

  meta_cxcl9$stats <- meta_gene$score_tbl %>%
    map(explore,
        variables = meta_cxcl9$variables,
        split_factor = 'cxcl9_group',
        what = 'table',
        pub_styled = TRUE)

# Testing for differences, Mann-Whitney test ------

  insert_msg('Testing')

  meta_cxcl9$test <- meta_gene$score_tbl %>%
    future_map(compare_variables,
               variables = meta_cxcl9$variables,
               split_factor = 'cxcl9_group',
               what = 'eff_size',
               types = 'wilcoxon_r',
               ci = FALSE,
               exact = FALSE,
               pub_styled = TRUE,
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))

# Violin plots -------

  insert_msg('Violin plots')

  meta_cxcl9$plots <-
    list(x = meta_gene$score_tbl,
         y = meta_cxcl9$test,
         z = globals$cohorts[names(meta_gene$score_tbl)]) %>%
    pmap(function(x, y, z) list(variable = meta_cxcl9$variables,
                                plot_title = paste(meta_cxcl9$labels,
                                                   z,
                                                   sep = ', '),
                                plot_subtitle = y$plot_cap) %>%
           pmap(plot_variable,
                x,
                split_factor = 'cxcl9_group',
                type = 'violin',
                y_lab = 'GSVA score',
                x_lab = 'CXCL9 strata',
                cust_theme = globals$common_theme) %>%
           map(~.x +
                 labs(tag = .x$labels$tag %>%
                        stri_replace_all(fixed = '\n', replacement = ', ') %>%
                        paste0('\n', .)) +
                 scale_fill_manual(values = c('steelblue', 'firebrick'))) %>%
           set_names(meta_cxcl9$variables))

# END -----

  plan('sequential')

  insert_tail()
