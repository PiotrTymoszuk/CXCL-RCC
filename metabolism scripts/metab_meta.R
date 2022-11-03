# Metabolic modeling with pooled rate estimates
# obtained by inverse variance method
# Pooled gene expression estimates were caclulated during
# the differential gene expression analysis


  insert_head()

# container ------

  meta_pool <- list()

# calculation of the pooled estimates ------

  meta_pool$cmm_expression <- dge_path$cmm_expression

# building a geneSBML model with the pooled gene expression estimates ------

  insert_msg('Building a geneSBML model')

  meta_pool$model <-
    build_geneSBML(x = set_names(meta_pool$cmm_expression$estimate,
                                 meta_pool$cmm_expression$entrez_id),
                   err = set_names(meta_pool$cmm_expression$error,
                                   meta_pool$cmm_expression$entrez_id),
                   database = Recon2D,
                   or_fun = 'mean',
                   and_fun = 'min',
                   x_default = 1,
                   err_method = 'mc',
                   n_iter = 1010,
                   ci_method = 'perc',
                   .parallel = TRUE)

# Reaction regulation for common regulated genes: Forest plot ------

  insert_msg('Common regulated reactions, Forest plot')

  ## regulation estimates

  meta_pool$cmm_react_reg <- meta_pool$model %>%
    components('regulation') %>%
    filter(react_id %in% meta$cmm_reactions,
           p_value < 0.05) %>%
    left_join(meta_plots$react_lex, by = 'react_id')

  ## Forest plot

  meta_pool$cmm_react_forest <-
    plot_forest(data = meta_pool$cmm_react_reg %>%
                  filter(class != 'other') %>%
                  mutate(fold_reg = log2(fold_reg),
                         lower_ci = log2(lower_ci),
                         upper_ci = log2(upper_ci)),
                regulation_variable = 'fold_reg',
                label_variable = 'react_id',
                p_variable = 'p_value',
                signif_level = 0.05,
                regulation_level = 0,
                lower_ci_variable = 'lower_ci',
                upper_ci_variable = 'upper_ci',
                plot_title = 'Common regulated reactions',
                plot_subtitle = 'Calculated for pooled gene regulation, CXCL9 high vs low',
                x_lab = expression('log'[2] * 'fold-regulation, CXCL9 high vs low'),
                cust_theme = globals$common_theme) +
    facet_grid(class ~ .,
               space = 'free',
               scales = 'free') +
    scale_color_manual(values = c(upregulated = 'firebrick',
                                  downregulated = 'steelblue',
                                  ns = 'gray60'),
                       labels = c(upregulated = 'activated',
                                  downregulated = 'inhibited',
                                  ns = 'ns'),
                       name = 'Regulation status') +
    guides(fill = 'none') +
    scale_y_discrete(labels = function(x) stri_replace(x, regex = '^R_', replacement = ''))

# END ------

  insert_tail()
