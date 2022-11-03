# Plotting of the results of metabolic modeling

  insert_head()

# container ------

  meta_plots <- list()

# globals ------

  insert_msg('Classification of the common reactions')

  meta_plots$react_lex <- read_excel('./input data/reaction_interest.xlsx')

  ## numbers of regulated reactions

  meta_plots$n_total <- meta$models %>%
    map(components, 'reactions') %>%
    map(length)

  meta_plots$n_reactions <- meta$models %>%
    map(components, 'regulation') %>%
    map(filter, p_value < 0.05) %>%
    map(mutate,
        reg_sign = ifelse(fold_reg >1, 'activated', 'inhibited'),
        reg_sign = factor(reg_sign, c('activated', 'inhibited'))) %>%
    map(count, reg_sign)

  ## labels with total and regulated reaction numbers

  meta_plots$n_labs <- meta_plots$n_reactions %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ') %>%
          paste(collapse = ', ')) %>%
    map2(meta_plots$n_total, .,
         paste, sep = ', ') %>%
    map(~paste('total reactions: n =', .x))

# Diagnostic plots: errors -------

  insert_msg('Plots of errors')

  meta_plots$errors <- meta$models %>%
    map(plot,
        type = 'errors',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts[names(meta$models)],
         ~.x +
           labs(title = .y,
                subtitle = 'Fold-regulation and regulation errors'))

# General regulation plots ------

  insert_msg('General regulation plots')

  meta_plots$regulation <- meta$models %>%
    map(plot,
        type = 'regulation',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts[names(meta$models)],
         ~.x +
           theme(panel.grid.major.x = element_blank()) +
           labs(title = .y,
                x = 'Reactions',
                y = expression('log'[2] * ' fold-regulation, CXCL9 high vs low'))) %>%
    map2(., meta_plots$n_labs,
         ~.x +
           labs(subtitle = .y))

# Top regulated reactions ------

  insert_msg('Top regulated reactions')

  meta_plots$top <- meta$models %>%
    map(plot,
        type = 'top',
        n_top = 10,
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts[names(meta$models)],
         ~.x +
           labs(title = .y,
                subtitle = 'General regulation of metabolic reactions'))

# Forest plots for the common regulated genes -------

  insert_msg('Forest plots for the common regulated genes')

  meta_plots$common <- meta$models %>%
    map(plot,
        type = 'forest',
        relevant.reactions = meta$cmm_reactions,
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohorts[names(meta$models)],
         ~.x +
           labs(title = .y,
                subtitle = 'General regulation of metabolic reactions'))

  ## adding the faceting information

  for (i in names(meta_plots$common)) {

    meta_plots$common[[i]]$data <-
      left_join(meta_plots$common[[i]]$data,
                meta_plots$react_lex,
                by = 'react_id')

  }

  meta_plots$common <- meta_plots$common %>%
    map(~.x +
          facet_grid(class ~ .,
                     space = 'free',
                     scales = 'free'))

# END -----

  insert_tail()
