# Regulation plots for FAOx, OxPhos and TCA

  insert_head()

# container ------

  meta_reg <- list()

# globals: FAOX, OxPhos and TCA reactions ------

  insert_msg('FAOx and OxPhos reactions')

  meta_reg$reacts$faox <- components(meta$models[[1]], 'reactions')

  meta_reg$reacts$faox <-
    meta_reg$reacts$faox[stri_detect(meta_reg$reacts$faox,
                                     regex = '(^R_FAOX)|(^R_ACCOAC)|(^R_PPA)')]

  meta_reg$reacts[c('oxphos',
                    'tca')] <- list(meta_hg$oxphos$reacts,
                                    meta_hg$tca$reacts)

  ## labels for the OxPhos reactions

  meta_reg$oxphos_labs <- c(R_CYOOm2 = 'Complex IV\nCYOOm2',
                            R_CYOOm3 = 'Complex IV\nCYOOm3',
                            R_NADH2_u10m = 'Complex I',
                            R_CYOR_u10m = 'Compex III',
                            R_SUCD1m = 'Complex II',
                            R_ATPS4m = 'Complex V')

  ## labels for the TCA reactions

  meta_reg$tca_labs <- meta_reg$reacts$tca %>%
    map(annotate_bigg, annotation_db = Recon2D) %>%
    unlist %>%
    stri_replace(regex = ',.*$', replacement = '') %>%
    stri_replace(regex = '\\s{1}\\(.*$', replacement = '') %>%
    stri_replace(fixed = '--', replacement = '-') %>%
    stri_replace(fixed = '/malate', replacement = '\nmalate') %>%
    stri_replace(fixed = 'dehydrogenase', replacement = '\ndehydrogenase') %>%
    tolower %>%
    stri_replace(fixed = 'coa', replacement = 'CoA') %>%
    set_names(meta_reg$reacts$tca)

# regulation estimates ------

  insert_msg('Reaction regulation estimates')

  meta_reg$regulation <- meta_reg$reacts %>%
    map(function(react) meta$models %>%
          map(components, 'regulation') %>%
          map(filter, react_id %in% react) %>%
          map(mutate,
              p_value = ifelse(is.na(p_value), 1, p_value),
              reg_sign = ifelse(p_value >= 0.05,
                                'ns',
                                ifelse(fold_reg > 1,
                                       'activated', 'inhibited'))))

  meta_reg$regulation$faox <-
    meta_reg$regulation$faox %>%
    map(mutate,
        fold_reg = log2(fold_reg),
        lower_ci = log2(lower_ci),
        upper_ci = log2(upper_ci),)

  meta_reg$regulation$oxphos <-
    meta_reg$regulation$oxphos %>%
    map(mutate,
        react_id = factor(react_id,
                          c('R_NADH2_u10m',
                            'R_SUCD1m',
                            'R_CYOR_u10m',
                            'R_CYOOm2',
                            'R_CYOOm3',
                            'R_ATPS4m')))

# Regulation plots for fatty acid oxidation ------

  insert_msg('FAOx plotting')

  ## numbers of activated and inhibited reactions

  meta_reg$faox$n_reacts <- meta_reg$regulation$faox %>%
    map(filter, reg_sign != 'ns') %>%
    map(count, reg_sign, .drop = FALSE) %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

  ## regulation plots

  meta_reg$faox$plots <-
    list(data = meta_reg$regulation$faox,
         plot_title = globals$cohorts[names(meta_reg$regulation$faox)],
         plot_subtitle = meta_reg$faox$n_reacts) %>%
    pmap(plot_sign,
         regulation_variable = 'fold_reg',
         p_variable = 'p_value',
         signif_level = 0.05,
         regulation_level = 0,
         x_lab = paste('FAOx reaction, total: n =',
                       nrow(meta_reg$regulation$faox[[1]])),
         y_lab = expression('log'[2] * ' fold-regulation vs CXCL9 low'),
         show_trend = FALSE,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(plot.tag = element_blank()) +
          scale_fill_manual(values = c(upregulated = 'firebrick',
                                       downregulated = 'steelblue',
                                       ns = 'gray60'),
                            labels = c(upregulated = 'activated',
                                       downregulated = 'inhibited',
                                       ns = 'ns'),
                            name = 'Regulation status'))

# Regulation plots for the OxPhos -----

  insert_msg('Regulation plots for OxPhos')

  meta_reg$oxphos$plots <-
    list(data = meta_reg$regulation$oxphos,
         plot_title = globals$cohorts[names(meta_reg$regulation$oxphos)]) %>%
    pmap(function(data, plot_title) data %>%
           ggplot(aes(x = fold_reg,
                       y = react_id,
                       color = reg_sign)) +
           geom_vline(xintercept = 1,
                      linetype = 'dashed') +
           geom_errorbarh(aes(xmin = lower_ci,
                              xmax = upper_ci),
                          height = 0) +
           geom_point(size = 2,
                      shape = 16) +
           scale_y_discrete(limits = rev(c('R_NADH2_u10m',
                                           'R_SUCD1m',
                                           'R_CYOR_u10m',
                                           'R_CYOOm2',
                                           'R_CYOOm3',
                                           'R_ATPS4m')),
                            labels = meta_reg$oxphos_labs) +
           scale_x_continuous(limits = c(0, 1.5)) +
           scale_color_manual(values = c(acitivated = 'firebrick',
                                         inhibited = 'steelblue',
                                         ns = 'gray60'),
                              name = 'Regulation status') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = plot_title,
                x = 'fold-regulation vs Collagen low'))

# Regulation plots for TCA ------

  insert_msg('Regulation plots for TCA')

  ## plotting table

  meta_reg$tca$plot_tbl <- meta_reg$regulation$tca %>%
    map2_dfr(., names(.), ~mutate(.x, cohort = .y)) %>%
    mutate(fold_reg = log2(fold_reg),
           significant = ifelse(reg_sign == 'ns', 'no', 'yes'),
           face = ifelse(reg_sign == 'ns', 'plain', 'bold'))

  ## bubble plot

  meta_reg$tca$plot <- meta_reg$tca$plot_tbl %>%
    ggplot(aes(x = cohort,
               y = react_id,
               fill = reg_sign,
               size = abs(fold_reg))) +
    geom_point(shape = 21) +
    geom_text(aes(label = signif(2^fold_reg, 2),
                  alpha = significant,
                  fontface = .data[['face']]),
              size = 2.5,
              hjust = 0.5,
              vjust = -1.4,
              show.legend = FALSE) +
    scale_y_discrete(limits = rev(meta_reg$reacts$tca),
                     labels = meta_reg$tca_labs) +
    scale_x_discrete(limits = c('tcga',
                                'cm10',
                                'cm25ev',
                                'cm25ni',
                                'gse73731',
                                'gse167093',
                                'emtab1980',
                                'reca'),
                     labels = globals$cohorts %>%
                       stri_replace(fixed = ' EVER', replacement = '\nEVER') %>%
                       stri_replace(fixed = ' NIVO', replacement = '\nNIVO') %>%
                       set_names(names(globals$cohorts))) +
    scale_fill_manual(values = c(activated = 'firebrick',
                                 inhibited = 'steelblue',
                                 ns = 'gray60'),
                      name = 'Regulation status') +
    scale_alpha_manual(values = c(no = 0.5,
                                  yes = 1)) +
    guides(size = 'none') +
    globals$common_theme +
    theme(axis.title = element_blank()) +
    labs(title = 'Regulation of TCA cycle reactions')

# END -----

  insert_tail()
