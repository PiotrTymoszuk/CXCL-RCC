# Hypergraphs for pathways of interest (common reactions)
# Rates based on the SMNL model fed with pooled gene expression estimates

  insert_head()

# container ------

  meta_hg <- list()

# globals ----

  insert_msg('Globals')

  meta_hg$exc_regex <-
    '(^o2_)|(^h_)|(^h2o_)|(^h2o2_)|(^nadph_)|(^nadp_)|(^nh4_)|(^accoa_)|(^coa_)|(^co2_)|(^nad_)|(^nadh_)|(^atp_)|(^adp_)|(^hco3)|(^fad_)|(^fadh2_)|(^amp_)|(^ppi_)|(^pi_)'

# Tryptophan metabolism reactions and metabolites of interest ----

  insert_msg('Tryptophan metabolism')

  ## reactions

  meta_hg$trp$reacts <- c('TRPHYDRO2',
                          '5HTRPDOX',
                          '5HLTDL',
                          'SRTNACT',
                          'SRTN23OX',
                          'SRTNMTX',
                          '5HOXINOXDA',
                          'ACSRTNMT',
                          'MELATNOX',
                          'MELATN23DOX',
                          'TRPO2',
                          'FKYNH',
                          'KYN3OX',
                          'LFORKYNHYD',
                          'KYN',
                          'KYNAKGAT')

  ## metabolites

  meta_hg$trp$metabs <- react_to_metab(meta_hg$trp$reacts,
                                       exc_regex = meta_hg$exc_regex) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## hypergraph

  meta_hg$trp$hg <- visualize(meta_pool$model,
                              rate_sep = ': ',
                              suffixes = 'none',
                              signif_type = 'raw',
                              fontsize = 0.6,
                              relevant.species = meta_hg$trp$metabs,
                              relevant.reactions = meta_hg$trp$reacts,
                              signif_digits = 3,
                              layoutType = 'dot',
                              plt.margins = c(20, 100, 20, 150))

# Glycine splitting -----

  insert_msg('Gylcine splitting')

  ## reactions

  meta_hg$gly$reacts <- components(meta_pool$model, 'reactions')

  meta_hg$gly$reacts <-
    meta_hg$gly$reacts[stri_detect(meta_hg$gly$reacts, regex = '^R_GCC.*')]

  ## metabolites

  meta_hg$gly$metabs <-
    react_to_metab(meta_hg$gly$reacts,
                   exc_regex = meta_hg$exc_regex) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ##hypergraph

  meta_hg$gly$hg <- visualize(meta_pool$model,
                              rate_sep = ': ',
                              suffixes = 'none',
                              signif_type = 'raw',
                              fontsize = 0.7,
                              relevant.species = meta_hg$gly$metabs,
                              relevant.reactions = meta_hg$gly$reacts,
                              signif_digits = 3,
                              layoutType = 'dot',
                              plt.margins = c(50, 100, 50, 150))

# Fatty acid oxidation -------

  insert_msg('Fatty acid oxidation')

  ## reactions, restricting to the significant ones

  meta_hg$faox$reacts <- components(meta_pool$model, 'regulation') %>%
    filter(stri_detect(react_id, regex = '(^R_FAOX)|(^R_ACCOAC)|(^R_PPA)'),
           p_value < 0.01) %>%
    .$react_id

  ## metabolites

  meta_hg$faox$metabs <-
    react_to_metab(meta_hg$faox$reacts,
                   exc_regex = meta_hg$exc_regex) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## Forest plot representation

  meta_hg$faox$reg_plot <-
    components(meta_pool$model, 'regulation') %>%
    filter(stri_detect(react_id, regex = '(^R_FAOX)|(^R_ACCOAC)')) %>%
    mutate(fold_reg = log(fold_reg),
           lower_ci = log(lower_ci),
           upper_ci = log(upper_ci))

  meta_hg$faox$reg_plot <- meta_hg$faox$reg_plot %>%
    plot_sign(regulation_variable = 'fold_reg',
              p_variable = 'p_value',
              signif_level = 0.05,
              regulation_level = 0,
              plot_title = 'Regulation of FAOx reactions, CXCL9 high vs low',
              plot_subtitle = paste('all FAOx reactions: n =',
                                    nrow(meta_hg$faox$reg_plot)),
              y_lab = expression('log'[2] * 'fold-regulation. CXCL9 high vs low'),
              x_lab = 'FAOx reaction',
              cust_theme = globals$common_theme)

  meta_hg$faox$reg_plot <-
    meta_hg$faox$reg_plot +
    labs(subtitle = paste(meta_hg$faox$reg_plot$labels$subtitle,
                          meta_hg$faox$reg_plot$labels$tag,
                          sep = ', ')) +
    theme(plot.tag = element_blank()) +
    scale_fill_manual(values = c('firebrick', 'steelblue', 'gray60'),
                      labels = c('activated', 'inhibited', 'ns'),
                      name = 'Regulation status')

  ## hypergraph

  meta_hg$faox$hg <- visualize(meta_pool$model,
                               rate_sep = ': ',
                               suffixes = 'none',
                               signif_type = 'raw',
                               fontsize = 0.5,
                               relevant.species = meta_hg$faox$metabs,
                               relevant.reactions = meta_hg$faox$reacts,
                               signif_digits = 3,
                               layoutType = 'dot',
                               plt.margins = c(50, 200, 50, 200))

# OxPhos ------

  insert_msg('OxPhos')

  ## reactions

  meta_hg$oxphos$reacts <- components(meta_pool$model, 'regulation') %>%
    filter(stri_detect(react_id,
                       regex = '(^R_NADH\\d+)|(^R_CYOOm)|(^R_SUCD\\d+)|(^R_CYOR)|(^R_ATPS4m)')) %>%
    .$react_id

  ## metabolites

  meta_hg$oxphos$metabs <-
    react_to_metab(meta_hg$oxphos$reacts,
                   exc_regex = '(^h_*)|(^h2o_*)|(^o2_*)|(^pi_)') %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## hypergraph

  meta_hg$oxphos$hg <- visualize(meta_pool$model,
                                 rate_sep = ': ',
                                 suffixes = 'none',
                                 signif_type = 'raw',
                                 fontsize = 0.7,
                                 relevant.species = meta_hg$oxphos$metabs,
                                 relevant.reactions = meta_hg$oxphos$reacts,
                                 signif_digits = 3,
                                 layoutType = 'dot',
                                 plt.margins = c(50, 50, 50, 50))

# Aminoacid transport -------

  insert_msg('Aminoacid transport')

  ## reactions

  meta_hg$aatx$reacts <-  meta_pool$cmm_react_reg %>%
    filter(class == 'AATx') %>%
    .$react_id

  ## metabolites

  meta_hg$aatx$metabs <- react_to_metab(meta_hg$aatx$reacts,
                                        annotation_db = reactions,
                                        exc_regex = meta_hg$exc_regex) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## hypergraph

  meta_hg$aatx$hg <- visualize(meta_pool$model,
                               rate_sep = ': ',
                               suffixes = 'none',
                               signif_type = 'raw',
                               fontsize = 0.5,
                               relevant.species = meta_hg$aatx$metabs,
                               relevant.reactions = meta_hg$aatx$reacts,
                               signif_digits = 3,
                               layoutType = 'dot',
                               plt.margins = c(50, 20, 50, 20))

# TCA ------

  insert_msg('TCA')

  ## reactions

  meta_hg$tca$reacts <- c('R_PDHm', 'R_PCm',
                          'R_CSm', 'R_ACONTm', 'R_ICDHxm',
                          'R_AKGDm', 'R_SUCOAS1m', 'R_SUCD1m',
                          'R_FUMm', 'R_MDHm', 'R_AKGMALtm')

  ## metabolites

  meta_hg$tca$metabs <-
    react_to_metab(meta_hg$tca$reacts,
                   annotation_db = reactions,
                   exc_regex = paste(meta_hg$exc_regex,
                                     '|(^akg_c)|(^mal_L_c)|(^gdp_m)|(^gtp_m)')) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## hypergraph

  meta_hg$tca$hg <- visualize(meta_pool$model,
                              rate_sep = ': ',
                              suffixes = 'none',
                              signif_type = 'raw',
                              fontsize = 0.5,
                              relevant.species = meta_hg$tca$metabs,
                              relevant.reactions = meta_hg$tca$reacts,
                              signif_digits = 3,
                              layoutType = 'circo',
                              plt.margins = c(50, 80, 50, 50))


# Glycolysis --------

  insert_msg('Glycolysis')

  ## reactions

  meta_hg$glycolysis$reacts <- c("R_HEX1",
                                 "R_PGI",
                                 "R_PFK",
                                 "R_FBA",
                                 "R_TPI",
                                 "R_GAPD",
                                 "R_PGK",
                                 "R_PGM",
                                 "R_ENO",
                                 "R_PYK",
                                 "R_G6PDH2r",
                                 "R_PGL",
                                 "R_GND",
                                 "R_RPE",
                                 "R_RPI",
                                 "R_TKT1")

  ## metabolites

  meta_hg$glycolysis$metabs <-
    react_to_metab(meta_hg$glycolysis$reacts,
                   annotation_db = reactions,
                   exc_regex = paste(meta_hg$exc_regex)) %>%
    .$metab_id %>%
    unlist %>%
    unique

  ## hypergraph

  meta_hg$glycolysis$hg <-
    visualize(meta_pool$model,
              rate_sep = ': ',
              suffixes = 'none',
              signif_type = 'raw',
              fontsize = 0.6,
              relevant.species = meta_hg$glycolysis$metabs,
              relevant.reactions = meta_hg$glycolysis$reacts,
              signif_digits = 3,
              layoutType = 'dot',
              plt.margins = c(50, 50, 50, 50))

# Saving the hypergraphs on the disc ------

  insert_msg('Saving the hypergraphs on the disc')

  ## tryptophan

  png(filename = './report/hypergraphs/trp.png',
      width = 3000,
      height = 2000,
      res = 600,
      pointsize = 5)

  plot(meta_hg$trp$hg)

  dev.off()

  ## glycine splitting

  png(filename = './report/hypergraphs/gly.png',
      width = 1500,
      height = 1200,
      res = 600,
      pointsize = 5)

  plot(meta_hg$gly$hg)

  dev.off()

  ## oxidative phosphorylation

  png(filename = './report/hypergraphs/oxphos.png',
      width = 2000,
      height = 1200,
      res = 600,
      pointsize = 5)

  plot(meta_hg$oxphos$hg)

  dev.off()

  ## aminoacid transporter

  png(filename = './report/hypergraphs/aatx.png',
      width = 4000,
      height = 1000,
      res = 600,
      pointsize = 5)

  plot(meta_hg$aatx$hg)

  dev.off()

  ## TCA

  png(filename = './report/hypergraphs/tca.png',
      width = 2400,
      height = 2000,
      res = 600,
      pointsize = 6)

  plot(meta_hg$tca$hg)

  dev.off()

  ## Glycolysis

  png(filename = './report/hypergraphs/glycolysis.png',
      width = 2000,
      height = 2000,
      res = 600,
      pointsize = 5)

  plot(meta_hg$glycolysis$hg)

  dev.off()

# END -----

  insert_tail()
