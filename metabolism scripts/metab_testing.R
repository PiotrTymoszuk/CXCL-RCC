# Testing for differentially regulated metabolic reactions, Recon 2 database.

  insert_head()

# container ------

  meta <- list()

# globals -------

  insert_msg('Globals setup')

  ## tables with regulation estimates and their errors

  meta$analysis_tbl <- dge$test_results %>%
    map(~.x$cxcl9) %>%
    map(mutate,
        error = abs(estimate/stat)) %>%
    map(~.x[c('gene_symbol', 'entrez_id', 'estimate', 'error')])

  ## bringing back to the linear scale
  ## errors transformed by the first-differential
  # equation log(2) * SEM(log2-Expression) * 2^log2-Expression

  meta$analysis_tbl <- meta$analysis_tbl %>%
    map(mutate,
        error = log(2) * error * 2^estimate,
        estimate = 2^estimate) %>%
    map(~filter(.x, complete.cases(.x)))

  ## loading the database

  data("Recon2")

# Construction of the geneSBML models -------

  insert_msg('Model construction')

  set.seed(1234)

  meta$models <- meta$analysis_tbl %>%
    map(~build_geneSBML(x = set_names(.x$estimate,
                                      .x$entrez_id),
                        err = set_names(.x$error,
                                        .x$entrez_id),
                        database = Recon2D,
                        or_fun = 'mean',
                        and_fun = 'min',
                        x_default = 1,
                        err_method = 'mc',
                        n_iter = 1010,
                        ci_method = 'perc',
                        .parallel = TRUE))

# Identification of significantly regulated reactions -----

  insert_msg('Significantly regulated reactions')

  meta$signif_reactions <- meta$models %>%
    map(components, 'regulation') %>%
    map(filter, p_value < 0.05) %>%
    map(~.$react_id)

# Significant reactions shared by at least 6 out of 8  study cohorts ------

  insert_msg('Common reactions')

  meta$cmm_sets <- utils::combn(names(meta$signif_reactions),
                                m = 6,
                                simplify = FALSE)

  meta$cmm_reactions <- meta$cmm_sets %>%
    map(~reduce(meta$signif_reactions[.x], intersect)) %>%
    reduce(union)

# END -----

  insert_tail()
