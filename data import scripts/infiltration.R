# Immune infiltration based on expression, QuanTIseq

  insert_head()

# container list ------

  deconv <- list()

# expression matrices, non-log scale ------

  insert_msg('Preparing expression matrices')

  deconv$annotation <- list(tcga = tcga$annotation$gene_symbol,
                            emtab1980 = emtab1980$annotation$gene_symbol,
                            gse73731 = gse73731$annotation$gene_symbol,
                            gse167093 = gse167093$annotation$gene_symbol,
                            reca = reca$annotation$gene_symbol,
                            cm10 = cm10$annotation$gene_symbol,
                            cm25ev = cm25ev$annotation$gene_symbol,
                            cm25ni = cm25ni$annotation$gene_symbol)

  deconv$expression <- list(tcga = tcga$expression %>%
                              filter(tissue_type == 'Tumor'),
                            emtab1980 = emtab1980$expression,
                            gse73731 = gse73731$expression,
                            gse167093 = gse167093$expression %>%
                              filter(tissue_type == 'Tumor'),
                            reca = reca$expression,
                            cm10 = cm10$expression,
                            cm25ev = cm25ev$expression,
                            cm25ni = cm25ni$expression)

  deconv$expression <- map2(deconv$expression,
                            deconv$annotation,
                            ~select(.x, patient_id, any_of(.y))) %>%
    map(column_to_rownames, 'patient_id') %>%
    map(t) %>%
    map(~2^.x)

# deconvolution -----

  insert_msg('Deconvolution')

  plan('multisession')

  deconv$results <- list(gene_expression = deconv$expression,
                         arrays = c(FALSE, TRUE,
                                    TRUE, TRUE,
                                    FALSE, FALSE,
                                    FALSE, FALSE)) %>%
    future_pmap(deconvolute,
                method = 'quantiseq')

  plan('sequential')

# saving the results ------

  insert_msg('saving the results')

  deconv <- deconv$results

  save(deconv, file = './input data/infiltration/quantiseq.RDa')

# END -----

  insert_tail()
