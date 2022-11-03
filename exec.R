# Executes the pipeline -------

  library(soucer)

  print(source_all(c('import.R',
                     'analysis.R',
                     'dge.R',
                     'metabolism.R',
                     'paper.R'),
                   crash = TRUE,
                   message = TRUE))

  save.image()
