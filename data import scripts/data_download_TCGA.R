# This script imports and clears the TCGA lung carcinoma mRNA seq and clinical data
# using TCGA Assembler (https://github.com/compgenome365/TCGA-Assembler-2/blob/master/TCGA-Assembler/)

  insert_head()

# data container -----

  import_tcga <- list()

# sourcing the assembler -----

  insert_msg('Launching the assembler')

  enter_directory('./TCGA Assembler')

  c('Module_A.R',
    'Module_B.R') %>%
    walk(source)

# reading the table with case IDs -----

  insert_msg('Reading the case ID table')

  import_tcga$patient_meta <- read_tsv('clinical.tsv')

  import_tcga$inputPatientIDs <- import_tcga$patient_meta$case_submitter_id

# downloading the clinical and expression information -----

  insert_msg('Downloading the data')

  import_tcga$clinical_path <- DownloadBiospecimenClinicalData(cancerType = 'KIRC',
                                                               saveFolderName = './KIRC data')

  import_tcga$expression_path <- DownloadRNASeqData(cancerType = 'KIRC',
                                                    saveFolderName = './KIRC data',
                                                    assayPlatform = 'gene.normalized_RNAseq',
                                                    inputPatientIDs = import_tcga$case_ids)

# END ----

  go_proj_directory()

  insert_tail()
