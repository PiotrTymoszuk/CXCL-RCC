# Clears the clinical ana exxpression data for the Checkmate010 and checkmate25
# cohorts

  insert_head()

# reading the clinical information and expression data ------

  insert_msg('Reading the data')

  cm10$clinical <-
    read_excel('./input data/checkmate/NIHMS1611472-supplement-Table_S1__Clinical_and_immune_phenotype_data_for_the_CheckMate_cohorts.xlsx',
               skip = 1)

  cm10$expression <-
    read_excel('./input data/checkmate/NIHMS1611472-supplement-Table_S4__RNA_expression__normalized_expression_matrix_.xlsx',
               sheet = 'RNA_Expression',
               skip = 1)

# Clearing the clinical information, mRNA available, primary tumor -------

  insert_msg('Cleaning the clinical information')

  cm10$clinical <- cm10$clinical %>%
    transmute(patient_id = SUBJID, ## patient's ID
              cohort = Cohort,
              arm = tolower(Arm),
              arm = factor(arm, c('everolimus', 'nivolumab')),
              sample_id = ifelse(RNA_ID == 'NA', NA, RNA_ID), ## ID of the RNAseq sample
              maf_tumor_id = ifelse(MAF_Tumor_ID == 'NA', NA, MAF_Tumor_ID), ## mutation analysis ID
              maf_normal_id = ifelse(MAF_Normal_ID == 'NA', NA, MAF_Normal_ID), ## mutation analysis ID, normal sample
              cnv_id = ifelse(CNV_ID == 'NA', NA, CNV_ID),
              cd8_if_id = ifelse(CD8_IF_ID == 'NA', NA, CD8_IF_ID), ## immunofluorescence sample ID
              sex = car::recode(Sex, "'M' = 'Male'; 'F' = 'Female'"),
              sex = factor(sex, c('Female', 'Male')),
              age = as.numeric(Age),
              MSKCC_risk = tolower(MSKCC),
              MSKCC_risk = factor(MSKCC_risk, c('favorable', 'intermediate', 'poor')),
              IMDC_risk = tolower(IMDC),
              IMDC_risk = factor(IMDC_risk, c('favorable', 'intermediate', 'poor')),
              sarc_diff = ifelse(Sarc == 0, 'no', 'yes'),
              sarc_diff = factor(sarc_diff, c('no', 'yes')),
              rhab_diff = ifelse(Rhab == 0, 'no', 'yes'),
              rhab_diff = factor(rhab_diff, c('no', 'yes')),
              sarc_rhab_diff = ifelse(Sarc_or_Rhab == 0, 'no', 'yes'),
              sarc_rhab_diff = factor(sarc_rhab_diff, c('no', 'yes')),
              neoadjuvant = tolower(Received_Prior_Therapy),
              neoadjuvant = factor(neoadjuvant, c('no', 'yes')),
              neoadjuvant_cycles = as.numeric(Number_of_Prior_Therapies),
              sample_therapy_start_days = as.numeric(Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy),
              biol_material = SampleType,
              tissue_type = tolower(Tumor_Sample_Primary_or_Metastasis),
              tissue_type = factor(tissue_type, c('primary', 'metastasis')),
              meta_site = tolower(Site_of_Metastasis),
              tumor_delta_vol = as.numeric(Tumor_Shrinkage),
              response_type = factor(ORR, c('PD', 'SD', 'PR', 'CRPR', 'CR')),
              benefit = factor(Benefit, c('NCB', 'ICB', 'CB')),
              rfs_days = as.numeric(PFS) * 30.437,
              relapse = as.numeric(PFS_CNSR),
              os_days = as.numeric(OS) * 30.437,
              death = as.numeric(OS_CNSR),
              purity = as.numeric(Purity),
              ploidy = as.numeric(Ploidy),
              TM_Area = as.numeric(TM_Area),
              TM_CD8 = as.numeric(TM_CD8),
              TM_CD8_Density = as.numeric(TM_CD8_Density),
              TC_Area = as.numeric(TC_Area),
              TC_CD8_Density = as.numeric(TC_CD8_Density),
              TM_TC_Ratio = as.numeric(TM_TC_Ratio),
              ImmunoPhenotype = as.numeric(ImmunoPhenotype),
              TM_CD8_percent = as.numeric(TM_CD8_PERCENT),
              tmb = as.numeric(TMB_Counts)) %>%
    filter(!is.na(sample_id),
           tissue_type == 'primary')

  ## splitting into cohorts

  cm25ev$clinical <- cm10$clinical %>%
    filter(cohort == 'CM-025',
           arm == 'everolimus')

  cm25ni$clinical <- cm10$clinical %>%
    filter(cohort == 'CM-025',
           arm == 'nivolumab')

  cm10$clinical <- cm10$clinical %>%
    filter(cohort == 'CM-010')

  cm10$clinical <- cm10$clinical %>%
    select(- IMDC_risk)

# Annotation table -------

  insert_msg('Annotation table')

  ## annotation table

  cm10$annotation <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = set_names(cm10$expression$gene_name,
                                                               cm10$expression$gene_name),
                                              column = 'ENTREZID',
                                              keytype = 'SYMBOL') %>%
    set_names(c('gene_symbol', 'entrez_id')) %>%
    filter(!is.na(entrez_id),
           !is.na(gene_symbol),
           !duplicated(entrez_id)) %>%
    as_tibble

  cm25ev$annotation <- cm10$annotation

  cm25ni$annotation <- cm10$annotation

# Clearing the expression table ------

  insert_msg('Clearig the expression table')

  cm10$expression <- cm10$expression %>%
    filter(gene_name %in% cm10$annotation$gene_symbol) %>%
    column_to_rownames('gene_name') %>%
    t %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

  cm25ev$expression <- cm10$expression %>%
    filter(sample_id %in% cm25ev$clinical$sample_id)

  cm25ni$expression <- cm10$expression %>%
    filter(sample_id %in% cm25ni$clinical$sample_id)

  cm10$expression <- cm10$expression %>%
    filter(sample_id %in% cm10$clinical$sample_id)

# Merging the expression data set with the clinical data -------

  insert_msg('Merging the clinical and expression data')

  cm10$expression <- left_join(cm10$clinical,
                               cm10$expression,
                               by = 'sample_id')

  cm25ev$expression <- left_join(cm25ev$clinical,
                                 cm25ev$expression,
                                 by = 'sample_id')

  cm25ni$expression <- left_join(cm25ni$clinical,
                                 cm25ni$expression,
                                 by = 'sample_id')

# END -----

  insert_tail()
