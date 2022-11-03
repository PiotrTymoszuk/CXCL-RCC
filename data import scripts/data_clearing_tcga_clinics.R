# This script clears clinical TCGA data

  insert_head()

# reading the relevant clinical data ------

  insert_msg('Rading the clinical text data')

  tcga_clinics <- c(drugs = './input data/TCGA/nationwidechildrens.org_clinical_drug_kirc.txt',
                    fup = './input data/TCGA/nationwidechildrens.org_clinical_follow_up_v1.0_kirc.txt',
                    clinical = './input data/TCGA/nationwidechildrens.org_clinical_patient_kirc.txt',
                    radiation = './input data/TCGA/nationwidechildrens.org_clinical_radiation_kirc.txt') %>%
    map(read_tsv)

# Clearing the general clinical data -----

  insert_msg('Clearing the general clinical data')

  tcga_clinics$clinical <- tcga_clinics$clinical[c(-1, -2), ] %>%
    transmute(patient_id = bcr_patient_barcode,
           histology = factor(car::recode(histologic_diagnosis,
                                          "'[Discrepancy]' = NA")),
           tumor_grade = factor(car::recode(tumor_grade, "'[Not Available]' = NA")),
           laterality = factor(laterality, c('Right', 'Left', 'Bilateral')),
           collection_type = ifelse(prospective_collection == 'YES',
                                    'prospective',
                                    ifelse(retrospective_collection == 'YES',
                                           'retrospective', NA)),
           collection_type = factor(collection_type),
           sex = factor(car::recode(gender,
                                    "'MALE' = 'Male'; 'FEMALE' = 'Female'")),
           race = factor(car::recode(race,
                                     "'WHITE' = 'White';
                                      'BLACK OR AFRICAN AMERICAN' = 'Black';
                                      'ASIAN' = 'Asian'"),
                         c('White', 'Black', 'Asian')),
           history_other_malignancy = factor(ifelse(history_other_malignancy == 'no',
                                                    'no', 'yes')),
           neoadjuvant = factor(ifelse(history_neoadjuvant_treatment == 'No',
                                       'no', 'yes')),
           age = as.numeric(age_at_initial_pathologic_diagnosis),
           pt_stage = factor(ajcc_tumor_pathologic_pt),
           pn_stage = factor(ajcc_nodes_pathologic_pn),
           pm_stage = factor(ifelse(ajcc_metastasis_pathologic_pm == '[Not Available]',
                                    NA, ajcc_metastasis_pathologic_pm)),
           p_stage = stri_replace(ajcc_pathologic_tumor_stage,
                                  fixed = 'Stage ',
                                  replacement = ''),
           p_stage = factor(p_stage, c('I', 'II', 'III', 'IV')),
           death = as.numeric(car::recode(vital_status,
                                          "'Alive' = 0; 'Dead' = 1")),
           last_contact_days_to = as.numeric(last_contact_days_to),
           death_days_to = as.numeric(death_days_to),
           os_days = ifelse(vital_status == 'Alive',
                            last_contact_days_to,
                            death_days_to),
           tumor_status = car::recode(tumor_status,
                                      "'TUMOR FREE' = 'no'; 'WITH TUMOR' = 'yes'"),
           tumor_status = factor(tumor_status, c('no', 'yes')),
           tumor_death = ifelse(vital_status == 'Dead' & tumor_status == 'yes',
                                '1', '0'),
           tumor_death = as.numeric(tumor_death))

# Clearing the drug data ----

  insert_msg('Clearing the drugs data')

  tcga_clinics$drugs <- tcga_clinics$drugs[c(-1, -2), ] %>%
    transmute(patient_id = bcr_patient_barcode,
              drug_type = car::recode(pharmaceutical_therapy_drug_name,
                                      "'[Unknown]' = NA;
                                   '[Not Available]' = NA;
                                   'Chemo, NOS' = NA;
                                   'Megace' = NA; ## no tumor specific-treatment
                                   '5 FU' = '5-FU';
                                   'Fluorouracil' = '5-FU';
                                   '5-fu' = '5-FU';
                                   'FU7' = '5-FU';
                                   '5-Fluorouracil' = '5-FU';
                                   'fluorouracil' = '5-FU';
                                   '5-fluorouracil' = '5-FU';
                                   '5-Fluorouracil?' = '5-FU';
                                   '5FU' = '5-FU';
                                   '5 FU' = '5-FU';
                                   'gemcitabine' = 'Gemcitabine';
                                   'Gemcitabine Injection' = 'Gemcitabine';
                                   'GEMCITABINE' = 'Gemcitabine';
                                   'Gemcitabine HCL' = 'Gemcitabine';
                                   'gemcitabine HCL' = 'Gemcitabine';
                                   'Gemcitibine' = 'Gemcitabine';
                                   'Gemzar' = 'Gemcitabine';
                                   'gemzar' = 'Gemcitabine';
                                   'CAPECITABINE' = 'Capecitabine';
                                   'Xeloda' = 'Capecitabine';
                                   'Camptosar' = 'Irinotecan';
                                   'Irinotecan Hydrochloride' = 'Irinotecan';
                                   'irinotecan' = 'Irinotecan';
                                   'Leucovorin Calcium' = 'FA';
                                   'Leucovorin calcium' = 'FA';
                                   'Leucovorin' = 'FA';
                                   'leucovorin' = 'FA';
                                   'folinic acid' = 'FA';
                                   'Folinic Acid' = 'FA';
                                   'Eloxatin' = 'Pt';
                                   'Cisplatin' = 'Pt';
                                   'oxaliplatin' = 'Pt';
                                   'cisplatin' = 'Pt';
                                   'Oxaliplatin' = 'Pt';
                                   'Carboplatin' = 'Pt';
                                   'Abraxane' = 'Tax';
                                   'Docetaxel' = 'Tax';
                                   'ABRAXANE' = 'Tax';
                                   'Tarceva' = 'other';
                                   'doxorubicin' = 'other';
                                   'cyclophosphamide' = 'other';
                                   'Dexamethasone' = 'other'"))

# Clearing the radiation data -----

  insert_msg('Clearing the radiation data')

  tcga_clinics$radiation <- tcga_clinics$radiation[c(-1, -2), ] %>%
    transmute(patient_id = bcr_patient_barcode,
              radiation_dose = ifelse(radiation_adjuvant_units == 'cGy',
                                      as.numeric(radiation_total_dose)/100,
                                      as.numeric(radiation_total_dose)),
              radiation_response = car::recode(treatment_best_response,
                                               "'Radiographic Progressive Disease' = 'Progression';
                                       'Stable Disease' = 'Stable';
                                       'Partial Response' = 'Partial response';
                                       'Complete Response' = 'Complete response'") %>%
                factor(c('Progression',
                         'Stable',
                         'Partial response',
                         'Complete response')),
              radiation_response_score = car::recode(treatment_best_response ,
                                                     "'Radiographic Progressive Disease' = 0;
                                             'Stable Disease' = 1;
                                             'Partial Response' = 2;
                                             'Complete Response' = 3") %>%
                as.numeric)

# Clearing the follow-up data -----

  insert_msg('Cleaning the follow-up data')

  tcga_clinics$fup <- tcga_clinics$fup[c(-1, -2), ] %>%
    transmute(patient_id = bcr_patient_barcode,
              relapse_days = as.numeric(new_tumor_event_dx_days_to),
              relapse_fup = ifelse(!is.na(relapse_days), 'yes', 'no'),
              relapse_fup = factor(relapse_fup),
              relapse_days = ifelse(is.na(relapse_days),
                                    as.numeric(last_contact_days_to),
                                    relapse_days))

  ## duplicate handling: the clinical record with
  ## the longest observation time is kept

  tcga_clinics$fup <- tcga_clinics$fup %>%
    dlply(.(patient_id),
          filter,
          relapse_days == max(relapse_days, na.rm = TRUE)) %>%
    map_dfr(~.x[1, ]) %>%
    as_tibble

# Merging the general clinical and fup data and clearing the relapse issue -----

  insert_msg('Merging and clearing the relapse issue')

  tcga_clinics$clinical <- tcga_clinics[c('clinical',
                                          'fup')] %>%
    reduce(left_join,
           by = 'patient_id')

  tcga_clinics$clinical <- tcga_clinics$clinical %>%
    mutate(relapse = as.numeric(relapse_fup) - 1,
           rfs_days = ifelse(relapse == 1,
                             ifelse(is.na(relapse_days),
                                    last_contact_days_to,
                                    relapse_days),
                             last_contact_days_to),
           os_days = ifelse(rfs_days > os_days,
                            rfs_days,
                            os_days),
           os_days = ifelse(is.na(os_days),
                            death_days_to,
                            os_days),
           rfs_days = ifelse(is.na(rfs_days),
                             last_contact_days_to,
                             rfs_days),
           os_days = ifelse(is.na(os_days),
                            last_contact_days_to,
                            os_days),
           rfs_days = ifelse(is.na(relapse_fup), NA, rfs_days),
           rfs_days = ifelse(rfs_days < 0, NA, rfs_days),
           os_days = ifelse(is.na(death), NA, os_days),
           os_days = ifelse(os_days < 0, NA, os_days),
           tumor_surv_days = os_days) %>%
    select(- relapse_fup,
           - relapse_days)

# END -----

  tcga_clinics <- tcga_clinics[c('clinical',
                                 'drugs',
                                 'radiation')]

  insert_tail()
