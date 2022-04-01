##############################################################################
# Name of file: 06_model.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: regression model
# Approximate run time: Unknown
##############################################################################

#load library for model
library(mgcv)
library(finalfit)

Location = '/conf/'

setwd(paste0(Location, 'EAVE/GPanalysis/analyses/BA.2-variant'))

output_dir = paste0("./output/TND")
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# Turn off scientific notation
options(scipen = 999)

# Does what it says
create_results_table <- function(model){
  #get ORs & 95%CIs
  result = round(exp(cbind("OR" = coef(model), confint.default(model, level=0.95), digits=3)), digits = 3) %>%
    as.data.frame() %>%
    select(-digits)
  
  names(result) = c('OR', 'lower', 'upper')
  
  result <- mutate(result, variable = rownames(result)) %>%
    select(variable, OR, lower, upper)
  
  result <- mutate( result, CI = paste0( '(', lower, '-', upper, ')') )  %>%
    select(-lower, -upper)
  
  names(result) = c('variable', 'OR', '95% CI')
  
  rownames(result) = 1:nrow(result)
  
  result
}

#load TND dataset
df_tnd <- readRDS('./data/df_tnd.rds')


####################### 0 Exploration ################################################

# Investigations on whether to use s gene positve/negative as a proxy for
# BA.2, B.1.1.529
# The issue is that s gene positive is also delta.
cutoff_date = as.Date('2022-01-06')

df_tnd_post = filter(df_tnd, date_ecoss_specimen >= cutoff_date)

# How good a proxy is s gene for BA.2, B.1.1.529
sgene_variant_table =  as.data.frame.matrix( table(df_tnd_post$sgene_classification, df_tnd_post$variant) )

write.csv(sgene_variant_table, paste0(output_dir, "/sgene_variant_table_tnd.csv"))

# How much of each variant do we have in the cutoff sequencing data
table(df_tnd_post %>% pull(variant)) 


##########################
# Study cohort descriptive
##########################

vacc_variant_table = as.data.frame( table(df_tnd$outcome, df_tnd$vacc_status_2))

write.csv(vacc_variant_table, paste0(output_dir, '/vacc_variant_table.csv'), row.names = FALSE)


#Population characteristics by test for Delta and Delta Plus datasets
dependent <- "outcome"
explanatory <- c("age", "subject_sex", "simd2020_sc_quintile", "EAVE_Smoke", "EAVE_BP",
                 "Q_DIAG_AF", "Q_DIAG_ASTHMA","Q_DIAG_BLOOD_CANCER","Q_DIAG_CCF","Q_DIAG_CEREBRALPALSY",        
                 "Q_DIAG_CHD","Q_DIAG_CIRRHOSIS","Q_DIAG_CONGEN_HD","Q_DIAG_COPD","Q_DIAG_DEMENTIA",         
                 "Q_DIAG_DIABETES_1", "Q_DIAG_DIABETES_2","Q_DIAG_EPILEPSY","Q_DIAG_FRACTURE",           
                 "Q_DIAG_NEURO","Q_DIAG_PARKINSONS","Q_DIAG_PULM_HYPER", "Q_DIAG_PULM_RARE",            
                 "Q_DIAG_PVD","Q_DIAG_RA_SLE", "Q_DIAG_RESP_CANCER","Q_DIAG_SEV_MENT_ILL",               
                 "Q_DIAG_SICKLE_CELL", "Q_DIAG_STROKE","Q_DIAG_VTE","Q_HOME_CAT",             
                 "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL", "n_risk_gps","bmi_impute","HB_residence")


study_cohort <- df_tnd %>% summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)

write.csv(study_cohort, paste0(output_dir, '/study_cohort.csv'), row.names = FALSE)

##########################
#BA.2 variant
##########################


outcome_table <- df_tnd %>% 
  filter(outcome %in% c("negative", "BA.2") ) %>%
  mutate(result = as.factor(result)) %>%
  summary_factorlist('outcome', explanatory, p=TRUE, add_dependent_label=TRUE)

write.csv(outcome_table, paste0(output_dir, '/tnd_outcome_table_BA2.csv'), row.names = FALSE)


# model 1
m1 <- gam(result ~ s(days) + s(age) + vacc_status + subject_sex + simd2020_sc_quintile + n_risk_gps
          , family=binomial, data=df_tnd, subset=outcome %in% c("negative", "BA.2"))

m1_result = create_results_table(m1)

write.csv(m1_result, paste0(output_dir, '/model_BA2_vs1.csv'), row.names = FALSE)

#model 2
m2 <- gam(result ~ s(days) + s(age) + vacc_status_2 + subject_sex + simd2020_sc_quintile + n_risk_gps
          , family=binomial, data=df_tnd, subset=outcome %in% c("negative", "BA.2"))

m2_result = create_results_table(m2)

write.csv(m2_result, paste0(output_dir, '/model_BA2_vs2.csv.csv'), row.names = FALSE)

##########################
# Omicron variant
##########################

outcome_table <- df_tnd %>% 
  filter(outcome %in% c("negative", "B.1.1.529") ) %>%
  mutate(result = as.factor(result)) %>%
  summary_factorlist('outcome', explanatory, p=TRUE, add_dependent_label=TRUE)

write.csv(outcome_table, paste0(output_dir, '/tnd_outcome_table_BA1.csv'), row.names = FALSE)

m3 <- gam(result ~ s(days) + s(age) + vacc_status + subject_sex + simd2020_sc_quintile + n_risk_gps
          , family=binomial, data=df_tnd, subset=outcome %in% c("negative", "B.1.1.529"))

m3_result = create_results_table(m3)

write.csv(m3_result, paste0(output_dir, '/model_B.1.1.529_vs.csv'), row.names = FALSE)

#model 4
m4 <- gam(result ~ s(days) + s(age) + vacc_status_2 + subject_sex + simd2020_sc_quintile + n_risk_gps
          , family=binomial, data=df_tnd, subset=outcome %in% c("negative", "B.1.1.529"))

m4_result = create_results_table(m4)

write.csv(m4_result, paste0(output_dir, '/model_B.1.1.529_vs2.csv'), row.names = FALSE)


###############################
#Pos. tests but not sequenced 
###############################

m5 <- gam(result ~ s(days) + s(age) + vacc_status + subject_sex + simd2020_sc_quintile + n_risk_gps,
          family=binomial, data=df_tnd, subset = result == 0 | result == 1 & 
            !(variant %in% c('B.1.1.529', 'BA.2')))

m5_result = create_results_table(m5)

write.csv(m5_result, paste0(output_dir, '/pos_no_seq_vs1.csv'), row.names = FALSE)


m6 <- gam(result ~ s(days) + s(age) + vacc_status_2 + subject_sex + simd2020_sc_quintile + n_risk_gps,
          family=binomial, data=df_tnd, subset = result == 0 | result == 1 & 
            !(variant %in% c('B.1.1.529', 'BA.2')))

m6_result = create_results_table(m6)

write.csv(m6_result, paste0(output_dir, '/pos_no_seq_vs2.csv'), row.names = FALSE)
