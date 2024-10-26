setwd('~/git/s4-heor/project/JNJ/')
source('create_clinical_table.R') 
library(gt)
library(tidyverse)
#data = readRDS('out/data_1005.4.RDS')
data = readRDS('out/data_1020.1.RDS')

# Merge recurrence type
x = data$RecurrenceType
res = ifelse(x %in% c('local', 'regional', 'Unknown'), 'local/regional', as.character(x))
res_f = factor(str_to_title(res), levels=c('Local/Regional', 'Distant', 'None'))
data$RecurrenceType = res_f
data %>% write_rds('out/data_1225.1.RDS')

create_clinical_table(data %>% mutate(fakeGroup='Total'), predictors=c('Gender', 'AgeGroup', 'Race', 'Ethnicity', 'SmokingStatus',
      'Histology', 'SurgeryType', 'RecurrenceType', 'BMI', 'Charlson', #'ECOG',
      'KRAS', 'EGFR', 'TP53', 
      'STK11', 'BRAF', 'MET', 'PIK3CA', 'ALK', 'MDM2'
      ), 
      group.variable='Stage',
      percentage.flag='columns',
      p.value.flag='Yes',
      clinical.table.title='Clinical Table', 
      clinical.table.file.name='out/clinical_table_1225_byStage.1.html')

create_clinical_table(data %>% mutate(fakeGroup='Total'), 
  predictors=c('Gender', 'AgeGroup', 'Race', 'Ethnicity', 'SmokingStatus',
               'Histology', 'SurgeryType', 'Stage', 'BMI', 'Charlson', #'ECOG',
               'KRAS', 'EGFR', 'TP53', 
               'STK11', 'BRAF', 'MET', 'PIK3CA', 'ALK', 'MDM2'
  ), 
  group.variable='RecurrenceType',
  percentage.flag='columns',
  p.value.flag='Yes',
  clinical.table.title='Clinical Table', 
  clinical.table.file.name='out/clinical_table_1225_byRecurrence.1.html')

# not run
create_clinical_table(data %>% mutate(fakeGroup='Total'), 
                      predictors=c('Gender', 'AgeGroup', 'Race', 'Ethnicity', 'SmokingStatus',
                         'Histology', 'SurgeryType', 
                         'StagePlus', 'RecurrenceType', # 'ECOG'
                         'BMI', # 'Charlson',
                         'KRAS', 'EGFR', 'TP53', 'MET'), 
                      group.variable='fakeGroup',
                      percentage.flag='columns',
                      p.value.flag='No',
                      clinical.table.title='Clinical Table', 
                      clinical.table.file.name='out/clinical_table_1006.2.html')