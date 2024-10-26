# nsclc_recurrence
Recurrence analyses for early stage non-small cell lung cancer patients

## Installation
```
install.packages(c('tidyverse', 'glue', 'sqldf'))
install.packages(c('gt','gtsummary')
install.packages(c('survminer','ggsurvfit', 'tidycmprsk'))
install.packages('devtools')
devtools::install_github("NightingaleHealth/ggforestplot")

```
## Data preparation and cleaning
* prepare_deid.R
* improve_biomarker.R
* improve_comorbidity.R
* prepare_table_extended.R

## baseline clinical table analyses
* create_clinical_table.R
* wrapup.clinical_table.R

## survival analyses
* surv_funcs.R
* wrapup.coxPH.R

## Competing risk analyses
* wrapup.recurrent_rate.R

## Logistic regression analyses
* feedback.R

