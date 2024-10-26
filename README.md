# nsclc_recurrence
Recurrence analyses for early stage non-small cell lung cancer patients

## Installation
```
install.packages(c('tidyverse', 'gt', 'survminer', 'sqldf'))
install.packages('devtools')
install.packages('gtsummary')
devtools::install_github("NightingaleHealth/ggforestplot")
```
## Data preparation and cleaning
* prepare_deid.R
* improve_biomarker.ipy
* improve_biomarker.R
* improve_comorbidity.R
* prepare_table_extended.R

## baseline clinical table analyses
* wrapup.clinical_table.R

## survival analyses
* wrapup.coxPH.R

## Logistic regression
* feedback.R

## Competing risk
* wrapup_draft.recurrent_rate.R
