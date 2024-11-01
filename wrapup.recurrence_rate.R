require(tidyverse) #dplyr, ggplot2, recode
require(glue)
require(survival) #Surv, overrided?

library(tidycmprsk) #cuminc, tbl_cuminc
library(ggsurvfit) ##ggcuminc

df = read.csv('out/patient_matrix_v1.01_filtered.csv')
data = readRDS('out/data_1005.4.RDS')

run_cuminc = function(data, strata=1, outcome=c('Recurrence', 'Death')) {
  fm = as.formula(glue('Surv(Months, Status) ~ {strata}'))
  ci = cuminc(fm, data=data)
  #browser()
  gg = ci %>%
    ggcuminc(outcome=outcome) + 
    labs( x = "Months after surgery", y='Cumulative incidence') + 
    add_confidence_interval() +
    add_risktable() +
    scale_x_continuous(breaks=seq(0,150,by=12))
  tbl = ci %>%
    tbl_cuminc(
      outcome=outcome,
      times = seq(12, 60, by=12),
      label_header = '**Year {time/12}**') #%>% transpose()
  list(gg=gg, tbl=tbl)
}


## competing risks: plus metastasis
data$end_type_plus = as.character(df$end_type)
data$end_type_plus[df$end_type=='Recurrence']='Local/Regional Recurrence'
data$end_type_plus[df$recurrence_type=='distant']='Metastasis'
data$Status = factor(data$end_type_plus, 
                     levels=c('LastFollowup', 'Local/Regional Recurrence', 'Metastasis','Death'))

## overall
res = run_cuminc(data, strata=1, outcome=c('Local/Regional Recurrence', 'Metastasis', 'Death'))
res$gg
res$tbl

## stratified by stage
res = run_cuminc(data, strata='Stage', outcome=c('Local/Regional Recurrence', 'Metastasis', 'Death'))
res$gg + ylim(0,0.5)
res$tbl
