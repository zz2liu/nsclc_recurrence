# run after analysis.r
setwd('~/git/heor/project/JNJ')
# conda activate rdev
#install.packages(c('ggsurvfit', 'tidycmprsk', 'tidyverse'))
require(tidyverse) #dplyr, ggplot2, recode
require(glue)
require(survival) #Surv, overrided?

library(tidycmprsk) #cuminc, tbl_cuminc
library(ggsurvfit) ##ggcuminc

df = read.csv('out/patient_matrix_v1.01_filtered.csv')
data = readRDS('out/data_1005.4.RDS')
#data$end_type = factor(df$end_type, 
#                       levels=c('LastFollowup', 'Recurrence', 'Death')) 

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
data$end_type_plus = as.character(data$end_type)
data$end_type_plus[data$RecurrenceType=='distant']='Metastasis'
data$Status = factor(data$end_type_plus, 
                     levels=c('LastFollowup', 'Recurrence', 'Metastasis','Death'))

## overall
res = run_cuminc(data, strata=1, outcome=c('Recurrence', 'Metastasis', 'Death'))
res$gg
res$tbl

## stratified by stage
res = run_cuminc(data, strata='Stage', outcome=c('Recurrence', 'Metastasis', 'Death'))
res$gg + ylim(0,0.5)
res$tbl





#### not used for draft below


## competing risks: overall
data$Status = data$end_type
res = run_cuminc(data, strata=1, outcome=c('Recurrence', 'Death'))
res$gg
res$tbl

res = run_cuminc(data, strata=1, outcome='Recurrence')
res$gg
res$tbl

ci = cuminc(Surv(Months, Status) ~ 1, data=data)
ci %>%
  ggcuminc(outcome=c('Recurrence','Death')) + 
  labs( x = "Months after surgery", y='Cumulative incidence of recurrence/death') + 
  add_confidence_interval() +
  add_risktable() +
  scale_x_continuous(breaks=seq(0,150,by=12))
ci %>%
  tbl_cuminc(
    outcome=c('Recurrence','Death'),
    times = seq(12, 60, by=12),
    label_header = '**Year {time/12}**') #%>% transpose()

## stratified by stage
ci = cuminc(Surv(Months, Status) ~ Stage, data=data)
ci %>%
  ggcuminc(outcome=c('Recurrence','Death')) + 
  labs( x = "Months after surgery", y='Cumulative incidence of recurrence/death') + 
  add_confidence_interval() +
  add_risktable() +
  scale_x_continuous(breaks=seq(0,150,by=12)) + #, limits=c(0,120)
  coord_cartesian(xlim = c(0, 120))
ci %>%
  tbl_cuminc(
    outcome=c('Recurrence','Death'),
    times = seq(12, 60, by=12), 
    label_header = "**Year {time/12}**")

## recurrence+death 
data = data %>%
  mutate(Status = recode(end_type, Death = 'Recurrence')) #%>% count(Status)
ci =cuminc(Surv(Months, Status) ~ Stage, data=data)   
ci %>%
  ggcuminc() + 
    labs( x = "Months after surgery", y='Cumulative incidence of recurrence+death') + 
    add_confidence_interval() +
    add_risktable() +
    scale_x_continuous(breaks=seq(0,150,by=12)) +
    coord_cartesian(xlim = c(0, 120))
ci %>%
  tbl_cuminc(
    times = seq(12, 60, by=12),
    label_header = '**Year {time/12}**')


## other analyses
# n year recurrence rate and tests
ci_stage %>%
  tbl_cuminc(
    times = 60, 
    label_header = "**{times/12}-year recurrence rate**") %>% 
  add_p()
ci_stage %>%
  tbl_cuminc(
    times = 12, 
    label_header = "**{times/12}-year recurrence rate**") %>% 
  add_p()

## recurrence only rate: death count as cencored
data = data %>% 
  mutate(Status = factor(if_else(end_type %in% c('Death', 'LastFollowup'), 'Cencored', end_type), #'Recurrence'
                               levels=c('Cencored', 'Recurrence')))
data %>% count(Status)
ci = cuminc(Surv(Months, 1) ~ Stage, data=data)
ci %>%
  tbl_cuminc(
    times = 60, 
    label_header = "**{times/12}-year recurrence rate (CI)**")

ci_stage =cuminc(Surv(Months, Status) ~ Stage, data=data)   
ci_stage %>%
  ggcuminc() + 
  labs( x = "Months after surgery", y='Cumulative incidence of recurrence (only)') + 
  add_confidence_interval() +
  add_risktable() +
  scale_x_continuous(breaks=seq(0,150,by=12))

summary(ci_stage)
# n year recurrence rate and tests
ci_stage %>%
  tbl_cuminc(
    #times = 60,
    times = seq(12, 60, by=12),
    label_header = '**Year {time/12}**')
  #add_p()