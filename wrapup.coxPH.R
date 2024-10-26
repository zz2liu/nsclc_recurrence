library(tidyverse)
library(glue)
library(gtsummary) #tbl_survfit, tbl_regression
library(survival) #coxph, survfit, Surv
library(survminer) #surv_fit, ggforest, ggsurvplot
#library(ggsurvfit)
source('surv_funcs.R')

#setwd('~/git/s4-heor/project/JNJ')
## extended data
df = read.csv('out/patient_matrix_v1.01_filtered.csv')
## recurrence type analysis
table(df$recurrence_type)
distant: 148
anyrecurrence: 333 (148   +   136   +    40    +    9 )


#data = readRDS('out/data_1020.1.RDS')
data = readRDS('out/data_1225.1.RDS')
stopifnot(all(data$person_id==df$person_id))

# genetic_factor = function(gene) {
#   factor(gene, levels=c('Unknown', 'S'))
# }
data_recurrence = data %>% 
  mutate(
     AgeGroup = as.factor(AgeGroup),
     Gender=as.factor(Gender),
     Ethnicity = factor(Ethnicity, levels=c('Not Hispanic or Latino', 'Hispanic or Latino', 'Unknown')),
     SmokingStatus = as.factor(SmokingStatus),
     Stage = as.factor(Stage),
     #ECOG = factor(ECOG, levels=c('0', '1', '>=2', 'Unknown')),
     BMI=as.factor(BMI),
     Charlson=as.factor(Charlson),
     EGFR=as.factor(EGFR),
     KRAS=as.factor(KRAS),
     MET=as.factor(MET),
     TP53=as.factor(TP53),
     Status = if_else(df$end_type=='LastFollowup', 0, 1), #Recurrence/meta/death
     Months = difftime(df$end_date,df$procedure_date, units='days')/30
     )
saveRDS(data_recurrence, 'out/data_recurrenc_1225.1.RDS')

median(data_recurrence$Months)

## Metastasis
#df = read.csv('out/patient_matrix_v1.01_filtered.csv')
#data = readRDS('out/data_1005.4.RDS')
data_metastasis = data_recurrence %>%
  mutate(Status =if_else(df$end_type=='Death'|.data$RecurrenceType=='Distant', 1, 0) )


# merge surgery type
# x = data$SurgeryType
# res = ifelse(x %in% c('Segmentectomy', 'Pneumonectomy'), 'Segmentectomy/Pneumonectomy', as.character(x))
# res_f = factor(str_to_title(res), levels=c('Local/Regional', 'Distant', 'None'))
# data$RecurrenceType = res_f
# data %>% write_rds('out/data_1225.1.RDS')

## CoxPH
coxph_plot = function(data, add_features="") {
  fit.coxph = coxph(Surv(Months, Status) ~ AgeGroup + Gender + Race + Ethnicity  +
                      SmokingStatus + Histology + StagePlus + SurgeryType +
                      BMI + Charlson + 
                      KRAS + EGFR + TP53 #+ STK11 #+ BRAF #+ MET #+ #+ BRAF  # + MDM2
                      , data = data)
  #fit.coxph %>%
  #  tbl_regression(exp=TRUE)
  ggforest(fit.coxph, data=data)
}

coxph_plot(data_recurrence |> 
             filter(StagePlus %in% c('IA', 'IB') & SurgeryType != 'Pneumonectomy') |>
             mutate(SurgeryType=factor(SurgeryType, levels=c('Lobectomy', 'Wedge Resection', 'Segmentectomy'))))
coxph_plot(data_metastasis |> 
             filter(StagePlus %in% c('IA', 'IB') & SurgeryType != 'Pneumonectomy') |>
             mutate(SurgeryType=factor(SurgeryType, levels=c('Lobectomy', 'Wedge Resection', 'Segmentectomy'))))

## KM function
## testing here
survdiff(Surv(time = Months, event = Status) ~ StagePlus, 
          data=data_metastasis_IAB) #p0.01
survdiff(Surv(time = Months, event = Status) ~ StagePlus, 
         data=data_metastasis %>% filter(StagePlus %in% c('IIA', 'IIB'))) 

## EGFR
table(data_recurrence$EGFR)
km_plot(data_recurrence, 'EGFR',
        out_prefix='out/RFS_by_EGFR', width=8)
survdiff(Surv(time = Months, event = Status) ~ EGFR, 
         data=data_recurrence %>% filter(EGFR %in% c('Significant', 'Negative/US'))) #p0.01

km_plot(data_metastasis, 'EGFR',
        out_prefix='out/240129_MFS_by_EGFR', xlim=c(0, 150), width=8)
survdiff(Surv(time = Months, event = Status) ~ EGFR, 
         data=data_metastasis %>% filter(EGFR %in% c('Significant', 'Negative/US'))) #p0.08
#p0.08

km_plot(data_recurrence, 'KRAS',
        out_prefix='out/RFS_by_KRAS', width=8)
survdiff(Surv(time = Months, event = Status) ~ KRAS, 
         data=data_recurrence %>% filter(KRAS %in% c('Significant', 'Negative/US'))) #p0.01

km_plot(data_metastasis, 'KRAS',
        out_prefix='out/MFS_by_KRAS', width=8)
survdiff(Surv(time = Months, event = Status) ~ KRAS, 
         data=data_metastasis %>% filter(KRAS %in% c('Significant', 'Negative/US'))) #p0.01

## Global recurrence
km_plot(data_recurrence, '1', pval=F, out_prefix='out/RFS_global')
## statify by stage
km_plot(data_recurrence, 'Stage',
        out_prefix='out/RFS_by_Stage', width=8)
ggsurvfit_plot(data_recurrence, 'Stage',
        out_prefix='out/RFS_by_Stage_ggsurvfit',
        width=7, height=7)
## statify by recurrence type
km_plot(data_recurrence, 'RecurrenceType',
        out_prefix='out/RFS_by_recurrence')
## stratify by surgery type
km_plot(data_recurrence, 'SurgeryType')

## filtered and stratified by surgery type
data_recurrence_surg = data_recurrence[
    data_recurrence$SurgeryType %in% c('Wedge Resection', 'Lobectomy'),]
km_plot(data_recurrence_surg, 'SurgeryType',
        out_prefix='out/RFS_wedge_lobe')

## limited and stratified by substage
data_recurrence_ss = data_recurrence %>% filter(!(StagePlus %in% c('I', 'II', 'III')))
km_plot(data_recurrence_ss, 'StagePlus',
        out_prefix='out/RFS_by_StagePlus', width=9)

### limited and stratified by 1A 1B only
data_recurrence_IAB = data_recurrence %>% filter(StagePlus %in% c('IA', 'IB'))
km_plot(data_recurrence_IAB, 'StagePlus',
        out_prefix='out/RFS_by_stageIAB', width=8)

# Metastasis Survival
## Global 
km_plot(data_metastasis, '1', pval=F, 
        out_prefix='out/MFS_global')
## statify by stage
km_plot(data_metastasis, 'Stage',
        out_prefix='out/MFS_by_Stage', width=8)
ggsurvfit_plot(data_metastasis, 'Stage',
               out_prefix='out/MFS_by_Stage_ggsurvfit',
               width=7, height=7)
## statify by metastasis type
km_plot(data_metastasis, 'RecurrenceType',
        out_prefix='out/MFS_by_RecurrenceType')  #issue here solved
## stratify by surgery type
km_plot(data_metastasis, 'SurgeryType')

## filtered and stratified by surgery type
data_metastasis_surg = data_metastasis[
  data_metastasis$SurgeryType %in% c('Wedge Resection', 'Lobectomy'),]
km_plot(data_metastasis_surg, 'SurgeryType',
        out_prefix='out/MFS_wedge_lobe')

## limited and stratified by substage
data_metastasis_ss = data_metastasis %>% filter(!(StagePlus %in% c('I', 'II', 'III')))
km_plot(data_metastasis_ss, 'StagePlus',
        out_prefix='out/MFS_by_stagePlus', width=9) #make sense at all??

### limited and stratified by 1A 1B only
data_metastasis_IAB = data_metastasis %>% filter(StagePlus %in% c('IA', 'IB'))
km_plot(data_metastasis_IAB, 'StagePlus',
        out_prefix='out/MFS_by_stageIAB', width=8)




###Old KM plots
## KM curves
fit = with(data, survfit(Surv(time = Months, event = Status) ~ 1))
ggsurvplot(fit, data=data
           , break.x.by=12
           , pval=T
           , risk.table='abs_pct', fontsize=3
           , surv.median.line='hv')

tbl_survfit(fit,
            probs = 0.5,
            label_header = "**Median survival (95% CI)**")
tbl_survfit(fit,
            times = seq(12, 60, by=12),
            label_header = '**Year {time/12}**')
            )

### KM curves, stratified by stage
fit = with(data, survfit(Surv(time = Months, event = Status) ~ Stage))
ggsurvplot(fit, data=data
               , break.x.by=12
               , pval=T
               , risk.table='abs_pct', fontsize=3
               , surv.median.line='hv'
               , xlim=c(0, 120))


tbl_survfit(fit,
            probs = 0.5,
            label_header = "**Median survival (95% CI)**"
)
tbl_survfit(fit,
            times = seq(12, 60, by=12),
            label_header = '**Year {time/12}**')
)
survdiff(Surv(time = Months, event = Status) ~ Stage, 
         data=data %>% filter(Stage %in% c('I', 'II')))
survdiff(Surv(time = Months, event = Status) ~ Stage, 
         data=data %>% filter(Stage %in% c('III', 'II')))

### stratified by substages
data_ = data %>% filter(!(StagePlus %in% c('I', 'II', 'III')))
fit = with(data_ , survfit(Surv(time = Months, event = Status) ~ StagePlus))
ggsurvplot(fit, data=data_
           , break.x.by=12
           , pval=T
           , risk.table='absolute', fontsize=3
           , surv.median.line='hv'
           , xlim=c(0, 120))
dim(data_)
dim(data)

tbl_survfit(fit,
            probs = 0.5,
            label_header = "**Median survival (95% CI)**"
)
tbl_survfit(fit,
            times = seq(12, 60, by=12),
            label_header = '**Year {time/12}**',
            title='Survival rate at each year')


### compare IA vs IB
fit = with(data %>% filter(StagePlus %in% c('IA', 'IB')), 
           survfit(Surv(time = Months, event = Status) ~ StagePlus))
ggsurvplot(fit, data=data
           , break.x.by=12
           , pval=T
           #, risk.table='absolute', fontsize=3
           , surv.median.line='hv')
survdiff(Surv(time = Months, event = Status) ~ StagePlus, 
         data=data %>% filter(StagePlus %in% c('IA', 'IB')))
survdiff(Surv(time = Months, event = Status) ~ StagePlus, 
         data=data %>% filter(StagePlus %in% c('IIA', 'IIB')))