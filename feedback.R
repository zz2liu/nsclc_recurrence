require(tidyverse) #dplyr, ggplot2, recode
require(glue)
library(gtsummary)

# dx to surgery time vs NGS tested or not
df = read.csv('out/patient_matrix_v1.01_filtered.csv')
data = readRDS('out/data_1225.1.RDS')
data_ = data %>% 
  inner_join(df, by=c('person_id')) |>
  select(c('person_id', 'AgeYear', 'Gender', 'EGFR', 'KRAS', 'diagnosis_date', 'procedure_date'))
data_ = data_ %>% mutate(dx2sx = as_date(procedure_date) - as_date(diagnosis_date))

boxplot(log2(as.numeric(data_$dx2sx)+1) ~ data_$EGFR)
 # not correlation at all
a = aov(log2(as.numeric(data_$dx2sx)+1) ~ data_$EGFR)
summary(a)
TukeyHSD(a)

boxplot(log2(as.numeric(data_$dx2sx)+1) ~ data_$KRAS)
# not correlation at all
a = aov(log2(as.numeric(data_$dx2sx)+1) ~ data_$KRAS)
summary(a)
TukeyHSD(a)

data_$dx2sx_days = as.numeric(data_$dx2sx)+1
data_$dx2sx_log2days = log2(data_$dx2sx_days)
# not log scale
g = ggplot(data_, aes(x=EGFR, y=dx2sx_days)) + #geom_boxplot()
   geom_violin() #draw_quantiles=c(0.25, 0.5, 0.75))
g + geom_boxplot(width=0.1) + scale_y_continuous(breaks=2^seq(0,12,by=2), trans='log2')
g
g +  scale_y_continuous(trans='log2')
g + scale_y_log2()
g + coord_trans(y='log2')
g +  scale_y_continuous(trans='log2') + coord_trans(y='log2')

g = ggplot(data_, aes(x=KRAS, y=dx2sx_days)) + 
  geom_boxplot()
g +  scale_y_continuous(name='Days from diagnosis to surgery', trans='log2', breaks=2^seq(0,10,by=2), limits=c(2^0, 2^10))

# pairwise changes in days
data_ |> drop_na(dx2sx_days) |> group_by(EGFR) |> summarize(med=median(dx2sx_days, drop_na=T))
data_ |> drop_na(dx2sx_days) |> group_by(KRAS) |> summarize(med=median(dx2sx_days, drop_na=T))
#data_ |> drop_na(dx2sx_log2days) |> group_by(KRAS) |> summarize(meddays=2^mean(dx2sx_log2days, drop_na=T))
# not correlation at all
a = aov(log2(as.numeric(data_$dx2sx)+1) ~ data_$EGFR)
summary(a)
TukeyHSD(a)


# multivariate logistic regression
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
dat = data_recurrence |> filter(RecurrenceType !='None') |> mutate(DistantRecurrence=if_else(RecurrenceType=='Distant', 1, 0))

fit = glm(DistantRecurrence ~ Gender + AgeGroup + Race + Ethnicity + SmokingStatus + 
            Histology + SurgeryType + BMI + Charlson + Stage, # + 
            #SurgeryYear+EGFR+KRAS+MET+TP53, 
          data=dat, family=binomial)
summary(fit)
tbl_regression(fit, exponentiate = T)

