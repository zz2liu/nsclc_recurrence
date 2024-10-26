#setwd('~/git/s4-heor/project/JNJ')
require(tidyverse) #dplyr, ggplot2
require(sqldf)

## patient_matrix
df = read.csv('out/patient_matrix_v1.01_filtered.csv')
pm = df[,c('person_id', 'procedure_date', 'end_type')]
pm$procedure_date = as.Date(pm$procedure_date)

## ECOG
# ref: https://oncologypro.esmo.org/oncology-in-practice/practice-tools/performance-scales,,,
mapping = read.csv('data/karnofsky_to_ecog.csv', skip=1)
df = read.csv('in/DeID Data/performance_scores.csv')
res = df[,c('person_id', 'reported_date', 'ecog')]
res$reported_date = as.Date(df$reported_date)

filtered = sqldf("select res.*
      from pm
      join res on pm.person_id=res.person_id and res.reported_date<=pm.procedure_date and ecog is not null")
ranked = sqldf("
  select person_id, ecog, reported_date, row_number() over (
    partition by person_id
    order by reported_date desc, ecog nulls last) as rn
  from filtered
  ")
nrow(df)
nrow(filtered)
nrow(pm_ecog)
pm_ecog = sqldf("select person_id, ecog, reported_date ecog_date from ranked where rn=1")
pm_ecog %>% count(ecog)
saveRDS(pm_ecog, 'out/pm_ecog.RDS')

## Charlson
df = read.csv('in/DeID Data/comorbidities_with_score.csv')
df$diagnosis_date = as.Date(df$diagnosis_date)
filtered = sqldf("select df.*, pm.procedure_date
      from pm
      join df on pm.person_id=df.person_id and df.diagnosis_date<=pm.procedure_date")
ranked = sqldf("
  select person_id, charlson_comorbidity_display_name
  , points, row_number() over (
    partition by person_id, charlson_comorbidity_display_name
    order by diagnosis_date desc) as rn
  from filtered
  ")
pm_charlson = sqldf("select person_id, sum(points) charlson 
  from ranked where rn=1 
  group by person_id")
nrow(df)
nrow(filtered)
nrow(pm_charlson)
pm_charlson %>% count(charlson)
saveRDS(pm_charlson, 'out/pm_charlson.RDS')

## Bimarkers
df = read.csv('in/DeID Data/genetics.csv')
### to add negative based on testing name
df %>% count(genetic_test_company,	genetic_test_panel) %>% write.csv('out/panel_name.csv')

### summary biomarker with recurrence state
long = df %>%
  inner_join(pm, by='person_id') %>%
  filter(is_clinically_significant=='S') %>%
  select(gene, person_id, end_type) %>%
  unique() %>%
  count(gene, end_type)
wide = long %>%
  pivot_wider(names_from='end_type', values_from='n', values_fill=0)
wide %>% write_csv('out/gene_recurrence.csv')

wide %>% mutate(NoRecurrence=LastFollowup+Death, Patients=NoRecurrence+Recurrence) %>%
  write_csv('out/gene_recurrence_.csv')

df %>% count(is_clinically_significant, sort=T)
filtered = sqldf("select df.*, pm.procedure_date
      from pm
      join df on pm.person_id=df.person_id and is_clinically_significant='S'")
nrow(filtered)
ranked = sqldf('select distinct person_id, gene from filtered')
nrow(ranked)
p_count = ranked %>% count(gene, sort=T)
#p_count = sqldf('select gene, count(*) patients from ranked group by gene order by patients desc ')

filtered_again = sqldf("select ranked.* from ranked join p_count using (gene) where n>=14.1")
filtered_again %>% count(gene, sort=T)
#sqldf("select count(distinct person_id) from filtered_again")
filtered_again %>% pull(person_id) %>% unique %>% length
### pivot by patient
filtered_again$status = 'S'
pm_biomarker = filtered_again %>% pivot_wider(names_from=gene, values_from=status, values_fill='Unknown')
saveRDS(pm_biomarker, 'out/pm_biomarker.RDS')

## BMI
df = read.csv('in/DeID Data/vitals.csv')
df %>% count(vital_description, sort=T)
processed = df %>% 
  filter(vital_description=='Body mass index (BMI) [Ratio]') %>% #nrow
  mutate(reported_date=as.Date(reported_date)) %>%
  mutate(bmi=as.numeric(vital_value)) %>%
  select(person_id, reported_date, bmi, vital_value, vital_unit)

processed %>% filter(is.na(bmi))
filtered = sqldf("select p.*, pm.procedure_date
      from pm
      join processed p on pm.person_id=p.person_id and reported_date<=procedure_date and bmi is not null")
nrow(filtered)
filtered %>% filter(bmi>20) %>% count(vital_unit)

ranked = sqldf("
  select *, row_number() over (
    partition by person_id
    order by reported_date desc, bmi) as rn
  from filtered
  ")
nrow(ranked)
pm_bmi = ranked %>% filter(rn==1) %>% select(person_id, bmi, reported_date) %>% rename(bmi_date=reported_date)
pm_bmi %>% nrow()
pm_bmi %>% filter(bmi<100) %>% pull(bmi) %>%  hist()
saveRDS(pm_bmi, 'out/pm_bmi.RDS')
