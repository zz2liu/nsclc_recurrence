#setwd('~/git/s4-heor/project/JNJ')
require(tidyverse) #dplyr, ggplot2
require(sqldf)

## patient_matrix
df = read.csv('out/patient_matrix_v1.01_filtered.csv')
pm = df[,c('person_id', 'procedure_date', 'end_type')]
pm$procedure_date = as.Date(pm$procedure_date)

## Bimarkers
df = read.csv('out/genetic_expanded.csv')
df = df %>% filter(is_clinically_significant!='UNKNOWN')
#keep only one status for a patient
#df$status = factor(df$is_clinically_significant, levels=c('panelNegative','US', 'S'))
mapping = data.frame(is_clinically_significant=c('panelNegative','US', 'S'), rank=c(1, 2, 3))
ranked = sqldf("select person_id, gene, max(rank) rank
      from df
      join mapping using (is_clinically_significant)
      group by person_id, gene")
genetic = sqldf("select person_id, gene, is_clinically_significant
      from ranked
      join mapping using (rank)
      ")

# sort(unique(df$is_clinically_significant))
# not other filtering needed here
#find the top significant genes, make a pivot table on gene

long = genetic %>%
  inner_join(pm, by='person_id') %>%
  select(person_id, gene, is_clinically_significant) %>%
  unique()

res = long %>% filter(is_clinically_significant=='S') %>% count(gene, sort=T)
top_gene = res %>% filter(n>14.1)
long %>% write_csv('out/genetic_expanded_res.csv')

#filtered = long %>% merge(top_gene, by=c('gene')) %>% select(-c('n')) %>% unique()
filtered = long[long$gene %in% top_gene$gene,]
##summary for each gene
summary = filtered %>% 
  count(gene,is_clinically_significant) %>%
  pivot_wider(names_from=is_clinically_significant, values_from=n) %>%
  arrange(desc(S))

filtered = filtered %>%
  mutate(ClinicalSignificance = if_else(is_clinically_significant=='S', 'Significant',
      if_else(is_clinically_significant %in% c('US', 'panelNegative'), 'Negative/US', 'Unknown')))
summary = filtered %>% 
  count(gene,ClinicalSignificance) %>%
  pivot_wider(names_from=ClinicalSignificance, values_from=n) %>%
  arrange(desc(Significant)) 

res = summary %>%
  mutate(Unknown=nrow(data)-rowSums(summary[,-1]))
gt(res, rowname_col='gene')
filtered %>% select(gene, ClinicalSignificance) %>% tbl_summary(by=ClinicalSignificance, sort='Significant')

pm_biomarker = filtered %>% 
  pivot_wider(id_cols=c('person_id'), names_from=gene, 
              values_from=ClinicalSignificance, values_fill='Unknown')
length(unique(pm_biomarker$person_id))
saveRDS(pm_biomarker, 'out/pm_biomarker_expanded.RDS')



## not run
### summary biomarker with recurrence state
long = df %>%
  inner_join(pm, by='person_id') %>%
  filter(is_clinically_significant=='S') %>%
  select(gene, person_id, end_type) %>%
  unique() %>%
  count(gene, end_type)
wide = long %>%
  pivot_wider(names_from='end_type', values_from='n', values_fill=0)
wide %>% write_csv('out/gene_recurrence_expanded.csv')

wide %>% mutate(NoRecurrence=LastFollowup+Death, Patients=NoRecurrence+Recurrence) %>%
  write_csv('out/gene_recurrence_expanded.csv')

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
pm_biomarker = readRDS('out/pm_biomarker.RDS')