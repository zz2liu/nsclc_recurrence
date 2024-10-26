setwd('~/git/s4-heor/project/JNJ')
require(tidyverse) #dplyr, ggplot2
require(sqldf)

## patient_matrix
df = read.csv('out/patient_matrix_v1.01_filtered.csv')
pm = df[,c('person_id', 'procedure_date', 'end_type')]
pm$procedure_date = as.Date(pm$procedure_date)

## icd
icd_history = read.csv('in/DeID Data/icd_history.csv')
icd_history = icd_history %>%
  select(c('person_id', 'icd_date', 'icd_code', 'icd_version')) %>%
  unique() %>%
  mutate(icd_date=as.Date(icd_date))
#icd before procedure for each person
icd_valid = sqldf("select person_id, icd_code, icd_version
                  from icd_history
                  join pm using (person_id)
                  where procedure_date>=icd_date")
icd_valid %>% write_csv('out/icd_valid.csv')

# python comorbidity_workflow.py convert -i ~/git/s4-heor/project/JNJ/out/icd_valid.csv -o ~/git/s4-heor/project/JNJ/out/icd_converted.csv
#normalize the icd code
icd_converted = read_csv('out/icd_converted.csv')
icd_converted %>%
  mutate(icd=str_replace(icd, '^(.{3})(\\d)', '\\1.\\2')) %>%
  write_csv('out/icd_converted_normal.csv')

# python comorbidity_workflow.py comorbidity -i ~/git/s4-heor/project/JNJ/out/icd_converted_normal.csv -o ~/git/s4-heor/project/JNJ/out/comorbidity.csv -s charlson
como = read_csv('out/comorbidity.csv')
como %>% count(wscore)
como %>% count(score)
pm_como = sqldf("select person_id, score charlson_score
                , case score when 0 then '0' when 1 then '1' else '2+'
                  end charlson
                from como")
pm_como %>% write_csv('out/pm_como.csv')
pm_como %>% saveRDS('out/pm_como.RDS')
