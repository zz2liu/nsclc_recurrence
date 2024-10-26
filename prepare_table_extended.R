# run after analysis.r
setwd('~/git/s4-heor/project/JNJ')
# conda activate rdev
# install.packages(c('tidyverse', 'gt', 'survminer', 'sqldf'))
require(tidyverse) #dplyr, ggplot2

df = read.csv('out/patient_matrix_v1.01_filtered.csv')
pm_biomarker = readRDS('out/pm_biomarker_expanded.RDS')
#pm_ecog = readRDS('out/pm_ecog.RDS')
pm_bmi = readRDS('out/pm_bmi.RDS')
pm_como = readRDS('out/pm_como.RDS')
df = df %>% 
  #left_join(pm_ecog, by='person_id') %>%
  left_join(pm_bmi, by='person_id') %>%
  left_join(pm_como, by='person_id') %>% #head()
  left_join(pm_biomarker, by='person_id') #%>% nrow()

# BMI: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjk_5fz99-BAxXElIkEHaAxCZYQFnoECA0QAw&url=https%3A%2F%2Fwww.cdc.gov%2Fobesity%2Fbasics%2Fadult-defining.html%23%3A~%3Atext%3DIf%2520your%2520BMI%2520is%2520less%2Cfalls%2520within%2520the%2520obesity%2520range.&usg=AOvVaw0GDXeVo190Z84mY-NvVH66&opi=89978449
data = readRDS('out/data_1005.1.RDS')
#df = sqldf("select df.* from df join data using (person_id)")
data = data %>% #count(ecog)
  mutate(BMI=if_else(df$bmi<18.5, 'Underweight <18.5',
                    if_else(df$bmi<25, 'Fit [18.5,25)',
                      if_else(df$bmi<30, 'Overweight [25,30)',
                        'Obese >30')), #else
                    missing='Unknown')) %>% #count(BMI)
  mutate(Charlson=replace_na(df$charlson, '0')) %>%
  mutate(KRAS=replace_na(df$KRAS, 'Unknown'), 
         EGFR=replace_na(df$EGFR, 'Unknown'), 
         TP53=replace_na(df$TP53, 'Unknown'),
         STK11=replace_na(df$STK11, 'Unknown'),
         BRAF=replace_na(df$BRAF, 'Unknown'),
         MET=replace_na(df$MET, 'Unknown'),
         PIK3CA=replace_na(df$PIK3CA, 'Unknown'),
         ALK=replace_na(df$ALK, 'Unknown'),
         MDM2=replace_na(df$MDM2, 'Unknown'))
saveRDS(data, 'out/data_1020.1.RDS')


## not run
data %>% count(Charlson)

## 5-year recurrence rate
data = data %>%
      mutate(Status = as.factor(as.integer(df$Event))) %>% #0 for censored, 1 for recurrence/dead
      mutate(Months = df$Months)


