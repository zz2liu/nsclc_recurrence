# Function responsible for creating clinical table (stratified by group.variable)
# over a series of predictor variables (predictors)
create_clinical_table = function(df, predictors, group.variable, percentage.flag, p.value.flag, clinical.table.title, 
         clinical.table.file.name, ...){
  
  for(i in 1:length(predictors)){
    # Patient Total by Group and Predictor Strata Levels
    print(i)
    clinical.summary = df %>%
      select(predictors, group.variable) %>%
      group_by(across({{group.variable}})) %>%
      count(.data[[predictors[i]]]) %>%
      mutate(Variable = predictors[i])
    
    clinical.summary.count.table = clinical.summary %>%
      select(-Variable) %>%
      spread(group.variable, n)
    
    clinical.summary.count.table[is.na(clinical.summary.count.table)] = 0
    
    # Produce Column Percentages based on Group Variable Stata Levels
    if(percentage.flag == "columns"){
      
      # Patient Totals by Group Variable Strata Levels
      clinical.subgroup.total = df %>%
        group_by(across({{group.variable}})) %>%
        summarise(Subgroup.Count = n()) %>%
        ungroup()
      
      clinical.summary = clinical.summary %>%
        inner_join(clinical.subgroup.total, by = group.variable) %>%
        mutate(Percentage = round(n / Subgroup.Count * 100, digits = 0)) %>%
        mutate(Count.String = paste(n, " (", Percentage, "%)", sep = "")) %>%
        select(-n) %>%
        rename(n = Count.String)
      
      # Produce Row Percentages based on Predictor Variable Stata Levels 
    }else{
      # Patient Totals by Predictor Variable Strata Levels
      clinical.subgroup.total = df %>%
        group_by(across(predictors[i])) %>%
        summarise(Subgroup.Count = n()) %>%
        ungroup()
      
      clinical.summary = clinical.summary %>%
        inner_join(clinical.subgroup.total, by = predictors[i]) %>%
        mutate(Percentage = round(n / Subgroup.Count * 100, digits = 0)) %>%
        mutate(Count.String = paste(n, " (", Percentage, "%)", sep = "")) %>%
        select(-n) %>%
        rename(n = Count.String)
    }
    if(p.value.flag == "Yes"){
      # Perform Chi-Square Test when sample size for all cells is sufficiently large
      if(!any(clinical.summary.count.table[, -1] == 0) & !(any(clinical.summary.count.table[, -1] < 20) & nrow(clinical.summary.count.table) == 2)){
        
        chi.squared.test = chisq.test(clinical.summary.count.table[, -1])
        
        chi.squared.test.p.value = signif(chi.squared.test$p.value, digits = 3)
        
        clinical.summary$Variable = paste(clinical.summary$Variable, " (p = ", 
                                          chi.squared.test.p.value, ")", sep = "")
        
        # Perform Fisher's Exact Test when sample size is insufficient for chi-squared test
      }else{
        #Quickfix: user Fisher Hybrid if ChiSq not appropriate.
        #if(!any(clinical.summary.count.table[, -1] == 0)) { #} & nrow(clinical.summary.count.table) == 2){
          fisher.test = fisher.test(clinical.summary.count.table[, -1], hybrid=T)
          
          fisher.test.p.value = signif(fisher.test$p.value, digits = 3)
          
          clinical.summary$Variable = paste(clinical.summary$Variable, " (p = ", 
                                            fisher.test.p.value, ")", sep = "")
        #}
      }
    }
    
    names(clinical.summary)[2] = "Value"
    
    clinical.summary = select(clinical.summary, group.variable, Variable, Value, n) 
    
    if(i == 1){
      clinical.summary.table = clinical.summary
    }else{
      clinical.summary.table = rbind.data.frame(clinical.summary.table,
                                                clinical.summary)
    }
  }
  
  #browser()
  ## QuickFix: remove group.variable from id_cols.
  ## QuickFix: fill missing with '0 (0%)'.
  clinical.summary.table = pivot_wider(clinical.summary.table, 
                                       id_cols = c("Variable", "Value"), 
                                       names_from = group.variable, 
                                       values_from = "n", values_fill="0 (0%)")
  
  # Calculate Total Column (Row Sums) and Output Resulting Table to .pdf
  if(percentage.flag == "columns"){
    group.variable.categories.indices = which(names(clinical.summary.table) != "Variable" & 
                                                names(clinical.summary.table) != "Value")
    
    clinical.summary.table.row.totals = rep(NA, nrow(clinical.summary.table))
    
    for(i in 1:nrow(clinical.summary.table)){
      # Calculate Row Total for Each Row in clinical summary table
      row.total = 0
      for(j in group.variable.categories.indices){
        predictor.level.count = unlist(strsplit(as.character(clinical.summary.table[i, j]), "(", fixed = TRUE))
        predictor.level.count = trimws(predictor.level.count[1])
        row.total = row.total + as.integer(predictor.level.count)
      }
      clinical.summary.table.row.totals[i] = row.total
    }
    
    clinical.summary.table.row.totals = data.frame(clinical.summary.table$Variable, clinical.summary.table$Value, 
                                                   clinical.summary.table.row.totals)
    
    names(clinical.summary.table.row.totals) = c("Variable", "Value", "Row.Total")
    
    # Totals by Variable in clinical Summary Table
    clinical.summary.table.variable.totals = clinical.summary.table.row.totals %>%
      group_by(Variable) %>%
      summarise(Variable.Total = sum(Row.Total))
    #browser()
    
    clinical.summary.table.row.totals = clinical.summary.table.row.totals %>%
      inner_join(clinical.summary.table.variable.totals, by = "Variable") %>%
      # Produce Percentages by Predictor Level within Variable
      mutate(Percentage = round(Row.Total / Variable.Total * 100, digits = 0)) %>%
      mutate(Total = paste(Row.Total, " (", Percentage, "%)", sep = "")) %>%
      select(Variable, Value, Total)
   
    clinical.summary.table = clinical.summary.table %>%
      inner_join(clinical.summary.table.row.totals, by = c("Variable", "Value")) %>%
      gt(rowname_col = "Value", groupname_col = "Variable") %>%
      tab_header(title = clinical.table.title)
  }else{
    
    clinical.summary.table = clinical.summary.table %>%
      gt(rowname_col = "Value", groupname_col = "Variable") %>%
      tab_header(title = clinical.table.title)
  }
  
  gtsave(clinical.summary.table, clinical.table.file.name)
  
  return(clinical.summary.table)
  
}
