swap_drugs_1_2 = function(df_) {

  for (dr_var in intersect(c('DrugName', 'Gnumber', 'drug_moa', 'Concentration'), colnames(df_))) {
    df_[[paste0('temp_', dr_var)]] <- df_[[paste0(dr_var, '_2')]]
    df_[[paste0(dr_var, '_2')]] <- df_[[dr_var]] 
    df_[[dr_var]] <- df_[[paste0('temp_', dr_var)]]
    df_[[paste0('temp_', dr_var)]] = NULL
  }
  df_
}