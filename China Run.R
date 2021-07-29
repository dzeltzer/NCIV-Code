
# rerunning table 3 w. the appropriate controls
# Based on pg. 2136 and czone_analysis_ipw_final.do
replicate_table_3 <- function(workfile_china, # dataset used to replicate table 3
                              col_2_controls, # groups of controls corresponding to
                              col_3_controls, # the columns in table 3
                              col_4_controls,
                              col_5_controls,
                              col_6_controls,
                              G, N #used for heteroscedasticity stds
                              ) {
  # Partial replication of table 3
  # eststo: ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw=d_tradeotch_pw_lag) t2 [aw=timepwt48], cluster(statefip) first
  tbl_3_col_2 <- ivreg(formula= get_formula_for_table_3(col_2_controls)
                       , weights= timepwt48
                       , data= workfile_china)
  print(get_stata_coef(G, N, tbl_3_col_2))
  
  tbl_3_col_3 <- ivreg(formula= get_formula_for_table_3(col_3_controls)
                       , weights= timepwt48
                       , data= workfile_china)
  print(get_stata_coef(G, N, tbl_3_col_3))
  
  tbl_3_col_4 <- ivreg(formula= get_formula_for_table_3(col_4_controls)
                       , weights= timepwt48
                       , data= workfile_china)
  print(get_stata_coef(G, N, tbl_3_col_4))
  
  tbl_3_col_5 <- ivreg(formula= get_formula_for_table_3(col_5_controls)
                       , weights= timepwt48
                       , data= workfile_china)
  print(get_stata_coef(G, N, tbl_3_col_5))
  
  tbl_3_col_6 <- ivreg(formula= get_formula_for_table_3(col_6_controls)
                       , weights= timepwt48
                       , data= workfile_china)
  print(get_stata_coef(G, N, tbl_3_col_6))
}

run_china_nciv <- function(data, # China data
                           instrument_form, # the form of the instrument, 
                           # for example: logarithmic: log(instrument2000), 
                           # linear: instrument2000
                           instrument, # the instrument for the NCIV test
                           col_2_controls, # groups of controls corresponding to
                           col_3_controls, # the columns in table 3
                           col_4_controls, 
                           col_5_controls,
                           col_6_controls,
                           weights, #the weights column name used for residualization 
                           variables_to_remove,
                           permutations,# the number of permutations in each 
                           # NCIV test
                           OOB, # use Out Of Bag RMSE instead of Cross Validation
                           # error
                           mtry_ratio, #the ratio of variables to check in 
                           # each iteration of RF algorithm
                           ntree, # the number of trees in the RF prediction
                           # algorithm used for NCIV test
                           title # the title for the output file and graph
                           )
  {
  rmse_rf_col_2 <- permutations.test.for.lm(data= data, instrument_form= instrument_form,
                                            instrument= instrument,
                                            controls= col_2_controls[2:length(col_2_controls)],
                                            weights=weights, variables_to_remove= variables_to_remove,
                                            title = paste(title,"col 2"), n_permutations= permutations,
                                            conditioned = T, OOB= OOB, saveplot=T,
                                            mtry_ratio=mtry_ratio, ntree=ntree)

  rmse_rf_col_3 <- permutations.test.for.lm(data= data, instrument_form= instrument_form,
                                            instrument= instrument,
                                            controls= col_3_controls[2:length(col_3_controls)],
                                            weights=weights, variables_to_remove= variables_to_remove,
                                            title = paste(title,"col 3"), n_permutations= permutations,
                                            conditioned = T, OOB= OOB, saveplot=T,
                                            mtry_ratio=mtry_ratio, ntree=ntree)

  rmse_rf_col_4 <- permutations.test.for.lm(data= data, instrument_form= instrument_form,
                                            instrument= instrument,
                                            controls= col_4_controls[2:length(col_4_controls)],
                                            weights=weights, variables_to_remove= variables_to_remove,
                                            title = paste(title,"col 4"), n_permutations= permutations,
                                            conditioned = T, OOB= OOB, saveplot=T,
                                            mtry_ratio=mtry_ratio, ntree=ntree)

  rmse_rf_col_5 <- permutations.test.for.lm(data= data, instrument_form= instrument_form,
                                            instrument= instrument,
                                            controls= col_5_controls[2:length(col_5_controls)],
                                            weights=weights, variables_to_remove= variables_to_remove,
                                            title = paste(title,"col 5"), n_permutations= permutations,
                                            conditioned = T, OOB= OOB, saveplot=T,
                                            mtry_ratio=mtry_ratio, ntree=ntree)
  
  rmse_rf_col_6 <- permutations.test.for.lm(data= data, instrument_form= instrument_form,
                                            instrument= instrument,
                                            controls= col_6_controls[2:length(col_6_controls)],
                                            weights=weights, variables_to_remove= variables_to_remove,
                                            title = paste(title,"col 6"), n_permutations= permutations,
                                            conditioned = T, OOB= OOB, saveplot=T,
                                            mtry_ratio=mtry_ratio, ntree=ntree)
  
  results <- rbind(c(rmse_rf_col_2, rmse_rf_col_3, rmse_rf_col_4,
                     rmse_rf_col_5, rmse_rf_col_6))
  write_csv(as.data.frame(results), sprintf("out//china_results_ntree%s.csv", ntree))
  }
  