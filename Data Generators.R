
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
                           specified_nco = NULL,
                           specified_for_sur_only = F
                           )
  {
  lin_p_vals_col_2 <- bonf.f.test.for.lm(data= data,
                                         instrument_form= "instrument2000",
                                         instrument= "instrument2000",
                                         controls= col_2_controls[2:length(col_2_controls)],
                                         weights="timepwt48",
                                         variables_to_remove= variables_to_remove,
                                         title= "Col 2",
                                         conditioned = T,
                                         specified_nco = specified_nco,
                                         specified_for_sur_only
                                         )
  
  lin_p_vals_col_3 <- bonf.f.test.for.lm(data= data, instrument_form= "instrument2000",
                                         instrument= "instrument2000",
                                         controls= col_3_controls[2:length(col_3_controls)],
                                         weights="timepwt48",
                                         variables_to_remove= variables_to_remove,
                                         title= "Col 3",
                                         conditioned = T,
                                         specified_nco = specified_nco,
                                         specified_for_sur_only)
  
  lin_p_vals_col_4 <- bonf.f.test.for.lm(data= data, instrument_form= "instrument2000",
                                         instrument= "instrument2000",
                                         controls= col_4_controls[2:length(col_4_controls)],
                                         weights="timepwt48",
                                         variables_to_remove= variables_to_remove,
                                         title= "Col 4",
                                         conditioned = T,
                                         specified_nco = specified_nco,
                                         specified_for_sur_only)
  
  lin_p_vals_col_5 <- bonf.f.test.for.lm(data= data, instrument_form= "instrument2000",
                                         instrument= "instrument2000",
                                         controls= col_5_controls[2:length(col_5_controls)],
                                         weights="timepwt48",
                                         variables_to_remove= variables_to_remove,
                                         title= "Col 5",
                                         conditioned = T,
                                         specified_nco = specified_nco,
                                         specified_for_sur_only)
  
  lin_p_vals_col_6 <- bonf.f.test.for.lm(data= data, instrument_form= "instrument2000",
                                         instrument= "instrument2000",
                                         controls= col_6_controls[2:length(col_6_controls)],
                                         weights="timepwt48",
                                         variables_to_remove= variables_to_remove,
                                         title= "Col 6",
                                         conditioned = T,
                                         specified_nco = specified_nco,
                                         specified_for_sur_only)
  
  results <- rbind(c(lin_p_vals_col_2, lin_p_vals_col_3, lin_p_vals_col_4,
                     lin_p_vals_col_5, lin_p_vals_col_6))
  return(results)
  }
  
