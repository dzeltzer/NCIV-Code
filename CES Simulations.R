# A simplifed version of CES function as describrd in 
# https://en.wikipedia.org/wiki/Constant_elasticity_of_substitution
# Q = (K^rho + L^rho)^(1/rho)
# where rho in [-inf, 1]
# when rho == inf -> CES(K, L) = MIN(K,L)
# when rho == 0 -> CES(K,L) = Const
# when rho == 1 -> CES(K, L) = K + L
run_ces_simulations <- function(n_value,# the number of observations in
                                # each iteration
                                number_of_all_ncs_value, # the number of 
                                # all variables in the output data frame
                                n_iterations, # the number iterations done
                                # in the simulations for each hyper parameters
                                # specification (number_of_good_ncs_values and 
                                # alpha_values). In each iteration we create a 
                                # simulated data frame, try different algorithms
                                # for NCIV test.
                                n_permutations, # the number of permutations
                                # in the the permutations tests (RF, Bagging)
                                # for each iteration
                                rejection_rate, # the rejection threshold rate that we 
                                # reject the null
                                number_of_good_ncs_values, # hyper parameter: vector of
                                # number of negative control variables (out 
                                # of number_of_all_ncs_value)
                                rho_values, # hyper parameters: a vector of rho - the value
                                # that determine the shape of CES function (desrtibed above)
                                single_nc_power, # the coefficient of the effect of
                                # the unmeasured confounder on the negative controls
                                ntree # the number of trees in the RF prediction
                                # algorithm used for NCIV test
                                ) {
  sg_y <- sg <- 1
  
  # Single lm p-val approach
  start_time <- Sys.time()
  results_bonf_raw <-  foreach (curr_rho = rho_values, 
                                .packages = c("dplyr","foreach"),
                                .export = c("get_min_p_val", "create_ces_nc", "ces_q","get_nc_col"),
                                .combine ="rbind") %do% {
                                  foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                           .packages = c("dplyr","foreach"),
                                           .export = c("get_min_p_val", "create_ces_nc", "ces_q","get_nc_col"),
                                           .combine ="cbind") %do% {
                                             mean(
                                               foreach (curr_iter = 1:n_iterations,
                                                        .packages = c("dplyr","foreach"),
                                                        .export = c("get_min_p_val", "create_ces_nc", "ces_q","get_nc_col"),
                                                        .combine ="cbind") %dopar% 
                                                 {
                                                   get_min_p_val(data= create_ces_nc(n= n_value, 
                                                                                     nc_power=single_nc_power*curr_n_number_of_good_ncs_value,
                                                                                     nc_power_split="uniform", 
                                                                                     number_of_first_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                     number_of_second_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                     number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value,
                                                                                     rho= curr_rho),
                                                                 number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                 number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value) < (rejection_rate/number_of_all_ncs_value)
                                                 }
                                             )
                                           }
                                }
  print(sprintf("Bonf: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) single_nc_power=%s rho_values=%s at %s",
                n_iterations, toString(number_of_good_ncs_values), number_of_all_ncs_value, single_nc_power,  toString(rho_values), Sys.time()))   
  print( Sys.time() - start_time)
  
  columns <- c("Scenario", "Algorithm", "n_value", "number_of_all_ncs_value", "n_iterations",
               "n_permutations", "number_of_good_ncs_value", "rho_value",
               "single_nc_power", "ntree", "rejection rate")
  
  results_bonf <- prepare_simulations_results(results_bonf_raw, columns,
                                              "Degrees", "Bonf", 
                                              single_nc_power, rho_values,
                                              number_of_good_ncs_values, n_value,
                                              number_of_all_ncs_value, n_iterations,
                                              n_permutations, ntree)
  
  print(results_bonf)  
  
  # F-test 
  start_time <- Sys.time()
  results_f_test_raw <- foreach (curr_rho = rho_values, 
                                 .packages = c("dplyr","foreach"),
                                 .export = c("get_f_p_val", "create_ces_nc","ces_q","get_nc_col"),
                                 .combine ="rbind") %do% {
                                   foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                            .packages = c("dplyr","foreach"),
                                            .export = c("get_f_p_val", "create_ces_nc","ces_q","get_nc_col"),
                                            .combine ="cbind") %do% {
                                              good_nc_number_iteration_results <-  foreach (curr_iter = 1:n_iterations,
                                                                                            .packages = c("dplyr","foreach"),
                                                                                            .export = c("get_f_p_val", "create_ces_nc","ces_q","get_nc_col"),
                                                                                            .combine ="cbind") %dopar% 
                                                {
                                                  get_f_p_val(data= create_ces_nc(n= n_value, 
                                                                                  nc_power=single_nc_power*curr_n_number_of_good_ncs_value,
                                                                                  nc_power_split="uniform", 
                                                                                  number_of_first_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                  number_of_second_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                  number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value,
                                                                                  rho= curr_rho),
                                                              number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                              number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value) < (rejection_rate)
                                                }
                                              mean(good_nc_number_iteration_results)
                                            }
                                 }
  print(sprintf("F-test: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) single_nc_power=%s rho_values=%s at %s",
                n_iterations, toString(number_of_good_ncs_values), number_of_all_ncs_value, single_nc_power,  toString(rho_values), Sys.time()))    
  print( Sys.time() - start_time)
  
  results_f_test <- prepare_simulations_results(results_f_test_raw, columns,
                                                "Degrees", "F test", 
                                                single_nc_power, rho_values,
                                                number_of_good_ncs_values, n_value,
                                                number_of_all_ncs_value, n_iterations,
                                                n_permutations, ntree)
  
  print(results_f_test)
  
  # SUR test
  start_time <- Sys.time()
  results_sur_raw <-   foreach (curr_rho = rho_values,
                                .packages = c("dplyr","foreach", "systemfit"),
                                .export = c("get_sur_p_val","create_ces_nc","ces_q","get_nc_col"),
                                .combine ="rbind") %do% {
                                  foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                           .packages = c("dplyr","foreach", "systemfit"),
                                           .export = c("get_sur_p_val","create_ces_nc","ces_q","get_nc_col"),
                                           .combine ="cbind") %do%
                                    {
                                      good_nc_number_iteration_results <- 
                                        foreach (curr_iter = 1:n_iterations,
                                                 .packages = c("dplyr","foreach", "systemfit"),
                                                 .export = c("get_sur_p_val","create_ces_nc","ces_q","get_nc_col"),
                                                 .combine ="cbind") %dopar% 
                                        {
                                          (get_sur_p_val(data= create_ces_nc(n= n_value, 
                                                                             nc_power=single_nc_power*curr_n_number_of_good_ncs_value,
                                                                             nc_power_split="uniform", 
                                                                             number_of_first_good_ncs = curr_n_number_of_good_ncs_value,
                                                                             number_of_second_good_ncs = curr_n_number_of_good_ncs_value,
                                                                             number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value,
                                                                             rho= curr_rho),
                                                         number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                         number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value) <= rejection_rate)
                                        }
                                      print(sprintf("SUR: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) single_nc_power=%s curr_rho=%s at %s",
                                                    n_iterations, curr_n_number_of_good_ncs_value, number_of_all_ncs_value, single_nc_power, curr_rho, Sys.time()))
                                      print(c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results)))
                                      mean(good_nc_number_iteration_results)
                                    }
                                }
  print( Sys.time() - start_time)
  
  results_sur <- prepare_simulations_results(results_sur_raw, columns,
                                             "Degrees", "SUR", 
                                             single_nc_power, rho_values,
                                             number_of_good_ncs_values, n_value,
                                             number_of_all_ncs_value, n_iterations,
                                             n_permutations, ntree)
  
  print(results_sur)  
  
  # RF test
  start_time <- Sys.time()
  results_rf_raw <-   foreach (curr_rho = rho_values,
                               .packages = c("dplyr","foreach"),
                               .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                               .combine ="rbind") %do% {
                                 foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                          .packages = c("dplyr","foreach"),
                                          .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                                          .combine ="cbind") %do%
                                   {
                                     good_nc_number_iteration_results <- 
                                       foreach (curr_iter = 1:n_iterations,
                                                .packages = c("dplyr","foreach"),
                                                .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                                                .combine ="cbind") %dopar% 
                                       {
                                         (permutations.test.for.lm(data= create_ces_nc(n= n_value, 
                                                                                       nc_power=single_nc_power*curr_n_number_of_good_ncs_value,
                                                                                       nc_power_split="uniform", 
                                                                                       number_of_first_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                       number_of_second_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                       number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value,
                                                                                       rho= curr_rho),
                                                                   instrument_form= "iv",
                                                                   instrument= "iv", 
                                                                   controls= c("t"),
                                                                   weights= NULL, variables_to_remove= c("h", "y", "iv"),
                                                                   title = "Simulation bad NC", n_permutations= n_permutations, 
                                                                   conditioned = T, OOB= T, saveplot=F,
                                                                   mtry_ratio= 1/3, ntree=ntree) <= rejection_rate)
                                       }
                                     print(sprintf("RF: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) single_nc_power=%s curr_rho=%s at %s",
                                                   n_iterations, curr_n_number_of_good_ncs_value, number_of_all_ncs_value, single_nc_power, curr_rho, Sys.time()))
                                     print(c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results)))
                                     mean(good_nc_number_iteration_results)
                                   }
                               }
  print( Sys.time() - start_time)
  
  results_rf <- prepare_simulations_results(results_rf_raw, columns,
                                            "Degree", "RF", 
                                            single_nc_power, rho_values,
                                            number_of_good_ncs_values, n_value,
                                            number_of_all_ncs_value, n_iterations,
                                            n_permutations, ntree)
  print(results_rf)  
  
  # Bagging
  start_time <- Sys.time()
  results_bag_raw <-   foreach (curr_rho = rho_values,
                                .packages = c("dplyr","foreach"),
                                .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                                .combine ="rbind") %do% {
                                  foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                           .packages = c("dplyr","foreach"),
                                           .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                                           .combine ="cbind") %do%
                                    {
                                      good_nc_number_iteration_results <- 
                                        foreach (curr_iter = 1:n_iterations,
                                                 .packages = c("dplyr","foreach"),
                                                 .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_ces_nc","ces_q","get_nc_col"),
                                                 .combine ="cbind") %dopar% 
                                        {
                                          (permutations.test.for.lm(data= create_ces_nc(n= n_value, 
                                                                                        nc_power=single_nc_power*curr_n_number_of_good_ncs_value,
                                                                                        nc_power_split="uniform", 
                                                                                        number_of_first_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                        number_of_second_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                        number_of_bad_ncs = number_of_all_ncs_value- 2*curr_n_number_of_good_ncs_value,
                                                                                        rho= curr_rho),
                                                                    instrument_form= "iv",
                                                                    instrument= "iv", 
                                                                    controls= c("t"),
                                                                    weights= NULL, variables_to_remove= c("h", "y", "iv"),
                                                                    title = "Simulation bad NC", n_permutations= n_permutations, 
                                                                    conditioned = T, OOB= T, saveplot=F,
                                                                    mtry_ratio= 1, ntree=ntree) <= rejection_rate)
                                        }
                                      print(sprintf("BAG: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) single_nc_power=%s curr_rho=%s at %s",
                                                    n_iterations, curr_n_number_of_good_ncs_value, number_of_all_ncs_value, single_nc_power, curr_rho, Sys.time()))
                                      print(c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results)))
                                      mean(good_nc_number_iteration_results)#c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results))
                                      
                                    }
                                }
  print( Sys.time() - start_time)
  results_bag <- prepare_simulations_results(results_bag_raw, columns,
                                             "Interactions", "Bagging", 
                                             single_nc_power, rho_values,
                                             number_of_good_ncs_values, n_value,
                                             number_of_all_ncs_value, n_iterations,
                                             n_permutations, ntree)  
  
  print(results_bag)
  
  results <- rbind(results_bonf,
                   results_f_test, 
                   results_sur,
                   results_rf,
                   results_bag)
  colnames(results) <- columns
  
  file_full_path <- file.path("out", paste("Simulations CES", gsub("[:]", " ",Sys.time())))
  write.csv(results, sprintf("%s.csv", file_full_path))
  
  return(results)  
}


