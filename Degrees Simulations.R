
run_degrees_simulations <- function(n_value, number_of_all_ncs_value, n_iterations,
                                         n_permutations, nc_power_value, rejection_rate,
                                         number_of_good_ncs_values, alpha_values,
                                         singal_nc_power, ntree) {
  sg_y <- sg <- 1
  
  # Single lm p-val approach
  start_time <- Sys.time()
  results_inter_bonf <-  foreach (curr_alpha = alpha_values, 
                                  .packages = c("dplyr","foreach"),
                                  .export = c("get_min_p_val", "create_degree_nc", "get_nc_col"),
                                  .combine ="rbind") %do% {
                                    foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                             .packages = c("dplyr","foreach"),
                                             .export = c("get_min_p_val", "create_degree_nc", "get_nc_col"),
                                             .combine ="cbind") %do% {
                                               mean(
                                                 foreach (curr_iter = 1:n_iterations,
                                                          .packages = c("dplyr","foreach"),
                                                          .export = c("get_min_p_val", "create_degree_nc", "get_nc_col"),
                                                          .combine ="cbind") %dopar% 
                                                   {
                                                     get_min_p_val(data= create_degree_nc(n= n_value,  nc_power=singal_nc_power*curr_n_number_of_good_ncs_value, nc_power_split="uniform", 
                                                                                          number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                          number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value,
                                                                                          alpha= curr_alpha),
                                                                   number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                   number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value) < (rejection_rate/number_of_all_ncs_value)
                                                   }
                                               )
                                             }
                                  }
  print(sprintf("Bonf: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) singal_nc_power=%s curr_alpha=%s at %s",
                n_iterations, toString(number_of_good_ncs_values), number_of_all_ncs_value, singal_nc_power,  toString(alpha_values), Sys.time()))   
  print( Sys.time() - start_time)
  print(results_inter_bonf)
  
  # F-test 
  start_time <- Sys.time()
  results_f_test <- foreach (curr_alpha = alpha_values, 
                             .packages = c("dplyr","foreach"),
                             .export = c("get_f_p_val", "create_degree_nc", "get_nc_col"),
                             .combine ="rbind") %do% {
                               foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                        .packages = c("dplyr","foreach"),
                                        .export = c("get_f_p_val", "create_degree_nc", "get_nc_col"),
                                        .combine ="cbind") %do% {
                                          good_nc_number_iteration_results <-  foreach (curr_iter = 1:n_iterations,
                                                                                        .packages = c("dplyr","foreach"),
                                                                                        .export = c("get_f_p_val", "create_degree_nc", "get_nc_col"),
                                                                                        .combine ="cbind") %dopar% 
                                            {
                                              get_f_p_val(data= create_degree_nc(n= n_value,  nc_power=singal_nc_power*curr_n_number_of_good_ncs_value, nc_power_split="uniform", 
                                                                                 number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                 number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value,
                                                                                 alpha= curr_alpha),
                                                          number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                          number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value) < (rejection_rate)
                                            }
                                          mean(good_nc_number_iteration_results)
                                        }
                             }
  print(sprintf("F-test: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) singal_nc_power=%s curr_alpha=%s at %s",
                n_iterations, toString(number_of_good_ncs_values), number_of_all_ncs_value, singal_nc_power,  toString(alpha_values), Sys.time()))    
  print( Sys.time() - start_time)
  print(results_f_test)
  
  # RF test
  start_time <- Sys.time()
  results_rf_interactions <-   foreach (curr_alpha = alpha_values,
                                        .packages = c("dplyr","foreach"),
                                        .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                        .combine ="rbind") %do% {
                                          foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                                   .packages = c("dplyr","foreach"),
                                                   .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                                   .combine ="cbind") %do%
                                            {
                                              good_nc_number_iteration_results <- 
                                                foreach (curr_iter = 1:n_iterations,
                                                         .packages = c("dplyr","foreach"),
                                                         .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                                         .combine ="cbind") %dopar% 
                                                {
                                                  (permutations.test.for.lm(data= create_degree_nc(n= n_value,  nc_power=singal_nc_power*curr_n_number_of_good_ncs_value, nc_power_split="uniform", 
                                                                                                   number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                                   number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value,
                                                                                                   alpha= curr_alpha),
                                                                            instrument_form= "iv",
                                                                            instrument= "iv", 
                                                                            controls= c("t"),
                                                                            weights= NULL, variables_to_remove= c("h", "y", "iv"),
                                                                            title = "Simulation bad NC", n_permutations= n_permutations, 
                                                                            conditioned = T, OOB= T, saveplot=F,
                                                                            mtry_ratio= 1/3, ntree=ntree) <= rejection_rate)
                                                }
                                              print(sprintf("RF: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) singal_nc_power=%s curr_alpha=%s at %s",
                                                            n_iterations, curr_n_number_of_good_ncs_value, number_of_all_ncs_value, singal_nc_power, curr_alpha, Sys.time()))
                                              print(c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results)))
                                              mean(good_nc_number_iteration_results)
                                            }
                                        }
  print( Sys.time() - start_time)
  print(results_rf_interactions)
  
  # Bagging
  start_time <- Sys.time()
  results_bag_interactions <-   foreach (curr_alpha = alpha_values,
                                         .packages = c("dplyr","foreach"),
                                         .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                         .combine ="rbind") %do% {
                                           foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
                                                    .packages = c("dplyr","foreach"),
                                                    .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                                    .combine ="cbind") %do%
                                             {
                                               good_nc_number_iteration_results <- 
                                                 foreach (curr_iter = 1:n_iterations,
                                                          .packages = c("dplyr","foreach"),
                                                          .export = c("permutations.test.for.lm", "run.rf.multiple.negative.controls","get_NC_matrix", "calculate.visualize.p.values","create_degree_nc", "get_nc_col"),
                                                          .combine ="cbind") %dopar% 
                                                 {
                                                   (permutations.test.for.lm(data= create_degree_nc(n= n_value,  nc_power=singal_nc_power*curr_n_number_of_good_ncs_value, nc_power_split="uniform", 
                                                                                                    number_of_good_ncs = curr_n_number_of_good_ncs_value,
                                                                                                    number_of_bad_ncs = number_of_all_ncs_value- curr_n_number_of_good_ncs_value,
                                                                                                    alpha= curr_alpha),
                                                                             instrument_form= "iv",
                                                                             instrument= "iv", 
                                                                             controls= c("t"),
                                                                             weights= NULL, variables_to_remove= c("h", "y", "iv"),
                                                                             title = "Simulation bad NC", n_permutations= n_permutations, 
                                                                             conditioned = T, OOB= T, saveplot=F,
                                                                             mtry_ratio= 1, ntree=ntree) <= rejection_rate)
                                                 }
                                               print(sprintf("BAG: Done %s iterations of curr_number_of_good_nc_value=%s from(%s) singal_nc_power=%s curr_alpha=%s at %s",
                                                             n_iterations, curr_n_number_of_good_ncs_value, number_of_all_ncs_value, singal_nc_power, curr_alpha, Sys.time()))
                                               print(c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results)))
                                               mean(good_nc_number_iteration_results)#c(mean(good_nc_number_iteration_results), sd(good_nc_number_iteration_results))
                                               
                                             }
                                         }
  print( Sys.time() - start_time)
  print(results_bag_interactions)
  
}


