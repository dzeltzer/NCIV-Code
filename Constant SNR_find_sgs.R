library(utils)
library(doParallel)

get_sgs_for_constant_SNR <- function(get_r2_func, hyper_parameter_vec, number_of_good_ncs_values,
                                     target_mse, n_observations_for_dis, number_of_all_ncs_value, single_nc_power, ntree) {
  start_time <- Sys.time()
  if (max(number_of_good_ncs_values) > number_of_all_ncs_value){
    abort("Cannot have value from number_of_good_ncs_values greater than number_of_all_ncs_value")
  }
  
  res_m40_3600 <-  foreach (hyper_parameter = hyper_parameter_vec,
                   .packages = c("dplyr","foreach"),
                   .export = c("get_r2_from_interactions_model", "get_r2_from_CES_model",
                               "create_interactions_multi_nc", "create_ces_nc", "get_nc_col",
                               "ces_q",
                               "get_model_mse", "prepare_z_NC", "get_NC_matrix",
                               "run.rf.multiple.negative.controls"),
                   .combine ="cbind") %dopar%
    {
      print(hyper_parameter)
      foreach (curr_n_number_of_good_ncs_value = number_of_good_ncs_values,
               .packages = c("dplyr","foreach"),
               .export = c("create_interactions_multi_nc"),
               .combine ="cbind") %dopar%
        {
          print(curr_n_number_of_good_ncs_value)
  
          print(sprintf("target mse: %f",target_mse))
  
          root_res <- uniroot(f= function(x) {
            get_r2_func(n_observations= n_observations_for_dis, 
                        number_of_all_ncs_value = number_of_all_ncs_value, # 
                        sg =  x, # sg =: the standard variation for the unmeasured confounder
                        curr_n_number_of_good_ncs_value = curr_n_number_of_good_ncs_value, # 
                        hyper_parameter= hyper_parameter, # curr_alpha /curr_rho
                        single_nc_power = single_nc_power, # single_nc_power = 
                        ntree =  ntree # 
                        ) - target_mse
          }, 
          interval = c(0, 60),
          tol= 0.1,
          extendInt="upX",
          trace=3
          )
          root_res[["hyper_parameter"]] <- hyper_parameter
          root_res[["number_of_good_ncs_values"]] <- curr_n_number_of_good_ncs_value
          print(root_res)
        }
    }
  
  constant_snr_coef_mat <- data.frame(matrix(unlist(as.data.frame(res_m40_3600)[1,]),
         nrow = length(hyper_parameter_vec),
         ncol = length(number_of_good_ncs_values)
         ),
         row.names = as.character(hyper_parameter_vec)
         ) %>% 
    `colnames<-`(number_of_good_ncs_values)
  
  print(Sys.time() - start_time)
  write.csv(constant_snr_coef_mat, file.path("out", "sgs_for_constant_mse", 
    paste0("inter_constant_snr_sgs_mse_m", number_of_all_ncs_value, "_target_mse",
           target_mse, Sys.Date(),
           ".csv")))
    return(constant_snr_coef_mat)
  
}



