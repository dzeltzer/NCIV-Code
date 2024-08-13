set.seed(12345)

# Loads the necessary libraries for the whole project. You can open the file 
# and see whether you need to install some libraries
source("Load Libraries.R")

# This file contains all the necessary functions for the NCIV approach
# sourcing the file will load the main functions:
# (1) permutations.test.for.felm
# (2) permutations.test.for.lm
# which perform a NCIV permutations test
# and the corresponding inner functions:
# run.rf.multiple.negative.controls (inner functions of )
# get_NC_matrix_felm
# get_NC_matrix
# calculate.visualize.p.values
# as well as:
# (3) get_stata_coef 
# which get heteroscedasticity consistent estimates for linear model
source("NCIV Main.R")

# Simulations
# this file loads all the necessary functions for the simulations: Data generations.
# sourcing the file will load the main functions:
# create_degree_nc()- creates a data frame scenario where the relation between 
# the unmeasured confounder and the IV is addition of linear and squared 
# create_interactions_multi_nc() - creates a data frame scenario where the 
# there are two unmeasured confounders that there interaction is correlated 
# with the IV
# create_ces_nc() - creates a data frame scenario where the relation between 
# the unmeasured confounders and the IV is via CES functions
# prepare_simulations_results() - convert data from foreach output to readable 
# simulations results 
source("Data Generators.R")

# this file loads benchmark tests algorithms
# get_oracle_p_val() - the theoretical best algorithm  
# get_min_p_val() - using Bonferroni correction 
# get_sur_p_val() - Using SUR (Seemingly Unrelated Regression)
# get_f_p_val() - Using F-test
source("Benchmark Algorithms.R")

#stopCluster(cl)
cl <- makeCluster(detectCores()*(3/4))
registerDoParallel(cl)

# Section 3 - Method
# Simulation 1 - advantage in multi-variable setting
# Y-axis: power (true negatives) 
# X-avis: number of vars (For now: with increasing SNR , in the 
# future with fixed signal to noise ratio) 
# t-test for the (apriori) strongest candidate
# t-test with Bonferroni correction for multiple hypothesis testing. 
# F-test
# SUR (explain why it’s not working)
# RF
# Simulation 2 - CES functions
# X-axis: elasticity of substitution
# Rest - the same


n_observations <- 100 # should be 1000 for level playing field, 100 to run fast
number_of_all_ncs_value <- 40 # 10 or 40
number_of_good_ncs_values <- c(2, 3, 5)# c(2, 5, 10, 20) #  c(2, 3, 5) number_of_good_ncs_values = c(2, 5, 10, 20) - used when number_of_all_ncs_value = 40
n_iterations <- 10
n_permutations <- 100
rejection_rate <- 0.05
ntree <- 100

# this file loads run_degrees_simulations() which runs simulations scenario
# for degrees
source("Degrees Simulations.R")



alpha_values = c(0, 0.25, 0.5, 0.75, 1)
single_nc_power <- 0.2 # 0.25 for M=10, 0.125 for M=40

# deg_not_null <- run_degrees_simulations(n_value= 50,
#                                         number_of_all_ncs_value = 10,
#                                         n_iterations = 2, 
#                                         n_permutations = 100,
#                                         rejection_rate = 0.05,
#                                         number_of_good_ncs_values = c(2, 3, 5),
#                                         alpha_values = c(0, 0.25, 0.5, 0.75, 1),
#                                         single_nc_power = 0.125,
#                                         ntree = 100,
#                                         is_null_scenario = F)
# 
# deg_null <- run_degrees_simulations(n_value= 1000,
#                                         number_of_all_ncs_value = 40,
#                                         n_iterations = 250, 
#                                         n_permutations = 100,
#                                         rejection_rate = 0.05,
#                                         number_of_good_ncs_values = c(2, 5, 10, 20),
#                                         alpha_values = c(0, 0.25, 0.5, 0.75, 1),
#                                         single_nc_power = 0.125,
#                                         ntree = 100,
#                                         is_null_scenario = T)



# this file loads run_interactions_simulations() which runs simulations scenario
# for interactions
source("Interactions Simulations.R")

n_observations <- 1000 # should be 1000 for level playing field, 100 to run fast
number_of_all_ncs_value <- 10 # 10 or 40
number_of_good_ncs_values <- c(2, 3, 5)# c(2, 5, 10, 20) #  c(2, 3, 5) number_of_good_ncs_values = c(2, 5, 10, 20) - used when number_of_all_ncs_value = 40

n_iterations <- 250
n_permutations <- 100
rejection_rate <- 0.05

alpha_values = c(0, 0.25, 0.5, 0.75, 1)
single_nc_power <- 0.5


# Increasing SNR

for (is_null_scenario in c(F,T)){
  
  inter_true_1 <- run_interactions_simulations(n_observations= n_observations, 
                                               number_of_all_ncs_value = number_of_all_ncs_value,
                                               n_iterations = n_iterations,
                                               n_permutations = n_permutations,
                                               rejection_rate = rejection_rate,
                                               number_of_good_ncs_values = number_of_good_ncs_values,
                                               alpha_values = alpha_values,
                                               single_nc_power = single_nc_power,
                                               ntree = ntree,
                                               is_null_scenario = is_null_scenario)
  
}


# this file loads run_ces_simulations() which runs simulations scenario
# with CES functions

source("CES Simulations.R")

n_observations <- 1000 # should be 1000 for level playing field, 100 to run fast
n_iterations <- 250

number_of_all_ncs_value <- 40 # 10 or 40
number_of_good_ncs_values <-  c(2, 5, 10, 20) #  c(2, 3, 5) number_of_good_ncs_values = c(2, 5, 10, 20) - used when number_of_all_ncs_value = 40
ntree <- 100

rho_values = c(-10, -1, 0, 0.1,0.5, 1) #c(-10,0.5, 1)#
single_nc_power <- 0.125
is_null_scenario <- T

for (is_null_scenario in c(T)){
  
ces_not_null_m40 <- run_ces_simulations(n_value= n_observations,
                                    number_of_all_ncs_value = number_of_all_ncs_value,
                                    n_iterations = n_iterations,
                                    n_permutations = n_permutations,
                                    rejection_rate = rejection_rate,
                                    number_of_good_ncs_values = number_of_good_ncs_values,
                                    rho_values = rho_values, # -10, -1,  0.1,0.5,
                                    single_nc_power = single_nc_power,
                                    sg = 1,
                                    constant_snr_sg_mat = NULL,
                                    ntree = ntree,
                                    is_null_scenario = T)

}


# Constant SNR

#source("Constant SNR Funcs.R")

# inter_target_mse <- 4177.824
# n_observations_for_dis <- 20000 # n_observations *10

# inter_constant_snr_sg_mat <- get_sgs_for_constant_SNR(get_mse_from_interactions_model, alpha_values, 
#                                                       number_of_good_ncs_values, inter_target_mse,
#                                                       n_observations_for_dis, number_of_all_ncs_value, 
#                                                       single_nc_power, ntree)

for (is_null_scenario in c(F)){
  
  is_null_scenario <- T
  alpha_values <- c(0, 0.25, 0.5, 0.75)

  inter_constant_snr_sg_mat <- read.csv(file.path("data", "sgs_mat", "SGS_Table_inter_SingleNC_1_RR_05.csv")) %>% 
    select(-1)
  
  n_iterations <- 250
  number_of_good_ncs_values <- c(2,4,7,10)
  number_of_all_ncs_value <- 20
  
  inter_constant_true <- run_interactions_simulations(n_observations= n_observations, 
                                                      number_of_all_ncs_value = number_of_all_ncs_value,
                                                      n_iterations = n_iterations,
                                                      n_permutations = n_permutations,
                                                      rejection_rate = rejection_rate,
                                                      number_of_good_ncs_values = number_of_good_ncs_values,
                                                      alpha_values = alpha_values,#alpha_values,#
                                                      single_nc_power = single_nc_power,
                                                      sg= NA,
                                                      constant_snr_sg_mat = inter_constant_snr_sg_mat,
                                                      ntree = ntree,
                                                      is_null_scenario = is_null_scenario,
                                                      saveplot = F)
}


# ces

number_of_all_ncs_value <- 40 # 10 or 40
number_of_good_ncs_values <-  c(2, 5, 10, 20) # c(2, 3, 5)#

target_mse_ces <- 2^2

# number_of_all_ncs_value <- 10 # 10 or 40
# number_of_good_ncs_values <- c(2, 3, 5)# c(2, 5, 10, 20) 
# 
# ces_constant_snr_sg_mat_m10 <- get_sgs_for_constant_SNR(get_r2_from_CES_model, rho_values, 
#                                                     number_of_good_ncs_values, target_mse_ces, 5*n_observations,
#                                                     number_of_all_ncs_value, single_nc_power, ntree)

ces_constant_snr_sg_mat_m40 <- get_sgs_for_constant_SNR(get_r2_from_CES_model, rho_values, 
                                                        number_of_good_ncs_values, target_mse, 5*n_observations,
                                                        number_of_all_ncs_value, single_nc_power, ntree)



ces_not_nullm_40 <- run_ces_simulations(n_value= n_observations,
                                        number_of_all_ncs_value = number_of_all_ncs_value,
                                        n_iterations = n_iterations,
                                        n_permutations = 100,
                                        rejection_rate = rejection_rate,
                                        number_of_good_ncs_values = number_of_good_ncs_values,
                                        rho_values = rho_values, # -10, -1,  0.1,0.5, 
                                        single_nc_power = single_nc_power,
                                        sg = 1,
                                        constant_snr_sg_mat = ces_constant_snr_sg_mat_m40,
                                        ntree = ntree,
                                        is_null_scenario = F)

# Section 4 - Examples

# The China Syndrome: Local Labor Market Effects of Import Competition in the
# United States (Autor, Dorn, and Hanson)

# This file loads all the necessary data functions, and variables groups for China article
# sourcing this will result in loading the main functions:
# (1) get_formula_for_table_3() 
# which creates the formula for replicating table 3 in the article
# as well as loading:
# (1) workfile_china
# basic data set used for replicating table 3 in the article
# (2) china_1990
# data set with only 1990 data with all available negative controls
# (3) china_1990_only_org_nc
# data set with only 1990 data with only the negative control used in the
# article - the outcome in 2000
# (4) col_2_controls - col_4_controls, location_controls, man_controls ...
# different groups of controls
source("China Init.R")

# This file loads all the necessary functions that runs replication of table 3,
# and NCIV test for China article
# sourcing this will result in loading the main functions:
# (1) replicate_table_3()
# (2) run_china_nciv() - runs 5 different controls specifications of NCIV
#     test on the IV (delta trade to other states)
source("China Run.R")


replicate_table_3(workfile_china, col_2_controls, col_3_controls,
                  col_4_controls, col_5_controls, col_6_controls, G, N)

# Autor - distribution of permutation test 
# Full specification - column 6

# bonf.f.test.for.lm(data= china_1990, instrument_form= "instrument2000",
#                          instrument= "instrument2000",
#                          controls= col_6_controls[2:length(col_6_controls)],
#                          weights="timepwt48",
#                          variables_to_remove= c("timepwt48", "instrument2000", "outcome2000"),
#                          conditioned = T)



# Table: Autor p values
# Rows: using only their negative control, using all negative controls
# Cols: different controls as in the paper

specified_nco_for_sur <- c(
  "reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
  "reg_escen", "reg_wscen", "reg_mount", "reg_pacif",
  "l_no_workers_totcbp", "l_sh_popedu_c", "l_sh_popfborn",
  "l_sh_empl_f", "l_sh_routine33", "l_task_outsource", "l_sh_empl",
  "l_sh_empl_mfg", "l_sh_empl_mfg_m", "l_sh_empl_mfg_f", "l_sh_empl_mfg_edu_nc",
  "l_sh_empl_mfg_edu_c", "l_sh_empl_mfg_age1634", "l_sh_empl_mfg_age3549", "l_sh_empl_mfg_age5064",
  "l_sh_empl_nmfg_edu_nc", "l_sh_empl_nmfg_edu_c", "l_sh_empl_nmfg_age1634", "l_sh_empl_nmfg_age3549",
  "l_sh_empl_nmfg_age5064", "l_sh_unempl", "l_sh_unempl_m", "l_sh_unempl_edu_nc",
  "l_sh_unempl_edu_c", "l_sh_unempl_age1634", "l_sh_unempl_age3549", "l_sh_unempl_age5064",
  "l_sh_ssadiswkrs", "l_avg_lnwkwage_mfg", "l_trans_totindiv_pc", "l_trans_totmed_pc",
  "l_trans_fedinc_pc", "l_trans_othinc_pc", "l_trans_unemp_pc", "l_trans_taaimp_pc",
  "l_trans_totedu_pc", "l_trans_ssaret_pc", "l_trans_ssadis_pc", "l_avg_hhincsum_pc_pw",
  "l_avg_hhincwage_pc_pw", "outcome1990"
)

# run NCIV test for the IV used in China article with all the available NC 
china_all_ncs_p_values <- run_china_nciv(data= china_1990, 
              instrument_form= "instrument2000",instrument= "instrument2000",
              col_2_controls= col_2_controls, col_3_controls= col_3_controls,
              col_4_controls= col_4_controls,col_5_controls= col_5_controls,
              col_6_controls= col_6_controls, weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000", "outcome2000"),
              specified_nco = specified_nco_for_sur,
              specified_for_sur_only = T
               )


# outcome1990 is the only NC for the original NC test used in the article in 
# our interpretation. Run NCIV test for when the only NC is the outcome1990
china_original_nc_only_p_values <- run_china_nciv(data= china_1990_only_org_nc, instrument_form= "instrument2000",
               instrument= "instrument2000", col_2_controls= col_2_controls,
               col_3_controls= col_3_controls, col_4_controls= col_4_controls,
               col_5_controls= col_5_controls, col_6_controls= col_6_controls,
               weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000"),
               specified_nco= c("outcome1990"))

china_p_values_table <- data.frame(rbind(china_all_ncs_p_values, china_original_nc_only_p_values))
colnames(china_p_values_table) <- paste("Column ", 2:6)
rownames(china_p_values_table) <- c("all_ncs", "original_nc_only")

china_p_values_table_path <- file.path("out", "china_p_values_table.csv")
write.csv(china_p_values_table, china_p_values_table_path)

# Figure variable importance 
# Y axis - when “y” is the label
# X axis - when IV is the label 
# Using specification of column 6 

draw_importance_graph(china_1990, y= "outcome1990",
                      instrument = "instrument2000",
                      controls= col_6_controls[2:length(col_6_controls)],
                      all_controls= all_controls,
                      weights = "timepwt48",
                      variables_to_remove = c("timepwt48", "outcome1990", "instrument2000"),
                      "Autor et. al (2014) - China's Import Shock - Variables Importance",
                      FALSE)


# "Using School Choice Lotteries to Test Measures of School Effectiveness" (Deming)

# This file loads all the necessary data and functions for school article
# sourcing the file will result in loading the main functions:
# (1) get_formula_for_table_1_lottery() 
# which creates the formula for replicating table 1 in the article
# (2) get_formula_for_iv_by_NC
# which residualize the iv given the NCs
# As well as loading cms data:
# (1) cms - the "raw" data, the output of the Stata do file 
# (2) curr_cms - only rows that have the lottery_fe (fixed effects), which is the
# unit of randomization for the lottery

source("School Init.R")

# This file consists the function that runs NCIV test for few potential IVs in Deming
# sourcing the file will result in loading the following function:
# (1) run_school_nciv() - which runs a permutations test (both lm and felm)
# (2) run_school_variable_importance() - which runs variable importance for 
# the controls specification decalred in the article
# As well as the inner function:
# (*) prepare_school_data
source("School Run.R")

# Deming - histogram of of permutation test
# A - his instrument
# B - our instrument (without home school)
# C- pure lottery
# Using the specification that he uses in the paper
run_school_nciv(data= curr_cms, ivs= c("lottery", "lott_VA","new_lott_VA"),
                                   controls_specifications = 3, #the specification 
                                   #in the article
                                   saveplot= T,
                                   iterations=1, permutations=1000, ntree= 1000,
                                   mtry_ratio= 1/3, OOB=T, randomize_lottery=F)



# run potential IVS for lottery article
#Tabe: Deming p-value
# Rows: instruments (his IV, raw lottery results, corrected IV)
# Cols: different controls (none, lottery fixed effects, home school) 

specified_nco_for_school <- c(
  "math_2001_imp", "math_2001_miss", "read_2001_imp", "read_2001_miss",
  "ch1_mod2mix_all_test", "ch2_mod2mix_all_test", "ch3_mod2mix_all_test", "hm_mod2mix_all_test"
)

school_p_values <- run_school_nciv(data= curr_cms, 
                                   ivs= c("lottery", "lott_VA","new_lott_VA"),
                controls_specifications = 1:4,
                specified_nco= specified_nco_for_school)

full_file_name <-  file_full_path <- file.path("out", "school_p_values.csv")
write.csv(school_p_values, full_file_name)

#run sanity check (replace lottery with Bernoulli RV)
# run_school_nciv(data= curr_cms, ivs= c("lottery"),
#                 controls_specifications = 1:4,
#                 saveplot= F, iterations= 1,
#                 permutations=20, ntree=20,
#                 mtry_ratio= 1/3, OOB=T, randomize_lottery=T)

#run_school_variable_importance(curr_cms, c("lottery", "lott_VA","new_lott_VA"))





