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

stopCluster(cl)
cl <- makeCluster(detectCores())
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


# Simulations
# this files file loads all the necessary functions for the simulations: Data generations, and benchmark
# tests algorithms
source("Data Generators.R")

source("Benchmark Algorithms.R")

source("Interactions Simulations.R")

run_interactions_simulations(n_value= 200, number_of_all_ncs_value = 10,
                             n_iterations = 5, n_permutations = 12,
                             rejection_rate = 0.05,
                             number_of_good_ncs_values = c(2, 3, 4),
                             alpha_values = c(0, 1),
                             single_nc_power = 0.5,
                             ntree = 10)

run_interactions_simulations(n_value= 200, number_of_all_ncs_value = 10,
                             n_iterations = 5, n_permutations = 12,
                             rejection_rate = 0.05,
                             number_of_good_ncs_values = c(2, 3, 4),
                             alpha_values = c(0, 1),
                             single_nc_power = 0,
                             ntree = 10)

source("Degrees Simulations.R")

run_degrees_simulations(n_value= 200, number_of_all_ncs_value = 10,
                        n_iterations = 5, n_permutations = 12,
                        rejection_rate = 0.05,
                        number_of_good_ncs_values = c(2, 3, 4),
                        alpha_values = c(0, 1),
                        single_nc_power = 0.5,
                        ntree = 10)

run_degrees_simulations(n_value= 200, number_of_all_ncs_value = 10,
                        n_iterations = 5, n_permutations = 12,
                        rejection_rate = 0.05,
                        number_of_good_ncs_values = c(2, 3, 4),
                        alpha_values = c(0, 1),
                        single_nc_power = 0,
                        ntree = 10)

source("CES Simulations.R")

run_ces_simulations(n_value= 200, number_of_all_ncs_value = 10,
                        n_iterations = 5, n_permutations = 12,
                        rejection_rate = 0.05,
                        number_of_good_ncs_values = c(2, 3, 4),
                        rho_values = c(-100, 0, 0.5),
                        single_nc_power = 0.5,
                        ntree = 10)

run_ces_simulations(n_value= 200, number_of_all_ncs_value = 10,
                    n_iterations = 5, n_permutations = 12,
                    rejection_rate = 0.05,
                    number_of_good_ncs_values = c(2, 3, 4),
                    rho_values = c(-100, 0, 0.5),
                    single_nc_power = 0,
                    ntree = 10)
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

permutations.test.for.lm(data= china_1990, instrument_form= "instrument2000",
                         instrument= "instrument2000",
                         controls= col_6_controls[2:length(col_6_controls)],
                         weights="timepwt48",
                         variables_to_remove= c("timepwt48", "instrument2000", "outcome2000"),
                         title = "China column 6", n_permutations= 100,
                         conditioned = T, OOB= T, saveplot=T,
                         mtry_ratio=1/3, ntree=100)

# Table: Autor p values
# Rows: using only their negative control, using all negative controls
# Cols: different controls as in the paper

# run NCIV test for the IV used in China article with all the available NC 
china_all_ncs_p_values <- run_china_nciv(data= china_1990, 
              instrument_form= "instrument2000",instrument= "instrument2000",
              col_2_controls= col_2_controls, col_3_controls= col_3_controls,
              col_4_controls= col_4_controls,col_5_controls= col_5_controls,
              col_6_controls= col_6_controls, weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000", "outcome2000"),
               permutations= 100, OOB= T, mtry_ratio= 1/3, ntree= 100, title= "China",
              saveplot= F)


# outcome1990 is the only NC for the original NC test used in the article in 
# our interpretation. Run NCIV test for when the only NC is the outcome1990
china_original_nc_only_p_values <- run_china_nciv(data= china_1990_only_org_nc, instrument_form= "instrument2000",
               instrument= "instrument2000", col_2_controls= col_2_controls,
               col_3_controls= col_3_controls, col_4_controls= col_4_controls,
               col_5_controls= col_5_controls, col_6_controls= col_6_controls,
               weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000"),
               permutations= 50, OOB= T, mtry_ratio= 1, ntree= 100, title= "China Orginal",
               saveplot= F)

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
                      "China column 6 spec")


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
                                   iterations=1, permutations=20, ntree= 100,
                                   mtry_ratio= 1/3, OOB=T, randomize_lottery=F)



# run potential IVS for lottery article
#Tabe: Deming p-value
# Rows: instruments (his IV, raw lottery results, corrected IV)
# Cols: different controls (none, lottery fixed effects, home school) 

school_p_values <- run_school_nciv(data= curr_cms, 
                                   ivs= c("lottery", "lott_VA","new_lott_VA"),
                controls_specifications = 1:4,
                saveplot= F,
                iterations=1, permutations=20, ntree= 100,
                mtry_ratio= 1/3, OOB=T, randomize_lottery=F)

full_file_name <-  file_full_path <- file.path("out", "school_p_values.csv")
write.csv(school_p_values, full_file_name)

#run sanity check (replace lottery with Bernoulli RV)
# run_school_nciv(data= curr_cms, ivs= c("lottery"),
#                 controls_specifications = 1:4,
#                 saveplot= F, iterations= 1,
#                 permutations=20, ntree=20,
#                 mtry_ratio= 1/3, OOB=T, randomize_lottery=T)

run_school_variable_importance(curr_cms, c("lottery", "lott_VA","new_lott_VA"))





