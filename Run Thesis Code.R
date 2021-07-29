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

cl <- makeCluster(detectCores())
registerDoParallel(cl)

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
# (*) run_school_nciv() - which runs a permutations test (both lm and felm)
source("School Run.R")

# run potential IVS for lottery article
run_school_nciv(data= curr_cms, ivs= c("lottery", "lott_VA","new_lott_VA"),
                iterations=1, permutations=5, ntree= 100,
                mtry_ratio= 1/3, OOB=T, randomize_lottery=F)

#run sanity check (replace lottery with Bernoulli RV)
# run_school_nciv(data= curr_cms, ivs= c("lottery"), iterations= 1,
#                 permutations=20, ntree=20,
#                 mtry_ratio= 1/3, OOB=T, randomize_lottery=T)

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

# run NCIV test for the IV used in China article with all the available NC 
run_china_nciv(data= china_1990, instrument_form= "instrument2000",
               instrument= "instrument2000", col_2_controls= col_2_controls,
               col_3_controls= col_3_controls, col_4_controls= col_4_controls,
               col_5_controls= col_5_controls, col_6_controls= col_6_controls,
               weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000", "outcome2000"),
               permutations= 100, OOB= T, mtry_ratio= 1/3, ntree= 100, title= "China")


# outcome1990 is the only NC for the original NC test used in the article in 
# our interpretation. Run NCIV test for when the only NC is the outcome1990
run_china_nciv(data= china_1990_only_org_nc, instrument_form= "instrument2000",
               instrument= "instrument2000", col_2_controls= col_2_controls,
               col_3_controls= col_3_controls, col_4_controls= col_4_controls,
               col_5_controls= col_5_controls, col_6_controls= col_6_controls,
               weights= "timepwt48", 
               variables_to_remove= c("timepwt48", "instrument2000"),
               permutations= 50, OOB= T, mtry_ratio= 1, ntree= 100, title= "China Orginal")

# Simulations
# this files file loads all the necessary functions for the simulations: Data generations, and benchmark
# tests algorithms
source("Data Generators.R")

source("Interactions Simulations.R")

run_interactions_simulations(n_value= 200, number_of_all_ncs_value = 10,
                             n_iterations = 2, n_permutations = 12,
                             nc_power_value=1,rejection_rate = 0.05,
                             number_of_good_ncs_values = c(2, 3),
                             alpha_values = c(0, 1),
                             singal_nc_power = 0.5,
                             ntree = 10)

source("Degrees Simulations.R")

run_degrees_simulations(n_value= 200, number_of_all_ncs_value = 10,
                        n_iterations = 2, n_permutations = 12,
                        nc_power_value=1,rejection_rate = 0.05,
                        number_of_good_ncs_values = c(2, 3),
                        alpha_values = c(0, 1),
                        singal_nc_power = 0.5,
                        ntree = 10)




