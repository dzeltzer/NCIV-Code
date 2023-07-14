library(systemfit)

# create a data frame with iv, treatment, 2 unmeasured confounders and negative 
# control variables where the unmeasured confounder affects the iv in 
# a combination of linear and interaction between the two unmeasured confounders
# the negative control variable are associated linearly with the unmeasured confounder
create_interactions_multi_nc <- function(n =100,# the number of observations in
                                         # each iteration
                                         sg =1, # the standard variation for the unmeasured confounder
                                         inst.effect= 0.2, # the coefficient of the effect of the 
                                         # instrumental variable on the treatment
                                         number_of_first_good_ncs= 5, # The number of negative control 
                                         # variable associated with the first unmeasured confounder
                                         number_of_second_good_ncs= 5, # The number of negative control 
                                         # variable associated with the second unmeasured confounder
                                         number_of_bad_ncs = 50, # The number of independent variables 
                                         # in the data frame which are not associated with either unmeasured confounder
                                         nc_power=1, # the sum of coefficient of the effect of the
                                         # unmeasured confounder on the negative control variables
                                         nc_power_split= "uniform", # the distribution of the effects 
                                         # of the unmeasured confounder on the negative control variables
                                         # can be uniform or exponential
                                         alpha= 0, # hyper parameter: the fraction of the linear effect of the unmeasured
                                         # confounder on the negative control variables. when alpha=1
                                         # there is only linear effect, when alpha=0 there is only
                                         # effect by the interactions of two unmeasured confounder
                                         is_null_scenario = F
) {
  sg_y <- sg
  h1 <- rnorm(n, 0, sg)
  h2 <- rnorm(n, 0, sg)
  iv <- (1-is_null_scenario)*((alpha/2)*h1 + (alpha/2)*h2 + (1-alpha)*h1*h2) + rnorm(n, 0, sg)
  t <- inst.effect*iv+ rnorm(n, 0, sg)
  y <- t + (alpha/2)*h1 + (alpha/2)*h2 + (1-alpha)*h1*h2 + rnorm(n, 0, sg_y)
  
  nc_power_for_single_nc  <- switch(   
    nc_power_split, 
    "uniform"=  rep(nc_power/number_of_first_good_ncs, number_of_first_good_ncs),
    "exp"= nc_power * c(1/(2^(1:(number_of_first_good_ncs-1))), 1/(2^(number_of_first_good_ncs-1)))
  )
  
  good_ncs1 <- controls.matrix <- foreach (pl = 1:number_of_first_good_ncs, .combine ="cbind")  %do% {
    get_nc_col(nc_power_for_single_nc[pl], h1, n, sg)
  }
  
  good_ncs2 <- controls.matrix <- foreach (pl = 1:number_of_second_good_ncs, .combine ="cbind")  %do% {
    get_nc_col(nc_power_for_single_nc[pl], h2, n, sg)
  }
  
  bad_ncs <- controls.matrix <- foreach (pl = 1:number_of_bad_ncs, .combine ="cbind")  %do% {
    rnorm(n, 0, sg)
  }
  colnames(good_ncs1) <- paste0("good_nc", 1:number_of_first_good_ncs)
  colnames(good_ncs2) <- paste0("good_nc", (number_of_first_good_ncs+1):(number_of_first_good_ncs+number_of_second_good_ncs))
  
  colnames(bad_ncs) <- paste0("bad_nc", 1:number_of_bad_ncs)
  
  return(data.frame(h1, h2,iv, t,y, good_ncs1,good_ncs2, bad_ncs))
}

# create a data frame with iv, treatment, unmeasured confounder and negative 
# control variables where the unmeasured confounder affects the iv in 
# a combination of linear and squared way.
# the negative control variable are associated linearly with the unmeasured confounder
create_degree_nc <- function(n =100, # the number of observations in
                             # each iteration
                             sg =1, # the standard variation for the unmeasured confounder
                             sg_y= 1, # the standard variation for the outcome variable
                             inst.effect= 0.2, # the coefficient of the effect of the 
                             # instrumental variable on the treatment
                             number_of_good_ncs= 5, # The number of negative control 
                             # variable associated with the unmeasured confounder
                             number_of_bad_ncs = 50, # The number of independent variables 
                             # in the data frame which are not associated with either unmeasured confounder
                             nc_power=1, # the sum of coefficient of the effect of the
                             # unmeasured confounder on the negative control variables
                             nc_power_split= "uniform", # the distribution of the effects 
                             # of the unmeasured confounder on the negative control variables
                             # can be uniform or exponential
                             alpha= 0, # hyper parameter: the fraction of the linear effect of the unmeasured
                             # confounder on the negative control variables. when alpha=1
                             # there is only linear effect, when alpha=0 there is only
                             # effect by the squared unmeasured confounder
                             is_null_scenario = F
                             ) {
 
  h <- rnorm(n, 0, sg)
  iv <- (1-is_null_scenario)*((alpha)*h + (1-alpha)*(h^2)) + rnorm(n, 0, sg)
  t <- inst.effect*iv+ rnorm(n, 0, sg)
  y <- t + (alpha)*h + (1-alpha)*(h^2) + rnorm(n, 0, sg_y)
  
  nc_power_for_single_nc  <- switch(   
    nc_power_split, 
    "uniform"=  rep(nc_power/number_of_good_ncs, number_of_good_ncs),
    "exp"= nc_power * c(1/(2^(1:(number_of_good_ncs-1))), 1/(2^(number_of_good_ncs-1)))
  )
  
  good_ncs <- controls.matrix <- foreach (pl = 1:number_of_good_ncs, .combine ="cbind")  %do% {
    get_nc_col(nc_power_for_single_nc[pl], h, n, sg)
  }
  colnames(good_ncs) <- paste0("good_nc", 1:number_of_good_ncs)
  
  if (number_of_bad_ncs > 0){
    bad_ncs <- controls.matrix <- foreach (pl = 1:number_of_bad_ncs, .combine ="cbind")  %do% {
      rnorm(n, 0, sg)
    }
    
    colnames(bad_ncs) <- paste0("bad_nc", 1:number_of_bad_ncs)
    
    return(data.frame(h,iv, t,y, good_ncs, bad_ncs))
  }
  return(data.frame(h,iv, t,y, good_ncs))
  
  
  
}

get_nc_col <- function(nc_power_for_single_nc, h, n, sg){
  return(nc_power_for_single_nc*h + rnorm(n,0, sg))
}

# create a data frame with iv, treatment, unmeasured confounder and negative 
# control variables where the unmeasured confounder affects the iv in 
# as an output of a simplified CES function
# the negative control variable are associated linearly with the unmeasured confounder
create_ces_nc <- function(n =100,# the number of observations in
                          # each iteration
                          sg =1, # the standard variation for the unmeasured confounder
                          sg_y= 1, # the standard variation for the outcome variable
                          inst.effect= 0.2, # the coefficient of the effect of the 
                          # instrumental variable on the treatment
                          number_of_first_good_ncs= 5, # The number of negative control 
                          # variable associated with the first unmeasured confounder
                          number_of_second_good_ncs= 5, # The number of negative control 
                          # variable associated with the second unmeasured confounder 
                          number_of_bad_ncs = 50, # The number of independent variables 
                          # in the data frame which are not associated with either unmeasured confounder
                          nc_power=1, # the sum of coefficient of the effect of the
                          # unmeasured confounder on the negative control variables
                          nc_power_split= "uniform", # the distribution of the effects 
                          # of the unmeasured confounder on the negative control variables
                          # can be uniform or exponential
                          rho= 0, # hyper parameters: the value
                          # that determine the shape of CES function (described above)
                          is_null_scenario = F
                          )
  {
  
  k <- runif(n, 0, 2) #(n, 0, sg)#
  l <- runif(n, 0, 2) #(n, 0, sg)#runif(n,-2,2)#
  iv <- (1-is_null_scenario)*ces_q(k, l, rho) + rnorm(n, 0, sg)
  t <- inst.effect*iv+ rnorm(n, 0, sg)
  y <- t + ces_q(k, l, rho) + rnorm(n, 0, sg_y)
  
  nc_power_for_single_nc  <- switch(   
    nc_power_split, 
    "uniform"=  rep(nc_power/number_of_first_good_ncs, number_of_first_good_ncs),
    "exp"= nc_power * c(1/(2^(1:(number_of_first_good_ncs-1))), 1/(2^(number_of_first_good_ncs-1)))
  )
  
  good_ncs1 <- controls.matrix <- foreach (pl = 1:number_of_first_good_ncs, .combine ="cbind")  %do% {
    get_nc_col(nc_power_for_single_nc[pl], k, n, sg)
  }
  
  good_ncs2 <- controls.matrix <- foreach (pl = 1:number_of_second_good_ncs, .combine ="cbind")  %do% {
    get_nc_col(nc_power_for_single_nc[pl], l, n, sg)
  }
  
  bad_ncs <- controls.matrix <- foreach (pl = 1:number_of_bad_ncs, .combine ="cbind")  %do% {
    rnorm(n, 0, sg)
  }
  colnames(good_ncs1) <- paste0("good_nc", 1:number_of_first_good_ncs)
  colnames(good_ncs2) <- paste0("good_nc", (number_of_first_good_ncs+1):(number_of_first_good_ncs+number_of_second_good_ncs))
  
  colnames(bad_ncs) <- paste0("bad_nc", 1:number_of_bad_ncs)
  
  return(data.frame(k, l ,iv, t,y, good_ncs1,good_ncs2, bad_ncs))
}

ces_q <- function(k, l, rho)
{
  if(rho == 0){
    return(1)
  }
  return((k^rho + l^rho)^(1/rho))
}

# convert 2-D p-values result matrix of 2-D (alpha * number_of_good_ncs) matrix
# to a readable long CSV format
prepare_simulations_results <- function(results_raw, #2-D (alpha * number_of_good_ncs) matrix
                                        # of the rejection fraction for the below specification
                                        columns, # The columns of the CSV output file
                                        scenario, # the scenario (Interactions, 
                                        # Degrees, CES ...) used 
                                        algo, # the prediction algorithm (Single lm p-val (Bonferroni),  
                                        # F-test, RF, CES ...) used 
                                        single_nc_power, # the coefficient of the effect of
                                        # the unmeasured confounder on the negative controls 
                                        alpha_values, # fraction of the linear effect of the unmeasured
                                        # confounder on the negative control variables. when alpha=1
                                        # there is only linear effect, when alpha=0 there is only
                                        # effect by the squared unmeasured confounder
                                        number_of_good_ncs_values, # number of negative control variables
                                        # (out of number_of_all_ncs_value)
                                        n_value, # the number of observations in
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
                                        ntree, # the number of trees in the RF prediction
                                        # algorithm used for NCIV test
                                        is_null_scenario,
                                        sg,
                                        constant_snr_sgs_mat = NULL # if needed - matrix of coefficients
                                        # for keeping constant SNR
                                        ) {
  

  colnames(results_raw) <- number_of_good_ncs_values
  results <- as.data.frame(cbind(alpha_values, results_raw)) %>% 
    pivot_longer(cols = matches("[0-9]+"), names_to = "number_of_good_ncs_values",
               values_to = "rejection_rate") %>%
    mutate(number_of_good_ncs_values  = as.integer(number_of_good_ncs_values ),
           scenario = scenario,
           algo = algo,
           n_value = n_value,
           number_of_all_ncs_value = number_of_all_ncs_value,
           n_iterations = n_iterations,
           n_permutations = n_permutations, 
           single_nc_power = single_nc_power,
           sg= sg,
           ntree = ntree,
           is_null_scenario= is_null_scenario) %>% 
    relocate(rejection_rate)
  
if (!(is.null(constant_snr_sgs_mat)))
  {
  sgs_power_mat <- as.data.frame(cbind(alpha_values, constant_snr_sgs_mat)) %>%  
    pivot_longer(
      cols = matches("[0-9]+"), names_to = "number_of_good_ncs_values",
      values_to = "sg") %>% 
    mutate(number_of_good_ncs_values  = as.integer(number_of_good_ncs_values ))
  
  results <- results %>% 
    select(-"sg") %>% 
    inner_join(sgs_power_mat, by= c("alpha_values", "number_of_good_ncs_values"))
  }

  return(results)
}
