library(systemfit)

create_degree_nc <- function(n =100, sg =1, sg_y= 1, inst.effect= 0.2, number_of_good_ncs= 5, 
                             number_of_bad_ncs = 50,
                             nc_power=1, nc_power_split= "uniform",
                             alpha= 0) {
  
  h <- rnorm(n, 0, sg)
  iv <- (alpha)*h + (1-alpha)*(h^2) + rnorm(n, 0, sg)
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

create_interactions_multi_nc <- function(n =100, sg =1, sg_y= 1, inst.effect= 0.2,
                                         number_of_first_good_ncs= 5, 
                                         number_of_second_good_ncs= 5, 
                                         number_of_bad_ncs = 50,
                                         nc_power=1, nc_power_split= "uniform",
                                         alpha= 0) {
  
  h1 <- rnorm(n, 0, sg)
  h2 <- rnorm(n, 0, sg)
  iv <- (alpha/2)*h1 + (alpha/2)*h2 + (1-alpha)*h1*h2 + rnorm(n, 0, sg)
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

get_nc_col <- function(nc_power_for_single_nc, h, n, sg){
  return(nc_power_for_single_nc*h + rnorm(n,0, sg))
}

create_ces_nc <- function(n =100, sg =1, sg_y= 1, inst.effect= 0.2,
                          number_of_first_good_ncs= 5, 
                          number_of_second_good_ncs= 5, 
                          number_of_bad_ncs = 50,
                          nc_power=1, nc_power_split= "uniform",
                          rho= 0) {
  
  k <- rnorm(n, 0, sg)#
  l <- rnorm(n, 0, sg)#runif(n,-2,2)#
  iv <- ces_q(k, l, 1)
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
  return((k^rho + l^rho)^(1/rho))
}

prepare_simulations_results <- function(results_raw, columns, scenario, algo,
                                        single_nc_power, alpha_values,
                                        number_of_good_ncs_values, n_value,
                                        number_of_all_ncs_value, n_iterations,
                                        n_permutations, ntree) {
  
  results <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(results) <- columns
  
  for (i in 1:length(alpha_values)){
    for (j in 1:length(number_of_good_ncs_values)){
      
      line <- c(scenario, algo, n_value, number_of_all_ncs_value, n_iterations,
                n_permutations, number_of_good_ncs_values[j],
                alpha_values[i], single_nc_power, ntree,
                results_raw[i,j])
      results <- rbind(results, line)
    }
  }
  return(as.matrix(results))
}