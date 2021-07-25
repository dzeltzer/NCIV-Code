library(systemfit)

create_multi_nc <- function(n =100, sg =1, sg_y= 1, inst.effect= 0.2, number_of_good_ncs= 5, 
                            number_of_bad_ncs = 50,
                            nc_power=1, nc_power_split= "uniform") {
  
  h <- rnorm(n, 0, sg)
  iv <- h + rnorm(n, 0, sg)
  t <- inst.effect*iv+ rnorm(n, 0, sg)
  y <- t + h + rnorm(n, 0, sg_y)
  
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



get_nc_col <- function(nc_power_for_single_nc, h, n, sg){
  return(nc_power_for_single_nc*h + rnorm(n,0, sg))
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


ces_q <- function(k, l, rho)
{
  return((k^rho + l^rho)^(1/rho))
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




get_oracle_p_val <- function(data, number_of_good_ncs, number_of_bad_ncs){
  controls <- paste0("good_nc", 1:number_of_good_ncs)
  s1 <- as.formula(sprintf("iv ~ %s", paste(controls, collapse="+")))
  nc_fit <- summary(lm(data = data, formula = s1))
  return(pf(nc_fit$fstatistic[1],nc_fit$fstatistic[2],nc_fit$fstatistic[3],lower.tail=FALSE))
  
}

get_f_p_val <- function(data, number_of_good_ncs, number_of_bad_ncs){
  controls <- paste0("good_nc", 1:number_of_good_ncs)
  if (number_of_bad_ncs > 0){
    controls <- c(controls,  paste0("bad_nc", 1:number_of_bad_ncs))
  }
  s1 <- lm("iv ~ t", data= data)
  data$iv_res <- s1$residuals
  s2 <- as.formula(sprintf("iv ~ %s", paste(controls, collapse="+")))
  nc_fit <- summary(lm(data = data, formula = s2))
  return(pf(nc_fit$fstatistic[1],nc_fit$fstatistic[2],nc_fit$fstatistic[3],lower.tail=FALSE))
  
}

get_sur_p_val <- function(data, number_of_good_ncs, number_of_bad_ncs){
  
  system_good <- lapply(paste0("good_nc", 1:number_of_good_ncs, " ~ iv"), as.formula)
  system_bad <- lapply(paste0("bad_nc", 1:number_of_bad_ncs, " ~ iv"), as.formula)
  
  fitsur <- systemfit( c(system_good, system_bad), method = "SUR", data=data)
  
  restrict <- paste0("eq", 1:(number_of_good_ncs+number_of_bad_ncs), "_iv=0")
  
  ft <- linearHypothesis( fitsur, restrict,  test = "F" )
  
  return(ft$`Pr(>F)`[2])
  
}

get_min_p_val <- function(data, number_of_good_ncs, number_of_bad_ncs){
  min_good <- min_bad <- 1
  if(number_of_good_ncs >0){
    min_good <- 
      min(foreach (curr_nc = 1:number_of_good_ncs, .combine ="cbind") %do% 
            {
              s1 <- lm("iv ~ t", data= data)
              data$iv_res <- s1$residuals
              s2 <- as.formula(sprintf("iv_res ~ good_nc%s", curr_nc))
              nc_fit <- summary(lm(data = data, formula = s2))
              pf(nc_fit$fstatistic[1],nc_fit$fstatistic[2],nc_fit$fstatistic[3],lower.tail=FALSE)
              
            })
  }
  if (number_of_bad_ncs >0){
    min_bad <- 
      min(foreach (curr_nc = 1:number_of_bad_ncs, .combine ="cbind") %do% 
            {
              
              s1 <- as.formula(sprintf("iv ~ bad_nc%s", curr_nc))
              nc_fit <- summary(lm(data = data, formula = s1))
              pf(nc_fit$fstatistic[1],nc_fit$fstatistic[2],nc_fit$fstatistic[3],lower.tail=FALSE)
              
            })
  }
  
  
  return(min(min_good, min_bad))
  
}

prepare_results <- function(scenario, number_of_good_ncs_values, alpha_values, nc_power, algo, data) {
  results <- data.frame(scenario= character(),
                        nc_power= character(),
                        alpha= character(), 
                        number_of_good_ncs= character(), 
                        algo= character(),
                        rejection_rate=character() ,
                        stringsAsFactors=FALSE) 
  for(number_of_good_ncs_ix in 1:length(number_of_good_ncs_values)){
    for (alpha_ix in 1:length(alpha_values)){
      line <-  data.frame(scenario, nc_power, alpha_values[alpha_ix], 
                          number_of_good_ncs_values[number_of_good_ncs_ix], algo,
                          data[alpha_ix, number_of_good_ncs_ix])
      results <- rbind(results, line)
    }
  }
  return(results)
}