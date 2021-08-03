
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
