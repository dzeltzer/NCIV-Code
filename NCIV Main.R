library(foreach)
library(doParallel)
library(car)
library(AER)
library(ggplot2)
library(randomForest)
library(tidyselect)
library(tidyverse)
library(rlang)
library(lfe)
library(stringr)


# run permutations test with fixed effects model. The function first residualize
# the iv and the controls according to the specified fixed effects model and 
# controls
# input: data, releavant formulas parameters for the permutations tests and for 
# the RF model
# output: permutations test p-value for the hypothesis that the IV can be
# predicted from the NCs
permutations.test.for.felm <- function(data, ins_lfe_formula, instrument,
                                       controls, fixed, variables_to_remove,
                                       title, n_permutations= 5, conditioned = T,
                                       OOB=F, saveplot=T,  mtry_ratio, ntree) {
  
  
  if (conditioned){
    felm_of_controls <-  felm(formula= ins_lfe_formula, data = data)
    z <- as.vector(felm_of_controls$residuals)
  }
  else{
    z <- pull(data, instrument)
  }
  NC <- get_NC_matrix_felm(select(data,-one_of(variables_to_remove)),  controls, fixed, T)
  
  controls_results <-run.rf.multiple.negative.controls(NC, z, F,
                                                                OOB,  mtry_ratio, ntree)
  
  permutations_results <- 
    foreach (ite = 1:n_permutations,.packages = c("lfe","dplyr","foreach"),
             .export = c("run.rf.multiple.negative.controls", "get_NC_matrix_felm", "rfcv"),
             .combine ="c")  %dopar% {
               run.rf.multiple.negative.controls(NC, z, T,
                                                          OOB,  mtry_ratio, ntree)
             }
  title <- sprintf("IV test: method= felm %s\nconditiond=%s, ntree=%s, OOB=%s", title, conditioned, ntree, OOB)
  
  return(calculate.visualize.p.values(controls_results, permutations_results ,title, OOB,
                                      saveplot))
}

# run single RF model on the residualiztion results
# input: data and parameters for the RF model
# output: MSE of the model IV = f(NCs)
run.rf.multiple.negative.controls <- function(NC, z, random=F, OOB,mtry_ratio=1/3, ntree){
  
  
  if (random){
    z <- sample(z)
  }
  
  if (!OOB){
    cv_result <- rfcv(NC, z, cv.fold=5,  step=ncol(NC), mtry=function(p) {max(1, floor(p/3))})
    return(cv_result$error.cv[1])
  }
  
  rf_model <- randomForest::randomForest(x= NC, y=z, ntree= ntree, mtry= mtry_ratio*ncol(NC))
  return(rf_model$mse[length(rf_model$mse)])
}

# residualize the NCs by the controls with fixed effects model
# input: data and the names of the control variables
# output: residuals of the NCs given the controls
get_NC_matrix_felm <- function(data, controls, fixed, conditioned = T) {
  assertthat::not_empty(data)
  assertthat::not_empty(controls)
  
  if (length(controls) == 1 & controls[1] == c("1")){
    return(data)
  }
  
  X <- data %>%
    # select(all_of(controls))
    select(one_of(controls, fixed))
  NC_raw <- data %>%
    select(-one_of(controls, fixed))
  
  if (!conditioned){
    return(NC_raw)
  }
  
  y <- as.vector(NC_raw[,1])
  names(y) <- "y"
  
  colnames(X)[1:(length(X))] <- paste0("x", 1:length(X))
  colnames(X)[length(X)] <- "fe"
  formula <- as.formula(sprintf("y~ %s| fe", paste(paste0("x", sep= 1:(length(X)-1)), collapse = " + ")))
  
  NC_res <- data.frame(felm(formula= formula, data= cbind(y, X))$residuals)  
  for (i in 2:ncol(NC_raw)) {
    y <- as.vector(NC_raw[,i])
    names(y) <- "y"
    colnames(X)[1:(length(X))] <- paste0("x", 1:length(X))
    colnames(X)[length(X)] <- "fe"
    NC_res <- cbind(NC_res, felm(formula= formula, data= cbind(y, X))$residuals) 
  }
  colnames(NC_res) <- colnames(NC_raw)
  return(NC_res)
}

# run permutations test with linear regression model. The function first residualize
# the iv and the controls according to the specified linear regression model and 
# controls
# input: data, releavant formulas parameters for the permutations tests and for 
# the RF model
# output: permutations test p-value for the hypothesis that the IV can be
# predicted from the NCs
permutations.test.for.lm <- function(data, instrument_form, instrument, controls,
                                     weights = NULL, variables_to_remove, title, n_permutations= 5, 
                                     conditioned = T, OOB= F, saveplot, mtry_ratio, ntree) {
  ins_formula <- as.formula(paste(instrument_form, "~", paste(controls, collapse="+")))
  lm_of_controls <-  lm(as.formula(ins_formula), data = data)
  if (!is.null(weights)){
    lm_of_controls <- rlang::eval_tidy(rlang::quo(
      lm(
        as.formula(ins_formula),
        data = data,
        weights = !!as.name(weights))
    ))
  }
  
  z <- lm_of_controls$residuals
  if (!conditioned){
    z <- pull(data, instrument)
  }
  
  NC <- get_NC_matrix(select(data,-one_of(variables_to_remove)),  controls, conditioned)
  
  controls_results <- run.rf.multiple.negative.controls(NC, z, F,
                                                               OOB, mtry_ratio, ntree) 
  
  permutations_results <- 
    foreach (ite = 1:n_permutations,.packages = c("dplyr","foreach"),
             .export = c("run.rf.multiple.negative.controls", "get_NC_matrix", "rfcv"),
             .combine ="c")  %dopar% {
               run.rf.multiple.negative.controls(NC, z, T,
                                                        OOB, mtry_ratio, ntree) 
             }
  title <- sprintf("IV test: method= lm %s\nconditiond=%s, OOB=%s, ntree=%s", title, conditioned, OOB, ntree)
  
  return(calculate.visualize.p.values(controls_results, permutations_results,
                                      title, OOB, saveplot))
}


# residualize the NCs by the controls with linear regression model
# input: data and the names of the control variables
# output: residuals of the NCs given the controls
get_NC_matrix <- function(data, controls, conditioned =T) {
  assertthat::not_empty(data)
  assertthat::not_empty(controls)
  
  if (length(controls) == 1 & controls[1] == c("1")){
    return(data)
  }
  
  
  X <- data %>%
    select(controls)
  NC_raw <- data %>%
    select(-one_of(controls))
  
  if (!conditioned){
    return(NC_raw)
  }
  
  y <- as.vector(NC_raw[,1])
  names(y) <- "y"
  NC_res <- data.frame(lm(y~., data= cbind(y, X))$residuals)  
  if (ncol(NC_raw) > 1){
    for (i in 2:ncol(NC_raw)) {
      y <- as.vector(NC_raw[,i])
      names(y) <- "y"
      NC_res <- cbind(NC_res, lm(y~., data= cbind(y, X))$residuals)
    }
  }
  
  
  colnames(NC_res) <- colnames(NC_raw)
  return(NC_res)
}


# present the results of the permutations test
# input: the MSE values for the results, and the corresponding permutations
# output: the p-value of the test
calculate.visualize.p.values <- function(controls_results, permutations_results, title, OOB, saveplot = F) {
  n <- length(permutations_results)
  if (saveplot){
    Xaxis <- "RMSE - CV"
    if (OOB){
      Xaxis <- "RMSE - OOB"
    }
    results <- data.frame(cbind(rep("permutaions", n), as.numeric(permutations_results)))
    binW <- (max(permutations_results)- min(permutations_results))/15
    if (binW == 0){
      binW = 0.00001
    }
    p <- ggplot(data= results)+
      geom_histogram(aes(x=permutations_results), binwidth = binW) +
      geom_vline(xintercept = controls_results, color = "blue") +
      xlab(Xaxis) +
      ggtitle(title)
    print(p)
    file_title <- stringr::str_trim(sprintf("out//%s",gsub("[.]|[:]|[\n]|[,]", " ",title)))
    ggsave(sprintf("%s.png", file_title), p)
    results <- data.frame(permutations_results, 
                          controls_results= c(controls_results, rep(NA, length(permutations_results)-1)))
    write.csv(results, sprintf("%s.csv", file_title))
  }
  
  
  p_values <- c(t.test(permutations_results, mu= controls_results, alternative = "greater")$p.value,
                sum(controls_results > permutations_results) /n
  )
  print(sprintf("P-values for %s:", title))
  print(sprintf("t-test: %s, permutations: %s", p_values[1], p_values[2]))
  return(p_values[2])
}

# get heteroscedasticity consistent estimates for linear model
# input: linear regression model, number of groups, observations 
# output: heteroscedasticity consistent estimates
get_stata_coef <- function(G, N, reg_model, hc_type= "HC2"){
  dfa <- (G/(G - 1)) * (N - 1)/reg_model$df.residual
  
  firm_c_vcov <- dfa * vcovHC(reg_model, type = hc_type, cluster = "group", adjust = T)
  return(coeftest(reg_model, vcov = firm_c_vcov))
}

# draw variable importance graph for the models: Y= f(NCs) and IV = f(NCs)
# given the controls using linear regression model
# input: data and names of Y, IV and controls
# output: graph
draw_importance_graph <- function(data, y, instrument, controls, all_controls,
                                  weights = NULL, variables_to_remove, title)  {
  
  ins_formula <- as.formula(paste(instrument, "~", paste(controls, collapse="+"))) 
  lm_of_controls <- rlang::eval_tidy(rlang::quo(
    lm(
      ins_formula,
      data = data,
      weights = !!as.name(weights))
  ))
  instrument_res <- lm_of_controls$residuals
  
  NC <- get_NC_matrix(select(data,-one_of(variables_to_remove)),  controls, T)
  
  model_ins <- randomForest::randomForest(NC, instrument_res)
  model_y <- randomForest::randomForest(NC, y)
  
  vi_scores_ins_raw <- importance(model_ins, type= 2, scale =T)
  vi_scores_ins <- vi_scores_ins_raw / sd(vi_scores_ins_raw)
  vi_scores_ins <- data.frame(variable=rownames(vi_scores_ins), vi_scores_ins)
  vi_scores_y_raw <- importance(model_y, type= 2, scale =T)
  vi_scores_y <- vi_scores_y_raw / sd(vi_scores_y_raw)
  vi_scores_y <- data.frame(variable=rownames(vi_scores_y), vi_scores_y)
  
  vi_scores <- vi_scores_ins %>%
    left_join (vi_scores_y,
               by = "variable")
  colnames(vi_scores)[2:3] <- c("ins_score", "y_score")
  
  
  p <- ggplot(vi_scores, aes(ins_score, y_score, label = variable,
                             color = variable %in% all_controls))+
    geom_point() +
    ggtitle(title)
  plotly::ggplotly(p)
}


# draw variable importance graph for the models: Y= f(NCs) and IV = f(NCs)
# given the controls using fixed effects model
# input: data and names of Y, IV and controls
# output: graph
draw_importance_graph_fle <- function(data, y, ins_lfe_formula,
                                      y_formula,
                                      controls, all_controls, fixed,
                                      variables_to_remove, title)  {
  
  felm_of_controls <-  felm(formula= ins_lfe_formula, data = data)
  instrument_res <- felm_of_controls$residuals
  
  y_res_felm <- felm(formula= y_formula, data = data)
  y_res <- y_res_felm$residuals
  
  NC <- get_NC_matrix_felm(select(data,-one_of(variables_to_remove)),  controls, fixed)
  
  model_ins <- randomForest::randomForest(NC, instrument_res)
  model_y <- randomForest::randomForest(NC, y_res)
  
  vi_scores_ins_raw <- importance(model_ins, type= 2, scale =T)
  vi_scores_ins <- vi_scores_ins_raw / sd(vi_scores_ins_raw)
  vi_scores_ins <- data.frame(variable=rownames(vi_scores_ins), vi_scores_ins)
  vi_scores_y_raw <- importance(model_y, type= 2, scale =T)
  vi_scores_y <- vi_scores_y_raw / sd(vi_scores_y_raw)
  vi_scores_y <- data.frame(variable=rownames(vi_scores_y), vi_scores_y)
  
  vi_scores <- vi_scores_ins %>%
    left_join (vi_scores_y,
               by = "variable")
  colnames(vi_scores)[2:3] <- c("ins_score", "y_score")
  
  
  p <- ggplot(vi_scores, aes(ins_score, y_score, label = variable,
                             color = variable %in% all_controls))+
    geom_point() +
    ggtitle(title)
  plotly::ggplotly(p)
}
