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
permutations.test.for.felm <- function(data,
                                       ins_lfe_formula, #the fixed effects formula
                                       # for the iv
                                       instrument, # the name of the instrument
                                       controls, # the names of the controls
                                       fixed, # the name of the fixed effects groups column
                                       variables_to_remove, # the names of the variables 
                                       title, #the title for the output file and graph
                                       n_permutations= 5, # the number of permutations in each 
                                       # NCIV test
                                       conditioned = T, #whether to residualize on the controls
                                       OOB=F, # use Out Of Bag RMSE instead of Cross Validation
                                       # error
                                       saveplot=T, # whether to save plot and CSV of the 
                                       # permutations test p
                                       mtry_ratio, #the ratio of variables to check in 
                                       # each iteration of RF algorithm
                                       ntree # the number of trees in the RF prediction
                                       # algorithm used for NCIV test
                                       ) {
  
  
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
run.rf.multiple.negative.controls <- function(NC, #matrix of the negative controls
                                              z, # the IV for prediction
                                              random=F, # whether to randomize the IV
                                              OOB,
                                              mtry_ratio=1/3,
                                              ntree){
  
  
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
get_NC_matrix_felm <- function(data, #matrix of the data (controls+NCs)
                               controls, # names of the controls
                               fixed, # the name of the fixed effects groups
                               conditioned = T
                               ) {
  assertthat::not_empty(data)
  assertthat::not_empty(controls)
  
  if (length(controls) == 1 & controls[1] == c("1")){
    return(data)
  }
  
  # filter in only controls+fe
  X <- data %>%
    select(one_of(controls, fixed))
  
  # filter in only NCs
  NC_raw <- data %>%
    select(-one_of(controls, fixed))
  
  if (!conditioned){
    return(NC_raw)
  }
  
  # for each NC: residualize the NC given the controls using fixed effects model
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


permutations.test.for.lm <- function(data, 
                                     instrument_form, # the form of the instrument, 
                                     # for example: logarithmic: log(instrument2000), 
                                     # linear: instrument2000
                                     instrument, # the name of the instrument
                                     controls,# the names of the controls
                                     weights = NULL, #the weights column name used for the residualization
                                     variables_to_remove, # the names of the variables not to used
                                     # as NCs
                                     title, #the title for the output file and graph
                                     n_permutations= 5, # the number of permutations in each 
                                     # NCIV test
                                     conditioned = T, #whether to residualize on the controls
                                     OOB= F, # use Out Of Bag RMSE instead of Cross Validation error
                                     saveplot, # whether to save plot and CSV of the 
                                     # permutations test
                                     mtry_ratio, # the ratio of variables to check in 
                                     # each iteration of RF algorithm
                                     ntree # the number of trees in the RF prediction
                                     # algorithm used for NCIV test
                                     ) {
  ins_formula <- as.formula(paste(instrument_form, "~", paste(controls, collapse="+")))
  lm_of_controls <-  lm(as.formula(ins_formula), data = data)
  weights_vec <- NULL
  if (!is.null(weights)){
    lm_of_controls <- rlang::eval_tidy(rlang::quo(
      lm(
        as.formula(ins_formula),
        data = data,
        weights = !!as.name(weights))
    ))
    weights_vec <- pull(data, weights)
    
  }
  
  z <- lm_of_controls$residuals
  
  if (!conditioned){
    z <- pull(data, instrument)
  }
  NC <- get_NC_matrix(select(data,-one_of(variables_to_remove)),  controls, weights_vec, conditioned)
  
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
get_NC_matrix <- function(data, #matrix of the data (controls+NCs)
                          controls, # names of the controls
                          weights_vec, #the weights vector used for residualization
                          conditioned =T #whether to residualize NCs on the controls
                          ) {
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
  
  # for each NC: residualize the NC given the controls using linear regression
  # model (with weights if declared)
  y <- as.vector(NC_raw[,1])
  names(y) <- "y"
  NC_res <- data.frame(lm(y~., weights = weights_vec, data= cbind(y, X))$residuals)  
  if (ncol(NC_raw) > 1){
    for (i in 2:ncol(NC_raw)) {
      y <- as.vector(NC_raw[,i])
      names(y) <- "y"
      NC_res <- cbind(NC_res, lm(y~., weights = weights_vec, data= cbind(y, X))$residuals)
    }
  }
  
  colnames(NC_res) <- colnames(NC_raw)
  return(NC_res)
}


# present the results of the permutations test
# input: the MSE values for the results, and the corresponding permutations
# output: the p-value of the test
calculate.visualize.p.values <- function(controls_results, #the MSE of predicting
                                         # the IV using NCs
                                         permutations_results, # a series of the MSE 
                                         # results of predicting randomized using NCs
                                         title, #the title for the output file and graph
                                         OOB, # use Out Of Bag RMSE instead of Cross Validation error
                                         saveplot = F # whether to save plot and CSV of the
                                         ) {
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
    file_title <- stringr::str_trim(sprintf("%s",gsub("[.]|[:]|[\n]|[,]", " ",title)))
    file_full_path <- file.path("out", file_title)
    ggsave(sprintf("%s.png", file_full_path), p)
    results <- data.frame(permutations_results, 
                          controls_results= c(controls_results, rep(NA, length(permutations_results)-1)))
    write.csv(results, sprintf("%s.csv", file_full_path))
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
get_stata_coef <- function(G, # number of groups
                           N, # observations
                           reg_model, #lm object
                           hc_type= "HC2" #heteroscedasticity type
                           ){
  dfa <- (G/(G - 1)) * (N - 1)/reg_model$df.residual
  
  firm_c_vcov <- dfa * vcovHC(reg_model, type = hc_type, cluster = "group", adjust = T)
  return(coeftest(reg_model, vcov = firm_c_vcov))
}

# draw variable importance graph for the models: Y= f(NCs) and IV = f(NCs)
# given the controls using linear regression model
# input: data and names of Y, IV and controls
# output: graph
draw_importance_graph <- function(data,
                                  y, #name of the outcome
                                  instrument, #name of the IV
                                  controls, # names of the controls used to residualize 
                                  # both the outcome (y) and the IV
                                  all_controls, #indicate which control variables are
                                  # used only for coloring 
                                  weights = NULL, #the weights used in the 
                                  # variables_to_remove
                                  variables_to_remove, # the names of the variables not to used
                                  # as NCs
                                  title # the name of the output file and graph
                                  )  {
  
  ins_formula <- as.formula(paste(instrument, "~", paste(controls, collapse="+"))) 
  ins_lm <-  lm(as.formula(ins_formula), data = data)
  y_formula <- as.formula(paste(y, "~", paste(controls, collapse="+"))) 
  y_lm <- lm(as.formula(y_formula), data = data)
  
  weights_vec <- NULL
  if (!is.null(weights)){
    ins_lm <- rlang::eval_tidy(rlang::quo(
      lm(
        ins_formula,
        data = data,
        weights = !!as.name(weights))
    ))
    
    y_lm <- rlang::eval_tidy(rlang::quo(
      lm(
        y_formula,
        data = data,
        weights = !!as.name(weights))
    ))
    
    weights_vec <- pull(data, weights)
  }
    
  instrument_res <- ins_lm$residuals
  y_res <- y_lm$residuals
  
  NC <- get_NC_matrix(select(data,-one_of(variables_to_remove)),  controls,
                      weights_vec, T)
  
  model_ins <- randomForest::randomForest(NC, instrument_res)
  model_y <- randomForest::randomForest(NC, y_res)
  
  vi_scores_ins_raw <- importance(model_ins, type= 2, scale =T)
  vi_scores_ins <- vi_scores_ins_raw / max(vi_scores_ins_raw)
  vi_scores_ins <- data.frame(variable=rownames(vi_scores_ins), vi_scores_ins)
  vi_scores_y_raw <- importance(model_y, type= 2, scale =T)
  vi_scores_y <- vi_scores_y_raw / max(vi_scores_y_raw)
  vi_scores_y <- data.frame(variable=rownames(vi_scores_y), vi_scores_y)
  
  vi_scores <- vi_scores_ins %>%
    left_join (vi_scores_y,
               by = "variable")
  colnames(vi_scores)[2:3] <- c("ins_score", "y_score")
  
  file_title <- paste("Variable Importance", title)
  
  p <- ggplot(vi_scores, aes(ins_score, y_score, label = variable,
                             color = variable %in% all_controls))+
    geom_point() +
    ggtitle(file_title)
  
  print(plotly::ggplotly(p))
  
  file_full_path <- file.path("out", file_title)
  ggsave(sprintf("%s.png", file_full_path), p)
  write.csv(vi_scores, sprintf("%s.csv", file_full_path))
}


# draw variable importance graph for the models: Y= f(NCs) and IV = f(NCs)
# given the controls using fixed effects model
# input: data and names of Y, IV and controls
# output: graph
draw_felm_importance_graph <- function(data,
                                      y, #name of the outcome
                                      ins_lfe_formula, #the fixed effects formula
                                      # for the iv
                                      y_formula, #the fixed effects formula
                                      # for the outcome
                                      controls, # names of the controls used to residualize 
                                      # both the outcome (y) and the IV
                                      all_controls, #indicate which control variables are
                                      # used only for coloring 
                                      fixed, # the name of the fixed effects groups column
                                      variables_to_remove, # the names of the variables not to used
                                      # as NCs
                                      title # the name of the output file and graph
                                      )  {
  
  felm_of_controls <-  felm(formula= ins_lfe_formula, data = data)
  instrument_res <- felm_of_controls$residuals
  
  y_res_felm <- felm(formula= y_formula, data = data)
  y_res <- y_res_felm$residuals
  
  NC <- get_NC_matrix_felm(select(data,-one_of(variables_to_remove)),  controls, fixed)
  
  model_ins <- randomForest::randomForest(NC, instrument_res)
  model_y <- randomForest::randomForest(NC, y_res)
  
  vi_scores_ins_raw <- importance(model_ins, type= 2, scale =T)
  vi_scores_ins <- vi_scores_ins_raw / max(vi_scores_ins_raw)
  vi_scores_ins <- data.frame(variable=rownames(vi_scores_ins), vi_scores_ins)
  vi_scores_y_raw <- importance(model_y, type= 2, scale =T)
  vi_scores_y <- vi_scores_y_raw / max(vi_scores_y_raw)
  vi_scores_y <- data.frame(variable=rownames(vi_scores_y), vi_scores_y)
  
  vi_scores <- vi_scores_ins %>%
    left_join (vi_scores_y,
               by = "variable")
  colnames(vi_scores)[2:3] <- c("ins_score", "y_score")
  
  
  file_title <- paste("Variable Importance", title)
  
  p <- ggplot(vi_scores, aes(ins_score, y_score, label = variable,
                             color = variable %in% all_controls))+
    geom_point() +
    ggtitle(file_title)
  
  print(plotly::ggplotly(p))
  
  file_full_path <- file.path("out", file_title)
  ggsave(sprintf("%s.png", file_full_path), p)
  write.csv(vi_scores, sprintf("%s.csv", file_full_path))
}
