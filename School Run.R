
run_school_nciv <- function(data, #cms school data
                            ivs, #the potential instrumental variables that
                            # need to be checked in NCIV permutations
                            controls_specifications = 1:4, #the controls specifications
                            # for the NCIV test
                            specified_nco,
                            randomize_lottery = F) {
  
  assertthat::assert_that(!(randomize_lottery &
                             !(length(ivs) == 1 & ivs[1] == "lottery")))
  
  
  models <-  c("2") #c("1,2")#
  prev_yrs <- c("all") #c("02", "2yr", "all")
  

    # There are 4 different models for the calculation of the 2nd stage,
    # VA= Value Added (of a school)
    # (1) Intercept for school only, no student prior tests: test_score = a*VA 
    # (2) add student prior scores: test_score = a*VA + b*prior_scores
    # (3) + (4) use demographic and other variables that we don't have
    for (model in models){
      # how many years the VA calculated is calculated on
      model.start.time <- Sys.time()
      for (prev_yr in prev_yrs)
      {
        print(sprintf("start NCIV tests on: model= %s prev_yr= %s", model, prev_yr))
        
        # the controls as specified in the article
        relevant_cols <- c("math_2002_imp", "read_2002_imp", "math_2002_imp_sq", "math_2002_imp_cub", 
                           "read_2002_imp_sq", "read_2002_imp_cub", "math_2002_miss", "read_2002_miss")
        
        
        data_reduced <- prepare_school_data(data, model, relevant_cols, 
                                            prev_yr, ivs, randomize_lottery)
        
        results <- matrix(data= NA, nrow = length(ivs), ncol = 4) #per request: Cols: different controls (none, lottery fixed effects, home school)
        iv_ix <- 0
        for (iv in ivs){
          iv_ix <- iv_ix +1
          iv.start.time <- Sys.time()
          print(sprintf("start processing for iv= %s", iv))
          title <- sprintf("iv=%s model= mod_%s_mix_%s", iv, model, prev_yr)
          

          #(Section 4.e) ii.	Cols: different controls (none (at all), none (specified model), lottery fixed effects, home school)
          # For each IV we present 4 different controls specifications for the NCIV test
          # note that a valid IV is an IV that doesn;t reject the null, therefore
          # generally, adding more controls can make it happen
          
          # (1) Linear model with no controls
          if (1 %in% controls_specifications){
            bonf.f.test.for.lm(data= data_reduced,
                                     instrument_form= iv,
                                     instrument= iv,
                                     controls= relevant_cols,
                                     weights=NULL,
                                     variables_to_remove= c("testz2003", "VA",ivs),
                                     title= title,
                                     conditioned = F,
                                      specified_nco = specified_nco
                                                )
          }
          

          # (2) Linear model with controls 
          if (2 %in% controls_specifications){
            bonf.f.test.for.lm(data= data_reduced,
                                      instrument_form= iv,
                                      instrument= iv,
                                      controls= relevant_cols,
                                      weights=NULL,
                                      variables_to_remove= c("testz2003", "VA",ivs),
                                      title= title,
                                      conditioned = T,
                                      specified_nco = specified_nco
            )
          }

          # (3) Fixed Effects Linear model with controls, the model used in the article
          if (3 %in% controls_specifications){
            bonf.f.test.for.felm(data= data_reduced,
                                      ins_lfe_formula = get_formula_for_iv_by_controls(iv, "lottery_FE", relevant_cols),
                                      instrument=iv,
                                      controls= relevant_cols, #c("1"),
                                      fixed= "lottery_FE",
                                      variables_to_remove= c("testz2003", "VA", ivs),
                                      title= title,
                                      conditioned = T,
                                      specified_nco = specified_nco)
          }

          
          # (4) Fixed Effects Linear model with controls + VA of the home school of the student
          if (4 %in% controls_specifications){
            bonf.f.test.for.felm(data= data_reduced,
                                      ins_lfe_formula = get_formula_for_iv_by_controls(iv, "lottery_FE", c(relevant_cols,sprintf("hm_mod%smix_%s_test", model, prev_yr))),
                                      instrument=iv,
                                      controls= c(relevant_cols,sprintf("hm_mod%smix_%s_test", model, prev_yr)), #c("1"),
                                      fixed= "lottery_FE",
                                      variables_to_remove= c("testz2003","VA", ivs),
                                      title= paste(title, "with hm"),
                                      conditioned = T,
                                     specified_nco = specified_nco)
          }

          print(sprintf("Done processing for iv= %s  Runtime is:", iv))
          print(Sys.time() - iv.start.time)
        }
        print(sprintf("Done working on: model= %s prev_yr= %s", model, prev_yr))
        print(Sys.time() - model.start.time)
        
      }
    
    
    colnames(results) <- c("lm: no controls", "lm: with controls", 
                           "felm with controls (as in the article",
                           "felm with controls + home school VAM")
    rownames(results) <- ivs
    print(results)
    return(results)
  }
}

# run RF variable importance for school data for desired ivs 
run_school_variable_importance <- function(data, #cms school data
                            ivs # the ivs that will be checked
                            )
 {
  models <-  c("2") #c("1,2")#
  prev_yrs <- c("all") #c("02", "2yr", "all")
  
  # There are 4 different models for the calculation of the 2nd stage,
  # VA= Value Added (of a school)
  # (1) Intercept for school only, no student prior tests: test_score = a*VA 
  # (2) add student prior scores: test_score = a*VA + b*prior_scores
  # (3) + (4) use demographic and other variables that we don't have
  for (model in models){
    # how many years the VA calculated is calculated on
    model.start.time <- Sys.time()
    for (prev_yr in prev_yrs)
    {
      print(sprintf("start variable importance on: model= %s prev_yr= %s", model, prev_yr))
      
      # the controls as specified in the article
      relevant_cols <- c("math_2002_imp", "read_2002_imp", "math_2002_imp_sq", "math_2002_imp_cub", 
                         "read_2002_imp_sq", "read_2002_imp_cub", "math_2002_miss", "read_2002_miss")
      
      data_reduced <- prepare_school_data(data, model, relevant_cols, 
                                          prev_yr, ivs, F)
      
      for (iv in ivs){
        iv.start.time <- Sys.time()
        print(sprintf("start processing for iv= %s", iv))
        title <- sprintf("iv=%s model= mod_%s_mix_%s", iv, model, prev_yr)
        
        # (3) Fixed Effects Linear model with controls, the model used in the article
        draw_felm_importance_graph(data_reduced,
                                   y= "testz2003",
                                   ins_lfe_formula = get_formula_for_iv_by_controls(iv, "lottery_FE", relevant_cols),
                                   y_formula= get_formula_for_iv_by_controls("testz2003", "lottery_FE", relevant_cols),
                                  controls= relevant_cols,
                                  all_controls= relevant_cols,
                                  fixed = "lottery_FE",
                                  variables_to_remove = c("testz2003","VA", ivs),
                                  paste("Deming (2014) - School Effectiveness - Variable Importance for", iv))
          
          print(sprintf("Done processing for iv= %s. Runtime is:", iv))
          print(Sys.time() - iv.start.time)
        }
      }
    print(sprintf("Done working on: model= %s prev_yr= %s", model, prev_yr))
    print(Sys.time() - model.start.time)
  }
}


prepare_school_data <- function(data, model, relevant_cols, prev_yr, ivs, randomize_lottery) {
  # gen VA=as_mod`m'`e'_`s'_`t'
  data$VA <- unlist(data[sprintf("as_mod%smix_%s_test", model, prev_yr)])
  
  data$lott_VA <- rep(-1, nrow(data))
  # gen lott_VA=`c'_mod`m'`e'_`s'_`t' if lottery==0
  data[data$lottery == 0,][,"lott_VA"] <- 
    data[data$lottery == 0,][,sprintf("hm_mod%smix_%s_test", model, prev_yr)]
  # replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
  data[data$lottery == 1,][,"lott_VA"] <- 
    data[data$lottery == 1,][,sprintf("ch1_mod%smix_%s_test", model, prev_yr)]
  
  
  all_controls <- c(#"math_1998_imp", "math_1998_miss",# "math_1998_imp_sq", "math_1998_imp_cub",       
                    #"read_1998_imp", "read_1998_miss",# "read_1998_imp_sq", "read_1998_imp_cub",       
                    #"math_1999_imp", "math_1999_miss",# "math_1999_imp_sq", "math_1999_imp_cub",       
                    #"read_1999_imp", "read_1999_miss",# "read_1999_imp_sq", "read_1999_imp_cub",       
                    #"math_2000_imp", "math_2000_miss",# "math_2000_imp_sq", "math_2000_imp_cub",       
                    #"read_2000_imp", "read_2000_miss",# "read_2000_imp_sq", "read_2000_imp_cub",       
                    "math_2001_imp", "math_2001_miss",# "math_2001_imp_sq", "math_2001_imp_cub",       
                    "read_2001_imp", "read_2001_miss"# , "read_2001_imp_sq", "read_2001_imp_cub"
  )
  
  iv_end_outcome = c("testz2003", "lott_VA", "VA","lottery_FE","lottery")
  
  #The VA of the school the student chose. For exmaple:ch2_mod1mix_02_test
  choice <- c(sprintf("ch1_mod%smix_%s_test", model, prev_yr),
              sprintf("ch2_mod%smix_%s_test", model, prev_yr),
              sprintf("ch3_mod%smix_%s_test", model, prev_yr),
              sprintf("hm_mod%smix_%s_test", model, prev_yr))
  
  # Create new_lott_VA: 0 when the student lost the lottery and VAM of the 
  # the 1st choice school had she won the lottery
  data$new_lott_VA <- rep(-1, nrow(data))
  data[data$lottery == 0,][,"new_lott_VA"] <- 0
  data[data$lottery == 1,][,"new_lott_VA"] <- 
    data[data$lottery == 1,][,sprintf("ch1_mod%smix_%s_test", model, prev_yr)]
  
  
  #leave only students with all the varibales 
  data_reduced <- data[,unique(c(relevant_cols, all_controls, iv_end_outcome, choice, ivs))]
  data_reduced <- data_reduced[complete.cases(data_reduced),]
  #problematic group, all chose the same school
  data_reduced <- data_reduced[which(data_reduced$lottery_FE != 14),]
  
  #used for sanity check
  if (randomize_lottery){
    prob <- mean(data_reduced$lottery)
    data_reduced_lot <- data_reduced
    data_reduced_lot$lottery <- rbernoulli(n= length(data_reduced_lot$lottery), p= prob)#rnorm(n=length(curr_cms_reduced_lot$lottery))#
    data_reduced <- data_reduced_lot
  }
  
  return(data_reduced)
}



