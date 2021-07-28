
run_school_nciv <- function(data, #cms school data
                            ivs, #the potential instrumental variables that
                            # need to be checked in NCIV permutations
                            iterations, # the number of iterations to run a full
                            # cycle of NCIV permutations for all IVs declared 
                            # just above (should be 1 in most cases)
                            permutations, # the number of permutations in each 
                            # NCIV test
                            ntree, # the number of trees in the RF prediction
                            # algorithm used for NCIV test
                            mtry_ratio, #the ratio of variables to check in 
                            # each iteration of RF algorithm
                            OOB, # use Out Of Bag RMSE instead of Cross Validation
                            # error
                            randomize_lottery #for sanity check, if True generate 
                            # Bernoulli RV instead of "lottery" variable
                            ) {
  models <-  c("2") #c("1,2")#
  prev_yrs <- c("all") #c("02", "2yr", "all")
  
  for (i in 1:iterations){
    # There are 4 different models for the calculation of the 2nd stage,
    # VA= Value Added (of a school)
    # (1) Intercept for school only, no student prior tests: test_score = a*VA 
    # (2) add student prior scores: test_score = a*VA + b*prior_scores
    # (3) + (4) use demographic and other variables that we don't have
    for (model in models){
      # how many years the VA calculated is calculated on
      for (prev_yr in prev_yrs)
      {
        print(sprintf("start working on: model= %s prev_yr= %s", model, prev_yr))
        
        # gen VA=as_mod`m'`e'_`s'_`t'
        data$VA <- unlist(data[sprintf("as_mod%smix_%s_test", model, prev_yr)])
        
        data$lott_VA <- rep(-1, nrow(data))
        # gen lott_VA=`c'_mod`m'`e'_`s'_`t' if lottery==0
        data[data$lottery == 0,][,"lott_VA"] <- 
          data[data$lottery == 0,][,sprintf("hm_mod%smix_%s_test", model, prev_yr)]
        # replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
        data[data$lottery == 1,][,"lott_VA"] <- 
          data[data$lottery == 1,][,sprintf("ch1_mod%smix_%s_test", model, prev_yr)]
        
        relevant_cols <- c("math_2002_imp", "read_2002_imp", "math_2002_imp_sq", "math_2002_imp_cub", 
                           "read_2002_imp_sq", "read_2002_imp_cub", "math_2002_miss", "read_2002_miss")
        
        
        all_controls <- c("math_1998_imp", "math_1998_miss",# "math_1998_imp_sq", "math_1998_imp_cub",       
                          "read_1998_imp", "read_1998_miss",# "read_1998_imp_sq", "read_1998_imp_cub",       
                          "math_1999_imp", "math_1999_miss",# "math_1999_imp_sq", "math_1999_imp_cub",       
                          "read_1999_imp", "read_1999_miss",# "read_1999_imp_sq", "read_1999_imp_cub",       
                          "math_2000_imp", "math_2000_miss",# "math_2000_imp_sq", "math_2000_imp_cub",       
                          "read_2000_imp", "read_2000_miss",# "read_2000_imp_sq", "read_2000_imp_cub",       
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
        
        results <- matrix(data= NA, nrow = length(ivs), ncol = 4) #per request: Cols: different controls (none, lottery fixed effects, home school)
        iv_ix <- 0
        for (iv in ivs){
          iv_ix <- iv_ix +1
          start.time <- Sys.time()
          print(sprintf("start processing for iv= %s permutations= %s", iv, permutations))
          title <- sprintf("iv=%s model= mod_%s_mix_%s", iv, model, prev_yr)
          
          #leave only students with all the varibales 
          curr_cms_reduced <- data[,unique(c(relevant_cols, all_controls, iv_end_outcome, choice, iv))]
          curr_cms_reduced <- curr_cms_reduced[complete.cases(curr_cms_reduced),]
          #problematic group, all chose the same school
          curr_cms_reduced <- curr_cms_reduced[which(curr_cms_reduced$lottery_FE != 14),]
          
          #used for sanity check
          if (randomize_lottery){
            prob <- mean(curr_cms_reduced$lottery)
            curr_cms_reduced_lot <- curr_cms_reduced
            curr_cms_reduced_lot$lottery <- rbernoulli(n= length(curr_cms_reduced_lot$lottery), p= prob)#rnorm(n=length(curr_cms_reduced_lot$lottery))#
            curr_cms_reduced <- curr_cms_reduced_lot
          }

          #(Section 4.e) ii.	Cols: different controls (none (at all), none (specified model), lottery fixed effects, home school)
          # For each IV we present 4 different controls models for the NCIV test
          # note that a valid IV is an IV that doesn;t reject the null, therefore
          # generally, adding more controls can make it happe
          # (1) Linear model with no controls
          results[iv_ix, 1] <- permutations.test.for.lm(data= curr_cms_reduced,
                                                        instrument_form = iv,
                                                        instrument= iv,
                                                        controls=relevant_cols,
                                                        weights = NULL,
                                                        variables_to_remove= c("testz2003", "lott_VA", "VA","lottery", iv),
                                                        title = title,
                                                        n_permutations= permutations,
                                                        conditioned = F,
                                                        OOB=OOB,
                                                        saveplot= T,
                                                        mtry_ratio= mtry_ratio,
                                                        ntree= ntree)
          
          # (2) Linear model with controls 
          results[iv_ix, 2] <- permutations.test.for.lm(data= curr_cms_reduced,
                                                        instrument_form = iv,
                                                        instrument= iv,
                                                        controls=relevant_cols,
                                                        weights = NULL,
                                                        variables_to_remove= c("testz2003", "lott_VA", "VA","lottery", iv),
                                                        title = title,
                                                        n_permutations= permutations,
                                                        conditioned = T,
                                                        OOB=OOB,
                                                        saveplot= T,
                                                        mtry_ratio= mtry_ratio,
                                                        ntree= ntree)
          
          # (3) Fixed Effects Linear model with controls, the model used in the article
          results[iv_ix, 3] <- permutations.test.for.felm(data= curr_cms_reduced,
                                                          ins_lfe_formula = get_formula_for_iv_by_NC(iv, "lottery_FE", relevant_cols),
                                                          instrument=iv,
                                                          controls= relevant_cols, #c("1"),
                                                          fixed= "lottery_FE",
                                                          variables_to_remove= c("testz2003", "lott_VA", "VA","lottery", iv),
                                                          title= title,
                                                          n_permutations= permutations,
                                                          conditioned = T,
                                                          OOB =OOB,
                                                          saveplot= T,
                                                          mtry_ratio= mtry_ratio,
                                                          ntree= ntree)
          
          # (4) Fixed Effects Linear model with controls + VA of the home school of the student
          results[iv_ix, 4] <- permutations.test.for.felm(data= curr_cms_reduced,
                                                          ins_lfe_formula = get_formula_for_iv_by_NC(iv, "lottery_FE", c(relevant_cols,sprintf("hm_mod%smix_%s_test", model, prev_yr))),
                                                          instrument=iv,
                                                          controls= c(relevant_cols,sprintf("hm_mod%smix_%s_test", model, prev_yr)), #c("1"),
                                                          fixed= "lottery_FE",
                                                          variables_to_remove= c("testz2003", "lott_VA", "VA","lottery", iv),
                                                          title= paste(title, "with hm"),
                                                          n_permutations= permutations,
                                                          conditioned = T,
                                                          OOB =OOB,
                                                          saveplot= T,
                                                          mtry_ratio= mtry_ratio,
                                                          ntree= ntree)
          
          
          print(sprintf("Done processing for iv= %s permutations= %s. Runtime is:", iv, permutations))
          print(Sys.time() - start.time)
        }
        print(sprintf("Done working on: model= %s prev_yr= %s", model, prev_yr))
        
      }
    }
    
    colnames(results) <- c("lm: no controls", "lm: with controls", 
                           "felm with controls (as in the article",
                           "felm with controls + home school VAM")
    rownames(results) <- ivs
    print(results)
    
    file_name <- sprintf("out//school_lottery_iv_%s_permutaions_%s_ntree_%s_mtry_ratio_%s_randomize_lottery_%s_iterations_%s.csv", 
            iv, permutations, ntree, mtry_ratio, randomize_lottery, i)
    write.csv(results, file_name)
  }
}


