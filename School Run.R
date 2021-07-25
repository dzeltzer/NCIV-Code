
run_school_nciv <- function(curr_cms, ivs, iterations, permutations, ntree,
                            mtry_ratio, OOB, randomize_lottery) {
  models <-  c("2") #c("1,2")#
  prev_yrs <- c("all") #c("02", "2yr", "all")
  
  for (i in 1:iterations){
    for (model in models){
      for (prev_yr in prev_yrs)
      {
        print(sprintf("start working on: model= %s prev_yr= %s", model, prev_yr))
        
        # gen VA=as_mod`m'`e'_`s'_`t'
        curr_cms$VA <- unlist(curr_cms[sprintf("as_mod%smix_%s_test", model, prev_yr)])
        
        curr_cms$lott_VA <- rep(-1, nrow(curr_cms))
        # gen lott_VA=`c'_mod`m'`e'_`s'_`t' if lottery==0
        curr_cms[curr_cms$lottery == 0,][,"lott_VA"] <- 
          curr_cms[curr_cms$lottery == 0,][,sprintf("hm_mod%smix_%s_test", model, prev_yr)]
        # replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
        curr_cms[curr_cms$lottery == 1,][,"lott_VA"] <- 
          curr_cms[curr_cms$lottery == 1,][,sprintf("ch1_mod%smix_%s_test", model, prev_yr)]
        
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
        #the VAM of the school the student chose. For exmaple:ch2_mod1mix_02_test
        choice <- c(sprintf("ch1_mod%smix_%s_test", model, prev_yr),
                    sprintf("ch2_mod%smix_%s_test", model, prev_yr),
                    sprintf("ch3_mod%smix_%s_test", model, prev_yr),
                    sprintf("hm_mod%smix_%s_test", model, prev_yr))
        
        #create new_lott_VA: 0 when the student lost the lottery and VAM of the 
        # the 1st choice school had she won the lottery
        curr_cms$new_lott_VA <- rep(-1, nrow(curr_cms))
        curr_cms[curr_cms$lottery == 0,][,"new_lott_VA"] <- 0
        curr_cms[curr_cms$lottery == 1,][,"new_lott_VA"] <- 
          curr_cms[curr_cms$lottery == 1,][,sprintf("ch1_mod%smix_%s_test", model, prev_yr)]
        
        #mod1 02
        # (Section 4.d) his instrument, our instrument (without home school) 
        
        results <- matrix(data= NA, nrow = length(ivs), ncol = 4) #per request: Cols: different controls (none, lottery fixed effects, home school)
        iv_ix <- 0
        for (iv in ivs){
          iv_ix <- iv_ix +1
          start.time <- Sys.time()
          print(sprintf("start processing for iv= %s permutations= %s", iv, permutations))
          title <- sprintf("iv=%s model= mod_%s_mix_%s", iv, model, prev_yr)
          
          curr_cms_reduced <- curr_cms[,unique(c(relevant_cols, all_controls, iv_end_outcome, choice, iv))]
          curr_cms_reduced <- curr_cms_reduced[complete.cases(curr_cms_reduced),]
          #problematic group, all chose the same school
          curr_cms_reduced <- curr_cms_reduced[which(curr_cms_reduced$lottery_FE != 14),]
          
          #sanity check
          if (randomize_lottery){
            prob <- mean(curr_cms_reduced$lottery)
            curr_cms_reduced_lot <- curr_cms_reduced
            curr_cms_reduced_lot$lottery <- rbernoulli(n= length(curr_cms_reduced_lot$lottery), p= prob)#rnorm(n=length(curr_cms_reduced_lot$lottery))#
            curr_cms_reduced <- curr_cms_reduced_lot
          }

          #(Section 4.e) ii.	Cols: different controls (none (at all), none (specified model), lottery fixed effects, home school)
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
