library(haven)

get_formula_for_table_1_lottery <- function(outcome, # the name of the outcome variable
                                            col_factors, # the name of the outcome variable
                                            fixed, # the name of the fixed effects groups
                                            endogen, # the name of the endogen treatment variable
                                            iv # the name of the instrumental variable
                                            ){
  return(as.formula(paste(outcome, "~ 1+ ", paste(col_factors, collapse=" + "),
                          "|", fixed,
                          "| (", endogen, "~", paste(c(iv, col_factors), collapse=" + "), ")",
                          "|", fixed))
  )
}

get_formula_for_iv_by_controls <- function(iv, # the name of the instrumental variable
                                           fixed,  # the name of the fixed effects groups
                                           controls # the name of the outcome variable
                                           ){
  return(as.formula(paste(iv, "~ ", paste(controls, collapse=" + "),
                          "|", fixed))
  ) 
}


cms_path <- file.path("data", "cms_added_columns.dta")
cms <- read_dta(cms_path) 
curr_cms <- dplyr::filter(cms, !is.na(cms$lottery_FE))
