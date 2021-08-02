library(haven)

get_formula_for_table_1_lottery <- function(outcome, col_factors, fixed, endogen, iv){
  return(as.formula(paste(outcome, "~ 1+ ", paste(col_factors, collapse=" + "),
                          "|", fixed,
                          "| (", endogen, "~", paste(c(iv, col_factors), collapse=" + "), ")",
                          "|", fixed))
  )
}

get_formula_for_iv_by_controls <- function(iv, fixed, controls){
  return(as.formula(paste(iv, "~ ", paste(controls, collapse=" + "),
                          "|", fixed))
  ) 
}


cms_path <- file.path("data", "cms_added_columns.dta")
cms <- read_dta(cms_path) 
curr_cms <- dplyr::filter(cms, !is.na(cms$lottery_FE))