library(haven)
library(dplyr)

workfile_china_raw_path <- "data\\workfile_china.dta"
workfile_china_raw <- read_dta(workfile_china_raw_path)

workfile_china <- workfile_china_raw  %>%
  left_join (
    workfile_china_raw %>%
      group_by(czone) %>%
      summarise(instrument1990= sum(d_tradeotch_pw_lag*(yr==1990 + 0)),
                exposure1990= sum(d_tradeusch_pw*(yr==1990 + 0)),
                outcome1990= sum(d_sh_empl_mfg*(yr==1990 + 0)),
                instrument2000= sum(d_tradeotch_pw_lag*(yr==2000 + 0)),
                exposure2000= sum(d_tradeusch_pw*(yr==2000 + 0)),
                outcome2000= sum(d_sh_empl_mfg*(yr==2000 + 0)),
      )
    , by= "czone")

#t2 should be removed because we check NCIV only on 1990 data
china_1990 <- workfile_china %>%
  filter(yr == 1990) %>%
  select(-c("l_tradeusch_pw", "l_tradeotch_pw",
            "exposure1990", "instrument1990", "exposure2000",  
            "czone", "yr", "t2", "statefip", "city"), 
         -starts_with("d"), -starts_with("relchg"), -starts_with("lnchg")) #Table 3 (pg. 2137)

assertthat::are_equal(ncol(workfile_china) , 197)
assertthat::are_equal(nrow(workfile_china) , 1444)

assertthat::are_equal(ncol(china_1990) , 67)
assertthat::are_equal(nrow(china_1990) , 722)


G <- length(unique(workfile_china$statefip))
N <- length(workfile_china$statefip)

# eststo: ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw=d_tradeotch_pw_lag) l_shind_manuf_cbp t2 [aw=timepwt48], cluster(statefip) first
col_2_controls <- c("t2", "l_shind_manuf_cbp")
census_controls <- c(#c("statefip",#	"city",
  "reg_midatl",	"reg_encen",	"reg_wncen",
  "reg_satl",	"reg_escen",	"reg_wscen",	"reg_mount",
  "reg_pacif")
col_3_controls <- c(col_2_controls, census_controls)
col_4_controls <- c(col_3_controls, "l_sh_popedu_c", "l_sh_popfborn", "l_sh_empl_f")
col_5_controls <-  c(col_3_controls, "l_sh_routine33", "l_task_outsource")
col_6_controls <-  c(col_4_controls, "l_sh_routine33", "l_task_outsource")
all_controls <- col_6_controls[2:length(col_6_controls)]

# Based on pg. 2135 second paragraph and czone_analysis_preperiod_final.do

#eststo: ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw_future=d_tradeotch_pw_lag_future) [aw=timepwt48] if yr==1970, cluster(statefip)
#eststo: ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw_future=d_tradeotch_pw_lag_future) [aw=timepwt48] if yr==1980, cluster(statefip)
#eststo: ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw_future=d_tradeotch_pw_lag_future) t1980 [aw=timepwt48] if yr>=1970 & yr<1990, cluster(statefip)

china_1990_only_org_nc <- china_1990 %>% select(c(col_6_controls[2:length(col_6_controls)], "timepwt48", "instrument2000", "outcome1990"))

assertthat::are_equal(ncol(china_1990_only_org_nc) , 17)
assertthat::are_equal(nrow(china_1990_only_org_nc) , 722)

get_formula_for_table_3 <- function(col_controls){
  return(as.formula(paste("d_sh_empl_mfg ~ d_tradeusch_pw +", paste(col_controls, collapse="+"),
                          paste("| . - d_tradeusch_pw + d_tradeotch_pw_lag"))))
}



location_controls <- c("factor(statefip)",	"city")
#mfg= manufacring
man_controls <- c("l_shind_manuf_cbp",
                  "l_sh_empl",	"l_sh_empl_mfg",	"l_sh_empl_mfg_m",	"l_sh_empl_mfg_f",
                  "l_sh_empl_mfg_edu_nc",	"l_sh_empl_mfg_edu_c",	"l_sh_empl_mfg_age1634",
                  "l_sh_empl_mfg_age3549", 	"l_sh_empl_mfg_age5064",	"l_sh_empl_nmfg",
                  "l_sh_empl_nmfg_m",	"l_sh_empl_nmfg_f",	"l_sh_empl_nmfg_edu_nc",
                  "l_sh_empl_nmfg_edu_c",	"l_sh_empl_nmfg_age1634",	"l_sh_empl_nmfg_age3549",
                  "l_sh_empl_nmfg_age5064")
#unemployment
uemp_controls <- c("l_sh_unempl",	"l_sh_unempl_m",	"l_sh_unempl_f",	"l_sh_unempl_edu_nc",	
                   "l_sh_unempl_edu_c",	"l_sh_unempl_age1634",	"l_sh_unempl_age3549",	"l_sh_unempl_age5064")
#Not-in-Labor-Force
nilf_controls <- c("l_sh_nilf",	"l_sh_nilf_m",	"l_sh_nilf_f",	"l_sh_nilf_edu_nc",
                   "l_sh_nilf_edu_c",	"l_sh_nilf_age1634",	"l_sh_nilf_age3549",	"l_sh_nilf_age5064")

trans_controls <- c("l_trans_totindiv_pc",	"l_trans_totmed_pc",	"l_trans_fedinc_pc",
                    "l_trans_othinc_pc",	"l_trans_unemp_pc",	"l_trans_taaimp_pc",
                    "l_trans_totedu_pc",	"l_trans_ssaret_pc",	"l_trans_ssadis_pc")

hh_controls <- c("l_avg_hhincsum_pc_pw",	"l_avg_hhincwage_pc_pw") #house hold

wage_controls <- c("l_sh_ssadiswkrs",	"l_avg_lnwkwage_mfg",	"l_avg_lnwkwage_nmfg")

demographics_controls <- c("l_popcount",	"l_no_workers_totcbp",	"l_shind_manuf_cbp",
                           "l_sh_popedu_c",	"l_sh_popfborn",	"l_sh_empl_f")

offshoring_controls <- c("l_sh_routine33",	"l_task_outsource")
