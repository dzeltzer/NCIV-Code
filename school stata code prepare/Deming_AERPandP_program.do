set more off
cd "C:\Users\itaiwalk\Desktop\Deming_AERPandP_datafolder\"
use cms_VAManalysis, clear

gen sample=future_grd>=4 & future_grd<=8 & miss_02==0
gen onmargin_sample=onmargin==1 & sample==1

*NOTES ON THE LOTTERY VARS*
*the variable lottery_FE is the unit of randomization*
*we condition on miss_02==0 in all specifications, because many students who previously in private schools applied, lost, and disappeared from CMS*
*We have no data for these kids*
*"enrolled" is enrolled in 1st choice*
*future_grd is the grade that a student expected be enrolled in the Fall of 2002 - lottery happened Spring 2002*
*onmargin is an indicator for students in priority groups where the probability of admission was greater than 0 and less than 1 - this is the analysis sample*

*Four different model variations - 1) different effect for every grade-year combo; 2) pool years, separate grades; 3) separate years, pool grades*
*Two different model approaches - 1) fixed effects; 2) average residual*
*three different outcomes - math and reading, average of the two*

*Create covariates*

*NOTE - I leave the code in for Model 3, which includes race and free/reduced price lunch eligibility*
*However, to protect student confidentiality, I cannot include demographic indicators in these data*
*I have commented out Model 3, but you can see the code and easily apply it to other data*
*You can also estimate Model 4, with twice-lagged scores, but without demographic data*

forvalues g=4(1)8 {
	gen year_`g'=.
	local f=`g'-1
	forvalues y=1998(1)2004 {
		local x=`y'-1
		replace year_`g'=`y' if grade`y'==`g' & grade`x'==`f'
	}
}
forvalues g=4(1)8 {
	egen schXcohortFE_`g'=group(school_`g' year_`g')
}
forvalues g=3(1)8 {
	egen test_`g'=rowmean(math_`g' read_`g')
}
forvalues y=3(1)7 {
	foreach x in math read {
		gen `x'_`y'_imp=`x'_`y'
		gen `x'_`y'_miss=`x'_`y'==.
		replace `x'_`y'_imp=0 if `x'_`y'_miss==1
		gen `x'_`y'_imp_sq=(`x'_`y'_imp)^2
		gen `x'_`y'_imp_cub=(`x'_`y'_imp)^3
	}
}
forvalues y=1998(1)2004 {
	foreach x in math read {
		gen `x'_`y'_imp=`x'z`y'
		gen `x'_`y'_miss=`x'z`y'==.
		replace `x'_`y'_imp=0 if `x'_`y'_miss==1
		gen `x'_`y'_imp_sq=(`x'_`y'_imp)^2
		gen `x'_`y'_imp_cub=(`x'_`y'_imp)^3
	}
}
forvalues y=1998(1)2004 {
	egen testz`y'=rowmean(mathz`y' readz`y')
}
*Need school by cohort FE*
forvalues y=1999(1)2004 {
	egen schXcohortFE_`y'=group(school`y' grade`y')
}

*Create peer prior test scores*
forvalues g=4(1)8 {
	local lag=`g'-1
	bysort school_`g': egen temp=mean(test_`lag') if onmargin_sample!=1
	bysort school_`g': egen avg_test`lag'=max(temp)
	drop temp
}

*Now estimate value-added models on pre-lottery data*
set more off
foreach t in /*math read*/ test {
	forvalues g=4(1)8 {
		local lag=`g'-1
	
		*Model 1 - levels only*
// 		xi: reg `t'_`g' if onmargin_sample!=1 & year_`g'==2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod1ar_02_`t'_`g'=mean(resid)
// 		drop resid
		
		xtmixed `t'_`g' || school_`g': if onmargin_sample!=1 & year_`g'==2002
		predict mod1mix_02_`t'_`g' if e(sample), reffects level(school_`g') 
		
// 		xtset school_`g' 
// 		xtreg `t'_`g' if onmargin_sample!=1 & year_`g'==2002, fe
// 		predict mod1FE_02_`t'_`g' if e(sample)==1, u
		
// 		xi: reg `t'_`g' if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod1ar_2yr_`t'_`g'=mean(resid)
// 		drop resid
		
		xi: xtmixed `t'_`g' i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
		predict mod1mix_2yr_`t'_`g' if e(sample), reffects level(school_`g') 
		
// 		xtset school_`g'
// 		xi: xtreg `t'_`g' i.year_`g' if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002, fe 
// 		predict mod1FE_2yr_`t'_`g' if e(sample)==1, u
	
// 		xi: reg `t'_`g' if onmargin_sample!=1 & year_`g'<=2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod1ar_all_`t'_`g'=mean(resid)
// 		drop resid
		
		xi: xtmixed `t'_`g' i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'<=2002
		predict mod1mix_all_`t'_`g' if e(sample), reffects level(school_`g')
	
// 		xtset school_`g'		
// 		xi: xtreg `t'_`g' i.year_`g' if onmargin_sample!=1 & year_`g'<=2002, fe 
// 		predict mod1FE_all_`t'_`g' if e(sample)==1, u
//	
// 		forvalues y=1998(1)2003 {
// 			xi: reg `t'_`g' i.year_`g' if onmargin_sample!=1 & year_`g'==`y', robust 
// 			predict resid if e(sample)==1, resid
// 			bysort school`y': egen mod1byG_`t'_`y'_`g'=mean(resid)
// 			drop resid
// 		}
		
		*Model 2 - gains (lagged scores only)
		
// 		xi: reg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'==2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod2ar_02_`t'_`g'=mean(resid)
// 		drop resid
		
		xi: xtmixed `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss || school_`g': if onmargin_sample!=1 & year_`g'==2002
		predict mod2mix_02_`t'_`g' if e(sample), reffects level(school_`g') 
//		
// 		xtset school_`g'
// 		xi: xtreg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'==2002, fe 
// 		predict mod2FE_02_`t'_`g' if e(sample)==1, u
//			
// 		xi: reg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod2ar_2yr_`t'_`g'=mean(resid)
// 		drop resid
		
		xi: xtmixed `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
		predict mod2mix_2yr_`t'_`g' if e(sample), reffects level(school_`g')

// 		xtset school_`g'
// 		xi: xtreg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002, fe 
// 		predict mod2FE_2yr_`t'_`g' if e(sample)==1, u
//	
// 		xi: reg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'<=2002
// 		predict resid if e(sample)==1, resid
// 		bysort school_`g': egen mod2ar_all_`t'_`g'=mean(resid)
// 		drop resid	
		
		xi: xtmixed `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'<=2002
		predict mod2mix_all_`t'_`g' if e(sample), reffects level(school_`g')

//		xtset school_`g'
// 		xi: xtreg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'<=2002, fe 
// 		predict mod2FE_all_`t'_`g' if e(sample)==1, u

// 		forvalues y=1998(1)2003 {
// 			xi: reg `t'_`g' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'==`y', robust 
// 			predict resid if e(sample)==1, resid
// 			bysort school`y': egen mod2byG_`t'_`y'_`g'=mean(resid)
// 			drop resid
// 		}
/*
		*Model 3 - demogs + prior scores + peer effects (standard VAM)
		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'==2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_02_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss || school_`g': if onmargin_sample!=1 & year_`g'==2002
		predict mod3mix_02_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'==2002, fe 
		predict mod3FE_02_`t'_`g' if e(sample)==1, u

		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_2yr_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
		predict mod3mix_2yr_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002, fe 
		predict mod3FE_2yr_`t'_`g' if e(sample)==1, u

		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss if onmargin_sample!=1 & year_`g'<=2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_all_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'<=2002
		predict mod3mix_all_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'<=2002, fe 
		predict mod3FE_all_`t'_`g' if e(sample)==1, u
			
		forvalues y=1998(1)2003 {
			xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss i.year_`g' if onmargin_sample!=1 & year_`g'==`y', robust 
			predict resid if e(sample)==1, resid
			bysort school`y': egen mod3byG_`t'_`y'_`g'=mean(resid)
			drop resid
		}*/
	}
}
/*
*Model 4 - add twice-lagged scores*

set more off
foreach t in math read test {
	forvalues g=5(1)8 {
		local lag=`g'-1
		local lag2=`g'-2
		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss if onmargin_sample!=1 & year_`g'==2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_02_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss || school_`g': if onmargin_sample!=1 & year_`g'==2002 
		predict mod4mix_02_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss if onmargin_sample!=1 & year_`g'==2002 , fe 
		predict mod4FE_02_`t'_`g' if e(sample)==1, u
	
		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_2yr_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002 
		predict mod4mix_2yr_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss i.year_`g' if onmargin_sample!=1 & year_`g'>=2001 & year_`g'<=2002 , fe 
		predict mod4FE_2yr_`t'_`g' if e(sample)==1, u

		xi: reg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss if onmargin_sample!=1 & year_`g'<=2002
		predict resid if e(sample)==1, resid
		bysort school_`g': egen mod3ar_all_`t'_`g'=mean(resid)
		drop resid
		xi: xtmixed `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss i.year_`g' || school_`g': || schXcohortFE_`g': if onmargin_sample!=1 & year_`g'<=2002 
		predict mod4mix_all_`t'_`g' if e(sample), reffects level(school_`g')
		xtset school_`g'
		xi: xtreg `t'_`g' male black hisp lunch avg_test`lag' math_`lag'_imp read_`lag'_imp math_`lag'_imp_sq read_`lag'_imp_sq math_`lag'_imp_cub read_`lag'_imp_cub math_`lag'_miss read_`lag'_miss math_`lag2'_imp read_`lag2'_imp math_`lag2'_imp_sq read_`lag2'_imp_sq math_`lag2'_imp_cub read_`lag2'_imp_cub math_`lag2'_miss read_`lag2'_miss i.year_`g' if onmargin_sample!=1 & year_`g'<=2002, fe 
		predict mod4FE_all_`t'_`g' if e(sample)==1, u

	}
}
*/
*Create VAM datasets*
forvalues g=4(1)8 {
	preserve
	drop if school_`g'==.	
	collapse (mean) mod*_`g', by(school_`g')
	rename school_`g' school
	save va_`g', replace	
	restore
}
*Now stack VAM datasets*
preserve
use va_4, clear
forvalues g=5(1)8 {
	merge 1:1 school using va_`g', keep(master using match) nogen
}
drop if school==0
drop if school>5999 & school<8000
drop if school>5999 & school!=8482
foreach var of varlist *_4 *_5 {
	replace `var'=. if school>5000
}
foreach var of varlist *_6 *_7 *_8 {
	replace `var'=. if school>4000 & school<5000
}
save va_schl, replace
restore

*restrict to analysis sample*
keep if sample==1
compress
drop mod*

*Merge on school VAMs*
gen school_as=schl_d20
forvalues c=1(1)3 {
	gen school_ch`c'=choice`c'_schl
}
gen school_hm=home0203

*My original code loops over 4 models, but if you run this you can only do 2*
*I comment out the original code and leave some that works on this dataset*

foreach schl in as ch1 ch2 ch3 hm {
	gen school=school_`schl'
	merge m:1 school using va_schl, keep(master match) nogen
	drop school
	/*forvalues m=1(1)3 {*/ forvalues m=1(1)2 {
		foreach e in /*FE ar*/ mix {
			foreach s in 02 2yr all {
				foreach t in /*math read*/ test {
					gen `schl'_mod`m'`e'_`s'_`t'=.
					*Add on leads, for 2004 test scores*
					gen `schl'lead_mod`m'`e'_`s'_`t'=.
					forvalues g=4(1)7 {
						local lead=`g'+1
						replace `schl'_mod`m'`e'_`s'_`t'=mod`m'`e'_`s'_`t'_`g' if future_grd==`g'
						replace `schl'lead_mod`m'`e'_`s'_`t'=mod`m'`e'_`s'_`t'_`lead' if future_grd==`g'
					}
					replace `schl'_mod`m'`e'_`s'_`t'=mod`m'`e'_`s'_`t'_8 if future_grd==8
				}
			}
		}
	}
/*	forvalues m=4(1)4 {
		foreach e in mix FE ar {
			foreach s in 02 2yr all {
				foreach t in /*math read*/ test {
					gen `schl'_mod`m'`e'_`s'_`t'=.
					forvalues g=5(1)8 {
						replace `schl'_mod`m'`e'_`s'_`t'=mod`m'`e'_`s'_`t'_`g' if future_grd==`g'
					}
				}
			}
		}
	}*/
// 	/*forvalues m=1(1)3 {*/ forvalues m=1(1)2 {
// 		foreach t in /*math read*/ test {
// 			forvalues y=1998(1)2003 {
// 				gen `schl'_mod`m'byG_`t'_`y'=.
// 				forvalues g=4(1)8 {
// 					replace `schl'_mod`m'byG_`t'_`y'=mod`m'byG_`t'_`y'_`g' if future_grd==`g'
// 				}
// 			}
// 		}
// 	}
drop mod*
}

*NOTE - In the paper I use the neighborhood school as the counterfactual school attended by lottery winners. 
*An alternative way to construct the counterfactual is to use choice sets plus information on the ex post probability of admission.*
*You can do this using the "margin" variables - I sketch out the procedure below*

*Each "margin" variable is the percentage chance that an applicant gets into that school, given the ranking of the choice and the applicant's priority group*
*One minus margin3 is the percentage chance that they will end up in their neighborhood school (unless choice3 is the neighborhood school)*

*For the analysis, we don't care about the margin1 variable*
*VA diff will be the difference between the first choice and the weighted average of what they would receive if they did not get first choice*
*This total is going to be VAMalt = (margin2)*(VAMch2) + (1-margin2)*(VAMch3)
*VAMch3 = (margin3)*VAMch3 + (1-margin3)*VAMhome
*So the whole expression is VAMalt = (margin2)*(VAMch2) + (1-margin2)*[(margin3*VAMch3) + (1-margin3)*VAMhome]

*Note that if the 2nd choice is a school to which a student is guaranteed access, margin2=1 and thus VAMalt=VAMch2*
*Alternatively margin2 and margin3 will be very low numbers, and VAMalt will be close to VAMhome*

*this procedure produces slightly more accurate estimates, but I used neighborhood schools in the paper for ease of explanation*

// /*forvalues m=1(1)3 {*/ forvalues m=1(1)2 {
// 	foreach e in /*FE ar*/ mix {
// 		foreach s in 02 2yr all {
// 			foreach t in /*math read*/ test {
// 				gen alt_mod`m'`e'_`s'_`t'=(margin2*ch2_mod`m'`e'_`s'_`t')+(((1-margin2)*(margin3*ch3_mod`m'`e'_`s'_`t'))+((1-margin3)*(hm_mod`m'`e'_`s'_`t')))
// 			}
// 		}
// 	}
// }
*Analysis Tables and Figures*

*Table A1 - descriptives*
egen grade_home=group(home0203 future_grd)
tabstat /*male black hisp lunch*/ mathz2002 readz2002 if testz2003!=. & as_mod1mix_02_test!=. & onmargin!=1, s(mean count)
tabstat /*male black hisp lunch*/ mathz2002 readz2002 if testz2003!=. & as_mod1mix_02_test!=. & onmargin==1, s(mean count)
areg onmargin /*male black hisp lunch*/ mathz2002 readz2002 if testz2003!=. & as_mod1mix_02_test!=., absorb(grade_home) robust


*Now create VAM that combines years according to correlation with 2003 - they are essentially weights - follows Chetty, Friedman and Rockoff (2013)*
/*foreach x in as ch1 hm {
	gen `x'_mod3byG_test_weight=.0275*`x'_mod3byG_test_1998+.1332*`x'_mod3byG_test_1999+.2383*`x'_mod3byG_test_2000+.2443*`x'_mod3byG_test_2001+.3567*`x'_mod3byG_test_2002
}
*/
*Table 1 - lottery is IV, VA in d20 school is endogenous variable, test score is outcome*
*IV=lottery 
set more off
xtset lottery_FE
foreach t in test /*math read*/{
	/*forvalues m=1(1)4 {*/ forvalues m=1(1)2 {
		foreach e in /*ar FE*/mix  {
			foreach s in 02 2yr all {
				foreach c in hm /*alt*/ {
					//t m e s c
					gen VA=as_mod`m'`e'_`s'_`t'
					bootstrap, reps(100) cluster(lottery_FE): xtivreg2 `t'z2003 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (VA=lottery), fe cluster(lottery_FE)
					test (VA=1)
					drop VA
				}
			}
		}
	}
}
*IV=lott_VA (originl)
set more off
xtset lottery_FE
foreach t in test /*math read*/{
	/*forvalues m=1(1)4 {*/ forvalues m=1(1)2 {
		foreach e in /*ar FE*/mix  {
			foreach s in 02 2yr all {
				foreach c in hm /*alt*/ {
					//t m e s c
					gen VA=as_mod`m'`e'_`s'_`t'
					gen lott_VA=`c'_mod`m'`e'_`s'_`t' if lottery==0
					replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
					bootstrap, reps(100) cluster(lottery_FE): xtivreg2 `t'z2003 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (VA=lott_VA), fe cluster(lottery_FE)
					test (VA=1)
					drop VA lott_VA
				}
			}
		}
	}
}
*IV= new_lott_VA
set more off
xtset lottery_FE
foreach t in test /*math read*/{
	/*forvalues m=1(1)4 {*/ forvalues m=1(1)2 {
		foreach e in /*ar FE*/mix  {
			foreach s in 02 2yr all {
				foreach c in hm /*alt*/ {
					//t m e s c
					gen VA=as_mod`m'`e'_`s'_`t'
					gen lott_VA=0
					replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
					bootstrap, reps(100) cluster(lottery_FE): xtivreg2 `t'z2003 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (VA=lott_VA), fe cluster(lottery_FE)
					test (VA=1)
					drop VA lott_VA
				}
			}
		}
	}
}
/*
gen VA=as_mod3byG_test_weight
gen lott_VA=hm_mod3byG_test_weight if lottery==0
replace lott_VA=ch1_mod3byG_test_weight if lottery==1
xi: xtivreg2 testz2003 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (VA=lott_VA), fe cluster(lottery_FE)
test (VA=1)
drop VA lott_VA
*/
*Now do 2004 outcomes - Table A2*

/*
foreach t in test /*math read*/ {
	/*forvalues m=1(1)3 {*/ forvalues m=1(1)2 {
		foreach e in ar mix {
			foreach s in 02 2yr all {
				gen VA=as_mod`m'`e'_`s'_`t'
				gen lott_VA=hm_mod`m'`e'_`s'_`t' if lottery==0
				replace lott_VA=ch1_mod`m'`e'_`s'_`t' if lottery==1
				
				gen leadVA=aslead_mod`m'`e'_`s'_`t'
				gen lott_leadVA=hmlead_mod`m'`e'_`s'_`t' if lottery==0
				replace lott_leadVA=ch1lead_mod`m'`e'_`s'_`t' if lottery==1

				bootstrap, reps(100) cluster(lottery_FE): xtivreg2 `t'z2004 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (leadVA=lott_leadVA) if future_grd!=8, fe cluster(lottery_FE)
				test (leadVA=1)
				bootstrap, reps(100) cluster(lottery_FE): xtivreg2 `t'z2004 /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (VA leadVA=lott_VA lott_leadVA) if future_grd!=8, fe cluster(lottery_FE)
				test (VA=1) (leadVA=1)
				test (VA+leadVA=1)
				drop VA lott_VA leadVA lott_leadVA
			}
		}
	}
}

*Table A3 - impact of winning the lottery on first stage, school Xs, achievement*

*gen black_hisp=black==1 | hisp==1
foreach x in /*black_hisp lunch*/ mathz2002 readz2002 {
	bysort schl_d20: egen pct`x'=mean(`x')
}
foreach x in enrolled enroll_home attend_magnet {
	xi: xtreg `x' lottery /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss, fe vce(cluster lottery_FE)
}
foreach x in /*pctblack_hisp pctlunch*/ pctmathz2002 pctreadz2002 mathz2003 readz2003 mathz2004 readz2004 {
	xi: xtivreg2 `x' /*male black hisp lunch*/ math_2002_imp read_2002_imp math_2002_imp_sq math_2002_imp_cub read_2002_imp_sq read_2002_imp_cub math_2002_miss read_2002_miss (enrolled=lottery), fe cluster(lottery_FE)
}
*/

