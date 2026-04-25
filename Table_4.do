************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates main paper, Table 4 (Pooled Treatment Effects of Night-Sleep and Nap Treatments)
* Last edited: 07 May 2021
************************************************
		clear all
		set more off
		set seed 1

		//Initiate anderson indices and fixed effects programs
		do "$s/anderson_index_program.do"
		cap reghdfe, compile
		
***************
*WORK OUTCOMES*
***************

	use "$d/typing_dataset.dta", clear

	egen treat_pool_diff = mean(treat_pool), by (pid)
	replace treat_pool=1 if treat_pool_diff!=0
	drop treat_pool_diff

	egen _treat_nap = max(treat_nap), by(pid)
	drop treat_nap
	rename  _treat_nap treat_nap

	drop if day_in_study==.
	drop if day_in_study==1

	rename productivity prod
	rename typing_time_hr typing

		*Standardize
		foreach var in prod typing earnings output {
			sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
			gen `var'_post_std = (`var' - r(mean)) / r(sd)

			bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
			bysort pid: egen `var'_pre = mean(`var'_pre_temp)
			sum `var'_pre if treat_pool == 0 & treat_nap == 0 & post_treatment==0
			gen `var'_pre_std = (`var'_pre - r(mean)) /r(sd)
			drop `var'_pre `var'_pre_temp
		}

		* Night Sleep
		foreach var in prod typing earnings output {

			cap drop baseline
			egen baseline 	= mean(`var'_pre_std), by(pid)

			reghdfe `var'_post_std treat_pool treat_nap extra_work baseline i.age fraction_high female long_day if post_treatment==1 & at_present_check==1, absorb(date day_in_study) cluster(pid)

			if "`var'" != "labor_index"{
				local pids_`var'_ns = e(N_clust)
				local obs_`var'_ns = e(N)

				lincom _b[treat_pool]
					local coef_`var'_ns = string(r(estimate),"%3.2f")
					local se_`var'_ns = string(r(se),"%3.2f")
					local tstat_`var'_ns = r(estimate)/r(se)
					local p_`var'_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_`var'_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

				lincom _b[baseline]
					local coef_`var'_base_ns = string(r(estimate),"%3.2f")
					local se_`var'_base_ns = string(r(se),"%3.2f")
					local p_`var'_base_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_`var'_base_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			else {
				local pids_avg_`var'_ns = e(N_clust)
				local obs_avg_`var'_ns = e(N)

				lincom _b[treat_pool]
					local coef_avg_`var'_ns = string(r(estimate),"%3.2f")
					local se_avg_`var'_ns = string(r(se),"%3.2f")
					local p_avg_`var'_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_avg_`var'_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))


				lincom _b[baseline]
					local coef_avg_`var'_base_ns = string(r(estimate),"%3.2f")
					local se_avg_`var'_base_ns = string(r(se),"%3.2f")
					local p_avg_`var'_base_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_avg_`var'_base_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		}

		* Naps
		foreach var in prod typing earnings output {

			cap drop baseline
			egen baseline 	= mean(`var'_pre_std), by(pid)

			reghdfe `var'_post_std treat_pool treat_nap baseline i.age fraction_high female long_day if post_treatment==1 & at_present_check==1 & day_in_study!=28, absorb(date day_in_study) cluster(pid)


				local pids_`var'_nap = e(N_clust)
				local obs_`var'_nap = e(N)

				lincom _b[treat_nap]
					local coef_`var'_nap = string(r(estimate),"%3.2f")
					local se_`var'_nap = string(r(se),"%3.2f")
					local tstat_`var'_nap = r(estimate)/r(se)
					local p_`var'_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_`var'_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

				lincom _b[baseline]
					local coef_`var'_base_nap = string(r(estimate),"%3.2f")
					local se_`var'_base_nap = string(r(se),"%3.2f")
					local p_`var'_base_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_`var'_base_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			
				* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_`var'_d_nap_ns = string(r(p), "%3.2f")
				local p_`var'_d_nap_ns1 = r(p)

				
			reghdfe `var'_post_std treat_pool treat_nap extra_work baseline i.age fraction_high female long_day if post_treatment==1 & at_present_check==1 & day_in_study!=28, absorb(date day_in_study) cluster(pid)

				local pids_`var'_nap = e(N_clust)
				local obs_`var'_nap = e(N)

				lincom _b[treat_nap]
					local coef_`var'_nap_break = string(r(estimate),"%3.2f")
					local se_`var'_nap_break = string(r(se),"%3.2f")
					local tstat_`var'_nap_break = r(estimate)/r(se)
					local p_`var'_nap_break = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_`var'_nap_break1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

				lincom _b[treat_nap] - _b[extra_work]
					local coef_`var'_nap_work = string(r(estimate),"%3.2f")
					local se_`var'_nap_work = string(r(se),"%3.2f")
					local tstat_`var'_nap_work = r(estimate)/r(se)
					local p_`var'_nap_work = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")

				lincom _b[baseline]
					local coef_`var'_base_nap = string(r(estimate),"%3.2f")
					local se_`var'_base_nap = string(r(se),"%3.2f")
					local p_`var'_base_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
		}
** Aggregate earnings in the participant level

* Residualizing the typing variables
describe earnings, varlist
	local type_vars `r(varlist)'
	foreach var in `type_vars' {
		reghdfe `var'  , absorb(day_in_study date) residuals
		cap drop `var'
		rename _reghdfe_resid `var'
	}

* Standardize residualized version
foreach var in  earnings {
			cap drop `var'_post_std
			cap drop `var'_pre
			cap drop `var'_pre_temp

			sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
			gen `var'_post_std = (`var' - r(mean)) / r(sd)

			bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
			bysort pid: egen `var'_pre = mean(`var'_pre_temp)
			sum `var'_pre if treat_pool == 0 & treat_nap == 0 & post_treatment==0

			gen `var'_pre_std = (`var'_pre - r(mean)) /r(sd)
			drop `var'_pre `var'_pre_temp
		}

* Generating different types of earnings for Nap comparisons
bys pid: egen earnings_post_std_break = mean(earnings_post_std) if extra_work==0 & post_treatment==1 & treat_nap==0 & day_in_study!=28 & at_present_check==1
bys pid: egen earnings_post_std_work = mean(earnings_post_std) if extra_work==1	& post_treatment==1 & treat_nap==0 & day_in_study!=28 & at_present_check==1

bys pid: egen earnings_break = mean(earnings) if extra_work==0  & post_treatment==1 & treat_nap==0 & day_in_study!=28 & at_present_check==1
bys pid: egen earnings_work = mean(earnings)  if extra_work==1	& post_treatment==1 & treat_nap==0 & day_in_study!=28 & at_present_check==1

keep if at_present_check==1

	collapse (mean) prod earnings* typing output (max) treat_pool treat_nap female age, by(pid post_treatment)

	keep if !mi(post_treatment)

	tempfile work
	save `work', replace


**************
* WELL-BEING *
**************

	use "$d/health_dataset.dta", clear

	* Psych Wellbeing
	* note: code will only work properly if you run it from preserve to restore in a single run
	preserve

	*residualize
	describe ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep, varlist
	local wb_vars `r(varlist)'
	foreach var in `wb_vars' {
		reghdfe `var'  , absorb(day_in_study date) residuals
		cap drop `var'
		rename _reghdfe_resid `var'
	}


	collapse (mean) ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep (max) treat_pool treat_nap female age, by(pid post_treatment)

	* Weighted index (Anderson 2008)

		* Step 1: define locals for the relevant arguments

			*Components
			local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep"

			*Covariates to be included linearly
			local covariates "female"

			*Treatment dummies
			local treatments "treat_nap treat_pool"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
			replace es_dep = -es_dep
			replace ds_g13c_stressed = -ds_g13c_stressed

		* Step 3: Just run this line! Note: this function will also generate variables necessary below
			anderson_index "`components'" "`treatments'" "`covariates'" "`factorcovariates'" "`subset'" `eststo'

			local pids_psych = e(N)
			local obs_psych = e(N)

			lincom _b[treat_pool]
				local coef_psych_ns = string(r(estimate),"%3.2f")
				local se_psych_ns = string(r(se),"%3.2f")
				local tstat_psych_ns = r(estimate)/r(se)
				local p_psych_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_psych_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap]
				local coef_psych_nap = string(r(estimate),"%3.2f")
				local se_psych_nap = string(r(se),"%3.2f")
				local tstat_psych_nap = r(estimate)/r(se)
				local p_psych_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_psych_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[summary_index_base]
				local coef_psych_base = string(r(estimate),"%3.2f")
				local se_psych_base = string(r(se),"%3.2f")
				local p_psych_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_psych_d_nap_ns = string(r(p), "%3.2f")
				local p_psych_d_nap_ns1 = r(p)


	tempfile psych
	save `psych', replace

	restore

	* Phys Wellbeing

	local health_vars "distance max_speed win_systolic_bp win_diastolic_bp"

	* Create standardized versions of postline variables
	foreach var in `health_vars' {
		sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
		gen `var'_post_std = (`var' - r(mean)) / r(sd)
	}

	*Create biking_task and win_bp
	egen biking_task = rowtotal(max_speed_post_std distance_post_std) if post_treatment == 1
	egen nomiss_post = rownonmiss(max_speed_post_std distance_post_std)
	replace biking_task = biking_task/nomiss_post if post_treatment == 1

	egen win_bp = rowtotal(win_systolic_bp_post_std win_diastolic_bp_post_std) if post_treatment == 1
	cap drop nomiss_post
	egen nomiss_post = rownonmiss(win_systolic_bp_post_std win_diastolic_bp_post_std)
	replace win_bp = win_bp/nomiss_post if post_treatment == 1

	* Create standardized versions of baseline variables
	foreach var in `health_vars'{
		bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
		bysort pid: egen `var'_pre = mean(`var'_pre_temp)
		sum `var'_pre if treat_pool == 0 & treat_nap == 0
		gen `var'_pre_std = (`var'_pre - r(mean)) / r(sd)
		drop `var'_pre `var'_pre_temp
	}

	* Create win_bp (summing blood pressure at pre-treatment)
	egen win_bp_pre = rowtotal(win_systolic_bp_pre_std win_diastolic_bp_pre_std) if post_treatment == 0
	cap drop nomiss_pre
	egen nomiss_pre = rownonmiss(win_systolic_bp_pre_std win_diastolic_bp_pre_std)
	replace win_bp_pre = win_bp_pre/nomiss_pre if post_treatment == 0

	replace win_bp = win_bp_pre if post_treatment == 0
	drop win_bp_pre

	* Adjust sign for different measures
	foreach var of varlist win_bp es_ill_week es_pain es_daily_act{
		replace `var' = -`var'
	}

	* Residualize
	local components "win_bp"
	foreach var in `components' {
		reghdfe `var' , absorb(day_in_study date) residuals
		cap drop `var'
		rename _reghdfe_resid `var'
	}
	local components "biking_task win_bp es_ill_week es_pain es_daily_act"
	collapse (mean) `components' (max) female age treat_pool treat_nap, by(pid post_treatment)

	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

		*Components
		local components "biking_task win_bp es_ill_week es_pain es_daily_act"

		*Covariates to be included linearly
		local controls "female"

		*Treatment dummies
		local treatments "treat_nap treat_pool"

		*Subset
		local subset "if post_treatment==1"

		*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
		local factorcovariates "age"

		*Eststo? Y = 1, N = 0
		local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_physical = e(N)
		local obs_physical = e(N)

		lincom _b[treat_pool]
			local coef_physical_ns = string(r(estimate),"%3.2f")
			local se_physical_ns = string(r(se),"%3.2f")
			local tstat_phys_ns = r(estimate)/r(se)
			local p_physical_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_physical_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_physical_nap = string(r(estimate),"%3.2f")
			local se_physical_nap = string(r(se),"%3.2f")
			local tstat_phys_nap = r(estimate)/r(se)
			local p_physical_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_physical_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
		lincom _b[summary_index_base]
			local coef_physical_base = string(r(estimate),"%3.2f")
			local se_physical_base = string(r(se),"%3.2f")
			local p_physical_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_physical_d_nap_ns = string(r(p), "%3.2f")
				local p_physical_d_nap_ns1 = r(p)


*******************
* WELLBEING INDEX *
*******************

	merge 1:1 pid post_treatment using `psych'
	
	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

		*Components
		local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act"

		*Covariates to be included linearly
		local controls "female"

		*Treatment dummies
		local treatments "treat_nap treat_pool"

		*Subset
		local subset "if post_treatment==1"

		*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
		local factorcovariates "age"

		*Eststo? Y = 1, N = 0
		local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_wellbeing = e(N)
		local obs_wellbeing = e(N)

		lincom _b[treat_pool]
			local coef_wellbeing_ns = string(r(estimate),"%3.2f")
			local se_wellbeing_ns = string(r(se),"%3.2f")
			local tstat_wellbeing_ns = r(estimate)/r(se)
			local p_wellbeing_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_wellbeing_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_wellbeing_nap = string(r(estimate),"%3.2f")
			local se_wellbeing_nap = string(r(se),"%3.2f")
			local tstat_wellbeing_nap = r(estimate)/r(se)
			local p_wellbeing_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_wellbeing_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_wellbeing_base = string(r(estimate),"%3.2f")
			local se_wellbeing_base = string(r(se),"%3.2f")
			local p_wellbeing_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
							
				lincom _b[treat_pool] - _b[treat_nap]
				local p_wellbeing_d_nap_ns = string(r(p), "%3.2f")
				local p_wellbeing_d_nap_ns1 = r(p)

	* Store Anderson index
	gen wellbeing_index       = summary_index
	gen wellbeing_index_base  = summary_index_base

	keep pid post_treatment ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act wellbeing_index wellbeing_index_base

	tempfile wellbeing
	save `wellbeing'
		
*******************************
* RISK AND SOCIAL PREFERENCES *
*******************************

use "$d/riskandsocial_dataset.dta", clear

	collapse (mean) risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent (max) female age treat_pool treat_nap, by(pid post_treatment)


	*** Risk Preferences

	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

			*Components
			local components "risk_switch_point loss_switch_point"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_nap treat_pool"

			*Subset
			local subset "if post_treatment==1"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			*Eststo? Y = 1, N = 0
			local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
	
	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_risk = e(N)
		local obs_risk = e(N)

		lincom _b[treat_pool]
			local coef_risk_ns = string(r(estimate),"%3.2f")
			local se_risk_ns = string(r(se),"%3.2f")
			local tstat_risk_ns = r(estimate)/r(se)
			local p_risk_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_risk_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_risk_nap = string(r(estimate),"%3.2f")
			local se_risk_nap = string(r(se),"%3.2f")
			local tstat_risk_nap = r(estimate)/r(se)
			local p_risk_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_risk_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_risk_base = string(r(estimate),"%3.2f")
			local se_risk_base = string(r(se), "%3.2f")
			local p_risk_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_risk_d_nap_ns = string(r(p), "%3.2f")
				local p_risk_d_nap_ns1 = r(p)


	*** Social Preferences

	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

		*Components
		local components "d_send u_send t_send u_avg_receive t_avg_amountsent"

		*Covariates to be included linearly
		local controls "female"

		*Treatment dummies
		local treatments "treat_nap treat_pool"

		*Subset
		local subset "if post_treatment==1"

		*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
		local factorcovariates "age"

		*Eststo? Y = 1, N = 0
		local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_social = e(N_clust)
		local obs_social = e(N)

		lincom _b[treat_pool]
			local coef_social_ns = string(r(estimate),"%3.2f")
			local se_social_ns = string(r(se),"%3.2f")
			local tstat_social_ns = r(estimate)/r(se)
			local p_social_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_social_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_social_nap = string(r(estimate),"%3.2f")
			local se_social_nap = string(r(se),"%3.2f")
			local tstat_social_nap = r(estimate)/r(se)
			local p_social_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_social_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_social_base = string(r(estimate),"%3.2f")
			local se_social_base = string(r(se),"%3.2f")
			local p_social_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_social_d_nap_ns = string(r(p), "%3.2f")
				local p_social_d_nap_ns1 = r(p)

	tempfile rs_pref
	save `rs_pref', replace

********************
* TIME PREFERENCES *
********************

	* Present Bias

		use "$d/pb_dataset.dta", clear

		keep if type == "censored"


		replace beta_post = 2 if beta_post > 2 & beta_post != .
		replace beta_base = 2 if beta_base > 2 & beta_base != .

		sum beta_base if convergence_code_post == 0 //only sample that's included in pb regressions
		replace beta_base = -9999999 if beta_base==.
		gen dummy_missing = (beta_base==-9999999)
		replace beta_base = r(mean) if beta_base==-9999999

		tempfile pb
		save `pb', replace

	* Savings

		use "$d/savings_dataset.dta", clear

		egen _treat_nap = max(treat_nap), by(pid)
		drop treat_nap
		rename  _treat_nap treat_nap

		egen treat_diff = mean(treat_pool), by(pid)
		replace treat_pool=1 if treat_diff!=0

		forval n=1/20{
			gen rewardrate`n'_dummy = (rewardrate`n' == float(0.02))
			replace rewardrate`n'_dummy = . if rewardrate`n' == .
		}

		egen fraction_high = rowtotal(rewardrate*_dummy)
		egen nomiss = rownonmiss(rewardrate*_dummy)
		replace fraction_high = fraction_high/nomiss

		* Residualized output
		local controls "female i.age risk_and_social_day default_amount rewardrate* max_pay_cog_tasks piece_rate_pb pid_interest"
		reghdfe deposits `controls', vce(cluster pid) absorb(surveyor date day_in_study) residuals
		cap drop deposits
		rename _reghdfe_resid deposits

		merge m:1 pid using `pb', gen(pb_merge)

		gen beta = beta_post
		replace beta = beta_base if post_treatment == 0

		collapse deposits beta dummy_missing (max) female age treat_pool treat_nap, by(pid post_treatment)


	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

		*Components
		local components "deposits beta"

		*Covariates to be included linearly
		local controls "female"

		*Treatment dummies
		local treatments "treat_nap treat_pool"

		*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
		local factorcovariates "i.age"

		*Subset
		local subset "if post_treatment==1"

		*Eststo? Y = 1, N = 0
		local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGH
		* signs are correct for time preferences already

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_timepref = e(N_clust)
		local obs_timepref = e(N)

		lincom _b[treat_pool]
			local coef_timepref_ns = string(r(estimate),"%3.2f")
			local se_timepref_ns = string(r(se),"%3.2f")
			local tstat_time_ns = r(estimate)/r(se)
			local p_timepref_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_timepref_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_timepref_nap = string(r(estimate),"%3.2f")
			local se_timepref_nap = string(r(se),"%3.2f")
			local tstat_time_nap = r(estimate)/r(se)
			local p_timepref_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_timepref_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_timepref_base = string(r(estimate),"%3.2f")
			local se_timepref_base = string(r(se),"%3.2f")
			local p_timepref_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_timepref_d_nap_ns = string(r(p), "%3.2f")
				local p_timepref_d_nap_ns1 = r(p)

*********************
* PREFERENCES INDEX *
*********************

	merge m:1 pid post_treatment using `rs_pref', gen(rs_merge)

	* Weighted (Anderson)


	* Step 1: define locals for the relevant arguments

			* Components
			local components "risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent deposits beta"

			* Covariates to be included linearly
			local controls "female dummy_missing"

			* Treatment dummies
			local treatments "treat_nap treat_pool"

			* Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			* Subset
			local subset "if post_treatment==1"

			* Eststo? Y = 1, N = 0
			local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
			* Signs have been changed already to make the average index!

	* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_pref = e(N_clust)
		local obs_pref = e(N)

		lincom _b[treat_pool]
			local coef_pref_ns = string(r(estimate),"%3.2f")
			local se_pref_ns = string(r(se),"%3.2f")
			local tstat_pref_ns = r(estimate)/r(se)
			local p_pref_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_pref_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_pref_nap = string(r(estimate),"%3.2f")
			local se_pref_nap = string(r(se),"%3.2f")
			local tstat_pref_nap = r(estimate)/r(se)
			local p_pref_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_pref_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))


		lincom _b[summary_index_base]
			local coef_pref_base = string(r(estimate),"%3.2f")
			local se_pref_base = string(r(se),"%3.2f")
			local p_pref_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_pref_d_nap_ns = string(r(p), "%3.2f")
				local p_pref_d_nap_ns1 = r(p)

	* Store Anderson index
	gen pref_post       = summary_index
	gen pref_pre  = summary_index_base

	keep pid post_treatment dummy_missing risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent deposits beta pref_post pref_pre

	tempfile preferences
	save `preferences', replace


**********************
* COGNITIVE FUNCTION *
**********************

	use "$d/cognitive_tasks_dataset.dta", clear

	merge 1:1 pid day_in_study using "$d/pvt_dataset.dta"

	egen treat_pool_diff = mean(treat_pool), by (pid)
	replace treat_pool=1 if treat_pool_diff!=0
	drop treat_pool_diff

	egen _treat_nap = max(treat_nap), by(pid)
	drop treat_nap
	rename  _treat_nap treat_nap

	egen max_pay_dummy = rowtotal(max_pay_hf_dummy max_pay_co_dummy)
	replace max_pay_dummy = 1 if max_pay_dummy == 2

	* Residualize
	local components "pv_perf hf_payment co_payment"
	foreach var in `components' {
		reghdfe `var', absorb(high_pay date day_in_study) residuals
		cap drop `var'
		rename _reghdfe_resid `var'
	}

	collapse `components' (max) female age treat_nap treat_pool, by(pid post_treatment)

	* Weighted index (Anderson 2008)

	* Step 1: define locals for the relevant arguments

			* Components
			local components "pv_perf hf_payment co_payment"

			* Covariates to be included linearly
			local controls "female"

			* Treatment dummies
			local treatments "treat_nap treat_pool"

			* Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			* Subset
			local subset "if post_treatment==1"

			* Eststo? Y = 1, N = 0
			local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		local pids_cogfunction = e(N_clust)
		local obs_cogfunction = e(N)

		lincom _b[treat_pool]
			local coef_cogfunction_ns = string(r(estimate),"%3.2f")
			local se_cogfunction_ns = string(r(se),"%3.2f")
			local tstat_cog_ns = r(estimate)/r(se)
			local p_cogfunction_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_cogfunction_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_cogfunction_nap = string(r(estimate),"%3.2f")
			local se_cogfunction_nap = string(r(se),"%3.2f")
			local tstat_cog_nap = r(estimate)/r(se)
			local p_cogfunction_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_cogfunction_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_cogfunction_base = string(r(estimate),"%3.2f")
			local se_cogfunction_base = string(r(se),"%3.2f")
			local p_cogfunction_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_cogfunction_d_nap_ns = string(r(p), "%3.2f")
				local p_cogfunction_d_nap_ns1 = r(p)

	tempfile cognitive
	save `cognitive', replace

*************
* ATTENTION *
*************

use "$d/salience_dataset.dta", clear

egen _treat_nap = max(treat_nap), by(pid)
drop treat_nap
rename  _treat_nap treat_nap

rename treat treat_pool
egen treat_pool_diff = mean(treat_pool), by (pid)
	replace treat_pool=1 if treat_pool_diff!=0
	drop treat_pool_diff

egen tag_pid = tag(pid) if post_treatment==1

* Residualized output
reghdfe output_section , absorb(day_in_study pid date) residuals
cap drop output_res
rename _reghdfe_resid output_res

* Generating average standardized output by high and salient sessions - post treatment
foreach var in  output_res  {

	bys pid: egen _avg_`var'_high_sal = mean(`var') if high==1 & salience==1 & post_treatment==1
	bys pid: egen avg_`var'_high_sal = mean(_avg_`var'_high_sal)

	bys pid: egen _avg_`var'_high_notsal = mean(`var') if high==1 & salience==0 & post_treatment==1
	bys pid: egen avg_`var'_high_notsal = mean(_avg_`var'_high_notsal)

	bys pid: egen _avg_`var'_low_sal = mean(`var') if high==0 & salience==1 & post_treatment==1
	bys pid: egen avg_`var'_low_sal = mean(_avg_`var'_low_sal)

	bys pid: egen _avg_`var'_low_notsal = mean(`var') if high==0 & salience==0 & post_treatment==1
	bys pid: egen avg_`var'_low_notsal = mean(_avg_`var'_low_notsal)

	gen attent_`var' = avg_`var'_high_sal - avg_`var'_high_notsal - (avg_`var'_low_sal - avg_`var'_low_notsal)

}

* Generating average standardized output by high and salient sessions - baseline

foreach var in output_res  {

	bys pid: egen _avg_`var'_high_sal_bas = mean(`var') if high==1 & salience==1 & post_treatment==0
	bys pid: egen avg_`var'_high_sal_bas = mean(_avg_`var'_high_sal_bas)

	bys pid: egen _avg_`var'_high_notsal_bas = mean(`var') if high==1 & salience==0 & post_treatment==0
	bys pid: egen avg_`var'_high_notsal_bas = mean(_avg_`var'_high_notsal_bas)

	bys pid: egen _avg_`var'_low_sal_bas = mean(`var') if high==0 & salience==1 & post_treatment==0
	bys pid: egen avg_`var'_low_sal_bas = mean(_avg_`var'_low_sal_bas)

	bys pid: egen _avg_`var'_low_notsal_bas = mean(`var') if high==0 & salience==0 & post_treatment==0
	bys pid: egen avg_`var'_low_notsal_bas = mean(_avg_`var'_low_notsal_bas)

	gen attent_`var'_bas = avg_`var'_high_sal_bas - avg_`var'_high_notsal_bas - (avg_`var'_low_sal_bas - avg_`var'_low_notsal_bas)

}

* standardize variables

foreach var in attent_output_res  {

	sum `var' if treat_pool==0 & treat_nap==0 & tag_pid == 1
	gen std_`var' = (`var' - `r(mean)')/`r(sd)'

	sum `var'_bas if treat_pool==0 & treat_nap==0 & tag_pid == 1
	gen std_`var'_base = (`var'_bas - `r(mean)')/`r(sd)'

}

	* Changing sign - IMPORTANT
	replace std_attent_output_res = -std_attent_output_res
	replace std_attent_output_res_base = -std_attent_output_res_base

	rename std_attent_output_res_base baseline

	* Collapsing

	collapse (mean) std_attent_output_res  baseline (max) treat_pool treat_nap age female, by(pid post_treatment)


	reghdfe std_attent_output_res treat_pool treat_nap baseline if post_treatment == 1, absorb(age female) vce(robust)

	rename baseline std_attent_output_base
	rename std_attent_output_res std_attent_output

	local pids_attention = e(N)
	local obs_attention = e(N)

	lincom _b[treat_pool]
		local coef_attention_ns = string(r(estimate),"%3.2f")
		local se_attention_ns = string(r(se),"%3.2f")
		local tstat_salience_ns = r(estimate)/r(se)
		local p_attention_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
		local p_attention_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

	lincom _b[treat_nap]
		local coef_attention_nap = string(r(estimate),"%3.2f")
		local se_attention_nap = string(r(se),"%3.2f")
		local tstat_salience_nap = r(estimate)/r(se)
		local p_attention_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
		local p_attention_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

	lincom _b[baseline]
		local coef_attention_base = string(r(estimate),"%3.2f")
		local se_attention_base = string(r(se),"%3.2f")
		local p_attention_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
		
		* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_attention_d_nap_ns = string(r(p), "%3.2f")
				local p_attention_d_nap_ns1 = r(p)
				
*******************
* COGNITION INDEX *
*******************

	merge 1:1 pid post_treatment using `cognitive', gen(cog_merge)

	replace std_attent_output = std_attent_output_base if post_treatment==0

	* Weighted index

	* Step 1: define locals for the relevant arguments

		*Components
		local components "pv_perf hf_payment co_payment std_attent_output"

		*Covariates to be included linearly
		local controls "female"

		*Treatment dummies
		local treatments "treat_nap treat_pool"

		*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
		local factorcovariates "age"

		*Subset
		local subset "if post_treatment==1"

		*Eststo? Y = 1, N = 0
		local eststo 0

	* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

	* Step 3: Just run this line!
		anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

		foreach var in pv_perf hf_payment co_payment std_attent_output {

		sum `var'_post_std if post_treatment==1
		sum `var'_pre_std if post_treatment==0

		}

		foreach var in female treat_pool treat_nap age {

		sum `var' if post_treatment==1, de

		}

		local pids_cogindex = e(N_clust)
		local obs_cogindex = e(N)

		lincom _b[treat_pool]
			local coef_cogindex_ns = string(r(estimate),"%3.2f")
			local se_cogindex_ns = string(r(se),"%3.2f")
			local tstat_cogindex_ns = r(estimate)/r(se)
			local p_cogindex_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_cogindex_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap]
			local coef_cogindex_nap = string(r(estimate),"%3.2f")
			local se_cogindex_nap = string(r(se),"%3.2f")
			local tstat_cogindex_nap = r(estimate)/r(se)
			local p_cogindex_nap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_cogindex_nap1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[summary_index_base]
			local coef_cogindex_base = string(r(estimate),"%3.2f")
			local se_cogindex_base = string(r(se),"%3.2f")
			local p_cogindex_base = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			
			* p-values
				
				lincom _b[treat_pool] - _b[treat_nap]
				local p_cogindex_d_nap_ns = string(r(p), "%3.2f")
				local p_cogindex_d_nap_ns1 = r(p)

	* Store Anderson index
	gen cogindex_post       = summary_index
	gen cogindex_pre  = summary_index_base


	keep pid post_treatment pv_perf hf_payment co_payment std_attent_output cogindex_post cogindex_pre

	keep if !mi(post_treatment)

	tempfile cognition
	save `cognition', replace

*****************
* OVERALL INDEX *
*****************

	use `work'
	keep if !mi(pid)
	merge 1:1 pid post_treatment using `wellbeing', gen(wellbeing_merge)
	merge 1:1 pid post_treatment using `preferences', gen(pref_merge)
	merge 1:1 pid post_treatment using `cognition', gen(cog_merge)
	
	save "$d/main_index_dataset.dta", replace 

	* Correcting baseline values of indices

		replace cogindex_post = cogindex_pre if post_treatment==0
		replace pref_post = pref_pre if post_treatment==0
		replace wellbeing_index = wellbeing_index_base if post_treatment==0

		replace earnings_post_std = earnings_pre_std if post_treatment==0
		replace earnings_post_std_break = earnings_pre_std if post_treatment==0
		replace earnings_post_std_work = earnings_pre_std if post_treatment==0

		* Earnings for Nap vs. Break
		gen _earnings_post_std_break = earnings_post_std if treat_nap==1
		replace _earnings_post_std_break = earnings_post_std_break if treat_nap==0

		* Earnings for Nap vs. Work
		gen _earnings_post_std_work = earnings_post_std if treat_nap==1
		replace _earnings_post_std_work = earnings_post_std_work if treat_nap==0

	eststo clear

		* Step 1: define locals for the relevant arguments

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_nap treat_pool"

				*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
				local factorcovariates "i.age"

				*Eststo? Y = 1, N = 0
				local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
				* Signs have been changed already to make the average index!

		* Step 3: Just run this line!
				local subset "if post_treatment==1"

				* 1. Nap vs Pooled
				local components "earnings_post_std pref_post wellbeing_index cogindex_post"
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" 1

				local pids_overall_napnonap = e(N_clust)
				local obs_overall_napnonap = e(N)

				lincom _b[treat_pool]
				local se_overall_ns = string(r(se),"%3.2f")
				local coef_overall_ns = string(r(estimate),"%3.2f")
				local p_overall_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")

				lincom _b[treat_nap]
				local coef_overall_nap_nonap = string(r(estimate),"%3.2f")
				local se_overall_nap_nonap = string(r(se),"%3.2f")
				local p_overall_nap_nonap = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				
				* p-val 
				
				* Get pval of differences
				lincom _b[treat_pool] - _b[treat_nap]
				local p_overall_d_nap_ns = string(r(p), "%3.2f")
			
				* Save the overall index
				gen overall_post  = summary_index
				gen overall_pre   = summary_index_base
				
				* 2. Nap vs Break
				cap drop _earnings
				gen _earnings = earnings_post_std_break if treat_nap==0
				replace _earnings = earnings_post_std if treat_nap==1

				local components "_earnings pref_post wellbeing_index cogindex_post"
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" 1

				local pids_avg_overall_napbreak = e(N_clust)
				local obs_avg_overall_napbreak = e(N)

				lincom _b[treat_nap]
				local coef_overall_nap_break = string(r(estimate),"%3.2f")
				local se_overall_nap_break = string(r(se),"%3.2f")
				local p_overall_nap_break = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")

				* 3. Nap vs Work
				cap drop _earnings
				gen _earnings = earnings_post_std_work if treat_nap==0
				replace _earnings = earnings_post_std if treat_nap==1

				local components "_earnings pref_post wellbeing_index cogindex_post"
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" 1

				local pids_overall_napwork = e(N_clust)
				local obs_overall_napwork = e(N)

				lincom _b[treat_nap]
				local coef_overall_nap_work = string(r(estimate),"%3.2f")
				local se_overall_nap_work = string(r(se),"%3.2f")
				local p_overall_nap_work = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")

save "$d/anderson_indices.dta", replace

clear
set obs 14
gen var_name = ""
gen p_val_ns1  = .
gen p_val_nap1 = .
gen p_val_d_nap_ns1 = .

local i = 1
foreach var in earnings prod typing output wellbeing physical psych cogindex cogfunction attention pref timepref social risk {
	replace var_name = "`var'" if _n == `i'
	foreach comp in ns1 nap1 d_nap_ns1 {
	di "`var'_`comp'"
		replace p_val_`comp' = `p_`var'_`comp''  if _n == `i'
	}
	local i = `i' + 1

}

save "$d/pvals_unadjusted.dta", replace 
export delimited "$d/pvals_unadjusted.csv", replace


/* Runs the randomization inference iterations and the p-value step-down
   procedure if RI is set as 1 in the master do-file */
if $RI == 1{

clear 
do "$sm/Table_4_RI"
clear

if "`c(os)'"=="MacOSX" | "`c(os)'"=="UNIX" {
    rsource using Scripts/Main_paper/Table_4_adjust_pval.R, rpath("/usr/local/bin/R") roptions(`"--vanilla"')
}
else {  // windows
    rsource using Scripts/Main_paper/Table_4_adjust_pval.R, rpath(`"c:\r\R-3.5.1\bin\Rterm.exe"') roptions(`"--vanilla"')  // change version number, if necessary
}

}

preserve 
//Importing adjusted p-values
import delimited "$d/pval_adj_table_pool.csv", clear 
gen p_val_adj_  = string(p_val_adj, "%03.2f")

local pval_w_mht_ns_earnings = p_val_adj_[1]
local pval_w_mht_ns_prod = p_val_adj_[13]
local pval_w_mht_ns_typing = p_val_adj_[14]
local pval_w_mht_ns_output = p_val_adj_[15]
local pval_w_mht_ns_wellbeing = p_val_adj_[2]
local pval_w_mht_ns_phys = p_val_adj_[22]
local pval_w_mht_ns_psych = p_val_adj_[23]

local pval_w_mht_nap_earnings = p_val_adj_[5]
local pval_w_mht_nap_prod = p_val_adj_[16]
local pval_w_mht_nap_typing = p_val_adj_[17]
local pval_w_mht_nap_output = p_val_adj_[18]
local pval_w_mht_nap_wellbeing = p_val_adj_[6]
local pval_w_mht_nap_phys = p_val_adj_[24]
local pval_w_mht_nap_psych = p_val_adj_[25]

local pval_w_mht_d_earnings = p_val_adj_[9]
local pval_w_mht_d_prod = p_val_adj_[19]  
local pval_w_mht_d_typing = p_val_adj_[20]
local pval_w_mht_d_output = p_val_adj_[21]
local pval_w_mht_d_wellbeing = p_val_adj_[10]
local pval_w_mht_d_phys = p_val_adj_[26]
local pval_w_mht_d_psych = p_val_adj_[27]

local pval_w_mht_ns_cogindex = p_val_adj_[3]
local pval_w_mht_ns_cog = p_val_adj_[28]
local pval_w_mht_ns_salience = p_val_adj_[29]
local pval_w_mht_ns_pref = p_val_adj_[4]
local pval_w_mht_ns_time = p_val_adj_[34]
local pval_w_mht_ns_social = p_val_adj_[35]
local pval_w_mht_ns_risk = p_val_adj_[36]

local pval_w_mht_nap_cogindex = p_val_adj_[7]
local pval_w_mht_nap_cog = p_val_adj_[30]
local pval_w_mht_nap_salience = p_val_adj_[31]
local pval_w_mht_nap_pref = p_val_adj_[8]
local pval_w_mht_nap_time = p_val_adj_[37]
local pval_w_mht_nap_social = p_val_adj_[38]
local pval_w_mht_nap_risk = p_val_adj_[39]

local pval_w_mht_d_cogindex = p_val_adj_[11]
local pval_w_mht_d_cog = p_val_adj_[32]
local pval_w_mht_d_salience = p_val_adj_[33]
local pval_w_mht_d_pref = p_val_adj_[12]
local pval_w_mht_d_time = p_val_adj_[40]
local pval_w_mht_d_social = p_val_adj_[41]
local pval_w_mht_d_risk = p_val_adj_[42]

restore  

* Outputting the results

	cd "$om/Tables"
		file open f using "Table4A_table_main_anderson_join1.tex", write replace
		file write f "\begin{tabular*}{1.7\textwidth}{l@{\extracolsep{\fill}}*{8}{c}}" _n ///
			"\toprule" _n ///
			"& \multicolumn{1}{c}{\textbf{OVERALL}} & \multicolumn{4}{c}{\textbf{WORK}} & \multicolumn{3}{c}{\textbf{WELL-BEING}}\\" _n ///
			///
			"\cmidrule(lr){2-2}\cmidrule(lr){3-6}\cmidrule(lr){7-9} & \textbf{Index} & \textbf{Earnings} & Productivity & Labor Supply & Output & \textbf{Index} & Physical & Mental \\" _n ///
			///
			" & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8)\\" _n ///
			"\midrule" _n ///
			///
			"\textbf{Night-Sleep Treatments} & \textbf{`coef_overall_ns'} & \textbf{`coef_earnings_ns'} & `coef_prod_ns' & `coef_typing_ns' & `coef_output_ns' & \textbf{`coef_wellbeing_ns'} & `coef_physical_ns' & `coef_psych_ns' \\" _n ///
			"& (`se_overall_ns')& (`se_earnings_ns') & (`se_prod_ns') & (`se_typing_ns') & (`se_output_ns')  & (`se_wellbeing_ns') & (`se_physical_ns') & (`se_psych_ns') \\" _n ///
			"& \{`p_overall_ns'\} & \{`p_earnings_ns'\} & \{`p_prod_ns'\} & \{`p_typing_ns'\} & \{`p_output_ns'\}  & \{`p_wellbeing_ns'\} & \{`p_physical_ns'\} & \{`p_psych_ns'\} \\" _n ///
			"& & [`pval_w_mht_ns_earnings'] & [`pval_w_mht_ns_prod'] & [`pval_w_mht_ns_typing'] & [`pval_w_mht_ns_output'] & [`pval_w_mht_ns_wellbeing'] & [`pval_w_mht_ns_phys'] & [`pval_w_mht_ns_psych'] \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Nap Treatment} & \textbf{`coef_overall_nap_nonap'} & \textbf{`coef_earnings_nap'} & `coef_prod_nap' & `coef_typing_nap' & `coef_output_nap'  & \textbf{`coef_wellbeing_nap'} & `coef_physical_nap' & `coef_psych_nap' \\" _n ///
			"& (`se_overall_nap_nonap')& (`se_earnings_nap') & (`se_prod_nap') & (`se_typing_nap') & (`se_output_nap')  & (`se_wellbeing_nap')  & (`se_physical_nap') &  (`se_psych_nap') \\" _n ///
			"& \{`p_overall_nap_nonap'\} & \{`p_earnings_nap'\} & \{`p_prod_nap'\} & \{`p_typing_nap'\} & \{`p_output_nap'\}  & \{`p_wellbeing_nap'\}  & \{`p_physical_nap'\} &  \{`p_psych_nap'\} \\" _n ///
			"& & [`pval_w_mht_nap_earnings'] & [`pval_w_mht_nap_prod'] & [`pval_w_mht_nap_typing'] & [`pval_w_mht_nap_output']  & [`pval_w_mht_nap_wellbeing'] & [`pval_w_mht_nap_phys'] & [`pval_w_mht_nap_psych'] \\" _n ///
			"\addlinespace" _n ///
			///
			"\midrule" _n ///
			"Participants     & `pids_overall_napnonap' & `pids_earnings_ns' & `pids_prod_ns' & `pids_typing_ns' & `pids_output_ns' & `pids_wellbeing' & `pids_physical' & `pids_psych' \\" _n ///
			///
			"Unadjusted \textit{p}-value NS vs. Nap & \{`p_overall_d_nap_ns'\} & \{`p_earnings_d_nap_ns'\} & \{`p_prod_d_nap_ns'\} & \{`p_typing_d_nap_ns'\} & \{`p_output_d_nap_ns'\} & \{`p_wellbeing_d_nap_ns'\} & \{`p_physical_d_nap_ns'\} & \{`p_psych_d_nap_ns'\} \\ " _n ///
			"FWER-corrected \textit{p}-value NS vs. Nap & & [`pval_w_mht_d_earnings'] & [`pval_w_mht_d_prod'] & [`pval_w_mht_d_typing'] & [`pval_w_mht_d_output']  & [`pval_w_mht_d_wellbeing'] & [`pval_w_mht_d_phys'] & [`pval_w_mht_d_psych'] \\" _n ///
			"\addlinespace" _n ///
			///
			"\end{tabular*}" _n 
			file close f
			
			file open f using "Table4B_table_main_anderson_join2.tex", write replace
		file write f "\begin{tabular*}{1.7\textwidth}{l@{\extracolsep{\fill}}*{7}{c}}" _n ///+
			"& \multicolumn{3}{c}{\textbf{COGNITION}} & \multicolumn{4}{c}{\textbf{PREFERENCES}}\\" _n ///
			///
			"\cmidrule(lr){2-4}\cmidrule(lr){5-8} & \textbf{Index} & Lab Tasks & Work Task & \textbf{Index} & Time & Social & Risk \\" _n ///
			///
			" & (9) & (10) & (11) & (12) & (13) & (14) & (15) \\" _n ///
			"\midrule" _n ///
			///
			"\textbf{Night-Sleep Treatments} & \textbf{`coef_cogindex_ns'} & `coef_cogfunction_ns' & `coef_attention_ns' & \textbf{`coef_pref_ns'} & `coef_timepref_ns' & `coef_social_ns' & `coef_risk_ns' \\" _n ///
			"& (`se_cogindex_ns') & (`se_cogfunction_ns') & (`se_attention_ns') & (`se_pref_ns') & (`se_timepref_ns') & (`se_social_ns') & (`se_risk_ns') \\" _n ///
			"& \{`p_cogindex_ns'\} & \{`p_cogfunction_ns'\} & \{`p_attention_ns'\} & \{`p_pref_ns'\} & \{`p_timepref_ns'\} & \{`p_social_ns'\} & \{`p_risk_ns'\} \\" _n ///
			"& [`pval_w_mht_ns_cogindex'] & [`pval_w_mht_ns_cog'] & [`pval_w_mht_ns_salience'] & [`pval_w_mht_ns_pref'] & [`pval_w_mht_ns_time'] & [`pval_w_mht_ns_social'] & [`pval_w_mht_ns_risk']\\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Nap Treatment} & \textbf{`coef_cogindex_nap'} & `coef_cogfunction_nap' & `coef_attention_nap' & \textbf{`coef_pref_nap'} & `coef_timepref_nap' & `coef_social_nap' & `coef_risk_nap' \\" _n ///
			"& (`se_cogindex_nap') & (`se_cogfunction_nap') & (`se_attention_nap') & (`se_pref_nap') & (`se_timepref_nap') & (`se_social_nap') & (`se_risk_nap') \\" _n ///
			"& \{`p_cogindex_nap'\} & \{`p_cogfunction_nap'\} & \{`p_attention_nap'\} & \{`p_pref_nap'\} & \{`p_timepref_nap'\} & \{`p_social_nap'\} & \{`p_risk_nap'\} \\" _n ///
			"& [`pval_w_mht_nap_cogindex'] & [`pval_w_mht_nap_cog'] & [`pval_w_mht_nap_salience'] & [`pval_w_mht_nap_pref'] & [`pval_w_mht_nap_time'] & [`pval_w_mht_nap_social'] & [`pval_w_mht_nap_risk'] \\" _n ///
			"\addlinespace" _n ///
			///
			"\midrule" _n ///
			"Participants     & `pids_cogindex' & `pids_cogfunction' & `pids_attention' & `pids_pref' & `pids_timepref' & `pids_social' & `pids_risk' \\" _n ///
			///
			"Unadjusted \textit{p}-value NS vs. Nap & \{`p_cogindex_d_nap_ns'\} & \{`p_cogfunction_d_nap_ns'\} & \{`p_attention_d_nap_ns'\} & \{`p_pref_d_nap_ns'\} & \{`p_timepref_d_nap_ns'\} & \{`p_social_d_nap_ns'\} & \{`p_risk_d_nap_ns'\} \\ " _n ///
			"FWER-corrected \textit{p}-value NS vs. Nap & [`pval_w_mht_d_cogindex'] & [`pval_w_mht_d_cog'] & [`pval_w_mht_d_salience'] & [`pval_w_mht_d_pref'] & [`pval_w_mht_d_time'] & [`pval_w_mht_d_social'] & [`pval_w_mht_d_risk'] \\" _n ///
			"\addlinespace" _n ///	
			///
			"\bottomrule" _n ///
			"\end{tabular*}" _n 
		file close f
		
		
