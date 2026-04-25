************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates main paper, Table 3 (Fully-Disaggregated Treatment Effects of Night-Sleep and Nap Treatments)
* Last edited: 07 May 2021
************************************************
		clear all
		set more off
		set seed 1

		//Initiate anderson indices and fixed effects programs
		do "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Scripts\anderson_index_program.do""
		cap reghdfe, compile

	*************************
	* LABOR SUPPLY OUTCOMES *
	*************************

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\typing_dataset.dta", clear

		egen _treat_nap = max(treat_nap), by(pid)
		drop treat_nap
		rename _treat_nap treat_nap

		drop if day_in_study==.
		drop if day_in_study==1

		rename productivity prod
		rename typing_time_hr typing

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i
		
		
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

			*Pooled treatment
			reghdfe `var'_post_std treat_pool_cell treat_nap_cell treat_int_cell baseline i.age fraction_high female long_day if post_treatment==1 & at_present_check==1, absorb(date day_in_study) cluster(pid)

			local pids_`var'_ns = e(N_clust)
			local obs_`var'_ns = e(N)

			lincom _b[treat_pool_cell]
				local coef_`var'_ns = string(r(estimate),"%3.2f")
				local se_`var'_ns = string(r(se),"%3.2f")
				local tstat_`var'_ns = r(estimate)/r(se)
				local p_`var'_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_`var'_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_`var'_nap2 = string(r(estimate),"%3.2f")
				local se_`var'_nap2 = string(r(se),"%3.2f")
				local tstat_`var'_nap2 = r(estimate)/r(se)
				local p_`var'_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_`var'_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_`var'_int = string(r(estimate),"%3.2f")
				local se_`var'_int = string(r(se),"%3.2f")
				local tstat_`var'_int = r(estimate)/r(se)
				local p_`var'_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_`var'_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_`var'_d_nap_ns = string(r(p), "%3.2f")
				local p_`var'_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_`var'_d_int_ns = string(r(p), "%3.2f")
				local p_`var'_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_`var'_d_int_nap = string(r(p), "%3.2f")
				local p_`var'_d_int_nap1 = r(p)

				
		*Devices and Devices + Incentives
			reghdfe `var'_post_std treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell baseline i.age fraction_high female long_day if post_treatment==1 & at_present_check==1, absorb(date day_in_study) cluster(pid)

		foreach j in s s_i int_s int_s_i nap{

			lincom _b[treat_`j'_cell]
				local coef_`var'_`j' = string(r(estimate),"%3.2f")
				local se_`var'_`j' = string(r(se),"%3.2f")
				local tstat_`var'_`j' = r(estimate)/r(se)
				local p_`var'_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_`var'_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		}

	* Aggregate earnings at the participant level

	* Residualizing the typing variables
		describe earnings, varlist
		local type_vars `r(varlist)'
		foreach var in `type_vars ' {
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

		keep if at_present_check==1

		collapse (mean) prod earnings* typing output (max) treat* female age, by(pid post_treatment)

		keep if !mi(post_treatment)

		tempfile work
		save `work', replace
		

	**************
	* WELL-BEING *
	**************

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\health_dataset.dta", clear
		
	* Psych Wellbeing

	preserve

	* Residualize

		describe ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep, varlist
		local wb_vars `r(varlist)'
		foreach var in `wb_vars' {
			reghdfe `var'  , absorb(day_in_study date) residuals
			cap drop `var'
			rename _reghdfe_resid `var'
		}

	* Change signs
	
		replace es_dep = -es_dep
		replace ds_g13c_stressed = -ds_g13c_stressed
		
		collapse (mean) ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep (max) treat_pool treat_s treat_s_i treat_nap  female age, by(pid post_treatment)

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i
		
		* Weighted index (Anderson 2008)

		*Pooled treatments
			* Step 1: define locals for the relevant arguments

				*Components
				local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep"

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

				*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
				local factorcovariates "age"

				*Subset
				local subset "if post_treatment==1"

				*Eststo? Y = 1, N = 0
				local eststo 0

			* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
				* Signs changed already

			* Step 3: Just run this line!
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

				local pids_psych = e(N)
				local obs_psych = e(N)

				lincom _b[treat_pool_cell]
					local coef_psych_ns = string(r(estimate),"%3.2f")
					local se_psych_ns = string(r(se),"%3.2f")
					local tstat_psych_ns = r(estimate)/r(se)
					local p_psych_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_psych_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
					
				lincom _b[treat_nap_cell]
					local coef_psych_nap2 = string(r(estimate),"%3.2f")
					local se_psych_nap2 = string(r(se),"%3.2f")
					local tstat_psych_nap2 = r(estimate)/r(se)
					local p_psych_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_psych_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

				lincom _b[treat_int_cell]
					local coef_psych_int = string(r(estimate),"%3.2f")
					local se_psych_int = string(r(se),"%3.2f")
					local tstat_psych_int = r(estimate)/r(se)
					local p_psych_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_psych_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))


			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_psych_d_nap_ns = string(r(p), "%3.2f")
				local p_psych_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_psych_d_int_ns = string(r(p), "%3.2f")
				local p_psych_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_psych_d_int_nap = string(r(p), "%3.2f")
				local p_psych_d_int_nap1 = r(p)

		*Devices and Devices + Incentives
			* Step 1: define locals for the relevant arguments

				*Components
				local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep"

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

				*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
				local factorcovariates "age"

				*Subset
				local subset "if post_treatment==1"

				*Eststo? Y = 1, N = 0
				local eststo 0

			* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
				* Signs changed already

			* Step 3: Just run this line!
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

				local pids_psych = e(N)
				local obs_psych = e(N)

				foreach j in s s_i int_s int_s_i nap{
					lincom _b[treat_`j'_cell]
						local coef_psych_`j' = string(r(estimate),"%3.2f")
						local se_psych_`j' = string(r(se),"%3.2f")
						local tstat_psych_`j' = r(estimate)/r(se)
						local p_psych_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
						local p_psych_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				}

d
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

	*Create win_bp (summing blood pressure at pre-treatment)
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

		collapse (mean) biking_task win_bp es_ill_week es_pain es_daily_act (max) female age treat_pool treat_s treat_s_i treat_nap , by(pid post_treatment)
		
		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i

	* Weighted index (Anderson 2008)

	*Pooled treatments
		* Step 1: define locals for the relevant arguments

			*Components
			local components "biking_task win_bp es_ill_week es_pain es_daily_act"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_physical_ns = string(r(estimate),"%3.2f")
				local se_physical_ns = string(r(se),"%3.2f")
				local tstat_physical_ns = r(estimate)/r(se)
				local p_physical_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_physical_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_physical_nap2 = string(r(estimate),"%3.2f")
				local se_physical_nap2 = string(r(se),"%3.2f")
				local tstat_physical_nap2 = r(estimate)/r(se)
				local p_physical_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_physical_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_physical_int = string(r(estimate),"%3.2f")
				local se_physical_int = string(r(se),"%3.2f")
				local tstat_physical_int = r(estimate)/r(se)
				local p_physical_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_physical_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_physical_d_nap_ns = string(r(p), "%3.2f")
				local p_physical_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_physical_d_int_ns = string(r(p), "%3.2f")
				local p_physical_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_physical_d_int_nap = string(r(p), "%3.2f")
				local p_physical_d_int_nap1 = r(p)


	*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

			*Components
			local components "biking_task win_bp es_ill_week es_pain es_daily_act"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_physical_`j' = string(r(estimate),"%3.2f")
					local se_physical_`j' = string(r(se),"%3.2f")
					local tstat_physical_`j' = r(estimate)/r(se)
					local p_physical_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_physical_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
	*******************
	* WELLBEING INDEX *
	*******************

		merge 1:1 pid post_treatment using `psych'
		
		

	* Weighted index (Anderson 2008)

	*Pooled treatments
		* Step 1: define locals for the relevant arguments

			*Components
			local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_wellbeing_ns = string(r(estimate),"%3.2f")
				local se_wellbeing_ns = string(r(se),"%3.2f")
				local tstat_wellbeing_ns = r(estimate)/r(se)
				local p_wellbeing_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_wellbeing_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_wellbeing_nap2 = string(r(estimate),"%3.2f")
				local se_wellbeing_nap2 = string(r(se),"%3.2f")
				local tstat_wellbeing_nap2 = r(estimate)/r(se)
				local p_wellbeing_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_wellbeing_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_wellbeing_int = string(r(estimate),"%3.2f")
				local se_wellbeing_int = string(r(se),"%3.2f")
				local tstat_wellbeing_int = r(estimate)/r(se)
				local p_wellbeing_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_wellbeing_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_wellbeing_d_nap_ns = string(r(p), "%3.2f")
				local p_wellbeing_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_wellbeing_d_int_ns = string(r(p), "%3.2f")
				local p_wellbeing_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_wellbeing_d_int_nap = string(r(p), "%3.2f")
				local p_wellbeing_d_int_nap1 = r(p)

	
	*Devices and Devices + incentives
		* Step 1: define locals for the relevant arguments

			*Components
			local components "ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_wellbeing_`j' = string(r(estimate),"%3.2f")
					local se_wellbeing_`j' = string(r(se),"%3.2f")
					local tstat_wellbeing_`j' = r(estimate)/r(se)
					local p_wellbeing_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_wellbeing_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
		* Store Anderson index
		gen wellbeing_index       = summary_index
		gen wellbeing_index_base  = summary_index_base

		keep pid post_treatment ds_a32_happiness ds_g2_life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act wellbeing_index wellbeing_index_base treat_s_i treat_int_s treat_int_s_i

		tempfile wellbeing
		save `wellbeing'

	*******************************
	* RISK AND SOCIAL PREFERENCES *
	*******************************

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\riskandsocial_dataset.dta", clear

		collapse (mean) risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent (max) female age treat_pool treat_s treat_s_i treat_nap , by(pid post_treatment)

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i
		
	* Risk Preferences

	* Weighted index (Anderson 2008)

	* Pooled treatments
		* Step 1: define locals for the relevant arguments

				*Components
				local components "risk_switch_point loss_switch_point"

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_risk_ns = string(r(estimate),"%3.2f")
				local se_risk_ns = string(r(se),"%3.2f")
				local tstat_risk_ns = r(estimate)/r(se)
				local p_risk_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_risk_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_risk_nap2 = string(r(estimate),"%3.2f")
				local se_risk_nap2 = string(r(se),"%3.2f")
				local tstat_risk_nap2 = r(estimate)/r(se)
				local p_risk_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_risk_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_int_cell]
				local coef_risk_int = string(r(estimate),"%3.2f")
				local se_risk_int = string(r(se),"%3.2f")
				local tstat_risk_int = r(estimate)/r(se)
				local p_risk_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_risk_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_risk_d_nap_ns = string(r(p), "%3.2f")
				local p_risk_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_risk_d_int_ns = string(r(p), "%3.2f")
				local p_risk_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_risk_d_int_nap = string(r(p), "%3.2f")		
				local p_risk_d_int_nap1 = r(p)		

	*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

				*Components
				local components "risk_switch_point loss_switch_point"

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_risk_`j' = string(r(estimate),"%3.2f")
					local se_risk_`j' = string(r(se),"%3.2f")
					local tstat_risk_`j' = r(estimate)/r(se)
					local p_risk_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_risk_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
	* Social Preferences

	* Weighted index (Anderson 2008)

	*Pooled treatments
		* Step 1: define locals for the relevant arguments

			*Components
			local components "d_send u_send t_send u_avg_receive t_avg_amountsent"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_social_ns = string(r(estimate),"%3.2f")
				local se_social_ns = string(r(se),"%3.2f")
				local tstat_social_ns = r(estimate)/r(se)
				local p_social_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_social_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap_cell]
				local coef_social_nap2 = string(r(estimate),"%3.2f")
				local se_social_nap2 = string(r(se),"%3.2f")
				local tstat_social_nap2 = r(estimate)/r(se)
				local p_social_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_social_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_social_int = string(r(estimate),"%3.2f")
				local se_social_int = string(r(se),"%3.2f")
				local tstat_social_int = r(estimate)/r(se)
				local p_social_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_social_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_social_d_nap_ns = string(r(p), "%3.2f")
				local p_social_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_social_d_int_ns = string(r(p), "%3.2f")
				local p_social_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_social_d_int_nap = string(r(p), "%3.2f")
				local p_social_d_int_nap1 = r(p)		
				
	*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

			*Components
			local components "d_send u_send t_send u_avg_receive t_avg_amountsent"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_social_`j' = string(r(estimate),"%3.2f")
					local se_social_`j' = string(r(se),"%3.2f")
					local tstat_social_`j' = r(estimate)/r(se)
					local p_social_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_social_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
		tempfile rs_pref
		save `rs_pref', replace

	********************
	* TIME PREFERENCES *
	********************

	* Present Bias

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\pb_dataset.dta", clear

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

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\savings_dataset.dta", clear

		egen _treat_nap = max(treat_nap), by(pid)
		drop treat_nap
		rename  _treat_nap treat_nap

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

		collapse deposits (max) female age treat_pool treat_s treat_s_i treat_nap, by(pid post_treatment)

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i
		
		merge m:1 pid using `pb', gen(pb_merge)

		gen beta = beta_post
		replace beta = beta_base if post_treatment == 0

	* Weighted index (Anderson 2008)

	*Pooled treatments
		* Step 1: define locals for the relevant arguments

			*Components
			local components "deposits beta"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "i.age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

		* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

			local pids_timepref = e(N_clust)
			local obs_timepref = e(N)

			lincom _b[treat_pool_cell]
				local coef_timepref_ns = string(r(estimate),"%3.2f")
				local se_timepref_ns = string(r(se),"%3.2f")
				local tstat_timepref_ns = r(estimate)/r(se)
				local p_timepref_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_timepref_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_timepref_nap2 = string(r(estimate),"%3.2f")
				local se_timepref_nap2 = string(r(se),"%3.2f")
				local tstat_timepref_nap2 = r(estimate)/r(se)
				local p_timepref_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_timepref_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_timepref_int = string(r(estimate),"%3.2f")
				local se_timepref_int = string(r(se),"%3.2f")
				local tstat_timepref_int = r(estimate)/r(se)
				local p_timepref_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_timepref_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_timepref_d_nap_ns = string(r(p), "%3.2f")
				local p_timepref_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_timepref_d_int_ns = string(r(p), "%3.2f")
				local p_timepref_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_timepref_d_int_nap = string(r(p), "%3.2f")		
				local p_timepref_d_int_nap1 = r(p)		
		
	*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

			*Components
			local components "deposits beta"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "i.age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGH

		* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

			local pids_timepref = e(N_clust)
			local obs_timepref = e(N)

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_timepref_`j' = string(r(estimate),"%3.2f")
					local se_timepref_`j' = string(r(se),"%3.2f")
					local tstat_timepref_`j' = r(estimate)/r(se)
					local p_timepref_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_timepref_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}

	*********************
	* PREFERENCES INDEX *
	*********************

		merge m:1 pid post_treatment using `rs_pref', gen(rs_merge)

	* Weighted (Anderson)

	*Pooled treatment
		* Step 1: define locals for the relevant arguments

				* Components
				local components "risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent deposits beta"

				* Covariates to be included linearly
				local controls "female dummy_missing"

				* Treatment dummies
				local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_pref_ns = string(r(estimate),"%3.2f")
				local se_pref_ns = string(r(se),"%3.2f")
				local tstat_pref_ns = r(estimate)/r(se)
				local p_pref_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_pref_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_pref_nap2 = string(r(estimate),"%3.2f")
				local se_pref_nap2 = string(r(se),"%3.2f")
				local tstat_pref_nap2 = r(estimate)/r(se)
				local p_pref_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_pref_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_pref_int = string(r(estimate),"%3.2f")
				local se_pref_int = string(r(se),"%3.2f")
				local tstat_pref_int = r(estimate)/r(se)
				local p_pref_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_pref_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_pref_d_nap_ns = string(r(p), "%3.2f")
				local p_pref_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_pref_d_int_ns = string(r(p), "%3.2f")
				local p_pref_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_pref_d_int_nap = string(r(p), "%3.2f")		
				local p_pref_d_int_nap1 = r(p)	

		*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

				* Components
				local components "risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent deposits beta"

				* Covariates to be included linearly
				local controls "female dummy_missing"

				* Treatment dummies
				local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_pref_`j' = string(r(estimate),"%3.2f")
					local se_pref_`j' = string(r(se),"%3.2f")
					local tstat_pref_`j' = r(estimate)/r(se)
					local p_pref_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_pref_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
	* Store Anderson index
		gen pref_post = summary_index
		gen pref_pre  = summary_index_base

		keep pid post_treatment dummy_missing risk_switch_point loss_switch_point d_send u_send t_send u_avg_receive t_avg_amountsent deposits beta pref_post pref_pre

		tempfile preferences
		save `preferences', replace

	**********************
	* COGNITIVE FUNCTION *
	**********************

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\cognitive_tasks_dataset.dta", clear

		merge 1:1 pid day_in_study using "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\pvt_dataset.dta"

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

		collapse pv_perf hf_payment co_payment (max) female age treat_nap treat_pool treat_s treat_s_i, by(pid post_treatment)

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i

	* Weighted index (Anderson 2008)

	* Pooled treatments
		* Step 1: define locals for the relevant arguments

				* Components
				local components "pv_perf hf_payment co_payment"

				* Covariates to be included linearly
				local controls "female"

				* Treatment dummies
				local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

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

			lincom _b[treat_pool_cell]
				local coef_cogfunction_ns = string(r(estimate),"%3.2f")
				local se_cogfunction_ns = string(r(se),"%3.2f")
				local tstat_cogfunction_ns = r(estimate)/r(se)
				local p_cogfunction_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogfunction_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_nap_cell]
				local coef_cogfunction_nap2 = string(r(estimate),"%3.2f")
				local se_cogfunction_nap2 = string(r(se),"%3.2f")
				local tstat_cogfunction_nap2 = r(estimate)/r(se)
				local p_cogfunction_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogfunction_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_cogfunction_int = string(r(estimate),"%3.2f")
				local se_cogfunction_int = string(r(se),"%3.2f")
				local tstat_cogfunction_int = r(estimate)/r(se)
				local p_cogfunction_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogfunction_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_cogfunction_d_nap_ns = string(r(p), "%3.2f")
				local p_cogfunction_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_cogfunction_d_int_ns = string(r(p), "%3.2f")
				local p_cogfunction_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_cogfunction_d_int_nap = string(r(p), "%3.2f")	
				local p_cogfunction_d_int_nap1 = r(p)
				
	*Devices and Devices + incentives
		* Step 1: define locals for the relevant arguments

				* Components
				local components "pv_perf hf_payment co_payment"

				* Covariates to be included linearly
				local controls "female"

				* Treatment dummies
				local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

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

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_cogfunction_`j' = string(r(estimate),"%3.2f")
					local se_cogfunction_`j' = string(r(se),"%3.2f")
					local tstat_cogfunction_`j' = r(estimate)/r(se)
					local p_cogfunction_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_cogfunction_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
		
	tempfile cognitive
	save `cognitive', replace

	*************
	* ATTENTION *
	*************

		use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\salience_dataset.dta", clear

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

	* Standardize variables

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

		collapse (mean) std_attent_output_res  baseline (max) treat_pool treat_s treat_s_i treat_nap age female, by(pid post_treatment)

		gen treat_int = treat_pool*treat_nap
		gen treat_int_s = treat_s*treat_nap
		gen treat_int_s_i = treat_s_i*treat_nap
		
		gen treat_pool_cell = treat_pool*(1-treat_nap)
		gen treat_s_cell = treat_s*(1-treat_nap)
		gen treat_s_i_cell = treat_s_i*(1-treat_nap)
		gen treat_nap_cell = treat_nap*(1-treat_pool)
		gen treat_int_cell = treat_int 
		gen treat_int_s_cell = treat_int_s
		gen treat_int_s_i_cell = treat_int_s_i
		
	*Pooled treatments
		reghdfe std_attent_output_res treat_pool_cell treat_nap_cell treat_int_cell baseline if post_treatment == 1, absorb(age female) vce(robust)

		local pids_attention = e(N)
		local obs_attention = e(N)

		lincom _b[treat_pool_cell]
			local coef_attention_ns = string(r(estimate),"%3.2f")
			local se_attention_ns = string(r(se),"%3.2f")
			local tstat_attention_ns = r(estimate)/r(se)
			local p_attention_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_attention_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_nap_cell]
			local coef_attention_nap2 = string(r(estimate),"%3.2f")
			local se_attention_nap2 = string(r(se),"%3.2f")
			local tstat_attention_nap2 = r(estimate)/r(se)
			local p_attention_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_attention_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

		lincom _b[treat_int_cell]
			local coef_attention_int = string(r(estimate),"%3.2f")
			local se_attention_int = string(r(se),"%3.2f")
			local tstat_attention_int = r(estimate)/r(se)
			local p_attention_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			local p_attention_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_attention_d_nap_ns = string(r(p), "%3.2f")
				local p_attention_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_attention_d_int_ns = string(r(p), "%3.2f")
				local p_attention_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_attention_d_int_nap = string(r(p), "%3.2f")	
				local p_attention_d_int_nap1 = r(p)
				
	*Devices and Devices + Incentives
		reghdfe std_attent_output_res treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell baseline if post_treatment == 1, absorb(age female) vce(robust)

		local pids_attention = e(N)
		local obs_attention = e(N)

		foreach j in s s_i int_s int_s_i nap{
			lincom _b[treat_`j'_cell]
				local coef_attention_`j' = string(r(estimate),"%3.2f")
				local se_attention_`j' = string(r(se),"%3.2f")
				local tstat_attention_`j' = r(estimate)/r(se)
				local p_attention_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_attention_`j'1 =  (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
		}

		rename baseline std_attent_output_base

	*******************
	* COGNITION INDEX *
	*******************

		merge 1:1 pid post_treatment using `cognitive', gen(cog_merge)

		replace std_attent_output_res = std_attent_output_base if post_treatment==0

	* Weighted index

	*Pooled treatments
		* Step 1: define locals for the relevant arguments

			*Components
			local components "pv_perf hf_payment co_payment std_attent_output_res"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

		* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

			foreach var in pv_perf hf_payment co_payment std_attent_output_res {

			sum `var'_post_std if post_treatment==1
			sum `var'_pre_std if post_treatment==0

			}

			foreach var in female treat_pool treat_nap age {

			sum `var' if post_treatment==1, de

			}

			local pids_cogindex = e(N_clust)
			local obs_cogindex = e(N)

			lincom _b[treat_pool_cell]
				local coef_cogindex_ns = string(r(estimate),"%3.2f")
				local se_cogindex_ns = string(r(se),"%3.2f")
				local tstat_cogindex_ns = r(estimate)/r(se)
				local p_cogindex_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogindex_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))


			lincom _b[treat_nap_cell]
				local coef_cogindex_nap2 = string(r(estimate),"%3.2f")
				local se_cogindex_nap2 = string(r(se),"%3.2f")
				local tstat_cogindex_nap2 = r(estimate)/r(se)
				local p_cogindex_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogindex_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			lincom _b[treat_int_cell]
				local coef_cogindex_int = string(r(estimate),"%3.2f")
				local se_cogindex_int = string(r(se),"%3.2f")
				local tstat_cogindex_int = r(estimate)/r(se)
				local p_cogindex_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_cogindex_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_cogindex_d_nap_ns = string(r(p), "%3.2f")
				local p_cogindex_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_cogindex_d_int_ns = string(r(p), "%3.2f")
				local p_cogindex_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_cogindex_d_int_nap = string(r(p), "%3.2f")	
				local p_cogindex_d_int_nap1 = r(p)
				
	*Devices and Devices + Incentives
		* Step 1: define locals for the relevant arguments

			*Components
			local components "pv_perf hf_payment co_payment std_attent_output_res"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

		* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

			local pids_cogindex = e(N_clust)
			local obs_cogindex = e(N)

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_cogindex_`j' = string(r(estimate),"%3.2f")
					local se_cogindex_`j' = string(r(se),"%3.2f")
					local tstat_cogindex_`j' = r(estimate)/r(se)
					local p_cogindex_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
					local p_cogindex_`j'1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}

	* Store Anderson index
		gen cogindex_post       = summary_index
		gen cogindex_pre  = summary_index_base

		rename std_attent_output_res std_attent_output
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
		
	* Correcting baseline values of indices

		replace cogindex_post = cogindex_pre if post_treatment==0
		replace pref_post = pref_pre if post_treatment==0
		replace wellbeing_index = wellbeing_index_base if post_treatment==0

		replace earnings_post_std = earnings_pre_std if post_treatment==0
		
		
	* Weighted (Anderson)
	eststo clear

		* Step 1: define locals for the relevant arguments

				*Covariates to be included linearly
				local controls "female"

				*Treatment dummies
				local treatments "treat_pool_cell treat_nap_cell treat_int_cell"

				*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
				local factorcovariates "i.age"

				*Eststo? Y = 1, N = 0
				local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
				* Signs have been changed already to make the average index!

		* Step 3: Just run this line!
				local subset "if post_treatment==1"

				local components "earnings_post_std pref_post wellbeing_index cogindex_post"
				anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" 1

			local pids_overall = e(N_clust)
			local obs_overall = e(N)

			lincom _b[treat_nap_cell]
				local coef_overall_nap2 = string(r(estimate),"%3.2f")
				local se_overall_nap2 = string(r(se),"%3.2f")
				local tstat_overall_nap2 = r(estimate)/r(se)
				local p_overall_nap2 = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_overall_nap21 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_pool_cell]
				local coef_overall_ns = string(r(estimate),"%3.2f")
				local se_overall_ns = string(r(se),"%3.2f")
				local tstat_overall_ns = r(estimate)/r(se)
				local p_overall_ns = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_overall_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_int_cell]
				local coef_overall_int = string(r(estimate),"%3.2f")
				local se_overall_int = string(r(se),"%3.2f")
				local tstat_overall_int = r(estimate)/r(se)
				local p_overall_int = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
				local p_overall_int1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			* Get pval of differences
			lincom _b[treat_nap_cell] - _b[treat_pool_cell]
				local p_overall_d_nap_ns = string(r(p), "%3.2f")
				local p_overall_d_nap_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_pool_cell]
				local p_overall_d_int_ns = string(r(p), "%3.2f")
				local p_overall_d_int_ns1 = r(p)
			lincom _b[treat_int_cell] - _b[treat_nap_cell]
				local p_overall_d_int_nap = string(r(p), "%3.2f")
				local p_overall_d_int_nap1 = r(p)
				
				
	*Devices and Devices + Incentives

		* Step 1: define locals for the relevant arguments

			*Components
			local components "earnings_post_std pref_post wellbeing_index cogindex_post"

			*Covariates to be included linearly
			local controls "female"

			*Treatment dummies
			local treatments "treat_s_cell treat_s_i_cell treat_nap_cell treat_int_s_cell treat_int_s_i_cell"

			*Covariates to be included inside "absorb" of reghdfe (if this is empty, add a factor var from controls)
			local factorcovariates "i.age"

			*Subset
			local subset "if post_treatment==1"

			*Eststo? Y = 1, N = 0
			local eststo 0

		* Step 2: REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT

		* Step 3: Just run this line!
			anderson_index "`components'" "`treatments'" "`controls'" "`factorcovariates'" "`subset'" `eststo'

			local pids_overall = e(N_clust)
			local obs_overall = e(N)

			foreach j in s s_i int_s int_s_i nap{
				lincom _b[treat_`j'_cell]
					local coef_overall_`j' = string(r(estimate),"%3.2f")
					local se_overall_`j' = string(r(se),"%3.2f")
					local tstat_overall_`j' = r(estimate)/r(se)
					local p_overall_`j' = string((2 * ttail(e(df_r), abs(r(estimate)/r(se)))), "%3.2f")
			}

	**********************
	* SIGNIFICANCE STARS *
	**********************

	* Anderson
		foreach var in overall prod typing earnings output wellbeing physical psych cogindex cogfunction attention pref timepref social risk {
			foreach i in ns s s_i nap nap2 int int_s int_s_i{
			di "`var'"
			di "`i'"
				 if `p_`var'_`i'' < 0.01 {
					local coef_`var'_`i' = "`coef_`var'_`i''***"
				}

				if `p_`var'_`i'' < 0.05 & `p_`var'_`i'' >= 0.01 {
					local coef_`var'_`i' = "`coef_`var'_`i''**"
				}

				 if `p_`var'_`i'' < 0.1 & `p_`var'_`i'' >= 0.05 {
					local coef_`var'_`i' = "`coef_`var'_`i''*"
				}
			}
		}
			

			 
**** Writing the table ****			
cd "C:\Users\Harshit\Downloads"
		file open f using "Table3A_table_main_int_separated_cell1_nomht.tex", write replace
			file write f "\begin{tabular*}{1.7\textwidth}{l@{\extracolsep{\fill}}*{8}{c}}" _n ///
			"\toprule" _n ///
			///
			"& \multicolumn{1}{c}{\textbf{OVERALL}} & \multicolumn{4}{c}{\textbf{WORK}} & \multicolumn{3}{c}{\textbf{WELL-BEING}}\\" _n ///
			///
			"\cmidrule(lr){2-2}\cmidrule(lr){3-6}\cmidrule(lr){7-9} & \textbf{Index} & \textbf{Earnings} & Productivity & Labor Supply & Output  & \textbf{Index} & Physical & Mental \\" _n ///
			///
			" & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\" _n ///
			"\midrule" _n ///
			///
			"\textbf{Devices+Encouragement Only} & `coef_overall_s' & `coef_earnings_s' & `coef_prod_s' & `coef_typing_s' & `coef_output_s'  & `coef_wellbeing_s' & `coef_physical_s' & `coef_psych_s' \\" _n ///
			"& (`se_overall_s') & (`se_earnings_s') & (`se_prod_s') & (`se_typing_s') & (`se_output_s')  & (`se_wellbeing_s') & (`se_physical_s') & (`se_psych_s') \\" _n /// 
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Incentives Only} & `coef_overall_s_i' & `coef_earnings_s_i' & `coef_prod_s_i' & `coef_typing_s_i' & `coef_output_s_i'  &  `coef_wellbeing_s_i' & `coef_physical_s_i' & `coef_psych_s_i' \\" _n ///
			" & (`se_overall_s_i') & (`se_earnings_s_i') & (`se_prod_s_i') & (`se_typing_s_i') & (`se_output_s_i')  & (`se_wellbeing_s_i') & (`se_physical_s_i') & (`se_psych_s_i') \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Nap Only} & `coef_overall_nap' & `coef_earnings_nap' & `coef_prod_nap' & `coef_typing_nap' & `coef_output_nap'  & `coef_wellbeing_nap' & `coef_physical_nap' & `coef_psych_nap' \\" _n ///
			"& (`se_overall_nap') & (`se_earnings_nap') & (`se_prod_nap') & (`se_typing_nap') & (`se_output_nap') & (`se_wellbeing_nap') & (`se_physical_nap') & (`se_psych_nap') \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Encouragement and Nap} & `coef_overall_int_s' & `coef_earnings_int_s' & `coef_prod_int_s' & `coef_typing_int_s' & `coef_output_int_s'  & `coef_wellbeing_int_s' & `coef_physical_int_s' & `coef_psych_int_s' \\" _n ///
			" & (`se_overall_int_s') & (`se_earnings_int_s') & (`se_prod_int_s') & (`se_typing_int_s') & (`se_output_int_s')  & (`se_wellbeing_int_s') & (`se_physical_int_s') & (`se_psych_int_s') \\" _n /// 
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Incentives and Nap} & `coef_overall_int_s_i' & `coef_earnings_int_s_i' & `coef_prod_int_s_i' & `coef_typing_int_s_i' & `coef_output_int_s_i'  & `coef_wellbeing_int_s_i' & `coef_physical_int_s_i' & `coef_psych_int_s_i' \\" _n /// 
			" & (`se_overall_int_s_i') & (`se_earnings_int_s_i') & (`se_prod_int_s_i') & (`se_typing_int_s_i') & (`se_output_int_s_i')  & (`se_wellbeing_int_s_i') & (`se_physical_int_s_i') & (`se_psych_int_s_i') \\" _n /// */
			"\addlinespace" _n ///
			///
			"\midrule" _n ///
			"Participants & `pids_overall' & `pids_earnings_ns'  & `pids_prod_ns' & `pids_typing_ns' & `pids_output_ns'  & `pids_wellbeing' & `pids_physical' & `pids_psych' \\" _n ///
			///
			"\end{tabular*}" _n
			file close f
			
		file open f using "Table3B_table_main_int_separated_cell2_nomht.tex", write replace
		file write f "\begin{tabular*}{1.7\textwidth}{l@{\extracolsep{\fill}}*{7}{c}}" _n ///+
			///
			"& \multicolumn{3}{c}{\textbf{COGNITION}} & \multicolumn{4}{c}{\textbf{PREFERENCES}} \\" _n ///
			///
			"\cmidrule(lr){2-4}\cmidrule(lr){5-8} & \textbf{Index} & Lab Tasks & Work Task & \textbf{Index} & Time & Social & Risk \\" _n ///
			///
			" & (9) & (10) & (11) & (12) & (13) & (14) & (15)\\" _n ///
			"\midrule" _n ///
			///
			"\textbf{Devices+Encouragement Only} & \textbf{`coef_cogindex_s'} & `coef_cogfunction_s' & `coef_attention_s' & \textbf{`coef_pref_s'} & `coef_timepref_s' & `coef_social_s' & `coef_risk_s' \\" _n ///
			"& (`se_cogindex_s') & (`se_cogfunction_s') & (`se_attention_s') & (`se_pref_s') & (`se_timepref_s') & (`se_social_s') & (`se_risk_s') \\" _n /// 
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Incentives Only} & `coef_cogindex_s_i' & `coef_cogfunction_s_i' & `coef_attention_s_i' & `coef_pref_s_i' & `coef_timepref_s_i' & `coef_social_s_i' & `coef_risk_s_i' \\" _n ///
			"& (`se_cogindex_s_i') & (`se_cogfunction_s_i') & (`se_attention_s_i') & (`se_pref_s_i') & (`se_timepref_s_i') & (`se_social_s_i') & (`se_risk_s_i') \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Nap Only} & `coef_cogindex_nap' & `coef_cogfunction_nap' & `coef_attention_nap' & `coef_pref_nap' & `coef_timepref_nap' & `coef_social_nap' & `coef_risk_nap' \\" _n ///
			"& (`se_cogindex_nap') & (`se_cogfunction_nap') & (`se_attention_nap') & (`se_pref_nap') & (`se_timepref_nap') & (`se_social_nap') & (`se_risk_nap') \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Encouragement and Nap} & `coef_cogindex_int_s' & `coef_cogfunction_int_s' & `coef_attention_int_s' & `coef_pref_int_s' & `coef_timepref_int_s' & `coef_social_int_s' & `coef_risk_int_s' \\" _n ///
			"& (`se_cogindex_int_s') & (`se_cogfunction_int_s') & (`se_attention_int_s') & (`se_pref_int_s') & (`se_timepref_int_s') & (`se_social_int_s') & (`se_risk_int_s') \\" _n ///
			"\addlinespace" _n ///
			///
			"\textbf{Devices+Incentives and Nap} & `coef_cogindex_int_s_i' & `coef_cogfunction_int_s_i' & `coef_attention_int_s_i' & `coef_pref_int_s_i' & `coef_timepref_int_s_i' & `coef_social_int_s_i' & `coef_risk_int_s_i' \\" _n ///
			"& (`se_cogindex_int_s_i') & (`se_cogfunction_int_s_i') & (`se_attention_int_s_i') & (`se_pref_int_s_i') & (`se_timepref_int_s_i') & (`se_social_int_s_i') & (`se_risk_int_s_i') \\" _n ///
			"\addlinespace" _n ///
			///
			"\midrule" _n ///
			"Participants     & `pids_cogindex' & `pids_cogfunction' & `pids_attention' & `pids_pref' & `pids_timepref' & `pids_social' & `pids_risk' \\" _n ///
			///
			"\bottomrule" _n ///
			"\end{tabular*}" _n
			file close f 
