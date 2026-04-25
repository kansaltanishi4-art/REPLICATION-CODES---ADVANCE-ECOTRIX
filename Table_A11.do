************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates Appendix Table 11 (Treatment Effects on Psychological and Physical Well-being)
* Last edited: 07 May 2021
************************************************

	clear all
	set more off
	set seed 1

	do "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Scripts\anderson_index_program.do"

	cap reghdfe, compile

use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\health_dataset.dta", clear

***
* Psych Wellbeing
***
	* note: code will only work properly if you run it from preserve to restore in a single run
	preserve
	
	*REMEMBER: ALWAYS MAKE SURE SIGNS ARE RIGHT
	replace es_dep = -es_dep
	replace ds_g13c_stressed = -ds_g13c_stressed

	*standardize
	rename ds_g2_life_possibility life_possibility
	foreach var in "ds_a32_happiness" "life_possibility" "ds_g1_satisfaction" "ds_g13c_stressed" "es_dep" {				
			sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
			gen `var'_post_std = (`var' - r(mean)) / r(sd)

			bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
			bysort pid: egen `var'_pre = mean(`var'_pre_temp)
			sum `var'_pre if treat_pool == 0 & treat_nap == 0 & post_treatment==0
			gen `var'_pre_std = (`var'_pre - r(mean)) /r(sd)
			drop `var'_pre `var'_pre_temp
			}

	tempfile psych
	save `psych', replace

	restore

***
* Phys Wellbeing
***
	
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

	*Standardize final vars
	foreach var in "biking_task" "win_bp" "es_ill_week" "es_pain" "es_daily_act" {				
		sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
		gen `var'_post_std = (`var' - r(mean)) / r(sd)

		bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
		bysort pid: egen `var'_pre = mean(`var'_pre_temp)
		sum `var'_pre if treat_pool == 0 & treat_nap == 0 & post_treatment==0
		gen `var'_pre_std = (`var'_pre - r(mean)) /r(sd)
		drop `var'_pre `var'_pre_temp
	}
	
	gen composite_health_baseline = (win_bp_pre_std + es_ill_week_pre_std + es_pain_pre_std + es_daily_act_pre_std) / 4
	
	tempfile phys
	save `phys'

***
* Regression - wellbeing
***
			
merge 1:1 pid day_in_study using `psych'

		*Components - day level
		local components "ds_a32_happiness life_possibility ds_g1_satisfaction ds_g13c_stressed win_bp" 
		
		foreach var in `components' {
			
			reghdfe `var'_post_std treat_nap treat_pool female `var'_pre_std if post_treatment==1, cluster(pid) absorb(age date day_in_study)

				//Key sleep coefs
				lincom _b[treat_pool]
				local coef_ns_`var' = string(r(estimate),"%3.2f")
				local se_ns_`var' = string(r(se),"%3.2f")
				local p_ns_`var' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				//Key nap coefs
				lincom _b[treat_nap]
				local coef_nap_`var' = string(r(estimate),"%3.2f")
				local se_nap_`var' = string(r(se),"%3.2f")
				local p_nap_`var' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			local obs_`var' = string(e(N))
			distinct(pid) if e(sample) == 1
			local pids_`var' = r(ndistinct)
		
			}
			
		*Components - pid level
		
		drop biking_task_pre_std
		rename composite_health_baseline biking_task_pre_std
		
		collapse (mean) es_dep* biking_task* es_ill_week* es_pain* es_daily_act* (max) female age treat_pool treat_nap, by(pid post_treatment)
		
		local components "es_dep biking_task es_ill_week es_pain es_daily_act" 
		
		foreach var in `components' {
			
			reghdfe `var'_post_std treat_nap treat_pool female `var'_pre_std if post_treatment==1, cluster(pid) absorb(age)

				//Key sleep coefs
				lincom _b[treat_pool]
				local coef_ns_`var' = string(r(estimate),"%3.2f")
				local se_ns_`var' = string(r(se),"%3.2f")
				local p_ns_`var' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				//Key nap coefs
				lincom _b[treat_nap]
				local coef_nap_`var' = string(r(estimate),"%3.2f")
				local se_nap_`var' = string(r(se),"%3.2f")
				local p_nap_`var' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			local obs_`var' = string(e(N))
			distinct(pid) if e(sample) == 1
			local pids_`var' = r(ndistinct)
		
			}
			
		local components "ds_a32_happiness life_possibility ds_g1_satisfaction ds_g13c_stressed es_dep biking_task win_bp es_ill_week es_pain es_daily_act" 
		foreach var in `components' {
			foreach i in ns nap {
				 if `p_`i'_`var'' < 0.01 {
					local coef_`i'_`var' = "`coef_`i'_`var''***"
				}
				
				if `p_`i'_`var'' < 0.05 & `p_`i'_`var'' >= 0.01 {
					local coef_`i'_`var' = "`coef_`i'_`var''**"
				}
				
				 if `p_`i'_`var'' < 0.1 & `p_`i'_`var'' >= 0.05 {
					local coef_`i'_`var' = "`coef_`i'_`var''*"
				}
			}

		}
	
	** Writing table **
		cd "C:\Users\Harshit\Downloads"
	file open f using "TableA11_wellbeing.tex", write replace
		file write f "\begin{tabular*}{1.15\textwidth}{l@{\extracolsep{\fill}}*{5}{c}}" _n ///
		"\toprule" _n ///
		"&\multicolumn{5}{c}{\textbf{Panel A: Standardized Psychological Well-being Components}} \\" _n ///
		"\cmidrule{2-6}" _n ///
"&\multicolumn{1}{c}{\Shortunderstack{Depression\\(1)}}&\multicolumn{1}{c}{\Shortunderstack{Happiness\\(2)}}&\multicolumn{1}{c}{\Shortunderstack{Life Possibility\\(3)}}&\multicolumn{1}{c}{\Shortunderstack{Life Satisfaction\\(4)}}&\multicolumn{1}{c}{\Shortunderstack{Stress\\(5)}} \\" _n ///
		"\midrule" _n ///
		"Night-Sleep Treatments & `coef_ns_es_dep' & `coef_ns_ds_a32_happiness' & `coef_ns_life_possibility' & `coef_ns_ds_g1_satisfaction' & `coef_ns_ds_g13c_stressed' \\" _n ///
		"& (`se_ns_es_dep') & (`se_ns_ds_a32_happiness') & (`se_ns_life_possibility') & (`se_ns_ds_g1_satisfaction') & (`se_ns_ds_g13c_stressed') \\" _n ///
		"Nap Treatment & `coef_nap_es_dep' & `coef_nap_ds_a32_happiness' & `coef_nap_life_possibility' & `coef_nap_ds_g1_satisfaction' & `coef_nap_ds_g13c_stressed' \\" _n ///
		"& (`se_nap_es_dep') & (`se_nap_ds_a32_happiness') & (`se_nap_life_possibility') & (`se_nap_ds_g1_satisfaction') & (`se_nap_ds_g13c_stressed') \\" _n ///
		"\midrule" _n ///
		"Participants &   `pids_es_dep' &   `pids_ds_a32_happiness' &   `pids_life_possibility' &   `pids_ds_g1_satisfaction' &   `pids_ds_g13c_stressed' \\" _n ///
		///
			"\bigskip \\" _n ///
		///
"&\multicolumn{5}{c}{\textbf{Panel B: Standardized Physical Well-being Components}} \\" _n ///
		"\cmidrule{2-6}" _n ///
"&\multicolumn{1}{c}{\Shortunderstack{Biking\\(6)}}&\multicolumn{1}{c}{\Shortunderstack{Illness\\(7)}}&\multicolumn{1}{c}{\Shortunderstack{Pain\\(8)}}&\multicolumn{1}{c}{\Shortunderstack{Daily Act.\\(9)}}&\multicolumn{1}{c}{\Shortunderstack{BP\\(10)}} \\" _n ///
		"\midrule" _n ///
		"Night-Sleep Treatments & `coef_ns_biking_task' & `coef_ns_es_ill_week' & `coef_ns_es_pain' & `coef_ns_es_daily_act' & `coef_ns_win_bp' \\" _n ///
		" & (`se_ns_biking_task') & (`se_ns_es_ill_week') & (`se_ns_es_pain') & (`se_ns_es_daily_act') & (`se_ns_win_bp') \\" _n ///
		"Nap Treatment &`coef_nap_biking_task' & `coef_nap_es_ill_week' & `coef_nap_es_pain' & `coef_nap_es_daily_act' & `coef_nap_win_bp' \\" _n ///
		"& (`se_nap_biking_task') & (`se_nap_es_ill_week') & (`se_nap_es_pain') & (`se_nap_es_daily_act') & (`se_nap_win_bp') \\" _n ///
		"\midrule" _n ///
		"Participants  &   `pids_biking_task'  &   `pids_es_ill_week' &   `pids_es_pain' &   `pids_es_daily_act' &   `pids_win_bp'  \\" _n ///
		"\bottomrule" _n ///
		"\end{tabular*}" _n
	file close f
