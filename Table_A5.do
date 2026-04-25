************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates Appendix Table 5 (Heterogeneous Treatment Effects for Sleep Outcome)
* Last edited: 07 May 2021
************************************************

clear all
set more off
set seed 1

	do "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Scripts\anderson_index_program.do"
	
	cap reghdfe, compile
	
**************
* Data setup *
**************

	use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\heterogeneity_dataset.dta", clear
	
	merge m:1 pid post_treatment using "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\anderson_indices.dta"
	drop _merge 
	
	* Add baseline naps variable 
	
	preserve 
	clear all 

	use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\baseline_cleaned.dta"
	keep pid c27 
	tempfile baseline_naps
	save `baseline_naps'

	restore 

	merge m:1 pid using "`baseline_naps'"
	keep if _merge==3
	drop _merge

	merge 1:1 pid day_in_study using "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\analysis_base.dta"
	keep if _merge==3
	drop _merge

	
	rename prod mean_prod
	rename typing mean_typing
	rename productivity prod 
	rename typing_time_hr typing
	rename earnings earnings
	
	rename ( sleep_night sleep_eff nap_time_mins c27 above_med_quality) ( ns_duration ns_quality nap_duration num_nap_bas above_median_quality)
	rename (cogindex_pre cogindex_post pref_pre  pref_post) (cognitive_index_base cognitive_index pref_index_base  pref_index)
	
	
	* Generate binary = 1 if the person naps 1 time or more at baseline
	gen nap_bas = 0 
	replace nap_bas = 1 if num_nap_bas!=0 
	
	*Above median age
	gen above_median_age = (age >= 3)
	
	* Fix inconsistencies in treatment assignment 
	egen treat_diff = mean(treat_pool), by (pid)
	replace treat_pool=1 if treat_diff!=0 
	
*************************
* NS Duration & Quality *
*************************

	foreach var in ns_duration ns_quality{ 

		*baseline
		cap drop _baseline
		cap drop baseline // necessary to divide in two lines because "baseline" may be ambiguous (and therefore drop will not work)
		gen _baseline	= `var' if post_treatment==0
		egen baseline 	= mean(_baseline), by(pid)
			
		local controls "female i.age"
					
		*benchmark regressions, no interaction
		reghdfe `var' treat_pool baseline `controls' if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
		
		local pids_`var' = string(e(N_clust))
		local obs_`var' = string(e(N))
		
		lincom _b[treat_pool]
		local coef_ns_`var'_1 = string(r(estimate),"%3.2f")
		local se_ns_`var'_1 = string(r(se),"%3.2f")
		local p_ns_`var'_1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
		
		*empty locals just to make loop run
		local coef_above_`var'_1 = 0.1
		local se_above_`var'_1 = 0.1
		local p_above_`var'_1 = 0.1
		
		local coef_int_`var'_1 = 0.1
		local se_int_`var'_1 = 0.1
		local p_int_`var'_1 = 0.1
		
		sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment == 1
		local mean_`var'_1 = string(r(mean), "%3.2f")
		local sd_`var'_1 = string(r(sd), "%3.2f")
			
		*interactions
		local i = 1
		foreach int in length eff nap_bas female{
			
			local i = `i' + 1
			
			*regressions
			if "`int'" == "nap_bas"{
				reghdfe `var' treat_pool baseline `controls' nap_bas 1.treat_pool#1.nap_bas if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			else if "`int'" == "female"{
				reghdfe `var' treat_pool baseline `controls' 1.treat_pool#1.female if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			else {
				reghdfe `var' treat_pool baseline `controls' i.above_median_`int' 1.treat_pool#1.above_median_`int' if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			
			*coefficient, se and p-val
				lincom _b[treat_pool]
				local coef_ns_`var'_`i' = string(r(estimate),"%3.2f")
				local se_ns_`var'_`i' = string(r(se),"%3.2f")
				local p_ns_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			if "`int'" == "female"{
				lincom _b[female]
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_pool#1.female]
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			
			else if "`int'" == "nap_bas"{
				lincom _b[nap_bas]
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_pool#1.nap_bas]
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			
			else{
				lincom _b[1.above_median_`int']
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_pool#1.above_median_`int']
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			}
		
		
		* Significance 

			foreach i in ns above int{
			forval j = 1/5 {
					if `p_`i'_`var'_`j'' < 0.01{
						 local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''***"
					}
					if `p_`i'_`var'_`j'' < 0.05 & `p_`i'_`var'_`j'' >= 0.01{
						 local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''**"
					}
					if `p_`i'_`var'_`j'' < 0.1 & `p_`i'_`var'_`j'' >= 0.05{
						 local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''*"
					} 
					
				}
				}
				}
				
**********************************************				
* First Stage Heterogeneity for Nap Duration *
**********************************************

foreach var in nap_duration {
			
		local controls "female i.age"
		
		* removed fraction_high long_day
					
		*benchmark regressions, no interaction
		reghdfe `var' treat_nap  `controls' if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
		
		local pids_`var' = string(e(N_clust))
		local obs_`var' = string(e(N))
		
		lincom _b[treat_nap]
		local coef_nap_`var'_1 = string(r(estimate),"%3.2f")
		local se_nap_`var'_1 = string(r(se),"%3.2f")
		local p_nap_`var'_1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
		local coef_int_`var'_1 = 0
		local se_int_`var'_1 = 0
		local p_int_`var'_1 = 0
		
		sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment == 1
		local mean_`var'_1 = string(r(mean), "%3.2f")
		local sd_`var'_1 = string(r(sd), "%3.2f")
		
		*interactions
		local i = 1
		foreach int in length eff nap_bas female{
			
			local i = `i' + 1
			
			*regressions
			if "`int'" == "nap_bas"{
				reghdfe `var' treat_nap  `controls' nap_bas 1.treat_nap#1.nap_bas if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			else if "`int'" == "female"{
				reghdfe `var' treat_nap  `controls' 1.treat_nap#1.female if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			else {
				reghdfe `var' treat_nap  `controls' i.above_median_`int' 1.treat_nap#1.above_median_`int' if post_treatment==1 & at_present_check==1, vce(cluster pid) absorb(date day_in_study)
			}
			
			*coefficient, se and p-val
			lincom _b[treat_nap]
			local coef_nap_`var'_`i' = string(r(estimate),"%3.2f")
			local se_nap_`var'_`i' = string(r(se),"%3.2f")
			local p_nap_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			if "`int'" == "female"{
				lincom _b[female]
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_nap#1.female]
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			
			else if "`int'" == "nap_bas"{
				lincom _b[nap_bas]
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_nap#1.nap_bas]
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			
			else{
				lincom _b[1.above_median_`int']
				local coef_above_`var'_`i' = string(r(estimate),"%3.2f")
				local se_above_`var'_`i' = string(r(se),"%3.2f")
				local p_above_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
				lincom _b[1.treat_nap#1.above_median_`int']
				local coef_int_`var'_`i' = string(r(estimate),"%3.2f")
				local se_int_`var'_`i' = string(r(se),"%3.2f")
				local p_int_`var'_`i' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			}
			
			
			*control mean and sd
			if "`int'" == "female"{
				sum `var' if treat_pool == 0 & treat_nap == 0 & female == 0 & post_treatment == 1
				local mean_`var'_`i' = string(r(mean), "%3.2f")
				local sd_`var'_`i' = string(r(sd), "%3.2f")
			}
			
			else if "`int'" == "nap_bas"{
				sum `var' if treat_pool == 0 & treat_nap == 0 & nap_bas == 0 & post_treatment == 1
				local mean_`var'_`i' = string(r(mean), "%3.2f")
				local sd_`var'_`i' = string(r(sd), "%3.2f")
			}
			
			else {
				sum `var' if treat_pool == 0 & treat_nap == 0 & above_median_`int' == 0 & post_treatment == 1
				local mean_`var'_`i' = string(r(mean), "%3.2f")
				local sd_`var'_`i' = string(r(sd), "%3.2f")
			}
		}
		
		* Significance 
		
			foreach i in nap int {
				forval j = 1/5{
					if `p_`i'_`var'_`j'' < 0.01{
						local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''***"
					}
					if `p_`i'_`var'_`j'' < 0.05 & `p_`i'_`var'_`j'' >= 0.01{
						local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''**"
					}
					if `p_`i'_`var'_`j'' < 0.1 & `p_`i'_`var'_`j'' >= 0.05{
						local coef_`i'_`var'_`j' = "`coef_`i'_`var'_`j''*"
					} 
					
				}
				}
				}

				
**** Exporting ****
cd "C:\Users\Harshit\Downloads"
	
	file open f using "TableA5_main_firststage_hte.tex", write replace
	file write f "\begin{tabular}{@{}c@{}}" _n ///
	"\begin{tabular}{l*{15}{c}}" _n ///
	"\toprule" _n ///
	"& \multicolumn{5}{c}{\textbf{Night Sleep Duration in Hours}} & \multicolumn{5}{c}{\textbf{Night Sleep Efficiency in \% }} & 	\multicolumn{5}{c}{\textbf{Nap Sleep Duration in Minutes}}\\" _n ///
	///
	"\cmidrule(lr){2-6}\cmidrule(lr){7-11}\cmidrule(lr){12-16} & & \makecell{X=Sleep \\ Duration} & \makecell{X=Sleep \\ Efficiency} & \makecell{X=Baseline \\ Naps} & \makecell{X=Female} & & \makecell{X=Sleep \\ Duration} & \makecell{X=Sleep \\ Efficiency} & \makecell{X=Baseline \\ Naps} & \makecell{X=Female} & & \makecell{X=Sleep \\ Duration} & \makecell{X=Sleep \\ Efficiency} & \makecell{X=Baseline \\ Naps} & \makecell{X=Female} \\" _n ///
	///
	"& (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15) \\" _n ///
	///
	"\midrule" _n ///
	/// 
	"Night-Sleep Treatments & `coef_ns_ns_duration_1' & `coef_ns_ns_duration_2' & `coef_ns_ns_duration_3' & `coef_ns_ns_duration_4' & `coef_ns_ns_duration_5' & `coef_ns_ns_quality_1' & `coef_ns_ns_quality_2' & `coef_ns_ns_quality_3' & `coef_ns_ns_quality_4' & `coef_ns_ns_quality_5' & & & & \\" _n /// 
	"& (`se_ns_ns_duration_1') & (`se_ns_ns_duration_2') & (`se_ns_ns_duration_3') & (`se_ns_ns_duration_4') & (`se_ns_ns_duration_5') & (`se_ns_ns_quality_1') & (`se_ns_ns_quality_2') & (`se_ns_ns_quality_3') & (`se_ns_ns_quality_4') & (`se_ns_ns_quality_5') & & & & \\" _n /// 
	"\addlinespace" _n ///
	///
	"X & & `coef_above_ns_duration_2' & `coef_above_ns_duration_3' & `coef_above_ns_duration_4' & `coef_above_ns_duration_5' & & `coef_above_ns_quality_2' & `coef_above_ns_quality_3' & `coef_above_ns_quality_4' & `coef_above_ns_quality_5' & & `coef_above_nap_duration_2' & `coef_above_nap_duration_3' & `coef_above_nap_duration_4' & `coef_above_nap_duration_5' \\" _n /// 
	"& & (`se_above_ns_duration_2') & (`se_above_ns_duration_3') & (`se_above_ns_duration_4') & (`se_above_ns_duration_5') & & (`se_above_ns_quality_2') & (`se_above_ns_quality_3') & (`se_above_ns_quality_4') & (`se_above_ns_quality_5') & & (`se_above_nap_duration_2') & (`se_above_nap_duration_3') & (`se_above_nap_duration_4') & (`se_above_nap_duration_5') \\" _n /// 
	"\addlinespace" _n ///
	///
	"Night-Sleep Treatments*X & & `coef_int_ns_duration_2' & `coef_int_ns_duration_3' & `coef_int_ns_duration_4' & `coef_int_ns_duration_5' & & `coef_int_ns_quality_2' & `coef_int_ns_quality_3' & `coef_int_ns_quality_4' & `coef_int_ns_quality_5' & & & & \\" _n /// 
	"& & (`se_int_ns_duration_2') & (`se_int_ns_duration_3') & (`se_int_ns_duration_4') & (`se_int_ns_duration_5') & & (`se_int_ns_quality_2') & (`se_int_ns_quality_3') & (`se_int_ns_quality_4') & (`se_int_ns_quality_5') & & & & \\" _n /// 
	"\addlinespace" _n ///
	///
	"Nap Treatment & & & & & & & & & & & `coef_nap_nap_duration_1' & `coef_nap_nap_duration_2' & `coef_nap_nap_duration_3' & `coef_nap_nap_duration_4' & `coef_nap_nap_duration_5'\\" _n /// 
	"& & & & & & & & & & & (`se_nap_nap_duration_1') & (`se_nap_nap_duration_2') & (`se_nap_nap_duration_3') & (`se_nap_nap_duration_4') & (`se_nap_nap_duration_5')\\" _n /// 
	"\addlinespace" _n ///
	///
	"Nap Treatment*X & & & & & & & & & & & & `coef_int_nap_duration_2' & `coef_int_nap_duration_3' & `coef_int_nap_duration_4' & `coef_int_nap_duration_5' \\" _n /// 
	"& & & & & & & & & & & & (`se_int_nap_duration_2') & (`se_int_nap_duration_3') & (`se_int_nap_duration_4') & (`se_int_nap_duration_5')\\" _n /// 
	"\addlinespace" _n ///
	///
	"\midrule" _n ///
	" DV Control Group Mean & `mean_ns_duration_1' & `mean_ns_duration_1' & `mean_ns_duration_1' & `mean_ns_duration_1' & `mean_ns_duration_1' & `mean_ns_quality_1' & `mean_ns_quality_1' & `mean_ns_quality_1' & `mean_ns_quality_1' & `mean_ns_quality_1' &&&& \\"  _n ///
	" Participants & `pids_ns_duration'& `pids_ns_duration' & `pids_ns_duration' & `pids_ns_duration' & `pids_ns_duration' & `pids_ns_quality' & `pids_ns_quality' & `pids_ns_quality' & `pids_ns_quality' & `pids_ns_quality' & `pids_nap_duration' & `pids_nap_duration' & `pids_nap_duration' & `pids_nap_duration' & `pids_nap_duration' \\"  _n ///
	///	
	"\bottomrule" _n ///
	"\end{tabular}" _n ///
	"\end{tabular}" _n
	file close f		
