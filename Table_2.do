************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates main paper, Table 2 (Treatment Effects on Sleep)
* Last edited: 07 May 2021
************************************************
use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\firststage_dataset.dta", clear

cap drop sleep_all
gen sleep_all= sleep_time
replace sleep_all = Sleep_Night if pid == 5591 //skipped a night in the data
replace sleep_all = Sleep_Night + nap_hr if pid == 5130 & day_in_study == 27 //raw actigraph data recorded nap_mins = half hour

	gen treat_int = treat_pool*treat_nap
	gen treat_int_s = treat_s*treat_nap
	gen treat_int_s_i = treat_s_i*treat_nap
	
	gen treat_pool_cell = treat_pool*(1-treat_nap)
	gen treat_s_cell = treat_s*(1-treat_nap)
	gen treat_s_i_cell = treat_s_i*(1-treat_nap)
	gen treat_np_cell = treat_nap*(1-treat_pool)
	gen treat_int_cell = treat_int 
	gen treat_int_s_cell = treat_int_s
	gen treat_int_s_i_cell = treat_int_s_i

	foreach type in "unpooled" "pooled" "cell" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}
		if "`type'" == "cell" {
		
			local treatments "treat_s_cell treat_s_i_cell treat_np_cell treat_int_s_cell treat_int_s_i_cell"
		}
		
		eststo clear
		
		preserve
		
		*Actigraph Sleep Time
		
			local controls "i.quart_age Female"

			capture drop _baseline baseline
			gen _baseline			= Sleep_Night if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
			
			reghdfe Sleep_Night `treatments' baseline `controls' if post_treatment==1, cluster(pid) absorb(date day_in_study)
			
			local pids_Sleep_Night = e(N_clust)
			local obs_Sleep_Night = e(N)
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_Sleep_Night_ns1 = string(r(estimate),"%3.2f")
				local se_Sleep_Night_ns1 = string(r(se),"%3.2f")
				local p_Sleep_Night_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_Sleep_Night_ns2 = string(r(estimate),"%3.2f")
				local se_Sleep_Night_ns2 = string(r(se),"%3.2f")
				local p_Sleep_Night_ns2 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_Sleep_Night_ns = string(r(estimate),"%3.2f")
				local se_Sleep_Night_ns = string(r(se),"%3.2f")
				local p_Sleep_Night_ns = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_Sleep_Night_nap = string(r(estimate),"%3.2f")
				local se_Sleep_Night_nap = string(r(se),"%3.2f")
				local p_Sleep_Night_nap = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[baseline]
				local coef_Sleep_Night_b = string(r(estimate),"%3.2f")
				local se_Sleep_Night_b = string(r(se),"%3.2f")
				local p_Sleep_Night_b = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "cell" {
			foreach j in s s_i int_s int_s_i np{

			lincom _b[treat_`j'_cell]
				local coef_Sleep_Night_`j' = string(r(estimate),"%3.2f")
				local se_Sleep_Night_`j' = string(r(se),"%3.2f")
				local tstat_Sleep_Night_`j' = r(estimate)/r(se)
				local p_Sleep_Night_`j' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			}
			* Get pval of differences
			lincom _b[treat_int_s_cell] - _b[treat_s_cell] - _b[treat_np_cell]
				local p_Sleep_Night_d_enc = string(r(p), "%3.2f")
			lincom _b[treat_int_s_i_cell] - _b[treat_s_i_cell] - _b[treat_np_cell]
				local p_Sleep_Night_d_inc = string(r(p), "%3.2f")
			}
				
			cap drop in_sleep
			gen in_sleep = e(sample)
	
		
		*Actigraph Time in Bed
			
			capture drop _baseline baseline
			gen _baseline			= act_inbed if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
					
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
			
			local pids_act_inbed = e(N_clust)
			local obs_act_inbed = e(N)
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_act_inbed_ns1 = string(r(estimate),"%3.2f")
				local se_act_inbed_ns1 = string(r(se),"%3.2f")
				local p_act_inbed_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_act_inbed_ns2 = string(r(estimate),"%3.2f")
				local se_act_inbed_ns2 = string(r(se),"%3.2f")
				local p_act_inbed_ns2 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_act_inbed_ns = string(r(estimate),"%3.2f")
				local se_act_inbed_ns = string(r(se),"%3.2f")
				local p_act_inbed_ns = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_act_inbed_nap = string(r(estimate),"%3.2f")
				local se_act_inbed_nap = string(r(se),"%3.2f")
				local p_act_inbed_nap = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[baseline]
				local coef_act_inbed_b = string(r(estimate),"%3.2f")
				local se_act_inbed_b = string(r(se),"%3.2f")
				local p_act_inbed_b = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "cell" {
			foreach j in s s_i int_s int_s_i np{

			lincom _b[treat_`j'_cell]
				local coef_act_inbed_`j' = string(r(estimate),"%3.2f")
				local se_act_inbed_`j' = string(r(se),"%3.2f")
				local tstat_act_inbed_`j' = r(estimate)/r(se)
				local p_act_inbed_`j' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			* Get pval of differences
			lincom _b[treat_int_s_cell] - _b[treat_s_cell] - _b[treat_np_cell]
				local p_act_inbed_d_enc = string(r(p), "%3.2f")
			lincom _b[treat_int_s_i_cell] - _b[treat_s_i_cell] - _b[treat_np_cell]
				local p_act_inbed_d_inc = string(r(p), "%3.2f")
			}
		
		*Actigraph sleep efficiency (Sleep in night/Time in bed)
		
		replace sleep_eff = sleep_eff*100
			
				capture drop _baseline baseline
				gen _baseline			= sleep_eff if post_treatment==0
				egen baseline 	= mean(_baseline), by(pid)
					
			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
				
			local pids_sleep_eff = e(N_clust)
			local obs_sleep_eff = e(N)	
				
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_sleep_eff_ns1 = string(r(estimate),"%3.2f")
				local se_sleep_eff_ns1 = string(r(se),"%3.2f")
				local p_sleep_eff_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_sleep_eff_ns2 = string(r(estimate),"%3.2f")
				local se_sleep_eff_ns2 = string(r(se),"%3.2f")
				local p_sleep_eff_ns2 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_sleep_eff_ns = string(r(estimate),"%3.2f")
				local se_sleep_eff_ns = string(r(se),"%3.2f")
				local p_sleep_eff_ns = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_sleep_eff_nap = string(r(estimate),"%3.2f")
				local se_sleep_eff_nap = string(r(se),"%3.2f")
				local p_sleep_eff_nap = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[baseline]
				local coef_sleep_eff_b = string(r(estimate),"%3.2f")
				local se_sleep_eff_b = string(r(se),"%3.2f")
				local p_sleep_eff_b = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "cell" {
				foreach j in s s_i int_s int_s_i np{
	
				lincom _b[treat_`j'_cell]
					local coef_sleep_eff_`j' = string(r(estimate),"%3.2f")
					local se_sleep_eff_`j' = string(r(se),"%3.2f")
					local tstat_sleep_eff_`j' = r(estimate)/r(se)
					local p_sleep_eff_`j' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				}
			* Get pval of differences
			lincom _b[treat_int_s_cell] - _b[treat_s_cell] - _b[treat_np_cell]
				local p_sleep_eff_d_enc = string(r(p), "%3.2f")
			lincom _b[treat_int_s_i_cell] - _b[treat_s_i_cell] - _b[treat_np_cell]
				local p_sleep_eff_d_inc = string(r(p), "%3.2f")
				}
				
		*Nap Sleep
		
		local controls "i.quart_age Female"
					
			reghdfe nap_hr `treatments' `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			local pids_nap_hr = e(N_clust)
			local obs_nap_hr = e(N)			
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_nap_hr_ns1 = string(r(estimate),"%3.2f")
				local se_nap_hr_ns1 = string(r(se),"%3.2f")
				local p_nap_hr_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_nap_hr_ns2 = string(r(estimate),"%3.2f")
				local se_nap_hr_ns2 = string(r(se),"%3.2f")
				local p_nap_hr_ns2 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_nap_hr_ns = string(r(estimate),"%3.2f")
				local se_nap_hr_ns = string(r(se),"%3.2f")
				local p_nap_hr_ns = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_nap_hr_nap = string(r(estimate),"%3.2f")
				local se_nap_hr_nap = string(r(se),"%3.2f")
				local p_nap_hr_nap = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))

			}
			if "`type'" == "cell" {
			foreach j in s s_i int_s int_s_i np{

			lincom _b[treat_`j'_cell]
				local coef_nap_hr_`j' = string(r(estimate),"%3.2f")
				local se_nap_hr_`j' = string(r(se),"%3.2f")
				local tstat_nap_hr_`j' = r(estimate)/r(se)
				local p_nap_hr_`j' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			* Get pval of differences
			lincom _b[treat_int_s_cell] - _b[treat_s_cell] - _b[treat_np_cell]
				local p_nap_hr_d_enc = string(r(p), "%3.2f")
			lincom _b[treat_int_s_i_cell] - _b[treat_s_i_cell] - _b[treat_np_cell]
				local p_nap_hr_d_inc = string(r(p), "%3.2f")
			}
					
			cap drop in_nap
			gen in_nap = e(sample)

		*24-Hour Sleep
			
				capture drop _baseline baseline
				gen _baseline			= sleep_all if post_treatment==0
				egen baseline 	= mean(_baseline), by(pid)
						
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			local pids_sleep_all = e(N_clust)
			local obs_sleep_all = e(N)		
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_sleep_all_ns1 = string(r(estimate),"%3.2f")
				local se_sleep_all_ns1 = string(r(se),"%3.2f")
				local p_sleep_all_ns1 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_sleep_all_ns2 = string(r(estimate),"%3.2f")
				local se_sleep_all_ns2 = string(r(se),"%3.2f")
				local p_sleep_all_ns2 = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_sleep_all_ns = string(r(estimate),"%3.2f")
				local se_sleep_all_ns = string(r(se),"%3.2f")
				local p_sleep_all_ns = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_sleep_all_nap = string(r(estimate),"%3.2f")
				local se_sleep_all_nap = string(r(se),"%3.2f")
				local p_sleep_all_nap = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			lincom _b[baseline]
				local coef_sleep_all_b = string(r(estimate),"%3.2f")
				local se_sleep_all_b = string(r(se),"%3.2f")
				local p_sleep_all_b = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "cell" {
			foreach j in s s_i int_s int_s_i np {

			lincom _b[treat_`j'_cell]
				local coef_sleep_all_`j' = string(r(estimate),"%3.2f")
				local se_sleep_all_`j' = string(r(se),"%3.2f")
				local tstat_sleep_all_`j' = r(estimate)/r(se)
				local p_sleep_all_`j' = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			* Get pval of differences
			lincom _b[treat_int_s_cell] - _b[treat_s_cell] - _b[treat_np_cell]
				local p_sleep_all_d_enc = string(r(p), "%3.2f")
			lincom _b[treat_int_s_i_cell] - _b[treat_s_i_cell] - _b[treat_np_cell]
				local p_sleep_all_d_inc = string(r(p), "%3.2f")
			}
					
			cap drop in_24
			gen in_24 = e(sample)
			
		lab var baseline "Baseline"

		restore
}


***
* SUMMARY STATS
***

preserve


foreach var in Sleep_Night sleep_eff act_inbed nap_hr sleep_all  { 

if "`var'" == "sleep_eff" replace `var' = `var'*100

sum `var' if treatment_group ==0 & post_treatment==1 & treat_nap==0

	local mean_`var' = string(r(mean), "%3.2f")
	local sd_`var' = string(r(sd), "%3.2f")
}
restore
***
* SIGNIFICANCE STARS
***

*Panel B
foreach var in Sleep_Night sleep_eff act_inbed nap_hr sleep_all {
	foreach i in ns ns1 ns2 nap s s_i int_s int_s_i np {
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
file close _all

			cd "C:\Users\Harshit\Downloads"
			file open f using "Table2_treatment_effects_on_sleep_cell.tex", write replace
				file write f "\begin{tabular*}{1.3\textwidth}{l@{\extracolsep{\fill}}*{5}{c}}" _n ///
			"\toprule \\" _n ///
			" & \textbf{Night Sleep} & \textbf{Time in Bed} & \textbf{Sleep Efficiency} & \textbf{Nap Sleep} & \textbf{24-Hr Sleep} \\" _n ///
			///
			" & (1) &(2) & (3) & (4) & (5) \\" _n ///
			"\midrule" _n ///
			"Devices+Encouragement Only & `coef_Sleep_Night_s' & `coef_act_inbed_s' & `coef_sleep_eff_s' & `coef_nap_hr_s' & `coef_sleep_all_s' \\" _n ///
			" & (`se_Sleep_Night_s') & (`se_act_inbed_s') & (`se_sleep_eff_s') & (`se_nap_hr_s') & (`se_sleep_all_s') \\" _n ///	
			"\addlinespace" _n ///
			///
			"Devices+Incentives Only & `coef_Sleep_Night_s_i' & `coef_act_inbed_s_i' & `coef_sleep_eff_s_i' & `coef_nap_hr_s_i' & `coef_sleep_all_s_i' \\" _n ///
			" & (`se_Sleep_Night_s_i') & (`se_act_inbed_s_i') & (`se_sleep_eff_s_i') & (`se_nap_hr_s_i') & (`se_sleep_all_s_i') \\" _n ///	
			"\addlinespace" _n ///
			///
			"Nap Only & `coef_Sleep_Night_np' & `coef_act_inbed_np' & `coef_sleep_eff_np' & `coef_nap_hr_np' & `coef_sleep_all_np' \\" _n ///
			" & (`se_Sleep_Night_np') & (`se_act_inbed_np') & (`se_sleep_eff_np') & (`se_nap_hr_np') & (`se_sleep_all_np') \\" _n ///	
			"\addlinespace" _n ///
			///
			"Devices+Encouragement and Nap & `coef_Sleep_Night_int_s' & `coef_act_inbed_int_s' & `coef_sleep_eff_int_s' & `coef_nap_hr_int_s' & `coef_sleep_all_int_s' \\" _n ///
			" & (`se_Sleep_Night_int_s') & (`se_act_inbed_int_s') & (`se_sleep_eff_int_s') & (`se_nap_hr_int_s') & (`se_sleep_all_int_s') \\" _n ///	
			"\addlinespace" _n ///
			///
			"Devices+Incentives and Nap & `coef_Sleep_Night_int_s_i' & `coef_act_inbed_int_s_i' & `coef_sleep_eff_int_s_i' & `coef_nap_hr_int_s_i' & `coef_sleep_all_int_s_i' \\" _n ///
			" & (`se_Sleep_Night_int_s_i') & (`se_act_inbed_int_s_i') & (`se_sleep_eff_int_s_i') & (`se_nap_hr_int_s_i') & (`se_sleep_all_int_s_i') \\" _n ///	
			"\addlinespace" _n ///
			///
			"\midrule" _n ///
			"Control Mean               & `mean_Sleep_Night' & `mean_act_inbed' & `mean_sleep_eff' & `mean_nap_hr' & `mean_sleep_all'  \\" _n ///
			"Control SD                & `sd_Sleep_Night' & `sd_act_inbed' & `sd_sleep_eff' & `sd_nap_hr' & `sd_sleep_all' \\" _n ///
			"Participant-nights               & `obs_Sleep_Night' & `obs_act_inbed' & `obs_sleep_eff' & `obs_nap_hr' & `obs_sleep_all' \\" _n ///
			"Participants               & `pids_Sleep_Night' & `pids_act_inbed' & `pids_sleep_eff' & `pids_nap_hr' & `pids_sleep_all' \\" _n ///
			"\bottomrule" _n ///
			"\end{tabular*}" _n
	file close f
