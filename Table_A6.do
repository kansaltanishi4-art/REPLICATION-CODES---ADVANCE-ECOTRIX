************************************************
* Sleep Project - Pedro Bessone, Gautam Rao, Heather Schofield, Frank Schilbach, and Mattie Toma
* Purpose: Replicates Appendix Table 6 (Treatment Effects on Sleep, Pooling Night-Sleep Treatments and Including Bounds)
* Last edited: 07 May 2021
************************************************


	************************************************	
	*0. Initial Setup
	************************************************

	clear all
	set more off
	set matsize 800
	
* Treatment Effects on Sleep Quantity

use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\firststage_dataset.dta", clear

cap drop sleep_all
gen sleep_all= sleep_time

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
			}
			cap drop in_sleep
			gen in_sleep = e(sample)
	
		
		*Actigraph Time in Bed
			
			capture drop _baseline baseline
			gen _baseline			= act_inbed if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
					
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
						
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
			}
		
		*Actigraph sleep efficiency (Sleep in night/Time in bed)
		
		replace sleep_eff = sleep_eff*100
			
				capture drop _baseline baseline
				gen _baseline			= sleep_eff if post_treatment==0
				egen baseline 	= mean(_baseline), by(pid)
					
			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
								
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
				}
				
		*Nap Sleep
		
		local controls "i.quart_age Female"
					
			reghdfe nap_hr `treatments' `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
						
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
			}
					
			cap drop in_nap
			gen in_nap = e(sample)

		*24-Hour Sleep
			
				capture drop _baseline baseline
				gen _baseline			= sleep_all if post_treatment==0
				egen baseline 	= mean(_baseline), by(pid)
						
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
						
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

	* LEE BOUNDS

	use "C:\Users\Harshit\Downloads\dataverse_files (1)\Replication_package_Economic_consequences_sleep\Replication_package_Economic_consequences_sleep\Datasets\firststage_dataset.dta", clear

	cap drop sleep_all
	gen sleep_all= sleep_time
	
	replace sleep_eff = sleep_eff*100
	
	// calculating attrition 
	
	* sleep night, night treatments
	gen treat_att = 1 if treat_pool==1 & Sleep_Night==. & post_treatment==1 
	replace treat_att=0 if treat_att==. & treat_pool==1 & post_treatment==1 
	// 6% attrition
	gen control_att=1 if treat_pool==0 & Sleep_Night==. & post_treatment==1 
	replace control_att=0 if control_att==. & treat_pool==0 & post_treatment==1 
	// 7.3% of attrition 
	
	* sleep night, nap treatment 
	gen treat_att_X= 1 if treat_nap==1 & Sleep_Night==. & post_treatment==1 
	replace treat_att_X=0 if treat_att==. & treat_nap==1 & post_treatment==1 
	// 7.4% attrition
	gen control_att_X=1 if treat_nap==0 & Sleep_Night==. & post_treatment==1 
	replace control_att_X=0 if control_att_X==. & treat_nap==0 & post_treatment==1 
	// 5.6% attrition 
			
	
	* nap, nap treatment
	gen nap_treat_att = 1 if treat_nap==1 & nap_hr==. & post_treatment==1 & day_in_study!=28 & at_present_check==1
	replace nap_treat_att=0 if nap_treat_att==. & treat_nap==1 & post_treatment==1 & day_in_study!=28 & at_present_check==1
	// 2.7% attrition
	gen nap_control_att=1 if treat_nap==0 & nap_hr==. & post_treatment==1 & day_in_study!=28 & at_present_check==1
	replace nap_control_att=0 if nap_control_att==. & treat_nap==0 & post_treatment==1 & day_in_study!=28 & at_present_check==1
	// 0% of attrition 
	
	// differential attrition = 2.7% ; control has no attrition by construction
		
	* constructing lee & manski 'samples', re running regs 

	
	foreach type in "unpooled" "pooled" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}

		eststo clear
		

	* Actigraph Sleep Time
	
	
			foreach bound in "lee"{
			foreach conf in "l" "u" {
			
			local controls "i.quart_age Female"
			
			capture drop _baseline baseline
			gen _baseline			= Sleep_Night if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
			
			
			if "`bound'" == "lee" {
			if "`conf'" == "l" {
			
			// for night sleep treatments
			
			preserve 
			
			// drop the highest 1.3% of treatment 
			
			egen rank_treat = rank(Sleep_Night) if treat_pool==1 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0.987, 1) & treat_pool==1 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
			
				reghdfe Sleep_Night `treatments' baseline `controls' if post_treatment==1, cluster(pid) absorb(date day_in_study)
			
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_sleep_s_lee_l = string(r(estimate),"%3.2f")
				local se_sleep_s_lee_l = string(r(se),"%3.2f")
				local p_sleep_s_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_s_lee_l_l = string(r(lb),"%3.2f")
				local coef_sleep_s_lee_l_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_sleep_s_i_lee_l = string(r(estimate),"%3.2f")
				local se_sleep_s_i_lee_l = string(r(se),"%3.2f")
				local p_sleep_s_i_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_s_i_lee_l_l = string(r(lb),"%3.2f")
				local coef_sleep_s_i_lee_l_u = string(r(ub),"%3.2f")
				
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_sleep_pool_lee_l = string(r(estimate),"%3.2f")
				local se_sleep_pool_lee_l = string(r(se),"%3.2f")
				local p_sleep_pool_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_pool_lee_l_l = string(r(lb),"%3.2f")
				local coef_sleep_pool_lee_l_u = string(r(ub),"%3.2f")
				
			}
			
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the lowest 1.8% of control
			
			egen rank = rank(Sleep_Night) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0, 0.018) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
				reghdfe Sleep_Night `treatments' baseline `controls' if post_treatment==1, cluster(pid) absorb(date day_in_study)
			
			if "`type'" == "pooled" {
			
				
			lincom _b[treat_nap]
				local coef_sleep_nap_lee_l = string(r(estimate),"%3.2f")
				local se_sleep_nap_lee_l = string(r(se),"%3.2f")
				local p_sleep_nap_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_nap_lee_l_l = string(r(lb),"%3.2f")
				local coef_sleep_nap_lee_l_u = string(r(ub),"%3.2f")
			
				
			}
			restore 
			}
			
			
			if "`conf'" == "u" {
			
			preserve 
			
			// drop lowest 1.3% of treatment 
			egen rank_treat = rank(Sleep_Night) if treat_pool==1 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0, 0.013) & treat_pool==1 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
			
				reghdfe Sleep_Night `treatments' baseline `controls' if post_treatment==1, cluster(pid) absorb(date day_in_study)
						
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_sleep_s_lee_u = string(r(estimate),"%3.2f")
				local se_sleep_s_lee_u = string(r(se),"%3.2f")
				local p_sleep_s_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_s_lee_u_l = string(r(lb),"%3.2f")
				local coef_sleep_s_lee_u_u = string(r(ub),"%3.2f")
				
				
			cap lincom _b[treat_s_i]
				local coef_sleep_s_i_lee_u = string(r(estimate),"%3.2f")
				local se_sleep_s_i_lee_u = string(r(se),"%3.2f")
				local p_sleep_s_i_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_s_i_lee_u_l = string(r(lb),"%3.2f")
				local coef_sleep_s_i_lee_u_u = string(r(ub),"%3.2f")
				
				
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_sleep_pool_lee_u = string(r(estimate),"%3.2f")
				local se_sleep_pool_lee_u = string(r(se),"%3.2f")
				local p_sleep_pool_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_pool_lee_u_l = string(r(lb),"%3.2f")
				local coef_sleep_pool_lee_u_u = string(r(ub),"%3.2f")
				
			}
			
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the highest 1.8% of control
			
			egen rank = rank(Sleep_Night) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0.982, 1) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
				reghdfe Sleep_Night `treatments' baseline `controls' if post_treatment==1, cluster(pid) absorb(date day_in_study)

			if "`type'" == "pooled" {
							
			lincom _b[treat_nap]
				local coef_sleep_nap_lee_u = string(r(estimate),"%3.2f")
				local se_sleep_nap_lee_u = string(r(se),"%3.2f")
				local p_sleep_nap_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_sleep_nap_lee_u_l = string(r(lb),"%3.2f")
				local coef_sleep_nap_lee_u_u = string(r(ub),"%3.2f")
			
				
			}
			restore 
			}
			
			}
			}
			}
			}
			
	foreach type in "unpooled" "pooled" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}

		eststo clear
	

		* Actigraph Time in Bed
		
		
			foreach bound in "lee" {
			foreach conf in "l" "u" {
			
			local controls "i.quart_age Female"
			
			
			capture drop _baseline baseline
			gen _baseline	= act_inbed if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
			
			cap drop in_sleep
			gen in_sleep=1 if Sleep_Night!=. & treat_pool!=. & treat_nap!=. & baseline!=. & post_treatment!=. & pid!=. 
		
			if "`bound'" == "lee" {
			if "`conf'" == "l" {
			
			preserve 
			
			// drop the highest 1.3% of treatment 
			egen rank_treat = rank(act_inbed) if treat_pool==1 & post_treatment==1 & in_sleep==1, unique 
			egen max_rank = max(rank_treat) 
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0.987, 1) & treat_pool==1 & post_treatment==1 & in_sleep==1
			drop rank_treat max_rank pct_treat
		
			// re run regs 
		
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
					
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_inbed_s_lee_l = string(r(estimate),"%3.2f")
				local se_inbed_s_lee_l = string(r(se),"%3.2f")
				local p_inbed_s_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_s_lee_l_l = string(r(lb),"%3.2f")
				local coef_inbed_s_lee_l_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_inbed_s_i_lee_l = string(r(estimate),"%3.2f")
				local se_inbed_s_i_lee_l = string(r(se),"%3.2f")
				local p_inbed_s_i_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_s_i_lee_l_l = string(r(lb),"%3.2f")
				local coef_inbed_s_i_lee_l_u = string(r(ub),"%3.2f")
				
				
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_inbed_pool_lee_l = string(r(estimate),"%3.2f")
				local se_inbed_pool_lee_l = string(r(se),"%3.2f")
				local p_inbed_pool_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_pool_lee_l_l = string(r(lb),"%3.2f")
				local coef_inbed_pool_lee_l_u = string(r(ub),"%3.2f")
								
			lincom _b[baseline]
				local coef_inbed_bsl_lee_l = string(r(estimate),"%3.2f")
				local se_inbed_bsl_lee_l = string(r(se),"%3.2f")
				local p_inbed_bsl_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the lowest 1.8% of control
			
			egen rank = rank(act_inbed) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0, 0.018) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
					
			
			if "`type'" == "pooled" {
				
			lincom _b[treat_nap]
				local coef_inbed_nap_lee_l = string(r(estimate),"%3.2f")
				local se_inbed_nap_lee_l = string(r(se),"%3.2f")
				local p_inbed_nap_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_nap_lee_l_l = string(r(lb),"%3.2f")
				local coef_inbed_nap_lee_l_u = string(r(ub),"%3.2f")
				
				
			}
			restore 
			} 
			
			
			if "`conf'" == "u" {
			
			preserve 
			
			// drop lowest 1.3% of treatment 
			
			egen rank_treat = rank(act_inbed) if treat_pool==1 & post_treatment==1 & in_sleep==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0, 0.013) & treat_pool==1 & post_treatment==1 & in_sleep==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
			
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
					
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_inbed_s_lee_u = string(r(estimate),"%3.2f")
				local se_inbed_s_lee_u = string(r(se),"%3.2f")
				local p_inbed_s_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_s_lee_u_l = string(r(lb),"%3.2f")
				local coef_inbed_s_lee_u_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_inbed_s_i_lee_u = string(r(estimate),"%3.2f")
				local se_inbed_s_i_lee_u = string(r(se),"%3.2f")
				local p_inbed_s_i_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_s_i_lee_u_l = string(r(lb),"%3.2f")
				local coef_inbed_s_i_lee_u_u = string(r(ub),"%3.2f")
				
				
			}
			
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_inbed_pool_lee_u = string(r(estimate),"%3.2f")
				local se_inbed_pool_lee_u = string(r(se),"%3.2f")
				local p_inbed_pool_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_pool_lee_u_l = string(r(lb),"%3.2f")
				local coef_inbed_pool_lee_u_u = string(r(ub),"%3.2f")
				
				
			lincom _b[baseline]
				local coef_inbed_bsl_lee_u = string(r(estimate),"%3.2f")
				local se_inbed_bsl_lee_u = string(r(se),"%3.2f")
				local p_inbed_bsl_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the highest 1.8% of control
			
			egen rank = rank(act_inbed) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0.982, 1) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
			reghdfe act_inbed `treatments' baseline `controls' if post_treatment==1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
					
			if "`type'" == "pooled" {
			lincom _b[treat_nap]
				local coef_inbed_nap_lee_u = string(r(estimate),"%3.2f")
				local se_inbed_nap_lee_u = string(r(se),"%3.2f")
				local p_inbed_nap_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_inbed_nap_lee_u_l = string(r(lb),"%3.2f")
				local coef_inbed_nap_lee_u_u = string(r(ub),"%3.2f")
				
			}
			restore 
			
			} 
			}
				
			
			}
			}  
			}
			
			

			
	
	foreach type in "unpooled" "pooled" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}

		eststo clear
		
			
		* Actigraph sleep efficiency (Sleep in night/Time in bed)
		
		
			foreach bound in "lee" {
			foreach conf in "l" "u" {
		
			
			local controls "i.quart_age Female"
			
			capture drop _baseline baseline
			gen _baseline			= sleep_eff if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
			
			cap drop in_sleep
			gen in_sleep=1 if Sleep_Night!=. & treat_pool!=. & treat_nap!=. & baseline!=. & post_treatment!=. & pid!=. 
			
			if "`bound'" == "lee" {
			if "`conf'" == "l" {
			
			preserve 
			
			// drop the highest 1.3% of treatment 
		
			egen rank_treat = rank(sleep_eff) if treat_pool==1 & post_treatment==1 & in_sleep==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0.987, 1) & treat_pool==1 & post_treatment==1 & in_sleep==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
		
			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
		
				
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_eff_s_lee_l = string(r(estimate),"%3.2f")
				local se_eff_s_lee_l = string(r(se),"%3.2f")
				local p_eff_s_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_s_lee_l_l = string(r(lb),"%3.2f")
				local coef_eff_s_lee_l_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_eff_s_i_lee_l = string(r(estimate),"%3.2f")
				local se_eff_s_i_lee_l = string(r(se),"%3.2f")
				local p_eff_s_i_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_s_i_lee_l_l = string(r(lb),"%3.2f")
				local coef_eff_s_i_lee_l_u = string(r(ub),"%3.2f")
				
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_eff_pool_lee_l = string(r(estimate),"%3.2f")
				local se_eff_pool_lee_l = string(r(se),"%3.2f")
				local p_eff_pool_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_pool_lee_l_l = string(r(lb),"%3.2f")
				local coef_eff_pool_lee_l_u = string(r(ub),"%3.2f")
				
				
			lincom _b[baseline]
				local coef_eff_bsl_lee_l = string(r(estimate),"%3.2f")
				local se_eff_bsl_lee_l = string(r(se),"%3.2f")
				local p_eff_bsl_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			// for nap treatment 
			
			preserve 
			
			// drop the lowest 1.8% of control
			
			egen rank = rank(sleep_eff) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0, 0.018) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
		
			if "`type'" == "pooled" {
			
				
			lincom _b[treat_nap]
				local coef_eff_nap_lee_l = string(r(estimate),"%3.2f")
				local se_eff_nap_lee_l = string(r(se),"%3.2f")
				local p_eff_nap_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_nap_lee_l_l = string(r(lb),"%3.2f")
				local coef_eff_nap_lee_l_u = string(r(ub),"%3.2f")
				
			lincom _b[baseline]
				local coef_eff_bsl_lee_l = string(r(estimate),"%3.2f")
				local se_eff_bsl_lee_l = string(r(se),"%3.2f")
				local p_eff_bsl_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			} 
			
			
			if "`conf'" == "u" {
			
			preserve 
			
			// drop lowest 1.3% of treatment 
			
			egen rank_treat = rank(sleep_eff) if treat_pool==1 & post_treatment==1 & in_sleep==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0, 0.013) & treat_pool==1 & post_treatment==1 & in_sleep==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 

			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
		
					
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_eff_s_lee_u = string(r(estimate),"%3.2f")
				local se_eff_s_lee_u = string(r(se),"%3.2f")
				local p_eff_s_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_s_lee_u_l = string(r(lb),"%3.2f")
				local coef_eff_s_lee_u_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_eff_s_i_lee_u = string(r(estimate),"%3.2f")
				local se_eff_s_i_lee_u = string(r(se),"%3.2f")
				local p_eff_s_i_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_s_i_lee_u_l = string(r(lb),"%3.2f")
				local coef_eff_s_i_lee_u_u = string(r(ub),"%3.2f")
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_eff_pool_lee_u = string(r(estimate),"%3.2f")
				local se_eff_pool_lee_u = string(r(se),"%3.2f")
				local p_eff_pool_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_pool_lee_u_l = string(r(lb),"%3.2f")
				local coef_eff_pool_lee_u_u = string(r(ub),"%3.2f")
				
				
			lincom _b[baseline]
				local coef_eff_bsl_lee_u = string(r(estimate),"%3.2f")
				local se_eff_bsl_lee_u = string(r(se),"%3.2f")
				local p_eff_bsl_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the highest 1.8% of control
			
			egen rank = rank(sleep_eff) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0.982, 1) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs
			
			reghdfe sleep_eff `treatments' baseline `controls' if post_treatment == 1 & in_sleep==1, cluster(pid) absorb(date day_in_study)
		
			
			if "`type'" == "pooled" {
				
			lincom _b[treat_nap]
				local coef_eff_nap_lee_u = string(r(estimate),"%3.2f")
				local se_eff_nap_lee_u = string(r(se),"%3.2f")
				local p_eff_nap_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_eff_nap_lee_u_l = string(r(lb),"%3.2f")
				local coef_eff_nap_lee_u_u = string(r(ub),"%3.2f")
				
			lincom _b[baseline]
				local coef_eff_bsl_lee_u = string(r(estimate),"%3.2f")
				local se_eff_bsl_lee_u = string(r(se),"%3.2f")
				local p_eff_bsl_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			} 
				
			}
			
			}
			} 
			}	
					
	foreach type in "unpooled" "pooled" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}

		eststo clear
			
					
		* Nap Sleep
		
			foreach bound in "lee"  {
			foreach conf in "l" "u" {
			
			local controls "i.quart_age Female"
			
			if "`bound'" == "lee" {
			if "`conf'" == "u" {
			
			preserve 
			
			egen rank_treat = rank(nap_hr) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0.973, 1) & treat_nap==0 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
		
			reghdfe nap_hr `treatments' `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
				
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_nap_s_lee_u = string(r(estimate),"%3.2f")
				local se_nap_s_lee_u = string(r(se),"%3.2f")
				local p_nap_s_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_s_lee_u_l = string(r(lb),"%3.2f")
				local coef_nap_s_lee_u_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_nap_s_i_lee_u = string(r(estimate),"%3.2f")
				local se_nap_s_i_lee_u = string(r(se),"%3.2f")
				local p_nap_s_i_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_s_i_lee_u_l = string(r(lb),"%3.2f")
				local coef_nap_s_i_lee_u_u = string(r(ub),"%3.2f")
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_nap_pool_lee_u = string(r(estimate),"%3.2f")
				local se_nap_pool_lee_u = string(r(se),"%3.2f")
				local p_nap_pool_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_pool_lee_u_l = string(r(lb),"%3.2f")
				local coef_nap_pool_lee_u_u = string(r(ub),"%3.2f")
				
			lincom _b[treat_nap]
				local coef_nap_nap_lee_u = string(r(estimate),"%3.2f")
				local se_nap_nap_lee_u = string(r(se),"%3.2f")
				local p_nap_nap_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_nap_lee_u_l = string(r(lb),"%3.2f")
				local coef_nap_nap_lee_u_u = string(r(ub),"%3.2f")

			}
			restore 
			}
		
			
			if "`conf'" == "l" {
			
			preserve 
			
			// drop lowest 2.7% of control
			
			egen rank_treat = rank(nap_hr) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0, 0.027) & treat_nap==0 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
			
			reghdfe nap_hr `treatments' `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
							
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_nap_s_lee_l = string(r(estimate),"%3.2f")
				local se_nap_s_lee_l = string(r(se),"%3.2f")
				local p_nap_s_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_s_lee_l_l = string(r(lb),"%3.2f")
				local coef_nap_s_lee_l_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_nap_s_i_lee_l = string(r(estimate),"%3.2f")
				local se_nap_s_i_lee_l = string(r(se),"%3.2f")
				local p_nap_s_i_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_s_i_lee_l_l = string(r(lb),"%3.2f")
				local coef_nap_s_i_lee_l_u = string(r(ub),"%3.2f")
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_nap_pool_lee_l = string(r(estimate),"%3.2f")
				local se_nap_pool_lee_l = string(r(se),"%3.2f")
				local p_nap_pool_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_pool_lee_l_l = string(r(lb),"%3.2f")
				local coef_nap_pool_lee_l_u = string(r(ub),"%3.2f")
				
			lincom _b[treat_nap]
				local coef_nap_nap_lee_l = string(r(estimate),"%3.2f")
				local se_nap_nap_lee_l = string(r(se),"%3.2f")
				local p_nap_nap_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_nap_nap_lee_l_l = string(r(lb),"%3.2f")
				local coef_nap_nap_lee_l_u = string(r(ub),"%3.2f")
			}
			restore 
			}
			
			}

			}
			} 
			}
			
			
			
		foreach type in "unpooled" "pooled" {
	
		if "`type'" == "pooled" {
			
			local treatments = "treat_pool treat_nap"
			
		}
		if "`type'" == "unpooled" {
		
			local treatments = "treat_s treat_s_i treat_nap"
		
		}

		eststo clear
			
			
		* 24-Hour Sleep
		
		
			foreach bound in "lee"  {
			foreach conf in "l" "u" {
			
			local controls "i.quart_age Female"
			
			capture drop _baseline baseline
			gen _baseline			= sleep_all if post_treatment==0
			egen baseline 	= mean(_baseline), by(pid)
						
			if "`bound'" == "lee" {
			if "`conf'" == "l" {
			
			preserve 
			
			// drop the highest 1.3% of treatment 
			
			egen rank_treat = rank(sleep_all) if treat_pool==1 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0.987, 1) & treat_pool==1 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
				
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_all_s_lee_l = string(r(estimate),"%3.2f")
				local se_all_s_lee_l = string(r(se),"%3.2f")
				local p_all_s_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_s_lee_l_l = string(r(lb),"%3.2f")
				local coef_all_s_lee_l_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_all_s_i_lee_l = string(r(estimate),"%3.2f")
				local se_all_s_i_lee_l = string(r(se),"%3.2f")
				local p_all_s_i_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_s_i_lee_l_l = string(r(lb),"%3.2f")
				local coef_all_s_i_lee_l_u = string(r(ub),"%3.2f")
				
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_all_pool_lee_l = string(r(estimate),"%3.2f")
				local se_all_pool_lee_l = string(r(se),"%3.2f")
				local p_all_pool_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_pool_lee_l_l = string(r(lb),"%3.2f")
				local coef_all_pool_lee_l_u = string(r(ub),"%3.2f")
				
			
			lincom _b[baseline]
				local coef_all_bsl_lee_l = string(r(estimate),"%3.2f")
				local se_all_bsl_lee_l = string(r(se),"%3.2f")
				local p_all_bsl_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			// for nap treatment 
			
			preserve 
			
			// drop the lowest 1.8% of control
			
			egen rank = rank(sleep_all) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0, 0.018) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			if "`type'" == "pooled" {
				
			lincom _b[treat_nap]
				local coef_all_nap_lee_l = string(r(estimate),"%3.2f")
				local se_all_nap_lee_l = string(r(se),"%3.2f")
				local p_all_nap_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_nap_lee_l_l = string(r(lb),"%3.2f")
				local coef_all_nap_lee_l_u = string(r(ub),"%3.2f")
			
			lincom _b[baseline]
				local coef_all_bsl_lee_l = string(r(estimate),"%3.2f")
				local se_all_bsl_lee_l = string(r(se),"%3.2f")
				local p_all_bsl_lee_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			}
		
			
			if "`conf'" == "u" {
			
			preserve 
			
			// drop lowest 1.3% of treatment 
			
			egen rank_treat = rank(sleep_all) if treat_pool==1 & post_treatment==1, unique 
			egen max_rank = max(rank_treat)
			gen pct_treat = rank_treat/max_rank
			drop if inrange(pct_treat, 0, 0.013) & treat_pool==1 & post_treatment==1
			drop rank_treat max_rank pct_treat
			
			// re run regs 
			
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_all_s_lee_u = string(r(estimate),"%3.2f")
				local se_all_s_lee_u = string(r(se),"%3.2f")
				local p_all_s_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_s_lee_u_l = string(r(lb),"%3.2f")
				local coef_all_s_lee_u_u = string(r(ub),"%3.2f")
				
			cap lincom _b[treat_s_i]
				local coef_all_s_i_lee_u = string(r(estimate),"%3.2f")
				local se_all_s_i_lee_u = string(r(se),"%3.2f")
				local p_all_s_i_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_s_i_lee_u_l = string(r(lb),"%3.2f")
				local coef_all_s_i_lee_u_u = string(r(ub),"%3.2f")
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_all_pool_lee_u = string(r(estimate),"%3.2f")
				local se_all_pool_lee_u = string(r(se),"%3.2f")
				local p_all_pool_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_pool_lee_u_l = string(r(lb),"%3.2f")
				local coef_all_pool_lee_u_u = string(r(ub),"%3.2f")
				
			
			lincom _b[baseline]
				local coef_all_bsl_lee_u = string(r(estimate),"%3.2f")
				local se_all_bsl_lee_u = string(r(se),"%3.2f")
				local p_all_bsl_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			// for nap treatment 
			
			preserve 
		
			// drop the highest 1.8% of control
			
			egen rank = rank(sleep_all) if treat_nap==0 & post_treatment==1, unique 
			egen max_rank = max(rank)
			gen pct = rank/max_rank
			drop if inrange(pct, 0.982, 1) & treat_nap==0 & post_treatment==1
			drop rank max_rank pct
			
			// re run regs 
			
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			
			if "`type'" == "pooled" {
			
				
			lincom _b[treat_nap]
				local coef_all_nap_lee_u = string(r(estimate),"%3.2f")
				local se_all_nap_lee_u = string(r(se),"%3.2f")
				local p_all_nap_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				local coef_all_nap_lee_u_l = string(r(lb),"%3.2f")
				local coef_all_nap_lee_u_u = string(r(ub),"%3.2f")
			
			lincom _b[baseline]
				local coef_all_bsl_lee_u = string(r(estimate),"%3.2f")
				local se_all_bsl_lee_u = string(r(se),"%3.2f")
				local p_all_bsl_lee_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			
			}
			
			}
			
			if "`bound'" == "manski" {
			
			if "`conf'" == "u" {
			
			preserve 
			
			egen min_sleep_all = min(sleep_all) if post_treatment==1
			egen max_sleep_all = max(sleep_all) if post_treatment==1
			
			replace sleep_all=max_sleep_all if treat_pool==1 & post_treatment==1 & sleep_all==. 
			replace sleep_all=min_sleep_all if treat_pool==0 & post_treatment==1 & sleep_all==. 
			
			// re run regs 
			
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
						
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_all_s_manski_u = string(r(estimate),"%3.2f")
				local se_all_s_manski_u = string(r(se),"%3.2f")
				local p_all_s_manski_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_all_s_i_manski_u = string(r(estimate),"%3.2f")
				local se_all_s_i_manski_u = string(r(se),"%3.2f")
				local p_all_s_i_manski_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_all_pool_manski_u = string(r(estimate),"%3.2f")
				local se_all_pool_manski_u = string(r(se),"%3.2f")
				local p_all_pool_manski_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_all_nap_manski_u = string(r(estimate),"%3.2f")
				local se_all_nap_manski_u = string(r(se),"%3.2f")
				local p_all_nap_manski_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			lincom _b[baseline]
				local coef_all_bsl_manski_u = string(r(estimate),"%3.2f")
				local se_all_bsl_manski_u = string(r(se),"%3.2f")
				local p_all_bsl_manski_u = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			}
			
			if "`conf'" == "l" {
			preserve 
			
			// replace treatment att with max and control att with minimum 
			
			egen min_sleep_all = min(sleep_all) if post_treatment==1
			egen max_sleep_all = max(sleep_all) if post_treatment==1
			
			replace sleep_all=min_sleep_all if treat_pool==1 & post_treatment==1 & sleep_all==. 
			replace sleep_all=max_sleep_all if treat_pool==0 & post_treatment==1 & sleep_all==. 
			// re run regs 
			
			reghdfe sleep_all `treatments' baseline `controls' if post_treatment==1 & day_in_study!=28, cluster(pid) absorb(date day_in_study)
			
			
			if "`type'" == "unpooled" {
			cap lincom _b[treat_s]
				local coef_all_s_manski_l = string(r(estimate),"%3.2f")
				local se_all_s_manski_l = string(r(se),"%3.2f")
				local p_all_s_manski_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			cap lincom _b[treat_s_i]
				local coef_all_s_i_manski_l = string(r(estimate),"%3.2f")
				local se_all_s_i_manski_l = string(r(se),"%3.2f")
				local p_all_s_i_manski_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			if "`type'" == "pooled" {
			cap lincom _b[treat_pool]
				local coef_all_pool_manski_l = string(r(estimate),"%3.2f")
				local se_all_pool_manski_l = string(r(se),"%3.2f")
				local p_all_pool_manski_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
				
			lincom _b[treat_nap]
				local coef_all_nap_manski_l = string(r(estimate),"%3.2f")
				local se_all_nap_manski_l = string(r(se),"%3.2f")
				local p_all_nap_manski_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			
			lincom _b[baseline]
				local coef_all_bsl_manski_l = string(r(estimate),"%3.2f")
				local se_all_bsl_manski_l = string(r(se),"%3.2f")
				local p_all_bsl_manski_l = (2 * ttail(e(df_r), abs(r(estimate)/r(se))))
			}
			restore 
			}
			}
			
			cap drop in_24
			gen in_24 = e(sample)
			
			}
			}
			}
	


***
* SUMMARY STATS
***

preserve


foreach var in Sleep_Night sleep_eff act_inbed nap_hr sleep_all  { 


sum `var' if treatment_group ==0 & post_treatment==1 & treat_nap==0

	local mean_`var' = string(r(mean), "%3.2f")
	local sd_`var' = string(r(sd), "%3.2f")
}
restore
	
	


*TABLE: LEEBOUNDS

cd "C:\Users\Harshit\Downloads"
file open f using "TableA6_first_stage_manual_att_bounds.tex", write replace
		file write f "\begin{tabular}{@{}c@{}}" _n ///
		"\begin{tabular}{l*{5}{c}}" _n ///
		"\toprule \\" _n ///
		" & \textbf{Night Sleep} & \textbf{Time in Bed} & \textbf{Sleep Efficiency} & \textbf{Nap Sleep} & \textbf{24-Hr Sleep} \\" _n ///
		///
		" & (1) & (2) & (3) & (4) & (5) \\" _n ///
		"\midrule" _n ///
		///
		"\textbf{Night-Sleep Treatments} & `coef_Sleep_Night_ns' & `coef_act_inbed_ns' & `coef_sleep_eff_ns' & `coef_nap_hr_ns' & `coef_sleep_all_ns' \\" _n ///
		" & (`se_Sleep_Night_ns') & (`se_act_inbed_ns') & (`se_sleep_eff_ns') & (`se_nap_hr_ns') & (`se_sleep_all_ns') \\" _n ///	
		"\hspace{.2cm} Lee Lower Bound & `coef_sleep_pool_lee_l' & `coef_inbed_pool_lee_l' & `coef_eff_pool_lee_l' & `coef_nap_pool_lee_l' & `coef_all_pool_lee_l' \\" _n ///
		"\hspace{.2cm} Lee Upper Bound & `coef_sleep_pool_lee_u' & `coef_inbed_pool_lee_u' & `coef_eff_pool_lee_u' & `coef_nap_pool_lee_u' & `coef_all_pool_lee_u' \\" _n ///
		"\hspace{.2cm} Confidence Interval & [`coef_sleep_pool_lee_l_l', `coef_sleep_pool_lee_u_u'] & [`coef_inbed_pool_lee_l_l', `coef_inbed_pool_lee_u_u'] & [`coef_eff_pool_lee_l_l', `coef_eff_pool_lee_u_u'] & [`coef_nap_pool_lee_l_l', `coef_nap_pool_lee_u_u'] & [`coef_all_pool_lee_l_l', `coef_all_pool_lee_u_u'] \\" _n ///
		"\addlinespace" _n ///
		///
		"\hspace{.1cm} \textbf{Devices+Encouragement} & `coef_Sleep_Night_ns1' & `coef_act_inbed_ns1' & `coef_sleep_eff_ns1' & `coef_nap_hr_ns1' & `coef_sleep_all_ns1' \\" _n ///
		" & (`se_Sleep_Night_ns1') & (`se_act_inbed_ns1') & (`se_sleep_eff_ns1') & (`se_nap_hr_ns1') & (`se_sleep_all_ns1') \\" _n ///	
		"\hspace{.2cm} Lee Lower Bound & `coef_sleep_s_lee_l' & `coef_inbed_s_lee_l' & `coef_eff_s_lee_l' & `coef_nap_s_lee_l' & `coef_all_s_lee_l' \\" _n ///
		"\hspace{.2cm} Lee Upper Bound & `coef_sleep_s_lee_u' & `coef_inbed_s_lee_u' & `coef_eff_s_lee_u' & `coef_nap_s_lee_u' & `coef_all_s_lee_u' \\" _n ///
		"\hspace{.2cm} Confidence Interval & [`coef_sleep_s_lee_l_l', `coef_sleep_s_lee_u_u'] & [`coef_inbed_s_lee_l_l', `coef_inbed_s_lee_u_u'] & [`coef_eff_s_lee_l_l', `coef_eff_s_lee_u_u'] & [`coef_nap_s_lee_l_l', `coef_nap_s_lee_u_u'] & [`coef_all_s_lee_l_l', `coef_all_s_lee_u_u'] \\" _n ///
		"\addlinespace" _n ///
		///
		"\hspace{.1cm} \textbf{Devices+Incentives} & `coef_Sleep_Night_ns2' & `coef_act_inbed_ns2' & `coef_sleep_eff_ns2' & `coef_nap_hr_ns2' & `coef_sleep_all_ns2' \\" _n ///
		" & (`se_Sleep_Night_ns2') & (`se_act_inbed_ns2') & (`se_sleep_eff_ns2') & (`se_nap_hr_ns2') & (`se_sleep_all_ns2') \\" _n ///	
		"\hspace{.2cm} Lee Lower Bound & `coef_sleep_s_i_lee_l' & `coef_inbed_s_i_lee_l' & `coef_eff_s_i_lee_l' & `coef_nap_s_i_lee_l' & `coef_all_s_i_lee_l' \\" _n ///
		"\hspace{.2cm} Lee Upper Bound & `coef_sleep_s_i_lee_u' & `coef_inbed_s_i_lee_u' & `coef_eff_s_i_lee_u' & `coef_nap_s_i_lee_u' & `coef_all_s_i_lee_u' \\" _n ///
		"\hspace{.2cm} Confidence Interval & [`coef_sleep_s_i_lee_l_l', `coef_sleep_s_i_lee_u_u'] & [`coef_inbed_s_i_lee_l_l', `coef_inbed_s_i_lee_u_u'] & [`coef_eff_s_i_lee_l_l', `coef_eff_s_i_lee_u_u'] & [`coef_nap_s_i_lee_l_l', `coef_nap_s_i_lee_u_u'] & [`coef_all_s_i_lee_l_l', `coef_all_s_i_lee_u_u'] \\" _n ///
		"\addlinespace" _n ///
		"\textbf{Nap Treatment} & `coef_Sleep_Night_nap' & `coef_act_inbed_nap' & `coef_sleep_eff_nap' & `coef_nap_hr_nap' & `coef_sleep_all_nap' \\" _n ///
		" & (`se_Sleep_Night_nap') & (`se_act_inbed_nap') & (`se_sleep_eff_nap') & (`se_nap_hr_nap') & (`se_sleep_all_nap') \\" _n ///	
		"\hspace{.2cm} Lee Lower Bound & `coef_sleep_nap_lee_l' & `coef_inbed_nap_lee_l' & `coef_eff_nap_lee_l' & `coef_nap_nap_lee_l' & `coef_all_nap_lee_l' \\" _n ///
		"\hspace{.2cm} Lee Upper Bound & `coef_sleep_nap_lee_u' & `coef_inbed_nap_lee_u' & `coef_eff_nap_lee_u' & `coef_nap_nap_lee_u' & `coef_all_nap_lee_u' \\" _n ///
		"\hspace{.2cm} Confidence Interval & [`coef_sleep_nap_lee_l_l', `coef_sleep_nap_lee_u_u'] & [`coef_inbed_nap_lee_l_l', `coef_inbed_nap_lee_u_u'] & [`coef_eff_nap_lee_l_l', `coef_eff_nap_lee_u_u'] & [`coef_nap_nap_lee_l_l', `coef_nap_nap_lee_u_u'] & [`coef_all_nap_lee_l_l', `coef_all_nap_lee_u_u'] \\" _n ///
		"\addlinespace" _n ///
		///
		"\midrule" _n ///			
		"Control Mean               & `mean_Sleep_Night' & `mean_act_inbed' & `mean_sleep_eff' & `mean_nap_hr' & `mean_sleep_all'  \\" _n ///
		"Control SD                & `sd_Sleep_Night' & `sd_act_inbed' & `sd_sleep_eff' & `sd_nap_hr' & `sd_sleep_all' \\" _n ///
		"\bottomrule" _n ///
		"\end{tabular}" _n ///
		"\end{tabular}" _n
	file close f
