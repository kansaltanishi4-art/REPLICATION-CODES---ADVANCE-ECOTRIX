************************
* Anderson Index program
************************

cap program drop anderson_index

program anderson_index

* Arguments: components and covariates
args components treatments controls factorcovariates subset eststo

{
cap drop weight*
cap drop summary_index
cap drop summary_index_base

di "Outcomes must have sign switched before running program!"

* Step 1: Standardize variables
	foreach var in `components' {
	di "`var'"
			
		* Create standardized versions of postline variable
		di "Post"
		cap drop `var'_post_std
				
		sum `var' if treat_pool == 0 & treat_nap == 0 & post_treatment==1
		gen `var'_post_std = (`var' - r(mean)) / r(sd)

		* Create standardized versions of baseline variable
		di "Base"
		cap drop `var'_pre_std
		
		bysort pid: egen `var'_pre_temp = mean(`var') if post_treatment == 0
		bysort pid: egen `var'_pre = mean(`var'_pre_temp)
		sum `var'_pre if treat_pool == 0 & treat_nap == 0 & post_treatment==0
		gen `var'_pre_std = (`var'_pre - r(mean)) /r(sd)
		drop `var'_pre `var'_pre_temp
	}

* Step 2: Run regressions with controls
		di "Got to step 2"
		
		local count = 1
		
		cap drop resid*
	
		* Index 1
		foreach var in `components' {
		
			mdesc `var'_pre_std
			* this if-else condition captures whether ALL observations for the baseline outcome are missing
			if r(total)!=r(miss) {
			
				*reghdfe `var' `covariates' average_index_baseline `subset', absorb(`factorcovariates') res(resid_`count')
				reghdfe `var'_post_std `covariates' `var'_pre_std `subset', absorb(`factorcovariates') res(resid_`count')
				
			}
			else {
			
				*reghdfe `var' `covariates' average_index_baseline `subset', absorb(`factorcovariates') res(resid_`count')
				reghdfe `var'_post_std `covariates' `subset', absorb(`factorcovariates') res(resid_`count')
			
			}
		
			* Standardize residuals
			sum resid_`count' if post_treatment==1
			gen resid_`count'_std = (resid_`count' - r(mean)) / r(sd)
			
			local count = `count' + 1
		
		}

* Step 3: Get efficient matrix
		di "Got to step 3"
	
		cor resid*_std, cov
		mat sigma = r(C)
		
		mat sigma_inv = inv(sigma)
		
		* create variables with matrix components
		local n: word count `components'		
		forval i = 1/`n'{
			forval j = 1/`n'{
				gen corr`i'_`j' = sigma_inv[`i',`j']
			}
		}
	
* Step 4: Calculate standardized index accounting for missings
	di "Got to step 4"
	
		local i = 1
		foreach var in `components' {
		
			forval j = 1/`n'{
				 replace corr`i'_`j' = . if `var'_post_std == . & post_treatment == 1
			}
			
			local i = `i' + 1		
		}	
		egen denominator = rowtotal(corr*) if post_treatment == 1
		cap drop nomiss
		egen nomiss = rownonmiss(corr*)
		replace denominator = . if nomiss == 0
		
		local i = 1
		foreach var in `components' {
		
			egen numerator_`i' = rowtotal(corr`i'_*) if post_treatment == 1
			gen weighted_`i' = `var'_post_std*(numerator_`i'/denominator) if post_treatment == 1
			gen weight_`i' = numerator_`i'/denominator
		
			local i = `i' + 1
		}
		
		egen summary_index = rowtotal(weighted*)
		replace summary_index = . if nomiss == 0
		
* Step 5: Calculate baseline standardized index accounting for missings
	di "Got to step 5"

		local i = 1
		foreach var in `components' {
		
			forval j = 1/`n'{
				 replace corr`i'_`j' = . if `var'_pre_std == . & post_treatment == 0
			}
			
			local i = `i' + 1		
		}	
		egen denominator_base = rowtotal(corr*) if post_treatment == 0
		*nomiss defined above
		replace denominator_base = . if nomiss == 0

		local i = 1
		foreach var in `components' {
		 
			egen numerator_base_`i' = rowtotal(corr`i'_*) if post_treatment == 0
			gen weighted_base_`i' = (numerator_base_`i'/denominator_base)*`var'_pre_std if post_treatment == 0
		
			local i = `i' + 1
		}
		
		egen _summary_index_base = rowtotal(weighted*) if post_treatment == 0
		replace _summary_index_base = . if nomiss == 0

		bys pid: egen summary_index_base = mean(_summary_index_base)
		drop _summary_index_base
	
* Step 6: Regress on standardized index
	
		if `eststo' == 1 {
			eststo: reghdfe summary_index `treatments' `controls' summary_index_base `subset', cluster(pid) absorb(`factorcovariates')
		}
		else{
			reghdfe summary_index `treatments' `controls' summary_index_base `subset', cluster(pid) absorb(`factorcovariates')
		}
		
drop corr* 
drop denominator
drop numerator_*
drop resid*
drop weighted_*
drop denominator_base
*drop weight*
*drop summary_index
*drop summary_index_base

}
end
