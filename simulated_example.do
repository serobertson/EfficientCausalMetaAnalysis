 

clear all

*Programs for estimators here
******Start of program here*******
program estimatorcalc_tx_constant, rclass
args a 

generate term1 = (S == 0) * (g_`a')

generate w_`a'= ((1-p)/(p*e_a`a'))
replace w_`a'=0 if w_`a'==.

generate term2 = (S != 0) *(A==`a') * w_`a' * (Y - g_`a')
replace term2=0 if term2==.

generate summand = term1 + term2 

summ summand
	local sumnum = r(sum)
	
count if S == 0
	local S0_size = r(N) 
	
count if S == 0 | S==1
	local total_size = r(N) 
	

*local estimate =  `sumnum' / `S0_size'
return local ESTIMATE_AIPW=`sumnum'/`S0_size'  /*EQUATION 10 in the draft*/


*IPW
generate term2_mod = (S != 0) *(A==`a') * w_`a' * (Y)
summ term2_mod
local term2_mod=r(sum)
return local ESTIMATE_IPW=`term2_mod'/`S0_size'  /*EQUATION 9 in the draft*/

*OM
summ term1
local term1=r(sum)
return local ESTIMATE_OM=`term1'/`S0_size' 


*Normalized IPW
generate normalization = (S != 0) * (A==`a') * w_`a'
	sum normalization
	local norm = 1 / r(sum) 
return local ESTIMATE_IPW2=`norm'*`term2_mod'

*Normalized DR
summ term2
local term2=r(sum)
return local ESTIMATE_AIPW2=`norm'*`term2' + `term1'/`S0_size' 




drop term*
drop summand
drop w*
drop normalization

end
******End of program here*******



*start simulation

set seed 12345678

*local N_runs = 10000
local N_runs = 1000

local scenario = 0

foreach total_size in  10000 /* 100000 */ {

*in_trial prevalence set to 20% and 50%
if `total_size' == 10000 {
	local trial_size_list "1000 2000 5000"
	}
	
*in_trial prevalence set to 2% and 5%
if `total_size' == 100000 {
	local trial_size_list "1000 2000 5000"
	}
	
foreach trial_size in `trial_size_list' {
 		
if `trial_size' == 5000 & `total_size' == 10000 {
	local intercept = 0
	}
	
if `trial_size' == 2000 & `total_size' == 10000 {
	local intercept =-2.036133 
	}
	
if `trial_size' == 1000 & `total_size' == 10000 {
	local intercept =-3.1289
	}
	
	
if `trial_size' == 5000 & `total_size' == 100000 {
	local intercept = -4.066406
	}
	
if `trial_size' == 2000 & `total_size' == 100000 {
	local intercept = -5.156250 
	}
	
if `trial_size' == 1000 & `total_size' == 100000 {
	local intercept =-5.9375
	}


forvalues scenario_trial_size = 1/2 {
	*set intercepts for trial prevalence, e.g. 1:1 trial, or 1:3
	
*1000
		if `scenario_trial_size' == 1 & `trial_size' == 1000 & `total_size' == 10000  {
			* ~ 1:1 trial size
				local trial_1_int = -1.388634
				local trial_2_int = 0.6610838
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 1000 & `total_size' == 10000 {
			* ~ 4:2:1 trial size
				local trial_1_int = -2.049594
				local trial_2_int = -0.9157119
	
			}
			
if `scenario_trial_size' == 1 & `trial_size' == 1000 &  `total_size' == 100000  {
			* ~ 1:1 trial size
				local trial_1_int = -1.829346
				local trial_2_int = 0.8972254
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 1000 & `total_size' == 100000 {
			*~ 4:2:1 trial size
				local trial_1_int = -2.477825
				local trial_2_int = -0.5427084
	
			}
		
*2000
		if `scenario_trial_size' == 1 & `trial_size' == 2000 & `total_size' == 10000  {
			* ~ 1:1 trial size
				local trial_1_int = -1.171908
				local trial_2_int = 0.53125
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 2000 & `total_size' == 10000 {
			* ~ 4:2:1 trial size
				local trial_1_int = -1.832412
				local trial_2_int = -0.9157119
	
			}
			
if `scenario_trial_size' == 1 & `trial_size' == 2000 &  `total_size' == 100000  {
			* ~ 1:1 trial size
				local trial_1_int = -1.729122
				local trial_2_int = 0.8636475
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 2000 & `total_size' == 100000 {
			*~ 4:2:1 trial size
				local trial_1_int = -2.373901
				local trial_2_int = -0.5703125
	
			}
*5000			
		if `scenario_trial_size' == 1 & `trial_size' == 5000 & `total_size' == 10000  {
			* ~ 1:1 trial size
				local trial_1_int = -0.7692247
				local trial_2_int = 0.234375
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 5000 & `total_size' == 10000 {
			* ~ 4:2:1 trial size
				local trial_1_int = -1.4375
				local trial_2_int = -1.218262
	
			}
			
		if `scenario_trial_size' == 1 & `trial_size' == 5000 & `total_size' == 100000   {
			* ~ 1:1 trial size
				local trial_1_int =  -1.562623
				local trial_2_int =  0.7578125
			}
			
		if `scenario_trial_size' == 2 & `trial_size' == 5000 & `total_size' == 100000  {
			*~ 4:2:1 trial size
				local trial_1_int = -2.217345
				local trial_2_int = -0.671875
	
			}
			
forvalues scenario_treatment_prev = 1/2 {
	*set intercepts for treatment prev, e.g. 1:1 , or 1:3
		
		if `scenario_treatment_prev' == 1 {
				local trial_1_tx_prev = 0.5
				local trial_2_tx_prev = 0.5
				local trial_3_tx_prev = 0.5
			}
			
		if `scenario_treatment_prev' == 2 {
			local trial_1_tx_prev = 1/2
			local trial_2_tx_prev = 1/3
			local trial_3_tx_prev = 2/3
	
			}	
			
local ++scenario	

		
tempname saved_results 



postfile `saved_results' double(total_size trial_size S1_size S2_size S3_size scenario_trial_size scenario_treatment_prev ///
								Y1_true Y0_true ///
								mu1_OM mu1_IPW1 mu1_IPW2 mu1_AIPW1 mu1_AIPW2 ///
								mu0_OM mu0_IPW1 mu0_IPW2 mu0_AIPW1 mu0_AIPW2 ///
								) using multi_sim_test_`scenario'_N_`total_size'.dta, replace 

								
								
qui forvalues i = 1/`N_runs'{

di in red "`i'"

clear

set obs `total_size'

matrix C = ( 1, 0.5, 0.5  ///
			\ 0.5, 	1, 0.5  ///
			\ 0.5, 0.5, 1 )
	
						
drawnorm X1 X2 X3, corr(C)

generate in_trial = runiform() < invlogit(`intercept'+ ln(2)* X1 + ln(2) * X2 + ln(2) * X3 )

gen xb1 = `trial_1_int' + ln(1.5)*X1 + ln(1.5)*X2  + ln(1.5)*X3  if in_trial == 1
gen xb2 = `trial_2_int' + ln(0.75)*X1 + ln(0.75)*X2 + ln(0.75)*X3  if in_trial == 1

gen p1 = 1        / ( 1 + exp(xb1) + exp(xb2) )  if in_trial == 1
gen p2 = exp(xb1) / ( 1 + exp(xb1) + exp(xb2) )  if in_trial == 1
gen p3 = exp(xb2) / ( 1 + exp(xb1) + exp(xb2) )  if in_trial == 1
 
gen UNIF = runiform()  if in_trial == 1

gen S = cond(UNIF < p1, 1, ///
		cond(UNIF < p1 + p2, 2, 3 ))  if in_trial == 1



replace S = 0 if in_trial != 1

tab S



*save trial sizes
count if S == 1
	local S1_size = r(N) 
count if S == 2
	local S2_size = r(N) 
count if S == 3
	local S3_size = r(N) 


* generate treatment in the trials
generate A = rbinomial(1, `trial_1_tx_prev') if S == 1
	replace A = rbinomial(1,  `trial_2_tx_prev') if S == 2
	replace A = rbinomial(1,  `trial_3_tx_prev') if S == 3



* generate outcome
generate Y1 = 0.5 - X1 - X2 -  X3  + rnormal()
generate Y0 = 1.5 + X1 + X2 + X3  + rnormal()

generate Y = A * Y1 + ( 1 - A ) * Y0
summ Y1

summ Y0



*******************************

*indicator for if target (S==0 , not in trials is the target*/ 
generate T=1 if S==0
replace T=0 if S!=0

*variable list for regression
local X X1 X2 X3


*Pr[IN_TRIAL==1|X]
*Pr[R=1|X]
generate IN_TRIAL=(1-T)
logit IN_TRIAL `X'
	predict p
	
	
*E[Y|X, IN_TRIAL=1, A=a]
forvalues arm = 0/1 {
regress Y `X' if IN_TRIAL==1 & A==`arm'
predict g_`arm'  /*g_a(X)*/
}

*Pr[A=a|X, IN_TRIAL=1]
logit A `X' if IN_TRIAL==1
predict e_a1

generate e_a0= 1- e_a1

*save true estimates
summ Y1 if S==0
	local Y1_true=r(mean)

summ Y0 if S==0
	local Y0_true=r(mean)

*Run program with constant tx

estimatorcalc_tx_constant 1 
return list 
	local mu1_OM = r(ESTIMATE_OM) 
	local mu1_IPW1 = r(ESTIMATE_IPW)
	local mu1_IPW2 = r(ESTIMATE_IPW2)
	local mu1_AIPW1 = r(ESTIMATE_AIPW)
	local mu1_AIPW2 = r(ESTIMATE_AIPW2)


estimatorcalc_tx_constant 0
return list 
    local mu0_OM = r(ESTIMATE_OM) 
	local mu0_IPW1 = r(ESTIMATE_IPW)
	local mu0_IPW2 = r(ESTIMATE_IPW2)
	local mu0_AIPW1 = r(ESTIMATE_AIPW)
	local mu0_AIPW2 = r(ESTIMATE_AIPW2)


/* store results */															
post `saved_results' 	(`total_size') (`trial_size') (`S1_size') (`S2_size') (`S3_size') (`scenario_trial_size') (`scenario_treatment_prev') ///
						(`Y1_true') (`Y0_true') ///
						(`mu1_OM') (`mu1_IPW1') (`mu1_IPW2') (`mu1_AIPW1') (`mu1_AIPW2') ///
						(`mu0_OM') (`mu0_IPW1') (`mu0_IPW2') (`mu0_AIPW1') (`mu0_AIPW2') 
} /* end of runs for a given scenario */
postclose `saved_results'

}
}
}
}	
