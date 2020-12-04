
	* Use real data.
	
	cls	
	capture clear
	set scheme lean3 /*or lean3*/
	cd "C:\Users\j.diong\Desktop\stats"
	*cd "C:\Users\r.herbert\Documents\Stata\Joanna Diong"
	import delimited using "subjects_data.csv", delimiters(",") varnames(1)

	/* Fit FP with robust SEs and extract predicted values. Show plot of fit.
	
	mfp: regress emgso activations, cluster(subject)
	fracplot, title("") xtitle("Activation (%MVC)") ytitle("EMG SO (V)") aspect(1) name(fplot, replace)
	  
	predict fracpred
	predict fracSE, stdp
	generate LCL = fracpred - invnormal(0.975) * fracSE
	generate UCL = fracpred + invnormal(0.975) * fracSE	
	di "`graphtext1'"
	sort activations
	twoway `graphtext1' ///
		(rarea LCL UCL activations, astyle(sunflowerlb)) ///
		(line fracpred activations, lpattern(solid) lcolor(black) lwidth(medthin)) ///
		, ///
		xtitle("Activation (%MVC)") ytitle("EMG SO (V)") legend(off) aspect(1) legend(off) name(rawAndFitted, replace)
	graph combine fplot raw rawAndFitted, cols(1) xcommon ycommon scale(0.6) */

	* Plot individual subject data with regression lines.
	
	egen subjectGroup = group(subject)
	sort activations
	forvalues subject = 1/7 {
		twoway(line emgso activations if subjectGroup == `subject', lpattern(solid)) ///
		(lfit emgso activations if subjectGroup == `subject', lpattern(solid) lcolor(red)) ///
		, ///
		legend(off) name(eachSubject`subject', replace) 
	}	
	graph combine eachSubject1 eachSubject2 eachSubject3 eachSubject4 eachSubject5 eachSubject6 eachSubject7, ///
		xcommon ycommon name(combined, replace)
	graph drop eachSubject1 eachSubject2 eachSubject3 eachSubject4 eachSubject5 eachSubject6 eachSubject7

	* Plot individual subject data on one plot.
	
	forvalues subject = 1/7 {
		local graphtext2 `graphtext2' (line emgso activations if subjectGroup == `subject')
	}	
	sort activations
	twoway `graphtext2' ///
		, ///
		xtitle("Activation (%MVC)") ytitle("EMG SO (V)") legend(off) aspect(1) legend(off) name(raw, replace)
		
		* Fit a linear mixed model, forced through the origin, with random slopes, and graph the fitted data.
		
	cls	
	mixed emgso activations, noconstant || subjectGroup: subjectGroup	
	quietly {
		predict mixedPred
		predict mixedSE, stdp
		generate mixedLCL = mixedPred - invnormal(0.975) * mixedSE
		generate mixedUCL = mixedPred + invnormal(0.975) * mixedSE	
		forvalues subject = 1/7 {
			local graphtext3 `graphtext3' (line emgso activations if subjectGroup == `subject')
		}	
		sort activations
		twoway `graphtext3' ///
			(rarea mixedLCL mixedUCL activations, astyle(sunflowerlb)) ///
			(line mixedPred activations, lpattern(solid)) ///
			, ///
			xtitle("Activation (%MVC)") ytitle("EMG SO (V)") legend(off) aspect(1) legend(off) name(mixed, replace)
	
		* List model predictions.
		
	preserve
		drop activations
		egen activations = seq() in 1/11, from(0)
		replace  activations = activations * 10
		predict newPred
		predict newSE, stdp
		generate newLCL = newPred - invnormal(0.975) * newSE
		generate newUCL = newPred + invnormal(0.975) * newSE	
		noisily list activations newPred newLCL newUCL in 1/11, noobs sep(0)
	restore
		
	}
		