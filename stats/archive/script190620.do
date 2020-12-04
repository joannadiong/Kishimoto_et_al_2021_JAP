/*
	* Set up.

	capture clear
	*set scheme lean1
	set obs 60

	* Make up some data.

	egen id = seq(), block(10)
	egen Act = seq(), from(1) to(10)
	replace Act = 0.1 * Act  // put actimVation on a 0-1 scale, to make coefficients manageable

	generate mV = .
	forvalues id = 1/10 {
		local a = 3 + runiform()
		local b  = 4 + runiform()
		replace mV = `a' * exp(`b' * Act) - 1 * rnormal() if id == `id'
	}
	replace mV = mV + 10 * rnormal()
	replace mV = 0 if mV < 0 
	list, sepby(id)

	* Fit a fractional polynomial model with robust SEs. Draw a simple graph.

	mfp: regress mV Act, cluster(id)
	fracplot, title("") xtitle("Activation") ytitle("mV") aspect(1)

	* Draw a fancier graph.
	
	mfp: regress mV Act, cluster(id)
	predict fracpred
	predict fracSE, stdp
	generate LCL = fracpred - invnormal(0.975) * fracSE
	generate UCL = fracpred + invnormal(0.975) * fracSE	
	sort Act
	forvalues id = 1/6 {
		local graphtext1 `graphtext1' (line mV Act if id == `id', lpattern(solid) lcolor(gs12))
	}	
	twoway `graphtext1' ///
		(rarea LCL UCL Act, astyle(sunflowerlb)) ///
		(line fracpred Act, lpattern(solid) lcolor(black) lwidth(medthin)) ///
		, ///
		xtitle("Activation") ytitle("mV") legend(off) aspect(1) legend(off) name(fpandlines, replace)	
*/


	* use real data
	capture clear
	set scheme lean1 /*or lean3*/

	cd "C:\Users\j.diong\Desktop\stats"
	import delimited using "subjects_data.csv", delimiters(",") varnames(1)

	mfp: regress emgso activations, cluster(subject)
	fracplot, title("") xtitle("Activation (%MVC)") ytitle("EMG SO (V)") aspect(1)
	  
	predict fracpred
	predict fracSE, stdp
	generate LCL = fracpred - invnormal(0.975) * fracSE
	generate UCL = fracpred + invnormal(0.975) * fracSE	
	sort activations

	foreach subject in subject { /*for loop needs fixing*/
		local graphtext1 `graphtext1' (line emgso activations if subject == `subject', lpattern(solid) lcolor(gs12))
	}	
	twoway `graphtext1' ///
		(rarea LCL UCL activations, astyle(sunflowerlb)) ///
		(line fracpred activations, lpattern(solid) lcolor(black) lwidth(medthin)) ///
		, ///
		xtitle("Activation (%MVC)") ytitle("EMG SO (V)") legend(off) aspect(1) legend(off) name(fpandlines, replace)	
