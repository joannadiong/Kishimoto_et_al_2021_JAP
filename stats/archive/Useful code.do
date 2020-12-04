
** COMPLIANCE ANALYSIS 2011

	** Exponential curves Torque, Angle in Control, SCI
	generate minangle = .
	label variable minangle "Minimum angle"
	generate maxangle = .
	label variable maxangle "Maximum angle"
	
	foreach subject of numlist 1/35 { /* Subjects: NOR = 1-20, SCI = 21-35 */
	  summarize angle if id == `subject' & !missing(angle)
	  local mina = r(min)
	  local maxa = r(max)
	  *egen anglerange`subject' = seq()if id == `subject', from(`mina') to(`maxa') // from() incorrectly specified
	  range anglerange`subject' `mina' `maxa' // generates columns angles within range
	  replace minangle = `mina' if id == `subject'
	  replace maxangle = `maxa' if id == `subject'
	  
	  quietly nl (torque = {A} * exp({B}*(angle-{C}))) if id == `subject' & !missing(angle), vce(cluster code)
	  quietly predict predictedtorque`subject' if id == `subject' & !missing(angle) // generates column torques
	  summarize predictedtorque`subject' if id == `subject' & !missing(predictedtorque`subject')
	  local mint = r(min)
	  local maxt = r(max)
	  *egen torquerange`subject' = seq()if id == `subject', from(`mint') to(`maxt')
	  range torquerange`subject' `mint' `maxt' // generates columns torques within range
	  *drop predictedtorque`subject'
	  
      matrix parameters = e(b) // subject-specific command execution 
	  matrix list parameters 
	  local A`subject' = parameters[1,1]
	  display "`A`subject''"
	  local B`subject' = parameters[1,2] 
	  display "`B`subject''" 
	  local C`subject' = parameters[1,3] 
	  display "`C`subject''" 
	  	        	  	  
	  sort anglerange`subject'
	  if `subject' <= 20 {
		local plotlinenor  "`plotlinenor' plot`subject'(lcolor(blue) lpattern(solid)) "
	  }		 
	  else {
		local plotlinesci "`plotlinesci' plot`subject'(lcolor(red) lpattern(solid)) "
	  }
	  *twoway function (torquerange`subject' = `A`subject'' * exp(`B`subject'' * (anglerange`subject' - `C`subject'')))
	} 


	** Remove outlier lfas using standardised residuals
	generate stdres = .
	foreach subject of numlist 1/35 {
      quietly regress lfas hipangle if id == `subject' & !missing(lfas)
	  quietly predict stdres`subject' if id == `subject' & !missing(lfas), rstandard 
	  quietly replace stdres = stdres`subject' if id == `subject' & !missing(lfas)
	  drop stdres`subject'
	}
		
	replace lfas = . if stdres >= 4 
		
	** Regress lfas on hip angle for each subject, obtain parameters
	generate slope = .
	label variable slope "Slope"
	generate intercept = .
	label variable intercept "Intercept"
	generate minhipangle = .
	label variable minhipangle "Min hip angle"
	generate maxhipangle = .
	label variable maxhipangle "Max hip angle"
	
	foreach subject of numlist 1/35 {
	  ** Extract slope and intercept
	  quietly regress lfas hipangle if id == `subject' & !missing(lfas)
	  local sl = _b[hipangle]
	  quietly replace slope = `sl' if id == `subject'
	  local in = _b[_cons]
	  quietly replace intercept = `in' if id == `subject'
	  
	  ** Extract min, max hip angle
	  quietly summarize hipangle if id == `subject' & !missing(hipangle)
	  local min = r(min)
	  quietly replace minhipangle = `min' if id == `subject'
	  local max = r(max)
	  quietly replace maxhipangle = `max' if id == `subject'
	}

	** Sensitivity analysis: T-test, regression with outlier removed (ID = SCI 4)
	ttest slope if id != 24, by(group)	
        regress slope sci if id != 24 // Report this result
	estimates store scils
	regress slope sci sex if id != 24
	estimates store scilssex
	display in green _n "Test if sex contributes to model ..."
	lrtest scilssex scils
	


** CONTRACTURE XSEC ANALYSIS 2010
  
  ** Unadjusted means, Gender and ls - adjusted differences, NOR and STR
  foreach vname of varlist lgfsl lgF1 alphag lnalphag ffsl lffsl ltfsl lgF100 lfF100 ltF100 chlg100 chlf100 chlt100 cat100 lgstrain100 lfstrain100 ltstrain100 {
    display _n _n _n in green "ANALYSIS OF: " "`vname'"
	local tableoptions row center format(%12.1fc)
	if "`vname'" == "lnalphag" local tableoptions row center format(%12.3fc)
	if "`vname'" == "ffsl" local tableoptions row center format(%12.0fc)
    table group if group == "NOR" | group == "STR",                                              ///
                contents(mean `vname' sd `vname' count `vname' min `vname' max `vname') `tableoptions'
    regress `vname' stroke ls if group == "NOR" | group == "STR"
    estimates store strls
    quietly regress `vname' stroke ls age if group == "NOR" | group == "STR"
    estimates store strlsage
    quietly regress `vname' stroke ls sex if group == "NOR" | group == "STR"
    estimates store strlssex
    if "`vname'" != "ffsl" quietly regress `vname' stroke ls ffsl if group == "NOR" | group == "STR"
    if "`vname'" != "ffsl" estimates store strlsfsl
    display in green _n "Test if gender contributes to model ..."
    lrtest strlssex strls
    display in green _n "Test if age contributes to model ..."
    lrtest strlsage strls
    if "`vname'" != "ffsl" display in green _n "Test if tension at fascicle slack length contributes to model ..."
    if "`vname'" != "ffsl" lrtest strlsfsl strls
	*There was a significant effect of gender for three variable, so run additional regressions with gender as a predictor.
	if "`vname'" == "lgF100" | "`vname'" == "ltfsl" | "`vname'" == "chlt100" regress `vname' stroke ls sex
  }

 ** Adjusted and unadjusted means, ls - adjusted differences, NOR and SCI
  foreach vname of varlist lgF1 alphag lnalphag ffsl lgfsl lffsl ltfsl lgF100 lfF100 ltF100 chlf100 chlt100 cat100 lgstrain100 lfstrain100 ltstrain100 {
    display _n _n _n in green "ANALYSIS OF: " "`vname'"
	local tableoptions row center format(%12.1fc)
	if "`vname'" == "lnalphag" local tableoptions row center format(%12.3fc)
	if "`vname'" == "ffsl" local tableoptions row center format(%12.0fc)
    table group if group == "NOR" | group == "SCI",                                              ///
                contents(mean `vname' sd `vname' count `vname' min `vname' max `vname') `tableoptions'
   	** Obtain group mean of ls
	quietly summarize ls
	local meanls = r(mean)
	local sdls = r(sd)
	** Linear regression for each vname adjusted for ls
	regress `vname' SCI ls if group == "NOR" | group == "SCI"
    estimates store scils
	** Create ls adjusted means for each vname
	matrix parameters = e(b) // variable-specific regression coefficients
	matrix list parameters
	  local b_SCI = parameters[1,1]
	  display "`b_SCI'"
	  local b_ls = parameters[1,2] 
	  display "`b_ls'" 
	  local _cons = parameters[1,3] 
	  display "`_cons'" 
	    quietly lincom `_cons' + (`b_ls' * `meanls') // ls adjusted NOR mean
	    local NOR = r(estimate)
	    display _n in green "Adjusted NOR mean of `vname' = `NOR'"
	    quietly lincom `_cons' + (`b_ls' * `meanls') + (`b_SCI') // ls adjusted SCI mean
	    local SCI = r(estimate)
	    display in green "Adjusted SCI mean of `vname' = `SCI'"
	** Continue likelihood ratio tests for age, sex effects for each vname
    quietly regress `vname' SCI ls age if group == "NOR" | group == "SCI"
    estimates store scilsage
    quietly regress `vname' SCI ls sex if group == "NOR" | group == "SCI"
    estimates store scilssex
    if "`vname'" != "ffsl" quietly regress `vname' SCI ls ffsl if group == "NOR" | group == "SCI"
    if "`vname'" != "ffsl" estimates store scilsfsl
    display in green _n "Test if gender contributes to model ..."
    lrtest scilssex scils
    display in green _n "Test if age contributes to model ..."
    lrtest scilsage scils
    if "`vname'" != "ffsl" display in green _n "Test if tension at fascicle slack length contributes to model ..."
    if "`vname'" != "ffsl" lrtest scilsfsl scils
  }


** CONTRACTURE XSEC GRAPHS (ROB'S CODE) 2010

    display in red _n _n "** NB: MUST RUN Contracture.do BEFORE RUNNING THIS DO FILE **" _n _n


    ********** GRAPH l-t CURVES ********** 

  cd "C:\Documents and Settings\rherbert\My Documents\Research projects\ARC DP 2007-9 (contracture)\Contracture analysis\Graph dump"

  set obs 100
  egen Tg = seq()
  replace lG = lG / 100
  replace ls = ls / 100
  local i = 1
  while `i' < 100 {
**local e0gi = e0g in `i'
    local alphagi = alphag in `i'
    local lGi = lG in `i'
    local groupi = group in `i'
    local codestringi = code in `i'
    local lsi = ls in `i'
    local maxTgi = maxTg in `i'
    if `alphagi' != . {
      local vname = "Lg" + "`codestringi'" + "`groupi'"
      local vnametomax = "Lg" + "`codestringi'" + "`groupi'" + "tomax"
      quietly generate `vname' =  (ln(Tg) - ln(1/`alphagi') + `alphagi'*`lGi') / `alphagi' * 100 / `lsi'
      quietly generate `vnametomax' = `vname' if Tg <= `maxTgi'
      if "`groupi'" == "NOR" local NORnames = "`NORnames'" + `" "' + "`vname'"
      if "`groupi'" == "NOR" local NORgraphs : display `"`NORgraphs'"' `"(line Tg "' `"`vname'"' `"tomax, lpattern(solid) lcolor(blue)) "'
      if "`groupi'" == "STR" local STRnames = "`STRnames'" + `" "' + "`vname'"
      if "`groupi'" == "STR" local STRgraphs : display `"`STRgraphs'"' `"(line Tg "' `"`vname'"' `"tomax, lpattern(solid) lcolor(red)) "'
      if "`groupi'" == "SCI" local SCInames = "`SCInames'" + `" "' + "`vname'"
      if "`groupi'" == "SCI" local SCIgraphs : display `"`SCIgraphs'"' `"(line Tg "' `"`vname'"' `"tomax, lpattern(solid) lcolor(red)) "'
    }
    local i = `i' + 1
  }

  replace lG = lG * 100
  replace ls = ls * 100

  egen lmeannorm = rowmean(`NORnames')
  egen lsdnorm = rowsd(`NORnames')
  generate llcinorm = lmeannorm - 1.96 * lsdnorm / sqrt(16)
  generate lucinorm = lmeannorm + 1.96 * lsdnorm / sqrt(16)
  egen lmeanstr = rowmean(`STRnames')
  egen lsdstr = rowsd(`STRnames')
  generate llcistr = lmeanstr - 1.96 * lsdstr / sqrt(10)
  generate lucistr = lmeanstr + 1.96 * lsdstr / sqrt(10)
  egen lmeansci = rowmean(`SCInames')
  egen lsdsci = rowsd(`SCInames')
  generate llcisci = lmeansci - 1.96 * lsdsci / sqrt(16)
  generate lucisci = lmeansci + 1.96 * lsdsci / sqrt(16)

    ** Graph length-tension curves for controls.

  twoway `NORgraphs',                                                                            ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1) label(1 "Control")) name(CONltcurves, replace)
  graph export CONltcurves.wmf, replace
  graph drop CONltcurves
  
    ** Graph length-tension curves for stroke patients and controls.

  twoway `NORgraphs' `STRgraphs',                                                                ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1 35) label(1 "Control") label(35 "Stroke")) name(STRCONltcurves, replace)
  graph export STRCONltcurves.wmf, replace
  graph drop STRCONltcurves

    ** Graph length-tension curves for SCI patients and controls.

  twoway `NORgraphs' `SCIgraphs',                                                                ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1 35) label(1 "Control") label(35 "SCI")) name(SCICONltcurves, replace)
  graph export SCICONltcurves.wmf, replace
  graph drop SCICONltcurves

    ** Graph summary data (95% CI about mean) for controls.

  twoway (line Tg llcinorm, lpattern(solid) lcolor(blue))                                        ///
         (line Tg lucinorm, lpattern(solid) lcolor(blue)),                                       ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1) label(1 "Control")) name(CONltcurvesCI, replace)
  graph export CONltcurvesCI.wmf, replace
  graph drop CONltcurvesCI

    ** Graph summary data (95% CI about mean) for stroke patients and controls.

  twoway (line Tg llcinorm, lpattern(solid) lcolor(blue))                                        ///
         (line Tg lucinorm, lpattern(solid) lcolor(blue))                                        ///
         (line Tg llcistr, lpattern(solid) lcolor(red))                                          ///
         (line Tg lucistr, lpattern(solid) lcolor(red)),                                         ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1 3) label(1 "Control") label(3 "Stroke")) name(STRCONltcurvesCI, replace)
  graph export STRCONltcurvesCI.wmf, replace
  graph drop STRCONltcurvesCI

    ** Graph summary data (95% CI about mean) for SCI patients and controls.

  twoway (line Tg llcinorm, lpattern(solid) lcolor(blue))                                        ///
         (line Tg lucinorm, lpattern(solid) lcolor(blue))                                        ///
         (line Tg llcisci, lpattern(solid) lcolor(red))                                          ///
         (line Tg lucisci, lpattern(solid) lcolor(red)),                                         ///
         xtitle("Length (% of leg length)") ysc(r(0 100)) xlabel(40(20)120, nogrid)              ///
         ytitle("Tension (N)") ysc(r(0 100)) ylabel(0(20)100, nogrid)                            ///
         scale(1.5) legend(ring(0) position(11) order(1 3) label(1 "Control") label(3 "SCI")) name(SCICONltcurvesCI, replace)
  graph export SCICONltcurvesCI.wmf, replace
  graph drop SCICONltcurvesCI


  ** GRAPHS TO PLOT INDIVIDUAL CURVES BY GROUP (Fascicle and tendon length-tension curves)

  foreach lname of varlist * {
    if substr("`lname'",1,2) == "lf" & (substr("`lname'",-4,3) == "nor" | substr("`lname'",-5,3) == "nor" | substr("`lname'",-6,3) == "nor") {
      local suffix = substr("`lname'",3,.)
      local nor : display `"`nor'"' `"(line lf"' `"`suffix'"' `" "' `"torqe"' `"`suffix'"' `", lpattern(solid) lcolor(blue)) "'
    }
    if substr("`lname'",1,2) == "lf" & (substr("`lname'",-4,3) == "sci" | substr("`lname'",-5,3) == "sci" | substr("`lname'",-6,3) == "sci") {
      local suffix = substr("`lname'",3,.)
      local sci : display `"`sci'"' `"(line lf"' `"`suffix'"' `" "' `"torqe"' `"`suffix'"' `", lpattern(solid) lcolor(red)) "'
    }
    if substr("`lname'",1,2) == "lf" & (substr("`lname'",-4,3) == "str" | substr("`lname'",-5,3) == "str" | substr("`lname'",-6,3) == "str") {
      local suffix = substr("`lname'",3,.)
      local str : display `"`str'"' `"(line lf"' `"`suffix'"' `" "' `"torqe"' `"`suffix'"' `", lpattern(solid) lcolor(red)) "'
    }
  }

    ** Control only.

  cd $pathofgraphdumpname
  twoway `nor',                                                                                  ///
         xtitle("Ankle torque (Nm)") xsc(r(-10 20)) xlabel(-10(5)20, nogrid)                     ///
         ytitle("Muscle fascicle length (cm)") ysc(r(0 10)) ylabel(0(2)10, nogrid)               ///
         scale(1.5) legend(ring(0) position(11) order(1) label(1 "Control")) name(lfvstorque_NOR, replace)
  graph export lfvstorque_NOR.wmf, replace

    ** SCI only.

  twoway `sci',                                                                                  ///
         xtitle("Ankle torque (Nm)") xsc(r(-10 20)) xlabel(-10(5)20, nogrid)                     ///
         ytitle("Muscle fascicle length (cm)") ysc(r(0 10)) ylabel(0(2)10, nogrid)               ///
         scale(1.5) legend(ring(0) position(11) order(1) label(1 "SCI")) name(lfvstorque_SCI, replace)
  graph export lfvstorque_sci.wmf, replace
  
    ** SCI and control.

  twoway `nor' `sci',                                                                            ///
         xtitle("Ankle torque (Nm)") xsc(r(-10 20)) xlabel(-10(5)20, nogrid)                     ///
         ytitle("Muscle fascicle length (cm)") ysc(r(0 10)) ylabel(0(2)10, nogrid)               ///
         scale(1.5) legend(ring(0) position(11) order(1 35) label(1 "Control")  label(35 "SCI")) name(lfvstorque_NORSCI, replace)
  graph export lfvstorque_NORSCI.wmf, replace

  
  ** Plot twoway function MTU length-tension curve
  
  twoway	(function y = (1/alphag) * (exp(alphag*(x-lG))) if id == 1, lpattern(solid) lcolor(blue)) ///
			(function y = (1/alphag) * (exp(alphag*(x-lG))) if id == 2, lpattern(solid) lcolor(blue)), ///
			ytitle("Tension (N)") ylabel(20(40)200) ///
			xtitle("Length (cm)") xlabel(0(10)40) ///
			legend(order(1 "Control" 21 "SCI") colfirst notextfirst position(10) ring(0))



** CONTRACTURE LONG ANALYSIS 2011

  
  ** Describe leg length, age and gender by group.

  foreach vname of varlist ls age height weight {
    display in green _n _n _n "ANALYSIS OF: " "`vname'"
    table group, contents(mean `vname' sd `vname' count `vname' min `vname' max `vname') row
  }
  tabulate sex group, column exact
 
  ** List missing values and outliers in variables of interest, set outliers to missing
  
  foreach var of varlist lgF1 lnalphag lffsl ltfsl lgF100 lfF100 ltF100 chlg100 chlf100 chlt100 {
    quietly summarize `var' if group == "Month3"
    quietly generate st`var' = (`var' - r(mean)) / r(sd) if group == "Month3"
    quietly summarize `var' if group == "Month12"
    quietly replace st`var'  = (`var'  - r(mean)) / r(sd) if group == "Month12"
    list code group `var' st`var' if abs(st`var') > 2 & !missing(st`var')
	drop st`var'
  }
  
  replace lnalphag = . if code == 8 & group == "Month3"
  replace ltfsl = . if code == 6 & group == "Month12"

    ** Reshape wide, view scatterplots, paired-ttest results  

  preserve
  
    drop date id month12 weight anklerom e0p alphap thetap e0d alphad thetad  lG ltF1 lfF1 minLg maxLg minFg maxFg minLf maxLf minFf maxFf a b
	reshape wide age sex ls lgF1 alphag lnalphag ffsl lgfsl lffsl ltfsl lgF100 lfF100 ltF100 chlg100 chlf100 chlt100 cat100 lgstrain100 lfstrain100 ltstrain100, /// 
	i(code) j(group) string
	compress
		
	** View scatterplots of Month12 on Month3, all outcomes
	
	cd `pathofgraphfilesname'
  	foreach var in "lgF1" "alphag" "lnalphag" "ffsl" "lgfsl" "lffsl" "ltfsl" "lgF100" "lfF100" "ltF100" "chlg100" "chlf100" "chlt100" "cat100" "lgstrain100" "lfstrain100" "ltstrain100" {
      summarize `var'Month12 `var'Month3
      twoway (scatter `var'Month12 `var'Month3) (lfit `var'Month12 `var'Month3, lpattern(solid)), ///
	    legend(off) ytitle("`var'Month12") xtitle("`var'Month3") name(`var')
      quietly graph export "`var'.png", replace
	}
    graph combine lgF1 alphag lnalphag ffsl lgfsl lffsl ltfsl lgF100 lfF100 ltF100 chlf100 chlt100 cat100 lgstrain100 lfstrain100 ltstrain100, /// 
	  ysize(18) xsize(20)
	quietly graph export "Combined scatterplots.png", replace
	graph drop Graph
	graph drop lgF1 alphag lnalphag ffsl lgfsl lffsl ltfsl lgF100 lfF100 ltF100 chlg100 chlf100 chlt100 cat100 lgstrain100 lfstrain100 ltstrain100
   	
	** Non-parametric tests of significance of Month12 on Month3, all outcomes
	
	cd `pathoflogfilesname'
	foreach var in "lgF1" "alphag" "lnalphag" "ffsl" "lgfsl" "lffsl" "ltfsl" "lgF100" "lfF100" "ltF100" "chlg100" "chlf100" "chlt100" "cat100" "lgstrain100" "lfstrain100" "ltstrain100" {
	  ** Non-parametric signrank test
	  *display _n _n _n in green "NON-PARAMETRIC SIGNIFICANCE TEST OF: " "`var'"
	  *signrank `var'Month12 = `var'Month3
	  
	  ** Paired ttest
	  display _n _n _n in green "PAIRED TTEST OF: " "`var'"
	  ttest `var'Month12 = `var'Month3 // unpaired ttest
	  
	  ** Create column of differences in Month12-Month3, linear regression
	  generate `var'diff = `var'Month12 - `var'Month3
	  label variable `var'diff "Difference of `var'"
	  *display _n _n _n in green "SHANK LENGTH - ADJUSTED REGRESSION OF: " "`var'"
	  *regress `var'diff lsMonth3 // but it is meaningless to adjust for covariates that are identical in both groups
	  *estimates store diffls
	}  
	
  *_estimates clear	
  restore



** INCIDENCE OF CONTRACTURE AFTER SCI ANALYSIS 2011

