
    ********** NOTES ********** 

      ** Author: Joanna Diong
      ** Date: 14 July 2019

      ** This do file analyses muscle activation (twitch interpolation) and EMG 
      ** data to obtain coefficients to predict EMG from activation in future in stroke. 
      **
      ** Data are in CSV format extracted from Spike2 files using Python.

    ********** PRELIMINARIES ********** 

  *version 13
  clear all
  clear matrix
  drop _all
  capture log close
   
      ** Settings.
  
  *ssc install mlt
  
  set more off
  set scheme lean3
  pause on

      ** Linux: Enter name of root directory and name of this do file. ** EDIT AS REQUIRED **

  local pathname = `""/home/joanna/Dropbox/Projects/activation_repo/stats/""'
  local dofilename = `""script.do""'

      ** Open a time- and date-stamped log file and copy a time- and date-stamped
      ** do file and data file to the log file directory.

  local pathandnameofdofile = `"""' + `pathname' + `dofilename' + `"""'
  local pathoflogfilesname = `"""' + `pathname' + "log_files/" + `"""'
  local pathofdatafilesname = `"""' + `pathname' + "data/" + `"""'
  local pathofgraphfilesname = `"""' + `pathname' + "graphs/" + `"""'
  local cdate = "c(current_date)"
  local ctime = "c(current_time)"
  local ctime = subinstr(`ctime',":","h",1)
  local ctime = subinstr("`ctime'",":","m",1)
  local logfilename = `"""' + "Log " + `cdate' + " " + "`ctime'" + "s.log" + `"""'
  local backupdofilename = `"""' + "Do " + `cdate' + " " + "`ctime'" + "s.txt" + `"""'
  cd `pathoflogfilesname'
  *log using `logfilename'
  *copy `pathandnameofdofile' `backupdofilename'
  cd `pathofdatafilesname'


    ********** IMPORT DATA ********** 

  cd `pathofdatafilesname'
  import delimited using "subjects_data.csv", delimiters(",") varnames(1)
  
  egen subjectGroup = group(subject)
  sort activations


    ********** GRAPH DATA ********** 
  
  cd `pathofgraphfilesname'

  foreach var of varlist torques emgso emgmg emglg {
/*
    * Plot individual subject data with regression lines.

    forvalues subject = 1/25 { 
      if `var' == torques {
        twoway (line activations `var' if subjectGroup == `subject', lpattern(solid)) ///
          (lfit activations `var' if subjectGroup == `subject', lpattern(solid) lcolor(red)) ///
          , ///
          xtitle("Torque [%MVC]") ytitle("Activation [%]") scale(0.8) legend(off) name(eachSubject`subject', replace)
      }
      if `var' == emgso {
        twoway (line `var' activations if subjectGroup == `subject', lpattern(solid)) ///
          (lfit `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(red)) ///
          , ///
          xtitle("Activation [%]") ytitle("EMG SO [mV]") scale(0.8) legend(off) name(eachSubject`subject', replace)
      }
      if `var' == emgmg {
        twoway (line `var' activations if subjectGroup == `subject', lpattern(solid)) ///
          (lfit `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(red)) ///
          , ///
          xtitle("Activation [%]") ytitle("EMG MG [mV]") scale(0.8) legend(off) name(eachSubject`subject', replace)
      }
      if `var' == emglg {
        twoway (line `var' activations if subjectGroup == `subject', lpattern(solid)) ///
          (lfit `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(red)) ///
          , ///
          xtitle("Activation [%]") ytitle("EMG LG [mV]") scale(0.8) legend(off) name(eachSubject`subject', replace)
      }
    }
    graph combine eachSubject1 eachSubject2 eachSubject3 eachSubject4 eachSubject5 ///
                  eachSubject6 eachSubject7 eachSubject8 eachSubject9 eachSubject10 ///
                  eachSubject11 eachSubject12 eachSubject13 eachSubject14 eachSubject15 ///
                  eachSubject16 eachSubject17 eachSubject18 eachSubject19 eachSubject20 ///
                  eachSubject21 eachSubject22 eachSubject23 eachSubject24 eachSubject25 ///
                  , ///
                  xcommon ycommon name(each_`var', replace)
    graph export each_`var'.tif, width(1200) height(900) replace
    graph drop eachSubject1 eachSubject2 eachSubject3 eachSubject4 eachSubject5 ///
               eachSubject6 eachSubject7 eachSubject8 eachSubject9 eachSubject10 ///
               eachSubject11 eachSubject12 eachSubject13 eachSubject14 eachSubject15 ///
               eachSubject16 eachSubject17 eachSubject18 eachSubject19 eachSubject20 ///
               eachSubject21 eachSubject22 eachSubject23 eachSubject24 eachSubject25
*/
    * Plot EMG-activation individual subject data on one plot.

    forvalues subject = 1/25 {
      local graphtext1 `graphtext1' (line `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
      local graphtext1_ `graphtext1_' (line activations `var' if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
    }
    sort activations

    if `var' == torques {
      twoway `graphtext1_' ///
        , ///
        xtitle("Torque [%MVC]") ytitle("Activation [%]") legend(off) aspect(1) name(all_`var', replace)
    }
    if `var' == emgso {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
        ytitle("EMG SO [mV]", size(small)) ylabel(,labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    if `var' == emgmg {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) ///
        ytitle("EMG MG [mV]", size(small)) ylabel(,labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    if `var' == emglg {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
        ytitle("EMG LG [mV]", size(small)) ylabel(,labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    graph export all_`var'.tif, width(1200) height(900) replace

    local graphtext1 = ""
    local graphtext1_ = ""

  }

    ********** ANALYSE DATA ********** 

  cls
  
    * Drop subject EMG data if trial EMG values are further than 4 SD from mean EMG 
  
  bys trial: egen avlnemgso = mean(lnemgso)
  bys trial: egen sdlnemgso = sd(lnemgso)
  
  bys trial: egen avlnemgmg = mean(lnemgmg)
  bys trial: egen sdlnemgmg = sd(lnemgmg)
  
  bys trial: egen avlnemglg = mean(lnemglg)
  bys trial: egen sdlnemglg = sd(lnemglg)
  
  forvalues subject = 1/25 {
    list subject lnemgso if subjectGroup == `subject' & (lnemgso  <= avlnemgso - 4 * sdlnemgso | lnemgso  >= avlnemgso + 4 * sdlnemgso)
    list subject lnemgmg if subjectGroup == `subject' & (lnemgmg  <= avlnemgmg - 4 * sdlnemgmg | lnemgmg  >= avlnemgmg + 4 * sdlnemgmg)
    list subject lnemglg if subjectGroup == `subject' & (lnemglg  <= avlnemglg - 4 * sdlnemglg | lnemglg  >= avlnemglg + 4 * sdlnemglg)
  }
  replace lnemgmg = . if subject == "sub17"
  sort subject trials
  
    * Analyse activations and logged EMG data
  
  foreach var of varlist torques lnemgso lnemgmg lnemglg {
    
    * Plot lnEMG-activation individual subject data on one plot.
    
    forvalues subject = 1/25 {
      local graphtext1 `graphtext1' (line `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
    }
    sort activations

    if `var' == lnemgso {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
        ytitle("EMG SO [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    if `var' == lnemgmg {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
        ytitle("EMG MG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    if `var' == lnemglg {
      twoway `graphtext1' ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
        ytitle("EMG LG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
        legend(off) aspect(1) name(all_`var', replace)
    }
    graph export all_`var'.tif, width(1200) height(900) replace

    local graphtext1 = ""
    

    * Fit a linear mixed model with random slopes to the logged EMG data, and graph the fitted data.
    
    display _n _n in green "ANALYSIS OF: " "`var'"

    xtmixed `var' activations || subjectGroup: activations, cov(unstructured) 
    estat icc
    mltrsq
    
    quietly {
      * get mean slope and 95% CI
      predict mixedPred
      predict mixedSE, stdp
      generate mixedLCL = mixedPred - invnormal(0.975) * mixedSE
      generate mixedUCL = mixedPred + invnormal(0.975) * mixedSE
      * get 95% prediction interval
      generate mixedLPL = mixedPred - invnormal(0.975) * mixedSE * sqrt(24) /* 1 subject dropped, values >4 SD from mean */
      generate mixedUPL = mixedPred + invnormal(0.975) * mixedSE * sqrt(24)
      
      * git fitted values
      predict emgFit, fitted
    }
    
      forvalues subject = 1/25 { 
          local graphtext2 `graphtext2' (line `var' activations if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
          local graphtext2_ `graphtext2_' (line activations `var' if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
      }
      sort activations

      if `var' == torques {
        twoway `graphtext2_', ///
          xtitle("Torque [%MVC]") ytitle("Activation [%]") legend(off) aspect(1) legend(off) name(mixed_`var', replace)
          graph export torques.tif, width(1200) height(900) replace
          graph export torques.eps, replace
      }
      if `var' == lnemgso {
        /*
        twoway `graphtext2' ///
        (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
              (line mixedPred activations, lpattern(solid)) ///
              , ///
              xtitle("Activation [%]") ytitle("EMG SO [ln(mV)]") legend(off) aspect(1) legend(off) name(mixed_`var', replace)
              graph export mixed_`var'.tif, width(1200) height(900) replace
        */
        twoway (rarea mixedLPL mixedUPL activations, astyle(sunflowerlb)) /// /* 95% PI */
              (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
              (line mixedPred activations, lpattern(solid)) ///
              , ///
              xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
              ytitle("EMG SO [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
              legend(off) aspect(1) name(mixed_`var', replace)
              graph export mixed_`var'.tif, width(1200) height(900) replace
	      graph export mixed_`var'.svg, replace
        
        preserve
        sort subjectGroup activations
        twoway (line emgFit activations, connect(ascending)), /// 
              xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
              ytitle("EMG SO [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
              legend(off) aspect(1) name(emgFit_`var', replace)
                graph export emgFit_`var'.tif, width(1200) height(900) replace
        restore
    
      }
      if `var' == lnemgmg {
        /*
        twoway `graphtext2' ///
        (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
              (line mixedPred activations, lpattern(solid)) ///
              , ///
              xtitle("Activation [%]") ytitle("EMG MG [ln(mV)]") legend(off) aspect(1) legend(off) name(mixed_`var', replace)
              graph export mixed_`var'.tif, width(1200) height(900) replace
        */
        twoway (rarea mixedLPL mixedUPL activations, astyle(sunflowerlb)) /// 
            (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
            (line mixedPred activations, lpattern(solid)) ///
            , ///
            xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
            ytitle("EMG MG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
            legend(off) aspect(1) name(mixed_`var', replace)
            graph export mixed_`var'.tif, width(1200) height(900) replace

        preserve
        sort subjectGroup activations
        twoway(line emgFit activations, connect(ascending)), /// 
            xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
            ytitle("EMG MG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
            legend(off) aspect(1) name(emgFit_`var', replace)
            graph export emgFit_`var'.tif, width(1200) height(900) replace
        restore
    
      }
      if `var' == lnemglg {
        /*
        twoway `graphtext2' ///
            (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
            (line mixedPred activations, lpattern(solid)) ///
            , ///
            xtitle("Activation [%]") ytitle("EMG LG [ln(mV)]") legend(off) aspect(1) legend(off) name(mixed_`var', replace)
            graph export mixed_`var'.tif, width(1200) height(900) replace
        */
        twoway (rarea mixedLPL mixedUPL activations, astyle(sunflowerlb)) /// 
            (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
            (line mixedPred activations, lpattern(solid)) ///
            , ///
            xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
            ytitle("EMG LG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
            legend(off) aspect(1) name(mixed_`var', replace)
            graph export mixed_`var'.tif, width(1200) height(900) replace

        preserve
        sort subjectGroup activations
        twoway(line emgFit activations, connect(ascending)), /// 
            xtitle("Activation [%]", size(small)) xlabel(, labsize(small)) /// 
            ytitle("EMG LG [ln(mV)]", size(small)) ylabel(, labsize(small)) /// 
            legend(off) aspect(1) name(emgFit_`var', replace)
            graph export emgFit_`var'.tif, width(1200) height(900) replace
        restore                         
    
      } 

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
    drop mixedPred mixedSE mixedLCL mixedUCL mixedLPL mixedUPL emgFit
    
    /*
    
    * Fit FP with robust SEs and extract predicted values. Show plot of fit.

    mfp: regress `var' activations, cluster(subject)
    quietly {
        predict fracpred
        predict fracSE, stdp
        generate LCL = fracpred - invnormal(0.975) * fracSE
        generate UCL = fracpred + invnormal(0.975) * fracSE
    }
    di "`graphtext2'"
    sort activations

    if `var' == torques {
      twoway `graphtext2_', ///
        xtitle("Torque (%MVC)") ytitle("Activation (%)") legend(off) aspect(1) name(frac_`var', replace)
        * graph export frac_`var'.tif, width(1200) height(900) replace
    }
    if `var' == emgso {
      twoway (rarea LCL UCL activations, astyle(sunflowerdb)) ///
        `graphtext2' ///
        (line fracpred activations, lpattern(solid) lcolor(black)) ///
        , ///
        xtitle("Activation (%)") ytitle("EMG SO (mV)") legend(off) aspect(1) name(frac_`var', replace)
        graph export frac_`var'.tif, width(1200) height(900) replace
    }
    if `var' == emgmg {
      twoway (rarea LCL UCL activations, astyle(sunflowerdb)) ///
        `graphtext2' ///
        (line fracpred activations, lpattern(solid) lcolor(black)) ///
        , ///
        xtitle("Activation (%)") ytitle("EMG MG (mV)") legend(off) aspect(1) name(frac_`var', replace)
        graph export frac_`var'.tif, width(1200) height(900) replace
    }
    if `var' == emglg {
      twoway (rarea LCL UCL activations, astyle(sunflowerdb)) ///
        `graphtext2' ///
        (line fracpred activations, lpattern(solid) lcolor(black)) ///
        , ///
        xtitle("Activation (%)") ytitle("EMG LG (mV)") legend(off) aspect(1) name(frac_`var', replace)
        graph export frac_`var'.tif, width(1200) height(900) replace
    }
    
    drop fracpred fracSE LCL UCL
    */
    
    local graphtext2 = ""
    local graphtext2_ = ""
    
  }

  graph combine all_emgso all_lnemgso mixed_lnemgso ///
        all_emgmg all_lnemgmg mixed_lnemgmg ///
        all_emglg all_lnemglg mixed_lnemglg ///
        , ///
        rows(3) ysize(12) xsize(12) name(all_mixed, replace) 
  graph export all_mixed.tif, width(4800) replace
  graph export all_mixed.eps, replace

  graph drop all_torques all_emgso all_emgmg all_emglg
  graph drop all_lnemgso all_lnemgmg all_lnemglg 
  *graph drop each_torques each_emgso each_emglg each_emgmg
  graph drop emgFit_lnemgso emgFit_lnemgmg emgFit_lnemglg 

  graph drop mixed_torques mixed_lnemgso mixed_lnemgmg mixed_lnemglg 
  *graph drop frac_torques frac_emgso frac_emgmg frac_emglg
  graph drop all_mixed


    ********** END ********** 


  *log close
  exit

  /*
  graph combine all_emgso all_emgmg all_emglg ///
    , ///
    rows(3) ysize(9) xsize(4) name(emg, replace) /* 8 x 600 dpi = 4800 */
    graph export emg.tif, width(2400) replace
    graph export emg.eps, replace

  graph combine mixed_emgso frac_emgso mixed_emgmg frac_emgmg mixed_emglg frac_emglg ///
    , ///
    rows(3) ysize(12) xsize(8) name(mixed_frac_all, replace) 
  graph export mixed_frac_all.tif, width(4800) replace
  graph export mixed_frac_all.eps, replace

  graph drop all_torques all_emgso all_emgmg all_emglg
  graph drop each_torques each_emgso each_emglg each_emgmg

  graph drop mixed_torques mixed_emgso mixed_emgmg mixed_emglg 
  graph drop frac_torques frac_emgso frac_emgmg frac_emglg
  graph drop mixed_frac_all
  graph drop emg
  */
