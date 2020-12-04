
    ********** NOTES ********** 

      ** Author: Joanna Diong
      ** Date: 14 July 2019

      ** This do file analyses muscle activation (twitch interpolation) and EMG 
      ** data to obtain coefficients to predict EMG from activation in future in stroke. 
      **
      ** Data are in CSV format extracted from Spike2 files using Python.

    ********** PRELIMINARIES ********** 

  version 16
  clear all
  clear matrix
  drop _all
  capture log close
   
      ** Settings.
   
  set more off
  set scheme lean3
  pause on

      ** Linux: Enter name of root directory and name of this do file. ** EDIT AS REQUIRED **

  local pathname = `""/home/joanna/Dropbox/Projects/activation/src/code/stats/""'
  local dofilename = `""norm.do""'

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
  
  forvalues subject = 1/25 { 
      local graphtext1 `graphtext1' (line emg_norm_mvc activations if subjectGroup == `subject', lpattern(solid) lcolor(gs12))
      local graphtext2 `graphtext2' (line emg_norm_mmax activations if subjectGroup == `subject', lpattern(solid) lcolor(gs12))
      local graphtext3 `graphtext3' (line emg_norm_mvc emg_norm_mmax if subjectGroup == `subject', lpattern(solid) lcolor(gs12))
  }
  sort activations
  
  twoway `graphtext1', xtitle("Activation [%]") ytitle("EMG [%MVC]") aspect(1) legend(off) name(norm_mvc)
  
  twoway `graphtext2', xtitle("Activation [%]") ytitle("EMG [%Mmax]") aspect(1) legend(off) name(norm_mmax)
  
  twoway `graphtext3', xtitle("EMG [%Mmax]") ytitle("EMG [%MVC]") aspect(1) legend(off) name(norm_corr)
  
  
  forvalues subject = 1/2 {
      regress emg_norm_mvc emg_norm_mmax if subjectGroup == `subject'
      matrix b = e(b)
      display b[1,1]
  }
  
  /*
  xtmixed emg_norm_mvc emg_norm_mmax || subjectGroup: emg_norm_mmax, cov(unstructured) 
  estat icc
  
  * model visualisation: for each subject
  predict emg_norm_mvc_Fit, fitted
  sort subjectGroup activations
  twoway (line emg_norm_mvc_Fit emg_norm_mmax, connect(ascending)), xtitle(EMG [%Mmax]) ytitle(EMG [*MVC]), ///
    graph export norm_pred.tif, width(1200) height(900) replace
  
  quietly {
      * get mean slope and 95% CI
      predict mixedPred
      predict mixedSE, stdp
      generate mixedLCL = mixedPred - invnormal(0.975) * mixedSE
      generate mixedUCL = mixedPred + invnormal(0.975) * mixedSE
      
      * get 95% prediction interval
      generate mixedLPL = mixedPred - invnormal(0.975) * mixedSE * sqrt(24) /* 1 subject dropped, values >4 SD from mean */
      generate mixedUPL = mixedPred + invnormal(0.975) * mixedSE * sqrt(24)
    }
    
  twoway (rarea mixedLPL mixedUPL emg_norm_mmax, astyle(sunflowerlb)) /// 
        (rarea mixedLCL mixedUCL emg_norm_mmax, astyle(sunflowerdb)) ///
        (line mixedPred emg_norm_mmax, lpattern(solid)) ///
        , ///
        xtitle("EMG [%Mmax]", size(medium)) xlabel(,labsize(medium)) /// 
        ytitle("EMG [%MVC]]", size(medium)) ylabel(,labsize(medium)) /// 
        legend(off) aspect(1) legend(off) name(norm_mvc_mmax, replace)
        
  graph combine norm_mvc norm_mmax norm_mvc_mmax, rows(1) ysize(4) xsize(12) name(all_norm, replace) 
  graph export all_norm.tif, width(1200) height(900) replace
  */  
  
  
    ********** END ********** 


  *log close
  exit
