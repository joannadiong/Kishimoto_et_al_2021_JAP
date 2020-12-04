  version 13
  clear all
  cls
  set scheme lean3
  
  * install mltrsq
  * findit mltrsq
  
  local pathname = `""/home/joanna/Dropbox/Projects/activation/src/code/stats/""'
  cd `pathname'
  
  import delimited using "subjects_data.csv", delimiters(",") varnames(1)
  
  egen subjectGroup = group(subject)
  sort activations
  
  * mixed c.lnemgso c.activations || subjectGroup: 
  xtmixed lnemgso activations || subjectGroup: activations, cov(unstructured) 
  estat icc
  mltrsq
  
  * model visualisation: for each subject
  predict emgFit, fitted
  sort subjectGroup activations
  twoway(line emgFit activations, connect(ascending)), xtitle(Activation [%]) ytitle(EMG SO [ln(mV)])
  
    
  * check our calculated CI are the same as margins CI
  margins, at(activations=(0(20)100)) 
  marginsplot, recast(line) recastci(rarea) aspect(1)
  graph export margins_lnemgso.tif, width(1200) height(900) replace
  
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
    
  forvalues subject = 1/25 { 
          local graphtext2 `graphtext2' (line lnemgso activations if subjectGroup == `subject', lpattern(solid) lcolor(gs10))
      }
  sort activations
  
  twoway (rarea mixedLPL mixedUPL activations, astyle(sunflowerlb)) /// /* 95% PI */
        `graphtext2' /// 
        (rarea mixedLCL mixedUCL activations, astyle(sunflowerdb)) ///
        (line mixedPred activations, lpattern(solid)) ///
        , ///
        xtitle("Activation [%]", size(small)) xlabel(,labsize(small)) /// 
        ytitle("EMG SO [ln(mV)]", size(small)) ylabel(,labsize(small)) /// 
        legend(off) aspect(1) legend(off) name(mixed_lnemgso, replace)
        graph export mixed_lnemgso.tif, width(1200) height(900) replace

