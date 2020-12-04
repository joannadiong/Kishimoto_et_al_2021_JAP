

    ********** NOTES ********** 

      ** Author: Joanna Diong
      ** Date: 1 January 2013

      ** This do file [describe what it does].
      **
      ** Data are [describe source and content of data file].
      **
      ** The code was checked by [name] on [date]


    ********** PRELIMINARIES ********** 

  version 13
  clear all
  clear matrix
  drop _all
  capture log close
   
      ** Settings.
   
  set more off
  set scheme lean3
  pause on
/*
*/
   
      ** Enter name of root directory and name of this do file. ** EDIT AS REQUIRED **

  local pathname = `""C:\Documents and Settings\Joanna Diong\My Documents\My Stata Files\""'
  local dofilename = `""Do file template.do""'

      ** Open a time- and date-stamped log file and copy a time- and date-stamped
      ** do file and data file to the log file directory.

  local pathandnameofdofile = `"""' + `pathname' + `dofilename' + `"""'
  local pathoflogfilesname = `"""' + `pathname' + "Log files\" + `"""'
  local pathofdatafilesname = `"""' + `pathname' + "Data\" + `"""'
  local pathofgraphfilesname = `"""' + `pathname' + "Graphs\" + `"""'
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


    ********** IMPORT DATA AND CLEAN ********** 

  cd `pathofdatafilesname'
  import excel "Hypothetical data file.xlsx", sheet("Sheet1") firstrow clear


    ********** RUN ANALYSIS ********** 

  summarize var1-var5


    ********** GRAPH DATA ********** 

  cd `pathofgraphfilesname'
  hist var2
  graph export "Example graph.tif", replace

  notes: best variable name convention: varname_group_subjectid
  
  foreach lname of varlist * {
    if substr("`lname'",1,5) == "alpha" & substr("`lname'",7,3) == "nor" {
      local suffix = substr("`lname'",6,.)
      local aa : display `"`aa'"' `"(line alpha"' `"`suffix'"' `" "' `"force"' `"`suffix'"' `", lpattern(solid) lcolor(blue)) "'
    }
    if substr("`lname'",1,5) == "alpha" & substr("`lname'",7,3) == "sci" {
      local suffix = substr("`lname'",6,.)
      local bb : display `"`bb'"' `"(line alpha"' `"`suffix'"' `" "' `"force"' `"`suffix'"' `", lpattern(solid) lcolor(red)) "'
    }
  }
  
  twoway `aa' `bb',                                           ///
    xtitle("Title in X (units)") xlabel(0(100)300) xsize(3)   ///
    ytitle("Title in Y (units)") ylabel(0(10)30) ysize(2)     ///
    legend(ring(0) position(2) order(1 "Control" 30 "SCI") colfirst notextfirst) text(100 -10 "a") scale(1.25) 
  *quietly graph export "aabb.tif", width(900) height(600) replace	
  *graph drop Graph
  graph save aabb, replace
  
  twoway `cc' `dd',                                           ///
    xtitle("Title in X (units)") xlabel(0(100)300) xsize(2)   ///
    ytitle("Title in Y (units)") ylabel(0(10)30) ysize(2)     ///
    legend(ring(0) position(2) order(1 "Control" 30 "SCI") colfirst notextfirst) text(100 -10 "b") scale(1.25)
  *quietly graph export "ccdd.tif", width(900) height(600) replace	
  *graph drop Graph
  graph save ccdd, replace
  graph combine aabb.gph ccdd.gph, col(1) xsize(3) ysize(4) //Edit xsize and ysize in twoway options
  graph export aabbccdd.tif, width(900) height(1200) replace  


    ********** END. ZZZZ  ********** 

  *log close
  exit
