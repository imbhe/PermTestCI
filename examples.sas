****************************************************************************************************
*                                                                                                  
* ¤Name:         	examples.sas                                                                      
*                                                                                                  
* ¤Description:  	Sample datasets and examples for PermTestCI macro.   
*                                                                                                  
* ¤Date:         	2021-12-22                                
*                                                                                                  
* ¤Author:       	Henrik Imberg                                                                                                                                  
*                                                                                                                                                                                                  
****************************************************************************************************;

options nomprint nomlogic nomacrogen nosymbolgen;

***  Sample datasets ***;

data pressure_wide;
	input SBPbefore SBPafter @@;
	diff = SBPafter - SBPbefore;
	datalines;
   120 128   124 131   130 131   118 127
   140 132   128 125   140 141   135 137
   126 118   130 132   126 129   127 135
   ;
run;

data pressure_long(drop = SBPbefore SBPafter);
	set pressure_wide;
	id = _n_;
	meas = 1;
	sbp = SBPbefore;
	output;
	meas = 2;
	sbp = SBPafter;
	output;
run;

data graze;
      length GrazeType $ 10;
      input GrazeType $ WtGain @@;
      datalines;
   controlled  45   controlled  62
   controlled  96   controlled 128
   controlled 120   controlled  99
   continuous  75   continuous  54
   continuous 112   continuous  69
   continuous 104   continuous  95
   continuous  53   continuous  21
   ;
run;

data testdata_large;
	grp = 1;
	do pid = 1 to 10000;
		x = RAND('NORMAL', 0, 1); 
		output;
	end;
	grp = 2;
	do pid = 1 to 10000;
		x = RAND('NORMAL', 0, 1); 
		output;
	end;
run;



*** Examples ***;

/* One sample */
%PermTestCI(	DATA = pressure_wide,
				TYPE = onesamplemean,
				VAR = diff);


/* Paired data, wide format. */
%PermTestCI(	DATA = pressure_wide,
				TYPE = pairedmeans,
				VAR1 = SBPbefore,
				VAR2 = SBPafter);


/* Paired data, long format. */
%PermTestCI(	DATA = pressure_long,
				TYPE = pairedmeans,
				VAR = sbp,
				VISITVAR = meas,
				VIS1 = 1,
				VIS2 = 2,
				IDVAR = id);


/* Compare above to paired t-test. */
ods graphics off;
proc ttest data = pressure_wide;
      paired SBPafter * SBPbefore;
run;


/* Compute signed ranks and compare permutation test to wilcoxon test */
data pressure_wide_ranks;
	set pressure_wide;
	absdiff = abs(diff);
run;
proc rank data = pressure_wide_ranks out = pressure_wide_ranks;
	var absdiff;
run;
data pressure_wide_ranks;
	set pressure_wide_ranks;
	sgnrank = absdiff * sign(diff);
run;
%PermTestCI(	DATA = pressure_wide_ranks,
				TYPE = onesamplemean,
				VAR = sgnrank);
proc univariate data = pressure_wide;
	ods select basicmeasures testsforlocation;
	var diff;
run;
/* 	p-value of PermTestCI on signed ranks agree with wilcoxon signed rank test in proc univariate
	on small data sets where all permutations may be computed. */


/* Two samples, small sample size */
%PermTestCI(	DATA = graze,
				TYPE = twosamplemeans,
				VAR = wtgain, 
				GROUP = grazetype);

/* Compare to t-test */
ods graphics off;
proc ttest data = graze;
	class grazetype;
	var wtgain;
run;

/* Compare to non-parametric permutation test and Hodges-Lehmann interval in proc npar1way. */
ods graphics off;
ods select datascorestest hodgeslehmann;
proc npar1way data = graze hl;
	class grazetype;
	var wtgain;
	exact scores = data;
run;
/* 	p-value of PermTestCI agree with exact test in proc npar1way 
	on small datasets where all permutations may be computed. */


/* Two samples, large sample size */
%PermTestCI(	DATA = testdata_large,
				TYPE = twosamplemeans,
				VAR = x, 
				GROUP = grp);
/* Finishes computation within a few seconds. */
