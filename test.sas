****************************************************************************************************
*                                                                                                  
* ¤Name:         	test.sas                                                                      
*                                                                                                  
* ¤Description:  	Check errors and warnings.   
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

data pressure_long(drop = SBPbefore SBPafter diff);
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


*** Check errors and warnings ***;
				
data test;
	set pressure_wide;
	diff = .;
run;
%PermTestCI(	DATA = test,
				TYPE = onesamplemean,
				VAR = diff);

data test;
	set pressure_wide;
	if _n_ < 6;
run;
%PermTestCI(	DATA = test,
				TYPE = onesamplemean,
				VAR = diff);

data test;
	set pressure_wide;
	if _n_ < 7;
run;
%PermTestCI(	DATA = test,
				TYPE = onesamplemean,
				VAR = diff);
		
data test;
	set pressure_wide;
	SBPbefore = .;
run;
%PermTestCI(	DATA = test,
				TYPE = pairedmeans,
				VAR1 = SBPbefore,
				VAR2 = SBPafter);

data test;
	set pressure_long;
	if meas = 1 then sbp = .;
run;
%PermTestCI(	DATA = test,
				TYPE = pairedmeans,
				VAR = sbp,
				VISITVAR = meas,
				VIS1 = 1,
				VIS2 = 2,
				IDVAR = id);

data test;
	set graze;
	if GrazeType eq "controlled" then WtGain = .;
run;
%PermTestCI(	DATA = test,
				TYPE = twosamplemeans,
				VAR = wtgain, 
				GROUP = grazetype);

data test;
	set graze;
	if GrazeType eq "controlled";
run;
%PermTestCI(	DATA = test,
				TYPE = twosamplemeans,
				VAR = wtgain, 
				GROUP = grazetype);

data test;
	set graze end = eof;
	if eof then GrazeType = "Hej";
	output;
run;
%PermTestCI(	DATA = test,
				TYPE = twosamplemeans,
				VAR = wtgain, 
				GROUP = grazetype);
