****************************************************************************************************
*                                                                                                  
* ¤Name:         	PermTestCI.sas                                                                      
*                                                                                                  
* ¤Description:  	SAS macro for computation of confidence intervals for the mean of one sample 
*					and for the mean difference between two independent or paired samples, 
*					by inversion of the Fisher-Pitman non-parametric permutation test. 
*                                                                                                  
* ¤Date:         	2021-12-22                                 
*                                                                                                  
* ¤Author:       	Henrik Imberg                                                                                                                                  
*                                                                                                                                                                                                  
****************************************************************************************************;

%macro PermTestCI(	DATA = , 
					WHERE = , 
					TYPE = , 
					VAR = , 
					GROUP = , 
					IDVAR = , 
					VISITVAR = , 
					VIS1 = ,
					VIS2 = ,
					VAR1 = , 
					VAR2 = , 
					OUTDATA = _PERMTESTCI_, 
					NPERM = ., 
					DECIMALS = 4,
					THETA0 = 0, 
					CLEVEL = 0.95, 
					SEED = 123456, 
					PRINT = Y);

/*
DATA 	Analysis dataset. 

WHERE	Where clause to subset a set of observations from the analysis data set. 
		Requires a valid data step where clause, e.g. WHERE = where ITT eq 1. 
		If nothing is specified, all observations of the analysis dataset are used. 

TYPE 	Type of confidence interval to the computed. 
		Valid types are ONESAMPLEMEAN, TWOSAMPLEMEANS and PAIREDMEANS, requesting a confidence interval for the mean of a single sample, for the difference between means of two independent samples, and for the difference between means of two paired samples, respectively. 
		For TYPE = PAIREDMEANS, the analysis data set may be structured either in wide (“ADSL”) or long (“ADVL”) format; see additional details about the IDVAR, VISITVAR, VIS1, VIS2, VAR1 and VAR2 parameters. 
		When TYPE = TWOSAMPLEMEANS, the difference between sample means is computed by subtracting the mean of the highest unformatted group from the mean of the lowest unformatted group, similar to the default option in the TTEST procedure. 
		When TYPE = PAIREDMEANS, the difference is computed by subtracting the mean of the first visit/measurement, as specified by VAR1 (data in wide “ADSL” format) or by the VIS1-value of VISITVAR (data in long “ADVL” format), from the mean of the second visit, as specified by VAR2 (data in wide “ADSL” format) or by the VIS2-value of VISITVAR (data in long “ADVL” format).

VAR 	Analysis variable for which a confidence interval will be computed. Required parameter for TYPE = ONSAMPLEMEAN and TWOSAMPLEMEANS. 
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in long (ADVL) format. 
		Invalid parameter for TYPE = PAIREDMEANS when the analysis data set is in wide (ADSL) format.

GROUP	Group variable. Required parameter for TYPE = TWOSAMPLEMEANS. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and PAIREDMEANS. 
		Only group variables with two distinct values are allowed. 
		The group variable may be either character or numeric. 

IDVAR 	ID variable. Each value of the ID-variable must uniquely determine a subject in the analysis data set. 
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in long (ADVL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in wide (ADSL) format. 
		Must be just in conjunction with the VISITVAR, VIS1 and VIS2 parameters. May not be used in conjunction with the VAR1 and VAR2 parameters. 

VISITVAR	Visit variable. The visit variable specifies visits/measurements for which paired differences will be computed. 
			Required parameter for TYPE = PAIREDMEANS when the analysis data set is in long (ADVL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in wide (ADSL) format. 
		Must be just in conjunction with the IDVAR, VIS1 and VIS2 parameters. 
		May not be used in conjunction with the VAR1 and VAR2 parameters. 
		Both numeric and character variables are accepted.

VIS1 	Specifies the value of VISITVAR for the first visit/measurement when TYPE = PAIREDMEANS and data is structured in long (ADVL) format.  
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in long (ADVL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in wide (ADSL) format. 
		Must be just in conjunction with the IDVAR, VISITVAR and VIS2 parameters. May not be used in conjunction with the VAR1 and VAR2 parameters. 

VIS2 	Specifies the value of VISITVAR for the second visit/measurement when TYPE = PAIREDMEANS and data is structured in long (ADVL) format.  
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in long (ADVL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in wide (ADSL) format. 
		Must be just in conjunction with the IDVAR, VISITVAR and VIS2 parameters. May not be used in conjunction with the VAR1 and VAR2 parameters.

VAR1 	Variable name of the first visit/measurement when TYPE = PAIREDMEANS and data is structured in wide (ADSL) format. 
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in wide (ADSL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in long (ADVL) format. 
		Must be just in conjunction with the VAR2 parameter. 
		May not be used in conjunction with the IDVAR, VISITVAR, VIS1 and VIS2 parameters. 

VAR2 	Variable name of the second visit/measurement when TYPE = PAIREDMEANS and data is structured in wide (ADSL) format. 
		Required parameter for TYPE = PAIREDMEANS when the analysis data set is in wide (ADSL) format. 
		Invalid parameter for TYPE = ONESAMPLEMEAN and TWOSAMPLEMEANS, and for TYPE = PAIREDMEANS when the analysis data set in in long (ADVL) format. 
		Must be just in conjunction with the VAR1 parameter. May not be used in conjunction with the IDVAR, VISITVAR, VIS1 and VIS2 parameters. 

OUTDATA		Output data set where results are stored. Default value is OUTDATA = _PermTestCI_. The following information is stored in the output data set:
		o MEAN: the one-sample mean. Only included in output data set if TYPE = ONESAMPLEMEAN.
		o DIFF: the difference between the two sample means. Only included in output data set if TYPE = TWOSAMPLEMEANS or PAIREDMEANS.
		o LCL: lower confidence limit of mean or mean difference.
		o UCL: upper confidence limit of mean or mean difference.
		o THETA0: value of the statistic (mean or difference between means) under the null hypothesis. By default, THETA0 is set to 0. Does not affect the confidence limits.
		o PVAL: two-sided p-value for testing the null hypothesis that the statistic (mean or difference between means) equals THETA0. 
		o CLEVEL: the confidence level of the confidence interval.
		o NPERM: number of permutations used for computation of the confidence interval.
		o METHOD: method used for computation of confidence interval. Equals “Approximate” when an approximate confidence interval using NPERM random permutations is computed and “Exact” when an exact confidence interval using all possible permutations is computed. 
		o NOTE: details on method and number of permutations used for computation of confidence interval.
		o WARN1: Contains a warning if the p-value computed from the permutation test has a 95% confidence interval that cover the value 1-CLEVEL. Empty string otherwise.
		o WARN2: Contains a warning if the confidence interval covers the real line. Empty string otherwise.
		o NPERM: Number of random permutations used for computing p-value and confidence intervals. 

NPERM	Number of random permutations used for computing p-value and confidence intervals. 
		Must be a positive integer or ‘.’ (missing). To use default number of permutations, leave unspecified or set NPERM = . , (missing). 
		Default number of permutations is
		 - 10 000 if N = 100
		 - 5000 if 100 < N = 200
		 - 2000 if 200 < N = 500
		 - 1000 otherwise,		
		where N is the number of subjects in the analysis data set. 
		If NPERM is greater than the total number of possible permutations, the macro uses all possible permutations rather than a random selection of permutations and changes NPERM accordingly. 

DECIMALS	Number of decimals used when printing results. Must be a positive integer. 	Default value is DECIMALS = 4.

THETA0 	Value of the statistic of interest (i.e. the mean if TYPE = ONESAMPLEMEAN and difference between means if TYPE = TWOSAMPLEMEANS or PAIREDMEANS) under the null hypothesis. Default value is THETA0 = 0. 

CLEVEL 	Desired confidence level. Must be strictly greater than 0 and strictly smaller than 1. Default value is CLEVEL = 0.95. 

SEED 	Seed for random number generator, which uses the SAS randgen routine. Must be a positive integer. Default value is SEED = 123456.

PRINT 	Should results be printed to open ODS output destination(s) (Y/N)? Default value is PRINT = Y. Results are printed to the current output destinations if PRINT = Y or YES. 

Note: all parameters that require character input are case insensitive.
*/

	%PUT NOTE: Macro PermTestCI starting execution.;

	/* Turn MLOGIC, SYMBOLGEN and MPRINT OFF. These are restored to current options upon completion of the macro */
	%LOCAL saveOptions;
	%LET saveOptions = %SYSFUNC(getoption(MLOGIC)) %SYSFUNC(getoption(SYMBOLGEN)) %SYSFUNC(getoption(MPRINT)) %SYSFUNC(getoption(notes)) linesize = %SYSFUNC(getoption(linesize));
	%LET _ptci_abort_ = 0;
	options nomlogic nosymbolgen nomprint linesize = max nonotes;


	/* Clear output data set */
	proc datasets nolist;
		delete &outdata.;
	quit;


	****** Computing time ******;

	%LET _ptci_sdtm_ = %SYSFUNC(datetime());


	****** Check input arguments ******;

	%ptci_argcheck;
	 

	****** Titles ******;

	data _ptci_oldtit_;
		set sashelp.VTitle;  
		if type = "T";
	run;
	proc sql noprint;
		select min(6, max(number)) into :titnum from  sashelp.vtitle;
	quit;

	data _null_;
		tit1 = "The PermTestCI macro";
		if lowcase("&type.") eq "twosamplemeans" then do;
			tit2 = "Non-parametric confidence interval for the mean difference between two independent samples";
		end;
		if lowcase("&type.") eq "onesamplemean" then do;
			tit2 = "Non-parametric confidence interval for the sample mean";
		end;
		if lowcase("&type.") eq "pairedmeans" then do;
			tit2 = "Non-parametric confidence interval for the mean difference between two paired samples";
		end;
		call symput("_ptci_tit1_", tit1);
		call symput("_ptci_tit2_", tit2);
	run;
	title%EVAL(&titnum. + 1) "&_ptci_tit1_.";
	title%EVAL(&titnum. + 2) "&_ptci_tit2_.";


	****** Prepare datasets ******;

	data _ptci_indata_;
		set &data.;
		&where.;
	run;	

	%LET _ptci_dsid_ = %SYSFUNC(open(_ptci_indata_));
	%LET _ptci_nobs_ = %SYSFUNC(attrn(&_ptci_dsid_, nobs));
	%LET _ptci_dsid_ = %SYSFUNC(close(&_ptci_dsid_));
	%IF &_ptci_nobs_. = 0 %THEN %DO;
		%ptci_term;
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Cannot find any records in the input data set &data. &where..;
		%ABORT;
	%END;



	****** ONESAMPLEMEAN ******;
	%IF %QLOWCASE(&type.) = onesamplemean %THEN %DO;

		%LET _ptci_msg_ = %BQUOTE(The message was generated for DATA = &data., VAR = &var.);
		%LET _ptci_mstat_ = mean;
		%LET _ptci_nullstatlab_ = Null statistic;

		%macro call_info;
			input_data = "&data.";
			input_var = "&var.";
		%mend call_info;

		%macro save_info;
			input_data input_var 
		%mend save_info;

		%macro labels;
			input_data = "Analysis dataset"
			input_var = "Analysis variable"
		%mend labels;

			
		%macro prepare;
			/* Read data */
			use _ptci_indata_(where = (not missing(&var.))); /* Load data with complete records */
			read all var {&var.} into x;
			N = nrow(x);

			/* Create headers and labels */
			header = "Statistics";
			mstatlab = "Sample mean";
			call symput("_ptci_header_", header);
			call symput("_ptci_mstatlab_", mstatlab);
		%mend prepare;

	%END; %* END ONESAMPLEMEAN;



	****** TWOSAMPLEMEANS ******;
	%IF %QLOWCASE(&type.) = twosamplemeans %THEN %DO;

		%LET _ptci_msg_ = %BQUOTE(The message was generated for DATA = &data., VAR = &var. and GROUP = &group.);
		%LET _ptci_mstat_ = diff;
		%LET _ptci_nullstatlab_ = Null difference;

		%macro call_info;
			input_data = "&data.";
			input_var = "&var.";
			input_group = "&group.";
		%mend call_info;

		%macro save_info;
			input_data input_var input_group
		%mend save_info;

		%macro labels;
			input_data = "Analysis dataset"
			input_var = "Analysis variable"
			input_group = "Group variable"
		%mend labels;

		%macro prepare;
			/* Read data */
			use _ptci_indata_(where = (not missing(&group.) and not missing(&var.))); /* Load data with complete cases. */ 
			read all var {&group.} into grp;
			read all var {&var.} into x;

			/* Re-arrange data */
			mingrp = min(grp);
			maxgrp = max(grp);
			x1 = x[loc(grp = min(grp))];
			x2 = x[loc(grp = max(grp))];
			x = x1 // x2;
			n1 = nrow(x1);
			n2 = nrow(x2);
			N = n1 + n2;
			grp = repeat(1, n1) // repeat(2, n2);

			/* Create headers and labels*/
			header = "Comparison of means";
			if type(mingrp) = "N" then 
				mstatlab = "Mean difference ("  + strip(char(mingrp))  + " - "  + strip(char(maxgrp)) + ")";
			else 
				mstatlab = "Mean difference ("  + strip(mingrp)  + " - "  + strip(maxgrp) + ")";
			call symput("_ptci_header_", header);
			call symput("_ptci_mstatlab_", mstatlab);

			/* Check that number of groups i correct */
			ngrps = ncol(unique(grp));
			call symputx("_ptci_ngrps_", ngrps);
			%IF %SYSEVALF(&_ptci_ngrps_. ^= 2) %THEN %DO;
				quit; /* Quit IML. */
				%ptci_term;
				%PUT %SYSFUNC(cat(ERR, OR)):Invalud input argument to permtestCI. Too many groups found in the GROUP = &group. variable. GROUP should be group variable with exactly two groups. Number of groups found is &_ptci_ngrps_.. &_ptci_msg_..;
				%ABORT;
			%END;	
		%mend prepare;

		%macro compute;
			/* Observed mean difference */
			&_ptci_mstat_. = mean(x1) - mean(x2);

			/* Permutation indexes. */
			if N < 500 then ncomb = comb(N, n1); /* Do not compute for large N */
			else ncomb = nperm + 1;

			if nperm >= ncomb then do;
				pix1 = allcomb(N, n1); /* Use combinations for faster execution. */
				nperm = ncomb;
				call symputx('nperm', nperm);
				method = "Exact";
				note = "Exact confidence interval computed using all possible " + strip(char(nperm)) + " permutations.";
			end; 
			else do;
				call randseed(&seed.);
				pix1 = rancomb(N, n1, &nperm. - 1);
				pix1 = 1:n1 // pix1;
				approx = 1;
				method = "Approximate";
				note = "Approximate confidence interval computed using " + strip(char(nperm)) + " random permutations.";
			end;

			/* Permuted data. Keep track of first group only. */
			perm_x1 = shape(x[pix1], nperm, n1);
			perm_grp1 = shape(grp[pix1], nperm, n1);		

			/* Reference distribution by permutations */
			sum = sum(x);
			perm_mean_grp1 = perm_x1[,1:n1][,:];
			perm_mean_grp2 = (sum - perm_mean_grp1 * n1) / n2;
			nulldist = perm_mean_grp1 - perm_mean_grp2;
			perm_n_grp1 = (perm_grp1[,1:n1] = 1)[,+]; /* Number of observations from group 1 in each permutation */
			relch = (n1 - perm_n_grp1) / n2 - perm_n_grp1 / n1; 
		%mend compute;

	%END; %* END TWOSAMPLEMEANS;



	****** PAIREDMEANS, long format ******;
	%IF %QLOWCASE(&type.) = pairedmeans and not %_isBlank_(&visitvar.) %THEN %DO;

		%PUT  %datatyp(&vis1.);
		%IF %datatyp(&vis1.) ^= NUMERIC %THEN %DO;
			%LET vis1 = %SYSFUNC(tranwrd(&vis1., %STR(%"), %str()));
			%LET vis2 = %SYSFUNC(tranwrd(&vis2., %STR(%"), %str()));
		%END;

		/* Prepare data */
		proc sort data = _ptci_indata_;
			by &idvar. &visitvar.;
		run;
		data _ptci_dat1_ _ptci_dat2_;
			set _ptci_indata_;
			%IF %datatyp(&vis1.) = NUMERIC %THEN %DO;
				if &visitvar. eq &vis1. then output _ptci_dat1_;
				if &visitvar. eq &vis2. then output _ptci_dat2_;
			%END;
			%ELSE %DO;
				if lowcase(&visitvar.) eq lowcase("&vis1.") then output _ptci_dat1_;
				if lowcase(&visitvar.) eq lowcase("&vis2.") then output _ptci_dat2_;
			%END;
		run;

		/* Check that visits exists in input data set. */
		%LET _ptci_dsid_ = %SYSFUNC(open(_ptci_dat1_));
		%LET _ptci_nobs1_ = %SYSFUNC(attrn(&_ptci_dsid_, nobs));
		%LET _ptci_dsid_ = %SYSFUNC(close(&_ptci_dsid_));
		%IF &_ptci_nobs1_. = 0 %THEN %DO;
			%ptci_term;
			%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument VIS1 = &vis1. found. Cannot find any records with &visitvar = &vis1. in the input data set &data. &where..;
			%ABORT;
		%END;

		%LET _ptci_dsid_ = %SYSFUNC(open(_ptci_dat2_));
		%LET _ptci_nobs2_ = %SYSFUNC(attrn(&_ptci_dsid_, nobs));
		%LET _ptci_dsid_ = %SYSFUNC(close(&_ptci_dsid_));
		%IF &_ptci_nobs2_. = 0 %THEN %DO;
			%ptci_term;
			%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument VIS2 = &vis2. found. Cannot find any records with &visitvar = &vis2. in the input data set &data. &where..;
			%ABORT;
		%END;

		/* Prepare data, continued */
		data _ptci_indata_(drop = &visitvar.);
			merge _ptci_dat1_(rename = (&var. = _var1_)) _ptci_dat2_(rename = (&var. = _var2_));
			by &idvar.;
			diff = _var2_ - _var1_;
		run;

		/* Check data */
		proc sql;
			create table _n_ as
			select count(diff) as n
			from _ptci_indata_(where = (not missing(diff)));
		quit;
		data _null_;
			set _n_;
			if n < 1 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI, TYPE = PAIREDMEANS. There are no valid observations.";
				call symput('_ptci_abort_', 1);
			end;
		run;
		proc sql noprint;
			drop table _n_;
		quit;	
		
		%IF &_ptci_abort_. %THEN %DO; 
			%ptci_term;
			%ABORT; 
		%END;


		%LET _ptci_msg_ = %BQUOTE(The message was generated for DATA = &data., VAR = &var., IDVAR = &idvar., VISITVAR = &visitvar., VIS1 = &vis1., VIS2 = &vis2.);
		%LET _ptci_mstat_ = diff;
		%LET _ptci_nullstatlab_ = Null difference;

		%macro call_info;
			input_data = "&data.";
			input_var = "&var.";
			input_visitvar = "&visitvar.";
			input_visits = "&vis1., &vis2.";
			%IF %datatyp(&visitvar.) = NUMERIC %THEN %DO;
				vis1 = &vis1.;
				vis2 = &vis2.;
			%END;
			%ELSE %DO;
				vis1 = "&vis1.";
				vis2 = "&vis2.";
			%END;
	
		%mend call_info;

		%macro save_info;
			input_data input_var input_visitvar input_visits
		%mend save_info;

		%macro labels;
			input_data = "Analysis dataset"
			input_var = "Analysis variable"
			input_visitvar = "Visit (measurement) variable"
			input_visits = "Visits (measurements)"
		%mend labels;
			
		%macro prepare;
			/* Read data */
			use _ptci_indata_(where = (not missing(diff)));
			read all var {diff} into x;
			N = nrow(x);		

			/* Create headers and labels */
			header = "Comparison of means";

			if type(vis1) = "N" then 
				mstatlab = "Mean difference ("  + strip(char(vis2))  + " - "  + strip(char(vis1)) + ")";
			else 
				mstatlab = "Mean difference ("  + strip(vis2)  + " - "  + strip(vis1) + ")";

			call symput("_ptci_header_", header);
			call symput("_ptci_mstatlab_", mstatlab);
		%mend prepare;
	
	%END; %* END PAIREDMEANS in long (advl) format;



	****** PAIREDMEANS, wide format ******;
	%IF %QLOWCASE(&type.) = pairedmeans and %_isBlank_(&visitvar.) %THEN %DO;

		/* Prepare data */
		data _ptci_indata_;
			set _ptci_indata_;
			diff = &var2. - &var1.;
		run;


		%LET _ptci_msg_ = %BQUOTE(The message was generated for DATA = &data., VAR1 = &var1., VAR2 = &var2.);
		%LET _ptci_mstat_ = diff;
		%LET _ptci_nullstatlab_ = Null difference;

		%macro call_info;
			input_data = "&data.";
			input_var = "&var1., &var2.";
		%mend call_info;

		%macro save_info;
			input_data input_var 
		%mend save_info;

		%macro labels;
			input_data = "Analysis dataset"
			input_var = "Analysis variables"
		%mend labels;

			
		%macro prepare;
			/* Read data */
			use _ptci_indata_(where = (not missing(diff))); /* Read complete records. */ 
			read all var {diff} into x;
			N = nrow(x);

			/* Create headers and labels */
			header = "Comparison of means";
			mstatlab = "Mean difference ("  + strip("&var2.")  + " - "  + strip("&var1.") + ")";
			call symput("_ptci_header_", header);
			call symput("_ptci_mstatlab_", mstatlab);
		%mend prepare;
	
	%END; %* END PAIREDMEANS in wide (adsl) format;



	****** ONESAMPLEMEAN and PAIREDMEANS, compute block. ******;
	%IF (%QLOWCASE(&type.) in pairedmeans onesamplemean) %THEN %DO;
		%macro compute;
			/* Module: create all cobminations of sign flippings. */
			/* From https://blogs.sas.com/content/iml/2011/01/05/creating-a-matrix-with-all-combinations-of-zeros-and-ones.html */
			start GetAllComb(n);
				rows = 2##n;
				x = j(rows, n);
				do j = 1 to n;
				PatLength = 2##(n-j); 
				PatRepl = 2##(j-1);   
				pattern = j(PatLength,1,1) // j(PatLength,1,-1);
				x[,j] = repeat(pattern, PatRepl);
				end;
				return( x );
			finish;

			/* Observed mean statistic */
			&_ptci_mstat_. = mean(x);

			/* Create sign flipping matrix. */
			if N < 100 then ncomb = 2 ** N; /* Do not compute for large N */
			else ncomb = nperm + 1;

			if nperm >= ncomb then do;		 
				sf = GetAllComb(N);
				nperm = ncomb;
				call symputx('nperm', nperm);
				method = "Exact";
				note = "Exact confidence interval computed using all possible " + strip(char(nperm)) + " permutations.";
			end;
			else do;
				sf = j(&nperm. - 1, N,.);
				call randseed(&seed.);
				call randgen(sf, 'BERNOULLI', 0.5);

				sf = 2 * sf - 1;
				sf = j(1, N, 1) // sf;
				approx = 1;
				method = "Approximate";
				note = "Approximate confidence interval computed using " + strip(char(nperm)) + " random permutations.";
			end;


			/* Reference distribution by permutations */
			nulldist = sf * x / N;
			mean_sgn = sf[,:];
			relch = - mean_sgn;
		%mend compute;
	%END; %* END Compute block.;



	ods proclabel 'PermTestCI';
	proc iml;
		/* For debugging: reset print; */ 

		/* Init */
		clevel = &clevel.;
		alpha = 1 - clevel;
		call symputx('_ptci_alpha_', alpha);
		theta0 = &theta0.;
		nperm = &nperm.;
		warn1 = " ";
		warn2 = " ";
		approx = 0;

		%call_info;

		%prepare;


		/* Compute default number of permutations if not supplied to macro. */
		if nperm < 0 then do;
			if N <= 100 then nperm = 10000; 
			else if N <= 200 then nperm = 5000; 
			else if N <= 500 then nperm = 2000; 
			else nperm = 1000; 
		end;
		call symputx('nperm', nperm);


		%compute;

			
		/* Translations */
		ix = loc(relch > -1);
		z = j(sum(relch <= -1), 1, .M);
		c = z // (&_ptci_mstat_. - nulldist[ix]) / (1 + relch[ix]);
		call sort(c);
		call sort(c); /* Need to sort twice, for some reason... */


		/* p-value */
		mstat0 = &_ptci_mstat_. - theta0;
		nulldist[1] = &_ptci_mstat_.; /* First record in nulldist contains observed statistic. */ 
		nulldist0 = nulldist + relch * theta0;
		pval = 2 * min(mean(nulldist0 <= (mstat0 + 0.000000000001)), mean(nulldist0 >= (mstat0 - 0.000000000001))); /* Ad hoc: add/subtract 1E-12 to account for rounding errors. */
		pval = min(1, pval);
		se = sqrt(pval * (1 - pval) / nperm);
		if approx & ( (pval > alpha & (pval - quantile('normal', 0.975) *se) < alpha) | (pval < alpha & (pval + quantile('normal', 0.975) * se) > alpha) ) then do;
			warn1 = "WARNING: In PermTestCI. Approximate p-value using &nperm. permutations has 95% confidence interval containing alpha = &_ptci_alpha_.. Increase the number of permutations of avoid this warning.";
		end;


		/* Confidence interval */
		pmin = 2 / nperm;
		clmin = 1 - 2 * pmin;
		npermmin = 2 / alpha;
		call symputx("_ptci_clmin_", clmin);
		call symputx("_ptci_npermmin_", npermmin);
		if pmin >= alpha then do;
			lcl = -1E300;
			ucl = 1E300;
			warn2 = "WARNING: In PermTestCI. Confidence interval can not be calculated at the desired confidence level. Confidence interval covers the real line.";
		end;
		else do;
			c = c[loc(c ^= .)];
			lix = ceil(nrow(c) * alpha / 2);
			uix = nrow(c) - lix + 1;
			lcl = c[lix - 1] + 0.000000000001; /* Small increment. */
			ucl = c[uix + 1] - 0.000000000001; /* Small reduction. */
			if N = 6 & lowcase("&type.") =  "onesamplemean" then do;
	 			warn2 = "WARNING: In PermTestCI. Confidence interval equals data range for TYPE = ONESAMPLEMEAN with N = 6 observations.";
			end;
			if N = 6 & lowcase("&type.") =  "pairedmeans" then do;
	 			warn2 = "WARNING: In PermTestCI. Confidence interval equals data range for TYPE = PAIREDMEANS with N = 6 observations.";
			end;
		end;	

		/* Write results to SAS data set &oudata. */
		create &outdata. var {%save_info N &_ptci_mstat_. lcl ucl clevel theta0 pval nperm method note warn1 warn2}; 
		append;      
		close &outdata.; 
	quit;


	/* Format ouput data set */
	proc datasets nolist;
		modify &outdata.;
			label 	&_ptci_mstat_. = "&_ptci_mstatlab_."
					lcl = "Lower confidence limit"
					ucl = "Upper confidence limit"
					theta0 = "&_ptci_nullstatlab_."
					pval = "p-value"
					clevel = "Confidence level"
					nperm = "Number of permutations used"
					N = "Number of observations used"
					method = "Method"
					note = "Note"
					%labels;
		format 	&_ptci_mstat_. 8.&decimals.
				lcl 8.&decimals.
				ucl 8.&decimals.
				theta0 8.&decimals.
				pval pvalue6.4;
	quit;	


	/* Print general information and results to output. */
	%IF (%QUPCASE(&print.) in Y YES) %THEN %DO;

		/* General information. */ 
		proc transpose data = &outdata. out = _info_;
			var %save_info clevel n nperm method;
		run;

		proc datasets nolist;
			modify _info_;
			label 	_name_ = "Name"
					_label_ = "Label"
					col1 = "Value";
		quit;

		data _info_;
			set _info_;
			col1 = strip(col1);
		run;

		ods proclabel 'PermTestCI';
		title%EVAL(&titnum. + 3) "General information";
		proc report data = _info_ noheader;
			column _label_ col1;
			define _label_ / width = 30;
			define col1 / width = 30 flow;
		run;


		/* Results. */
		ods proclabel 'PermTestCI';
		title%EVAL(&titnum. + 3) "&_ptci_header_.";
		proc print data = &outdata. label noobs;
			var &_ptci_mstat_. lcl ucl theta0 pval;
		run;

		proc datasets nolist;
			delete _info_;
		quit;
	%END;

	/* Drop unnecessary variables from output data set. */
	data &outdata.;
		set &outdata.;
		drop N %save_info;
	run;

	/* Print warnings to log and output. */
	proc transpose data = &outdata. out = _ptci_warn_;
		var warn1 warn2;
	run;

	proc datasets nolist;
		modify _ptci_warn_;
		label 	_name_ = "Name"
				col1 = "Warning";
	quit;

	title%EVAL(&titnum. + 3);
    data _null_ ;
		set _ptci_warn_;
		where col1 ne " "; 
		%IF (%QUPCASE(&print.) in Y YES) %THEN %DO;
		    file print;
		   	put col1;
		%END;
		file log;
		put col1;
    run;

	/* Time elapsed during execution. */
	%LET _ptci_edtm_ = %SYSFUNC(datetime());
	%LET _ptci_runtm_ = %SYSFUNC(putn(&_ptci_edtm_ - &_ptci_sdtm_, 12.4));


	/* Clean-up */
	%ptci_term;


	/* Print completion info to log. */
	data _null_ ;
	    file log;
	   	put "NOTE: Macro permtestCI completed in &_ptci_runtm_. seconds. &_ptci_msg_. using &nperm. permutations.";
    run;


%mend PermTestCI;


%macro _isBlank_(param);   
	%sysevalf(%superq(param)=,boolean) 
%mend _isBlank_; 


%macro ptci_argcheck;
/* Check input arguments to PermTestCI macro. */

	options minoperator;
	%LET _ptci_abort_ = 0;


	/* Check that required input parameter are supplied. */
	%IF %_isBlank_(&data.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)): Invalid input argument to PermTestCI. DATA argument is missing. Please specify an analysis data set.;
		%ABORT;
	%END;

	%IF %_isBlank_(&type.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. TYPE parameter is missing. Valid arguments to TYPE are ONESAMPLEMEAN, TWOSAMPLEMEANS and PAIREDMEANS.;
		%ABORT;
	%END;

	%IF %_isBlank_(&nperm.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input NPERM = &nperm.. Please specify a positive integer. Omit from macro call to use default value.;
		%ABORT;
	%END;

	%IF %_isBlank_(&outdata.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. OUTDATA parameter is missing. Please specify an output data set. Omit from macro call to use default value. Default value is _PERMTESTCI_.;
		%ABORT;
	%END;

	%IF %_isBlank_(&decimals.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input DECIMALS = &decimals.. Please specify a positive integer. Omit from macro call to use default value. Default value is 4.;
		%ABORT;
	%END;

	%IF %_isBlank_(&theta0.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input THETA0 = &theta0.. Please specify a numeric value. Omit from macro call to use default value. Default value is 0.;
		%ABORT;
	%END;

	%IF %_isBlank_(&clevel.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input CLEVEL = &clevel.. Please specify a numeric value between 0 and 1. Omit from macro call to use default value.. Default value is 0.95.;
		%ABORT;
	%END;

	%IF %_isBlank_(&seed.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input SEED = &seed.. Please specify a positive integer. Omit from macro call to use default value. Default value is 123456.;
		%ABORT;
	%END;


	/* Check that data types are correct. */
	%IF %datatyp(&data.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input DATA = &data.. Please specify an input data set.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&data.) = _null_ %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input DATA = &data.. Please specify an input data set.;
		%ABORT;
	%END;

	%IF %datatyp(&where.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input WHERE = &where.. Please specify a valid data step where clause.;
		%ABORT;
	%END;

	%IF %datatyp(&type.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input TYPE = &type.. Valid types are ONESAMPLEMEAN, TWOSAMPLEMEANS and PAIREDMEANS.;
		%ABORT;
	%END;

	%IF %datatyp(&var.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input VAR = &var.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&group.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input GROUP = &group.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&idvar.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input IDVAR = &idvar.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&visitvar.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input VISITVAR = &visitvar.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&var1.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input VAR1 = &var1.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&var2.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input VAR2 = &var2.. Please specify a valid variable name.;
		%ABORT;
	%END;

	%IF %datatyp(&outdata.) ^= CHAR %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input OUTDATA = &outdata.. Please specify a valid dataset name.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&outdata.) = _null_ %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input OUTDATA = &outdata.. Please specify an output data set.;
		%ABORT;
	%END;

	%IF &nperm. ne . and (%datatyp(&nperm.) ^= NUMERIC or %SYSFUNC(round(&nperm.)) ne &nperm. or &nperm. < 0) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input NPERM = &nperm.. Please specify a positive integer.;
		%ABORT;
	%END;

	%IF %datatyp(&decimals.) ^= NUMERIC or %SYSFUNC(round(&decimals.)) ne &decimals. or &decimals. < 0 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input DECIMALS = &decimals.. Please specify a positive integer.;
		%ABORT;
	%END;

	%IF %datatyp(&theta0.) ne NUMERIC %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input THETA0 = &theta0.. Please specify a numeric value.;
		%ABORT;
	%END;

	%IF %datatyp(&clevel.) ne NUMERIC or &clevel. <= 0 or &clevel. >= 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input CLEVEL = &clevel.. Please specify a numeric value between 0 and 1.;
		%ABORT;
	%END;

	%IF %datatyp(&seed.) ^= NUMERIC or %SYSFUNC(round(&seed.)) ne &seed. or &seed. < 0 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input SEED = &seed.. Please specify a positive integer.;
		%ABORT;
	%END;



	/* Check that correct number of parameters are supplied. */
	%IF %SYSFUNC(countw(&data., "!$%&()*+,-/;<^| ")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input parameter DATA = &data. found. Too many analysis data sets specified. Please specify a single analysis data set.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw(&type.)) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid input parameter TYPE = &type. found. Too many types specified. Please specify a single test type.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&var.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid parameter VAR = &var. found. Too many analysis variables specified. Please specify a single analysis variable.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&group.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid parameter GROUP = &group. found. Too many group variables specified. Please specify a single group variable.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&idvar.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument IDVAR = &idvar. found. Too many id-variables specified. Please specify a single id-variable.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&visitvar.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument VISITVAR = &visitvar. found. Too many visit-variables specified. Please specify a single visit variable.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&var1.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument VAR1 = &var1. found. Too many analysis variables specified. Please specify a single analysis variable.;
		%ABORT;
	%END;

	%IF %SYSFUNC(countw("&var2.")) > 1 %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument VAR2 = &var2. found. Too many analysis variables specified. Please specify a single analysis variable.;
		%ABORT;
	%END;


	/* Check that analysis dataset and variables exists. */
	%IF not %SYSFUNC(exist(&data.)) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument DATA = &data. found. Input data set &data. does not exist.;
		%ABORT;
	%END;

	%IF not (%QLOWCASE(&type.) in onesamplemean twosamplemeans pairedmeans) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. TYPE = &type. is not recognized. Valid types are ONESAMPLEMEAN, TWOSAMPLEMEANS and PAIREDMEANS.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = twosamplemeans %THEN %DO;
		data _null_;
			dsid = open("&data.");
			if varnum(dsid, "&var.") = 0 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VAR = &var. found. Variable &var. does not exist in analysis data set &data..";
				call symput('_ptci_abort_', 1);
			end;
			if varnum(dsid, "&group.") = 0 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument GROUP = &group. found. Variable &group. does not exist in analysis data set &data..";
				call symput('_ptci_abort_', 1);
			end;
		run;
	%END;

	%IF %QLOWCASE(&type.) = onesamplemean %THEN %DO;
	data _null_;
		dsid = open("&data.");
		if varnum(dsid, "&var.") = 0 then do;
			put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VAR = &var. found. Variable &var. does not exist in analysis data set &data..";
			call symput('_ptci_abort_', 1);
		end;
	run;
	%END;

	%IF %QLOWCASE(&type.) = pairedmeans %THEN %DO;
		%IF %_isBlank_(&visitvar.) %THEN %DO;
			data _null_;
				dsid = open("&data.");
				if varnum(dsid, "&var1.") = 0 then do;
					put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VAR1 = &var1. found. Variable &var1. does not exist in analysis data set &data..";
					call symput('_ptci_abort_', 1);				
				end;
				if varnum(dsid, "&var2.") = 0 then do;
					put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VAR2 = &var2. found. Variable &var2. does not exist in analysis data set &data..";
					call symput('_ptci_abort_', 1);
				end;
			run;
		%END;
		%ELSE %DO;
			data _null_;
				dsid = open("&data.");
				if varnum(dsid, "&var.") = 0 then do;
					put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VAR = &var. found. Variable &var. does not exist in analysis data set &data..";
					call symput('_ptci_abort_', 1);
				end;
				if varnum(dsid, "&visitvar.") = 0 then do;
					put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument VISITVAR = &visitvar.. found. Variable &visitvar. does not exist in analysis data set &data..";
					call symput('_ptci_abort_', 1);
				end;
				if varnum(dsid, "&idvar.") = 0 then do;
					put "ERR" "OR: Invalid input argument to PermTestCI. Invalid argument IDVAR = &idvar.. found. Variable &idvar. does not exist in analysis data set &data..";
					call symput('_ptci_abort_', 1);
				end;
			run;
		%END;	
	%END;

	%IF &_ptci_abort_. %THEN %DO; 		
		options &saveOptions.; 
		%ABORT; 
	%END;


	/* Check that certain combinations of input paramters are correctly used. */
	%IF (%QLOWCASE(&type.) in onesamplemean twosamplemeans) and %_isBlank_(&var.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. but no analysis variable specified in the VAR argument.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = twosamplemeans and %_isBlank_(&group.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. but GROUP variable is missing. A GROUP variable is required for TYPE = &type..;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = twosamplemeans and not %_isBlank_(&visitvar.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. and VISITVAR = &visitvar.. TYPE = &type. and VISITVAR are incompatible. TYPE = TWOSAMPLEMEANS requires a GROUP variable. VISITVAR is only used for paired samples, specified by TYPE = PAIREDMEANS.;
		%ABORT;
	%END;

	%IF (%QLOWCASE(&type.) in onesamplemean twosamplemeans) and (not %_isBlank_(&var1.) or not %_isBlank_(&var2.)) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VAR1 = &var1. and VAR2 = &var2.. TYPE = &type. is not compatible with VAR1 and VAR2. TYPE = &type. requires a single analysis variable, specified by the VAR argument. VAR1 and VAR2 may be used with TYPE = PAIREDMEANS for paired data in wide (adsl) format.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = onesamplemean and not %_isBlank_(&group.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. and GROUP = &group.. TYPE = &type. is not compatible with a GROUP variable. A GROUP variable is only allowed for TYPE = TWOSAMPLEMEANS.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = onesamplemean and not %_isBlank_(&visitvar.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. and VISITVAR = &visitvar.. TYPE = &type. and VISITVAR are incompatible. VISITVAR is only used for paired samples, specified by TYPE = PAIREDMEANS.;
		%ABORT;
	%END;	

	%IF %QLOWCASE(&type.) = pairedmeans and not %_isBlank_(&group.) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type. and GROUP = &group.. TYPE = &type. is not compatible with a GROUP variable. A VISITVAR may be specified with TYPE = &type. for data in long (advl) format. A GROUP variable is only allowed for TYPE = TWOSAMPLEMEANS.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = pairedmeans and (%_isBlank_(&visitvar.) and (%_isBlank_(&var1.) or %_isBlank_(&var2.))) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., VAR = &var., VAR1 = &var1., VAR2 = &var2.. TYPE = &type. requires VISITVAR and VAR to be specified for data in long (advl) format, or no VAR1 and VAR2 to be specified for data in wide (adsl) format.;
		%ABORT;
	%END;	

	%IF %QLOWCASE(&type.) = pairedmeans and (not %_isBlank_(&visitvar.) and  (not %_isBlank_(&var1.) or not %_isBlank_(&var2.))) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., VAR1 = &var1., VAR2 = &var2.. VISITVAR may only be used in combination with the VAR argument for data in long (advl) format. The VAR1 and VAR2 arguments are used for data in wide (adsl) format.;
		%ABORT;
	%END;	

	%IF %QLOWCASE(&type.) = pairedmeans and (not %_isBlank_(&visitvar.) and %_isBlank_(&var.)) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., VAR = &var.. TYPE = &type. requires an analysis variable to be speficied in the VAR argument when a VISITVAR is supplied.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = pairedmeans and (not %_isBlank_(&visitvar.) and %_isBlank_(&idvar.)) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., IDVAR = &idvar.. TYPE = &type. requires an ID variable to be speficied in the IDVAR argument when a VISITVAR is supplied.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = pairedmeans and (not %_isBlank_(&visitvar.) and (%_isBlank_(&vis1.) or %_isBlank_(&vis2.))) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., VIS1 = &vis1., VIS2 = &vis2.. TYPE = &type. two unique visit values to be specified in the VIS1 and VIS2 arguments when a VISITVAR is supplied.;
		%ABORT;
	%END;

	%IF %QLOWCASE(&type.) = pairedmeans and (%_isBlank_(&visitvar.) and not %_isBlank_(&idvar.)) %THEN %DO;
		options &saveOptions.; 
		%PUT %SYSFUNC(cat(ERR, OR)):Invalid input arguments to PermTestCI. TYPE = &type., VISITVAR = &visitvar., IDVAR = &idvar.. TYPE = &type. requires an visit variable to be specified in VISITVAR argumen when an IDVAR supplied.;
		%ABORT;
	%END;


	/* Check where clause. */
	%IF not %_isBlank_(&where.) %THEN %DO;
		%IF %SYSFUNC(length(%BQUOTE(&where.))) < 5 %THEN %DO;
			options &saveOptions.; 
			%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument WHERE = &where. found. Argument specified to WHERE must be a valid data step where clause.;
			%ABORT;
		%END;
		%ELSE %IF %LOWCASE(%SUBSTR("&where"., 2, 6)) ^= where %THEN %DO;
			options &saveOptions.; 
			%PUT %SYSFUNC(cat(ERR, OR)):Invalid input argument to PermTestCI. Invalid argument WHERE = &where. found. Argument specified to WHERE must be a valid data step where clause.;
			%ABORT;
		%END;
	%END;

	/* Check that valid data is supplied. */
	%IF %QLOWCASE(&type.) = twosamplemeans %THEN %DO;
		proc sql;
			create table _n_ as
			select 	count(n) as ngrps
			from 	(	select count(&var.)	 as n
						from &data.(where = (not missing(&var.)))
						&where.
						group by &group.);
		quit;
		data _null_;
			set _n_;
			if ngrps < 2 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI, TYPE = TWOSAMPLEMEANS. The GROUP variable does not have two levels with valid data.";
				call symput('_ptci_abort_', 1);
			end;
			else if ngrps > 2 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI, TYPE = TWOSAMPLEMEANS. The GROUP variable has more than two levels.";
				call symput('_ptci_abort_', 1);
			end;
		run;
		proc sql noprint;
			drop table _n_;
		quit;
	%END;
	%IF %QLOWCASE(&type.) = onesamplemean %THEN %DO;
		proc sql;
			create table _n_ as
			select count(&var.)	as n
			from &data.(where = (not missing(&var.)))
			&where.;
		quit;
		data _null_;
			set _n_;
			if n < 1 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI, TYPE = ONESAMPLEMEAN. There are no valid observations.";
				call symput('_ptci_abort_', 1);
			end;
		run;
		proc sql noprint;
			drop table _n_;
		quit;	
	%END;
	%IF %QLOWCASE(&type.) = pairedmeans and %_isBlank_(&visitvar.) %THEN %DO;
		proc sql;
			create table _n_ as
			select count(&var1.)	as n
			from &data.(where = (not missing(&var1.) and not missing(&var2.)))
			&where.;
		quit;
		data _null_;
			set _n_;
			if n < 1 then do;
				put "ERR" "OR: Invalid input argument to PermTestCI, TYPE = PAIREDMEANS. There are no valid observations.";
				call symput('_ptci_abort_', 1);
			end;
		run;
		proc sql noprint;
			drop table _n_;
		quit;	
	%END;
	
	%IF &_ptci_abort_. %THEN %DO; 		
		options &saveOptions.; 
		%ABORT; 
	%END;

%mend ptci_argcheck;



%macro ptci_term;
/* Reset titles and remove temporary datasets. */

	%macro settitle(num, text);
		title&num. &text.;
	%mend settitle;
	
	title;
	data _null_;
		set _ptci_oldtit_;
		call execute('%settitle('||number||','||text||')');
	run;

	proc datasets nolist;
		delete _ptci_indata_ _ptci_warn_ _ptci_oldtit_ %IF %SYSFUNC(exist(_ptci_dat1_)) %THEN %DO; _ptci_dat1_ _ptci_dat2_ %END;;
	quit;

	options &saveOptions.; 

%mend ptci_term;
