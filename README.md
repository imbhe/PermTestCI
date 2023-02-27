# PermTestCI
SAS macro for computation of confidence intervals for the mean of one sample and for the mean difference between two independent or paired samples, by inversion of the Fisher-Pitman non-parametric permutation test.

#### Requirements
SAS/IML software.

#### Note of performance
To reduce computation time, the number of permutations by default ranges from 1,000 (n>500) to 10,000 (n≤100), where n is the number of observations in the dataset. The number of permutations can be set manually using the NPERM argument if a larger number of permutations is desired. 

It is strongly recommended to have the macrogen, mlogicmprint, and symbolgen options turned OFF when running the macro. Printing macro execution information to the log may otherwise result in a substantial increase of computation time and is recommended for debugging purposes only.

## Running the macro with default options

#### One sample
```
%PermTestCI2(DATA = my_data, /* Analysis dataset */
             TYPE = ONESAMPLEMEAN, /* Confidence interval for mean */
             VAR = my_var /* Analysis variable*/
);
``` 

#### Two independent samples
```
%PermTestCI2(DATA = my_data, /* Analysis dataset */
             TYPE = TWOSAMPLEMEANS, /* CI for mean difference*/
             VAR = my_var, /* Analysis variable*/
             GROUP = my_grp /* Group variable */
);
``` 

#### Two paired samples, data in wide (ADSL) format
```
%PermTestCI2(DATA = my_adsl, /* Analysis dataset, wide format */
             TYPE = PAIREDMEANS, /* CI for mean difference*/
             VAR1 = my_var_vis1, /* Analysis variable, first visit */
             VAR2 = my_var_vis2, /* Analysis variable, second visit */
);
``` 

#### Two paired samples, data in long (ADVL) format
```
%PermTestCI2(DATA = my_advl, /* Analysis dataset, long format */
             TYPE = PAIREDMEANS, /* CI for mean difference*/
             VAR = my_var_vis1, /* Analysis variable */
             VISITVAR = my_visit_var, /* Visit variable */
             IDVAR = my_id_var, /* ID variable */
             VIS1 = val1, /* Value of VISITVAR for first visit */
             VIS2 = val2 /* Value of VISITVAR for second visit */
);
``` 
Additional examples may be found in examples.sas.

#### Results
The macro computes a 95% confidence interval for the mean of one sample or for the mean difference between two independent or paired samples, as specified by the TYPE parameter. 

When TYPE = TWOSAMPLEMEANS, the difference between sample means is computed by subtracting the mean of the highest unformatted group from the mean of the lowest unformatted group, similar to the default option in the TTEST procedure. 

When TYPE = PAIREDMEANS, the difference is computed by subtracting the mean of the first visit/measurement, as specified by VAR1 (data in wide “ADSL” format) or by the VIS1-value of VISITVAR (data in long “ADVL” format), from the mean of the second visit, as specified by VAR2 (data in wide “ADSL” format) or by the VIS2-value of VISITVAR (data in long “ADVL” format).

Output is printed to all open ODS destinations. Additionally, results are stored in the work.\_permtestci\_ dataset, containing the following variables:
* MEAN: the one-sample mean. 
Only included in output data set if TYPE = ONESAMPLEMEAN.
* DIFF: the difference between the two sample means. 
Only included in output data set if TYPE = TWOSAMPLEMEANS or PAIREDMEANS.
*	LCL: lower confidence limit of mean or mean difference.
*	UCL: upper confidence limit of mean or mean difference.
*	THETA0: value of the statistic (mean or difference between means) under the null hypothesis. By default, THETA0 is set to 0. Affect the p-value of the test but does not affect the confidence limits.
*	PVAL: two-sided p-value for testing the null hypothesis that the statistic (mean or difference between means) equals THETA0. 
*	CLEVEL: the confidence level of the confidence interval.
*	NPERM: number of permutations used for computation of the confidence interval.
*	METHOD: method used for computation of confidence interval. Equals “Approximate” when an approximate confidence interval using NPERM random permutations is computed and “Exact” when an exact confidence interval using all possible permutations is computed. 
*	NOTE: details on method and number of permutations used for computation of confidence interval.
*	WARN1: Contains a warning if the p-value computed from the permutation test has a 95% confidence interval that covers the value 1-CLEVEL, i.e., does not lie strictly below the desired significance level. Empty string otherwise.
*	WARN2: Contains a warning if the confidence interval covers the real line or equals data range (TYPE = ONESAMPLEMEANS or PAIREDMEANS with N = 6 observations). Empty string otherwise.
*	NPERM: The number of permutations used for computation of the confidence interval. By default chosen as
    - 10 000 if N ≤ 100
    -	5000 if 100 < N ≤ 200
    -	2000 if 200 < N ≤ 500
    -	1000 otherwise,    
 where N is the number of subjects in the analysis data set. If this number is greater than the total number of possible permutations, the macro uses all possible permutations rather than a random selection of permutations. 

The macro issues a warning if:
*	The p-value computed from the permutation test has a 95% confidence interval that covers the value 1-CLEVEL, i.e., does not lie strictly below the desired significance level. In this case, results are ambiguous and it is recommended to increase the number of permutations. 
*	The confidence interval covers the real line. This occurs if the number of permutations is insufficient for computation of a confidence interval with the desired confidence level. In this case, it is recommended to increase the number of permutations. If the maximum number of permutations already is reached, it is recommended to reduce the confidence level. 

Warnings are printed to the log and to all open ODS destinations. Warnings are also stored in the WARN1 and WARN2 variables of the output data set. 

Upon completion, the time elapsed during execution of the macro is printed to the log. 
