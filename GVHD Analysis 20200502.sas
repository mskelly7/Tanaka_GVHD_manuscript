*TITLE: Evaluating associations between antibiotic exposures and graft-versus- host 
disease in hematopoietic stem cell transplant recipients

Biostatistician: Rebecca Young, rty5@duke.edu
PI: Matthew Kelly
First Author: John Tanaka

Secondary Outcomes:
Secondary outcomes will include:  acute GVHD mortality, chronic GVHD diagnosis,
and chronic GVHD mortality. 
	**Two exploratory secondary outcomes were added (not in SAP). These are acute GVHD in skin only, and acute GVHD 
	in liver or gut only (no skin involvement)

Secondary hypotheses: receipt of an antibiotic regimen with an anaerobic spectrum of activity 
(piperacillin- tazobactam OR carbapenem) will be associated with:
-	acute GVHD diagnosis involving the gut and/or liver (with or without skin involvement)
-	1-year acute GVHD mortality
-	chronic GVHD diagnosis
-	5-year chronic GVHD mortality


Statistical Analysis Method: Our primary analysis will use Cox proportional hazards regression to compare the hazard 
of liver or gut acute GVHD in patients who received an anaerobic antibiotic regimen (group 1) to patients who received 
a non-anaerobic antibiotic regimen only (group 2). This will be an adjusted model including the following covariates: 
year of HSCT, age, sex, race, HSCT indication, HSCT donor/source, HLA matching, preparative regimen intensity, and GVHD 
prophylaxis regimen. All covariates will be assessed at the time of the HSCT, and therefore will not be time-varying. 
Unadjusted Cox proportional hazard models will be run as a sensitivity analysis. 

PART ONE: PRIMARY HYPOTHESIS

Primary outcomes:
A.) Acute GVHD of gut or liver only

PART TWO: Secondary Hypothesis - the following will happen for each of the secondary outcomes

a.) Extended Cox Model to check proportional hazards assumption of primary predictor
b.) Extended Cox Model to check PH in covariates
c.) Primary Model: Adjusted Cox Model
d.) Sensitivity Analysis: Unadjusted Cox Model
e.) Sensitivity Analysis: Adjusted Cox Model in adults only
f.) Sensitivity Analysis: Adjusted Cox Model in children only


Secondary Outcomes:

2.) Acute GVHD of skin only
3.) Acute GVHD of gut or liver, regardless of skin
4.) Acute GVHD Mortality
5.) Chronic GVHD
6.) Chronic GVHD Mortality
7.) Acute GVHD of skin (any skin)
8.) Any non-GVHD mortality
*;

*Loading data*;

libname raw 'U:\Rebecca\John\Raw Data';
libname ana 'U:\Rebecca\John\Analytic Data';

*Setting covariates*;

%let catcovariates = female race_white hememalignancy  myeloablative gvhdprophygroup hsct_donor hsct_source;
%let contcovariates = age_dot matches bmt_year;

*Loading analytic data file*;
data a;
set ana.gvhd_analytic_20191023;
run;

**************************************
PRIMARY MODEL
*************************************;

*********************************************************
*1.) Acute GVHD of gut or liver (regardless of skin)
*********************************************************;

*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_agvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_agvhd;
race_whitet = race_white * lof_agvhd;
hememalignancyt = hememalignancy * lof_agvhd;
transplanttypet = transplanttype * lof_agvhd;
myeloablativet = myeloablative * lof_agvhd;
age_dott = age_dot * lof_agvhd;
matchest = matches * lof_agvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the main hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only';
run;

proc freq data=a; tables aGVHD_livergut*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of aGVHD in liver or gut only';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

proc freq data=a; tables aGVHD_livergut*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut , kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergut*aGVHD_livergut(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of aGVHD in liver or gut , kids only';
where kid = 1;
run;

proc freq data=a; tables aGVHD_livergut*abx_class; where kid = 1; run;


proc freq data=a; tables aGVHD_livergutonly*abx_class; where kid = 1; run;


*********************************************************
*2.) Acute GVHD of skin only
*********************************************************;



*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_agvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;


proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_agvhd;
race_whitet = race_white * lof_agvhd;
hememalignancyt = hememalignancy * lof_agvhd;
transplanttypet = transplanttype * lof_agvhd;
myeloablativet = myeloablative * lof_agvhd;
age_dott = age_dot * lof_agvhd;
matchest = matches * lof_agvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in skin only';
run;

proc freq data=a; tables aGVHD_skinonly*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of aGVHD in skin only';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in skin only, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of aGVHD in skin only, adults only';
where kid = 0;
run;

proc freq data=a; tables aGVHD_skinonly*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in skin only, kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_skinonly*aGVHD_skinonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of aGVHD in skin only, kids only';
where kid = 1;
run;

proc freq data=a; tables aGVHD_skinonly*abx_class; where kid = 1; run;

proc freq data=a; tables died_agvhd*abx_class; where kid=0; run;


*********************************************************
*4. Acute GVHD mortality
*********************************************************;



*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_agvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;


proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_agvhd;
race_whitet = race_white * lof_agvhd;
hememalignancyt = hememalignancy * lof_agvhd;
transplanttypet = transplanttype * lof_agvhd;
myeloablativet = myeloablative * lof_agvhd;
age_dott = age_dot * lof_agvhd;
matchest = matches * lof_agvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD mortality';
run;

proc freq data=a; tables died_agvhd*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of aGVHD mortality';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD mortality, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of aGVHD mortality, adults only';
where kid = 0;
run;

proc freq data=a; tables died_agvhd*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD mortality, kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_mortality*died_agvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of aGVHD mortality, kids only';
where kid = 1;
run;

proc freq data=a; tables died_agvhd*abx_class; where kid=1; run;



*********************************************************
*5. Chronic GVHD incidence
*********************************************************;



*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd*cgvhd(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'chronic GVHD Proportionality Test';
abx_classt = abx_class * lof_cgvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;


proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd*cgvhd(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'chronic GVHD Proportionality Test';
femalet = female * lof_cgvhd;
race_whitet = race_white * lof_cgvhd;
hememalignancyt = hememalignancy * lof_cgvhd;
transplanttypet = transplanttype * lof_cgvhd;
myeloablativet = myeloablative * lof_cgvhd;
age_dott = age_dot * lof_cgvhd;
matchest = matches * lof_cgvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd*cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd';
run;

proc freq data=a; tables cgvhd*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd*cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of cgvhd ';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd*cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd , adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd*cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of cgvhd , adults only';
where kid = 0;
run;

proc freq data=a; tables cgvhd*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd*cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd , kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd*cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of cgvhd , kids only';
where kid = 1;
run;

proc freq data=a; tables cgvhd*abx_class; where kid=1; run;



*********************************************************
*6. Chronic GVHD mortality
*********************************************************;



*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_cgvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;


proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_cgvhd;
race_whitet = race_white * lof_cgvhd;
hememalignancyt = hememalignancy * lof_cgvhd;
transplanttypet = transplanttype * lof_cgvhd;
myeloablativet = myeloablative * lof_cgvhd;
age_dott = age_dot * lof_cgvhd;
matchest = matches * lof_cgvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd mortality';
run;

proc freq data=a; tables died_cgvhd*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of cgvhd mortality';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd mortality, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of cgvhd mortality, adults only';
where kid = 0;
run;

proc freq data=a; tables died_cgvhd*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of cgvhd mortality, kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_cgvhd_mortality*died_cgvhd(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of cgvhd mortality, kids only';
where kid = 1;
run;

proc freq data=a; tables died_cgvhd*abx_class; where kid=1; run;

proc freq data=a; tables abx_class * kid / chisq; run;



*********************************************************
*7. Acute GVHD Skin (any Skin GVHD)
*********************************************************;



*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhdskin*agvhd_skin(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_cgvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;


proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhdskin*agvhd_skin(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_cgvhd;
race_whitet = race_white * lof_cgvhd;
hememalignancyt = hememalignancy * lof_cgvhd;
transplanttypet = transplanttype * lof_cgvhd;
myeloablativet = myeloablative * lof_cgvhd;
age_dott = age_dot * lof_cgvhd;
matchest = matches * lof_cgvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhdskin*agvhd_skin(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of agvhd skin';
run;

proc freq data=a; tables agvhd_skin*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhdskin*agvhd_skin(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of agvhd skin';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhdskin*agvhd_skin(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of avhd skin, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhdskin*agvhd_skin(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of agvhd skin, adults only';
where kid = 0;
run;

proc freq data=a; tables agvhd_skin*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhdskin*agvhd_skin(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of agvhd skin, kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhdskin*agvhd_skin(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of agvhd skin, kids only';
where kid = 1;
run;

proc freq data=a; tables agvhd_skin*abx_class; where kid=1; run;


*********************************************************
*1.) Acute GVHD of gut or liver only (NO  skin)
*********************************************************;

*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
abx_classt = abx_class * lof_agvhd;
proportionality_test :test abx_classt;
run;

* B.) Extended Cox Model to check PH in covariates*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class &catcovariates &contcovariates  femalet race_whitet hememalignancyt transplanttypet myeloablativet 
 age_dott matchest;
hazardratio abx_class /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
femalet = female * lof_agvhd;
race_whitet = race_white * lof_agvhd;
hememalignancyt = hememalignancy * lof_agvhd;
transplanttypet = transplanttype * lof_agvhd;
myeloablativet = myeloablative * lof_agvhd;
age_dott = age_dot * lof_agvhd;
matchest = matches * lof_agvhd;
proportionality_test :test femalet ;
run;

*!!! C.)Adjusted Cox PH model of the main hypothesis !!!*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only';
run;

proc freq data=a; tables aGVHD_livergutonly*abx_class;  run;


*D.) Sensitivity Analysis: Unadjusted Cox Model *;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Unadjusted, Incidence of aGVHD in liver or gut only';
run;


*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 unadjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

proc freq data=a; tables aGVHD_livergutonly*abx_class; where kid=0; run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class &catcovariates &contcovariates;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut , kids only';
where kid=1;
run;
proc phreg data=a ;
class abx_class (ref = '2')   ;
model lof_agvhd_livergutonly*aGVHD_livergutonly(0) = abx_class  ;
hazardratio abx_class /cl=wald diff=ref;
title 'Aim 1 undjusted, Incidence of aGVHD in liver or gut , kids only';
where kid = 1;
run;

proc freq data=a; tables lof_agvhd_livergutonly*abx_class; where kid = 1; run;




***************************************
Part two: anaerobic antibiotic duration

Repeating primary analysis, where the predictor is the category of anaerobic antibiotic duration

*;




*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class durationcat (ref = '2') &catcovariates  ;
model lof_agvhd2to4*aGVHD_2to4(0) = durationcat &catcovariates &contcovariates  durationcatt;
hazardratio durationcat /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
durationcatt = durationcat * lof_agvhd;
proportionality_test :test durationcatt;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class durationcat (ref = '0') &catcovariates  ;
model lof_agvhd2to4*aGVHD_2to4(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only';
run;

*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class durationcat (ref = '0') &catcovariates  ;
model lof_agvhd2to4*aGVHD_2to4(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class durationcat (ref = '2') &catcovariates  ;
model lof_agvhd2to4*aGVHD_2to4(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, kids only';
where kid=1;
run;




***************************************
Part two: anaerobic antibiotic duration

Repeating primary analysis, where the predictor is the category of anaerobic antibiotic duration. 
Repeating for any liver and or gut GVHD

*;




*	A.) Extended Cox Model to check proportional hazards assumption of the anaerobic antibiotic exposure variable.;

proc phreg data=a ;
class durationcat (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = durationcat &catcovariates &contcovariates  durationcatt;
hazardratio durationcat /cl=wald diff=ref;
title 'Acute GVHD Proportionality Test';
durationcatt = durationcat * lof_agvhd;
proportionality_test :test durationcatt;
run;

*!!! C.)Adjusted Cox PH model of the secondary hypothesis !!!*;

proc phreg data=a ;
class durationcat (ref = '0') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only';
run;
proc freq data=a; tables durationcat; run;
*e.) Sensitivity Analysis: Adjusted Cox Model in adults only*;
proc phreg data=a ;
class durationcat (ref = '0') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, adults only';
where kid = 0;
run;

*f.) Sensitivity Analysis: Adjusted Cox Model in children only*;

proc phreg data=a ;
class durationcat (ref = '2') &catcovariates  ;
model lof_agvhd_livergut*aGVHD_livergut(0) = durationcat &catcovariates &contcovariates;
hazardratio durationcat /cl=wald diff=ref;
title 'Aim 1 Adjusted, Incidence of aGVHD in liver or gut only, kids only';
where kid=1;
run;


************************************************
PART THREE: MORTALITY FROM NON-GVHD CAUSES

*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_died_notGVHD*died_notGVHD(0) = abx_class &catcovariates &contcovariates  abx_classt;
hazardratio abx_class /cl=wald diff=ref;
title 'Mortality from causes other than GVHD - Proportionality Test';
abx_classt = abx_class * lof_died_notGVHD;
proportionality_test :test abx_classt;
run; *proportional hazard, p=0.37*;
proc phreg data=a ;
class abx_class (ref = '2') &catcovariates  ;
model lof_died_notGVHD*died_notGVHD(0) = abx_class &catcovariates &contcovariates ;
hazardratio abx_class /cl=wald diff=ref;
title 'Mortality from causes other than GVHD';
run;
proc freq data=a; tables died_notGVHD*abx_class; run;
proc phreg data=a ;
class abx_class (ref = '2')    ;
model lof_died_notGVHD*died_notGVHD(0) = abx_class;
hazardratio abx_class /cl=wald diff=ref;
title 'Mortality from causes other than GVHD';
run;

%let catcovariates = female race_white hememalignancy  myeloablative gvhdprophygroup hsct_donor hsct_source;
%let contcovariates = age_dot matches bmt_year;

proc phreg data=a ;
class abx_class (ref = '2') female race_white hememalignancy  myeloablative gvhdprophygroup hsct_donor hsct_source age_dot matches bmt_year;
model lof_died_notGVHD*died_notGVHD(0) = abx_class age_dot matches bmt_year female race_white hememalignancy  myeloablative gvhdprophygroup hsct_donor hsct_source;
hazardratio abx_class /cl=wald diff=ref;
title 'Mortality from causes other than GVHD';
run;
