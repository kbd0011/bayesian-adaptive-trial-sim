/* ============================================================================
   tcga_phreg.sas
   Purpose  : Parallel SAS implementation of the TCGA-BRCA Cox PH regression
              from R/09_tcga_cox.R. PROC PHREG fits the same model
              (HR_status + age_decade, ties=efron) and produces a Schoenfeld
              residual diagnostic via the ASSESS PH statement.
   Inputs   : sas/data/tcga_brca.csv  (written by R/07_tcga_data.R)
   Outputs  : sas/output/tcga_phreg.rtf
   ============================================================================ */

%let csv     = sas/data/tcga_brca.csv;
%let outdir  = sas/output;
options dlcreatedir;
libname out "&outdir";

proc import datafile="&csv"
    out=tcga
    dbms=csv
    replace;
    getnames=yes;
    guessingrows=max;
run;

ods graphics on;
ods rtf file="&outdir/tcga_phreg.rtf";

proc phreg data=tcga plots=survival;
   class hr_status (ref="HR-") / param=ref;
   model times*status(0) = hr_status age_decade / ties=efron risklimits;
   /* Resampling-based test of the proportional-hazards assumption per
      Lin, Wei, and Ying (1993). Should flag hr_status PH violation,
      matching R cox.zph result (p~=0.013). */
   assess ph / resample seed=20260513;
   title "TCGA-BRCA Cox PH: hr_status + age_decade";
run;

ods rtf close;
title;
