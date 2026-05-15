/* ============================================================================
   tcga_lifetest.sas
   Purpose  : Parallel SAS implementation of the TCGA-BRCA Kaplan-Meier
              analysis from R/08_tcga_km.R. PROC LIFETEST produces KM
              curves stratified by hormone-receptor status with a log-rank
              test. Cross-checks numerical agreement with R survival::survfit.
   Inputs   : sas/data/tcga_brca.csv  (written by R/07_tcga_data.R)
   Outputs  : sas/output/tcga_lifetest.rtf
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
ods rtf file="&outdir/tcga_lifetest.rtf";

proc lifetest data=tcga
              plots=(survival(atrisk cb=hw) hazard logsurv)
              method=km;
   time times*status(0);
   strata hr_status / test=logrank;
   title "TCGA-BRCA Overall Survival by Hormone-Receptor Status (KM)";
run;

ods rtf close;
title;
