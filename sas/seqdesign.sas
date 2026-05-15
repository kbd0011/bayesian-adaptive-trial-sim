/* ============================================================================
   seqdesign.sas
   Purpose  : Cross-validate the R/{rpact} group-sequential design in SAS via
              PROC SEQDESIGN. Should produce the same z-boundaries as
              outputs/tables/design_boundaries.csv from R/01_design_params.R.
   Inputs   : None (parameters set in code, matching config.yml).
   Outputs  : sas/output/seqdesign_boundary.rtf, seqdesign_boundary.csv
   Run via  : SAS OnDemand for Academics.
   ============================================================================ */

%let outdir = sas/output;
options dlcreatedir;
libname out "&outdir";

ods graphics on;
ods rtf file="&outdir/seqdesign_boundary.rtf";

/* ----------------------------------------------------------------------------
   Group-sequential design: O'Brien-Fleming alpha + beta spending, k=2
   stages with interim at 50% information.
   alpha = 0.025 one-sided upper rejection, beta = 0.20 (power = 0.80).
   ---------------------------------------------------------------------------- */
proc seqdesign altref=0.70
               errspend
               plots=(asn boundary errspend);
   OneSidedOBF: design nstages=2
                       method=errfuncobf
                       alt=upper
                       alpha=0.025
                       beta=0.20
                       info=cum(0.5 1.0);
   samplesize model=twosamplesurvival(
       nullhazard = 0.025         /* monthly control hazard = 0.30/12 */
       hazard     = 0.0175        /* nullhazard x HR(0.70) = 0.025*0.70 */
       accrate    = 12            /* enrollments per month */
       acctime    = 10            /* months of accrual */
       totaltime  = 24);          /* total study duration in months */
   ods output Boundary=boundary OutSeqDesign=design_out;
run;

proc print data=boundary; run;

/* Export boundary table to CSV for direct comparison vs R/{rpact} output. */
proc export data=boundary
    outfile="&outdir/seqdesign_boundary.csv"
    dbms=csv replace;
run;

ods rtf close;
