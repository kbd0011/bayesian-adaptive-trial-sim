.PHONY: all sims report tcga sas-data sensitivity test clean help

R := Rscript

help:
	@echo "make sims        - run all simulations (fixed + adaptive, all scenarios)"
	@echo "make sensitivity - sweep the futility threshold (4 thresholds x 6 scenarios x 300 sims)"
	@echo "make tcga        - run TCGA-BRCA survival analyses (also writes sas/data/tcga_brca.csv)"
	@echo "make sas-data    - just regenerate sas/data/tcga_brca.csv (subset of tcga)"
	@echo "make report      - render Quarto report"
	@echo "make test        - run testthat suite"
	@echo "make all         - sims + tcga + report (excludes sensitivity by default)"
	@echo "make clean       - remove generated outputs"

sims:
	$(R) R/04_run_all_sims.R

sensitivity:
	$(R) R/11_futility_sensitivity.R

tcga: sas-data
	$(R) R/08_tcga_km.R
	$(R) R/09_tcga_cox.R
	$(R) R/10_tcga_bayes_aft.R

sas-data:
	$(R) R/07_tcga_data.R

report:
	quarto render report/index.qmd
	quarto render report/sap_section.qmd

test:
	$(R) tests/testthat.R

all: sims tcga report

clean:
	rm -rf outputs/sims/*.rds outputs/figures/*.pdf outputs/figures/*.png \
	       outputs/tables/*.csv sas/data/*.csv report/_output/ \
	       report/*.html report/*_files/
