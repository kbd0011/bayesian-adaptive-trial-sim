# -----------------------------------------------------------------------------
# Purpose : Pull TCGA-BRCA clinical data, derive a tidy survival dataset for
#           overall-survival analyses stratified by hormone-receptor status.
# Inputs  : RTCGA.clinical::BRCA.clinical (bundled with package)
# Outputs : data/tcga_brca_clean.rds  (tibble: times, status, hr_status,
#                                       er, pr, her2, age, age_decade, barcode)
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))
suppressPackageStartupMessages({
  library(RTCGA)
  library(RTCGA.clinical)
})

tic_step("TCGA-BRCA data prep")

data(BRCA.clinical, package = "RTCGA.clinical")
cli::cli_alert_info("BRCA.clinical raw: {nrow(BRCA.clinical)} subjects x {ncol(BRCA.clinical)} cols")

# --- Extract survival + receptor/age columns ---------------------------------
raw <- suppressWarnings(
  RTCGA::survivalTCGA(
    BRCA.clinical,
    extract.cols = c(
      "patient.breast_carcinoma_estrogen_receptor_status",
      "patient.breast_carcinoma_progesterone_receptor_status",
      "patient.lab_proc_her2_neu_immunohistochemistry_receptor_status",
      "patient.age_at_initial_pathologic_diagnosis"
    )
  )
)

# --- Tidy column names + types -----------------------------------------------
tcga <- tibble::as_tibble(raw) |>
  dplyr::rename(
    barcode = bcr_patient_barcode,
    status  = patient.vital_status,
    er      = patient.breast_carcinoma_estrogen_receptor_status,
    pr      = patient.breast_carcinoma_progesterone_receptor_status,
    her2    = patient.lab_proc_her2_neu_immunohistochemistry_receptor_status,
    age     = patient.age_at_initial_pathologic_diagnosis
  ) |>
  dplyr::mutate(
    status = as.integer(status),
    age    = suppressWarnings(as.numeric(age)),
    er     = tolower(as.character(er)),
    pr     = tolower(as.character(pr)),
    her2   = tolower(as.character(her2))
  )

# --- Derive HR_status --------------------------------------------------------
# HR+ if ER or PR is positive; HR- if both negative; otherwise NA (drop).
tcga <- tcga |>
  dplyr::mutate(
    hr_status = dplyr::case_when(
      er == "positive" | pr == "positive"      ~ "HR+",
      er == "negative" & pr == "negative"      ~ "HR-",
      TRUE                                      ~ NA_character_
    )
  )

# --- Filter to analyzable rows ----------------------------------------------
n0 <- nrow(tcga)
tcga <- tcga |>
  dplyr::filter(
    !is.na(times), times > 0,
    !is.na(status),
    !is.na(hr_status),
    !is.na(age)
  ) |>
  dplyr::mutate(age_decade = floor(age / 10))

dropped <- n0 - nrow(tcga)
cli::cli_alert_info("Dropped {dropped} rows for missing times/status/HR/age")
cli::cli_alert_info("Final analyzable: {nrow(tcga)} subjects, {sum(tcga$status == 1)} events")

# --- Sanity prints -----------------------------------------------------------
cli::cli_h2("HR_status distribution")
print(table(tcga$hr_status, useNA = "always"))
cli::cli_h2("Median follow-up (days), all subjects")
print(stats::median(tcga$times))
cli::cli_h2("Age distribution")
print(summary(tcga$age))

# --- Save --------------------------------------------------------------------
out_path <- file.path(PATH_DATA, "tcga_brca_clean.rds")
saveRDS(tcga, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}")

# Also export a CSV for the SAS PROC LIFETEST/PROC PHREG parallel implementation
sas_csv_dir <- here::here("sas", "data")
fs::dir_create(sas_csv_dir)
sas_csv_path <- file.path(sas_csv_dir, "tcga_brca.csv")
utils::write.csv(tcga, sas_csv_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {sas_csv_path}}")

toc_step()
