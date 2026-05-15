# Run each test file in a separate Rscript subprocess. This avoids state
# accumulation across test files that has been observed to break rstan's
# stanc parser (causing "parser failed badly" / V8 / C-stack errors) when
# multiple files load the heavy R/00_setup.R + rstan stack in sequence.
#
# Run via: Rscript tests/testthat.R
suppressPackageStartupMessages(library(testthat))

files <- list.files(here::here("tests"),
                    pattern = "^test-.*\\.R$",
                    full.names = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")
# `stop_on_failure = TRUE` causes testthat::test_file to throw if any test
# fails -> Rscript exits non-zero -> system2 returns non-zero status, which
# this runner then aggregates and propagates back to the OS exit code.
results <- vapply(files, function(f) {
  cat("\n==== Running", basename(f), "====\n")
  status <- system2(
    rscript,
    args = c("-e", shQuote(sprintf(
      'suppressPackageStartupMessages(library(testthat)); testthat::test_file("%s", stop_on_failure = TRUE)',
      f
    )))
  )
  status
}, integer(1))

if (any(results != 0)) {
  cat("\nFailures in:", paste(basename(files[results != 0]), collapse = ", "), "\n")
  quit(status = 1)
}
cat("\nAll test files passed.\n")
