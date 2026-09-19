#!/usr/bin/env Rscript
# scripts/verify_stored_results.R
# Re-run the simulation pipeline from scratch and check that it reproduces the
# stored tables in results/tables/ (T1, T2, T6, T7).
#
# Usage (from the repository root): Rscript scripts/verify_stored_results.R
#
# The stored CSVs are restored byte-for-byte afterwards, so the working tree is
# left unchanged. Exit status is 0 if every table matches within `tol`.

tol    <- 1e-10
tables <- c("T1_theoretical", "T2_monte_carlo",
            "T6_misspecification", "T7_bimodal_test")
dir    <- file.path("results", "tables")
paths  <- file.path(dir, paste0(tables, ".csv"))

stored_bytes <- lapply(paths, function(p) readBin(p, "raw", file.size(p)))
stored       <- lapply(paths, read.csv)

# Fresh run: no Monte Carlo cache.
unlink(file.path("results", "mc_cache"), recursive = TRUE)
PROJECT_ROOT <- getwd()
source(file.path("scripts", "run_all.R"))

ok <- TRUE
cat("\n=== Verification against stored results (tol =", tol, ") ===\n")
for (i in seq_along(tables)) {
  new <- read.csv(paths[i])
  old <- stored[[i]]
  same_shape <- identical(dim(new), dim(old)) && identical(names(new), names(old))
  num  <- vapply(old, is.numeric, logical(1))
  diff <- if (same_shape && any(num)) {
    max(abs(as.matrix(new[num]) - as.matrix(old[num])), na.rm = TRUE)
  } else NA_real_
  chr_same <- same_shape && identical(new[!num], old[!num])
  pass <- same_shape && chr_same && !is.na(diff) && diff <= tol
  ok <- ok && pass
  cat(sprintf("  %-22s %s  max|diff| = %s\n", tables[i],
              if (pass) "OK  " else "FAIL", format(diff, digits = 3)))
  writeBin(stored_bytes[[i]], paths[i])   # restore stored file
}

unlink(file.path("results", "mc_cache"), recursive = TRUE)
unlink(file.path("results", "figures"), recursive = TRUE)
unlink(file.path("results", "tables", "biaxial_diagnostics.rds"))

cat(if (ok) "\nAll stored tables reproduced.\n" else "\nMISMATCH: see FAIL rows above.\n")
quit(status = if (ok) 0 else 1)
