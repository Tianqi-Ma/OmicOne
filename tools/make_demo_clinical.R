# Write a small synthetic cohort so the Clinical & survival step can be tested
# without real patient data.
#
#   Rscript tools/make_demo_clinical.R            # -> demo_clinical.csv
#   Rscript tools/make_demo_clinical.R out.csv
#
# The cohort is built with a real effect in it, so you can tell a working
# analysis from a broken one:
#   - stage III dies faster than stage I  -> log-rank p should be small
#   - risk_score is continuous            -> exercises the median / tertile /
#                                            optimal cutpoint controls
#   - age is noise                        -> its Cox HR should sit near 1
# About 40% of patients are censored (42% with seed 2024: 69 events in 120), so
# the Kaplan-Meier curves get censoring ticks and the median of at least one
# arm is reached.
#
# This cohort uses its own P001... ids, so it cannot be joined to a MAF. For the
# WES pipeline, the clinical table that ships with maftools is the demo:
#   system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
# It is loaded automatically with the TCGA LAML demo MAF (Import MAF -> Demo
# data): 200 patients keyed by Tumor_Sample_Barcode, with FAB_classification,
# days_to_last_followup and Overall_Survival_Status. Seven of the 200 have no
# MAF record and 12 have no finite follow-up, which exercises the survival
# step's "no MAF record = sequenced WT" option (188 patients analysed with it,
# 182 without).

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "demo_clinical.csv"

set.seed(2024)
n <- 120

stage <- sample(c("I", "II", "III"), n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
# mean survival in months by stage, spread far enough apart that three clearly
# separated Kaplan-Meier curves are the visible "it worked" signal
scale <- c(I = 70, II = 30, III = 10)[stage]
risk  <- round(stats::runif(n, 0, 1), 3)
# risk_score is prognostic too, on a multiplicative scale so it stays
# independent of stage: a high scorer survives about a third as long as a low
# scorer, which n = 120 can detect without washing stage out
true_time <- stats::rexp(n, rate = 1 / (scale * exp(-1.4 * (risk - 0.5))))

# administrative censoring at 60 months, plus random loss to follow-up
cens_at  <- pmin(60, stats::rexp(n, rate = 1 / 70))
observed <- pmin(true_time, cens_at)
event    <- as.integer(true_time <= cens_at)

clin <- data.frame(
  sample     = sprintf("P%03d", seq_len(n)),
  os_months  = round(observed, 1),
  os_status  = ifelse(event == 1, "Dead", "Alive"),
  stage      = stage,
  age        = round(stats::rnorm(n, 62, 11)),
  sex        = sample(c("F", "M"), n, replace = TRUE),
  risk_score = risk,
  stringsAsFactors = FALSE
)

utils::write.csv(clin, out, row.names = FALSE)

cat(sprintf("Wrote %s: %d patients, %d events (%.0f%% censored).\n",
            out, nrow(clin), sum(event), 100 * mean(event == 0)))
cat("\nIn the app -> Clinical & survival, map the columns like this:\n")
cat("  Sample id        sample\n")
cat("  Follow-up time   os_months     (unit: Months)\n")
cat("  Outcome          os_status\n")
cat("  Cox covariates   stage, age, sex, risk_score\n")

cat("\nWith seed 2024 and n = 120 you should get, near enough:\n")
cat("  Stratify by stage         log-rank p ~ 5e-05\n")
cat("                            median OS  ~ 30 / 21 / 9 months for I / II / III\n")
cat("  Stratify by risk_score    median   cutpoint  p ~ 0.03\n")
cat("                            tertile  cutpoint  p ~ 0.01\n")
cat("                            optimal  cutpoint  p ~ 0.002\n")
cat("  Cox: stage III HR ~ 3.5 (p ~ 6e-05), risk_score HR ~ 2.9 (p ~ 0.009),\n")
cat("       age and sex both non-significant.\n")

cat("\nNote what the three cutpoints do to the same data: p drops from 0.03 to\n")
cat("0.002 purely by choosing the split that maximises separation. That is the\n")
cat("'optimal is exploratory' warning in the app, made concrete.\n")
