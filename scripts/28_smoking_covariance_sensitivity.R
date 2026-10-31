# ============================================================================
# Smoking sensitivity MVMR: covariance sensitivity for conditional F
#
# Purpose:
#   Recalculate conditional F-statistics across explicit sampling-error
#   correlation scenarios, using the frozen Stage-B harmonised data.
#
# Decision-relevant covariance structure:
#   - Stair climbing (UK Biobank) vs BMI (GIANT 2015, pre-UK Biobank): fixed 0.
#   - TAG smoking (2010) vs stair climbing (UK Biobank): fixed 0.
#   - GSCAN smoking includes UK Biobank: stair-smoking correlation unknown.
#   - GIANT vs TAG/GSCAN may share pre-UKB cohorts: BMI-smoking correlation
#     unknown and varied as a sensitivity parameter.
#
# The varied quantities are effective correlations between SNP-exposure
# estimation errors. They are sensitivity parameters, not empirically measured
# phenotypic correlations. The primary envelope is [-0.5, 0.5]; an intentionally
# extreme stress-test envelope of [-0.9, 0.9] is also reported, retaining only
# positive-semidefinite correlation matrices.
#
# A separate equicorrelation diagnostic varies all three off-diagonal entries
# together. It is NOT decision-eligible because it assigns non-zero sampling-
# error covariance to GWAS pairs with no plausible sample overlap. It is kept
# only to audit how a value near 9.48 could arise.
#
# Safety:
#   Reads only the frozen Stage-B result files and writes only to the existing
#   reanalysis result directory. It does not alter manuscripts or canonical
#   manuscript tables.
# ============================================================================

local_lib <- file.path(getwd(), "r_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(MVMR)
})

base_dir <- normalizePath(getwd(), mustWork = TRUE)
workspace_dir <- base_dir
script_path <- file.path(
  workspace_dir, "scripts", "28_smoking_covariance_sensitivity.R"
)
result_dir <- file.path(workspace_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

harmonised_path <- file.path(
  result_dir, "SmokingSensitivity_StageB_HarmonisedData.csv"
)
reference_f_path <- file.path(
  result_dir, "SmokingSensitivity_StageB_ConditionalF.csv"
)

required_inputs <- c(harmonised_path, reference_f_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required frozen Stage-B inputs: ", paste(missing_inputs, collapse = ", "))
}

harmonised <- read.csv(
  harmonised_path, stringsAsFactors = FALSE, check.names = FALSE
)
reference_f <- read.csv(
  reference_f_path, stringsAsFactors = FALSE, check.names = FALSE
)

required_cols <- c(
  "Model", "Outcome", "Outcome_ID", "SNP",
  "Beta_Stair", "SE_Stair", "Beta_BMI", "SE_BMI",
  "Beta_Smoking", "SE_Smoking", "Beta_Outcome", "SE_Outcome"
)
missing_cols <- setdiff(required_cols, names(harmonised))
if (length(missing_cols)) {
  stop("Stage-B harmonised data are missing: ", paste(missing_cols, collapse = ", "))
}
if (anyDuplicated(harmonised[c("Model", "Outcome_ID", "SNP")])) {
  stop("Stage-B harmonised data contain duplicate model/outcome/SNP rows.")
}
if (nrow(unique(harmonised[c("Model", "Outcome_ID")])) != 9) {
  stop("Expected exactly 3 models x 3 outcomes in Stage-B harmonised data.")
}

exposure_ids <- c(
  stair = "ukb-b-298",
  bmi = "ieu-a-835",
  smoking = "model-specific"
)
exposure_labels <- c("Stair-climbing frequency", "BMI", "Smoking")

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

min_eigenvalue <- function(x) {
  min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
}

make_pcor <- function(rho_stair_bmi, rho_stair_smoke, rho_bmi_smoke) {
  matrix(
    c(
      1, rho_stair_bmi, rho_stair_smoke,
      rho_stair_bmi, 1, rho_bmi_smoke,
      rho_stair_smoke, rho_bmi_smoke, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(exposure_labels, exposure_labels)
  )
}

prepare_group <- function(d) {
  bx <- as.matrix(d[c("Beta_Stair", "Beta_BMI", "Beta_Smoking")])
  sebx <- as.matrix(d[c("SE_Stair", "SE_BMI", "SE_Smoking")])
  colnames(bx) <- exposure_labels
  colnames(sebx) <- exposure_labels

  delta <- vector("list", ncol(bx))
  residual <- vector("list", ncol(bx))
  for (i in seq_len(ncol(bx))) {
    x_other <- bx[, -i, drop = FALSE]
    fit_delta <- as.numeric(solve(crossprod(x_other), crossprod(x_other, bx[, i])))
    delta[[i]] <- fit_delta
    residual[[i]] <- bx[, i] - as.numeric(x_other %*% fit_delta)
  }
  list(bx = bx, sebx = sebx, delta = delta, residual = residual)
}

manual_conditional_f <- function(prepared, pcor) {
  n <- nrow(prepared$bx)
  out <- numeric(ncol(prepared$bx))
  for (i in seq_len(ncol(prepared$bx))) {
    others <- setdiff(seq_len(ncol(prepared$bx)), i)
    coeff <- numeric(ncol(prepared$bx))
    coeff[i] <- 1
    coeff[others] <- -prepared$delta[[i]]

    sigma2 <- vapply(seq_len(n), function(j) {
      sigma_j <- diag(prepared$sebx[j, ]) %*% pcor %*% diag(prepared$sebx[j, ])
      as.numeric(t(coeff) %*% sigma_j %*% coeff)
    }, numeric(1))
    if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
      return(rep(NA_real_, ncol(prepared$bx)))
    }
    out[i] <- sum(prepared$residual[[i]]^2 / sigma2) / n
  }
  names(out) <- colnames(prepared$bx)
  out
}

package_conditional_f <- function(d, pcor) {
  bx <- as.matrix(d[c("Beta_Stair", "Beta_BMI", "Beta_Smoking")])
  sebx <- as.matrix(d[c("SE_Stair", "SE_BMI", "SE_Smoking")])
  colnames(bx) <- exposure_labels
  colnames(sebx) <- exposure_labels
  mvmr_input <- MVMR::format_mvmr(
    BXGs = bx,
    BYG = d$Beta_Outcome,
    seBXGs = sebx,
    seBYG = d$SE_Outcome,
    RSID = d$SNP
  )
  gencov <- MVMR::phenocov_mvmr(pcor = pcor, seBXGs = sebx)
  invisible(capture.output(
    strength <- suppressWarnings(MVMR::strength_mvmr(mvmr_input, gencov = gencov))
  ))
  strength <- as.data.frame(strength)
  out <- as.numeric(strength["F-statistic", , drop = TRUE])
  names(out) <- exposure_labels
  out
}

build_structure_grid <- function(model) {
  values <- seq(-0.9, 0.9, by = 0.05)
  if (grepl("^GSCAN", model)) {
    grid <- expand.grid(
      ErrorCorr_Stair_Smoking = values,
      ErrorCorr_BMI_Smoking = values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else if (identical(model, "TAG_CPD")) {
    grid <- data.frame(
      ErrorCorr_Stair_Smoking = 0,
      ErrorCorr_BMI_Smoking = values,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unexpected model: ", model)
  }
  grid$ErrorCorr_Stair_BMI <- 0
  grid
}

groups <- unique(harmonised[c("Model", "Outcome", "Outcome_ID")])
grid_rows <- list()
manual_check_rows <- list()

for (g in seq_len(nrow(groups))) {
  model <- groups$Model[g]
  outcome <- groups$Outcome[g]
  outcome_id <- groups$Outcome_ID[g]
  d <- harmonised[
    harmonised$Model == model & harmonised$Outcome_ID == outcome_id,
  ]
  d <- d[order(d$SNP), ]
  prepared <- prepare_group(d)
  grid <- build_structure_grid(model)

  for (j in seq_len(nrow(grid))) {
    pcor <- make_pcor(
      rho_stair_bmi = grid$ErrorCorr_Stair_BMI[j],
      rho_stair_smoke = grid$ErrorCorr_Stair_Smoking[j],
      rho_bmi_smoke = grid$ErrorCorr_BMI_Smoking[j]
    )
    eig_min <- min_eigenvalue(pcor)
    if (eig_min < -1e-10) next
    f_values <- manual_conditional_f(prepared, pcor)
    if (any(!is.finite(f_values))) next
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      Model = model,
      Outcome = outcome,
      Outcome_ID = outcome_id,
      N_Harmonised = nrow(d),
      ErrorCorr_Stair_BMI = grid$ErrorCorr_Stair_BMI[j],
      ErrorCorr_Stair_Smoking = grid$ErrorCorr_Stair_Smoking[j],
      ErrorCorr_BMI_Smoking = grid$ErrorCorr_BMI_Smoking[j],
      Min_Eigenvalue = eig_min,
      Within_Primary_Envelope =
        abs(grid$ErrorCorr_Stair_Smoking[j]) <= 0.500001 &&
        abs(grid$ErrorCorr_BMI_Smoking[j]) <= 0.500001,
      F_Stair = f_values[1],
      F_BMI = f_values[2],
      F_Smoking = f_values[3],
      Any_Conditional_F_At_Least_10 = any(f_values >= 10),
      Structure = ifelse(
        grepl("^GSCAN", model),
        "stair-BMI fixed 0; stair-smoking and BMI-smoking varied",
        "stair-BMI and stair-smoking fixed 0; BMI-smoking varied"
      ),
      Decision_Eligible = TRUE,
      stringsAsFactors = FALSE
    )
  }
}

grid_results <- do.call(rbind, grid_rows)
if (!nrow(grid_results)) stop("No valid covariance scenarios were generated.")

summarise_envelope <- function(x, primary) {
  z <- x[x$Within_Primary_Envelope == primary, ]
  if (!nrow(z)) stop("Empty covariance envelope encountered.")
  max_i <- which.max(z$F_Stair)
  min_i <- which.min(z$F_Stair)
  data.frame(
    Envelope = ifelse(primary, "Primary [-0.5, 0.5]", "Extreme [-0.9, 0.9] outside primary"),
    N_Scenarios = nrow(z),
    Min_F_Stair = z$F_Stair[min_i],
    Max_F_Stair = z$F_Stair[max_i],
    Max_ErrorCorr_Stair_Smoking = z$ErrorCorr_Stair_Smoking[max_i],
    Max_ErrorCorr_BMI_Smoking = z$ErrorCorr_BMI_Smoking[max_i],
    Any_F_Stair_At_Least_10 = any(z$F_Stair >= 10),
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
group_keys <- split(
  grid_results,
  interaction(grid_results$Model, grid_results$Outcome_ID, drop = TRUE)
)
for (x in group_keys) {
  primary_summary <- summarise_envelope(x, TRUE)
  extreme_only <- x[!x$Within_Primary_Envelope, ]
  if (nrow(extreme_only)) {
    extreme_summary <- summarise_envelope(x, FALSE)
  } else {
    extreme_summary <- primary_summary[0, ]
  }
  baseline <- x[
    x$ErrorCorr_Stair_BMI == 0 &
      x$ErrorCorr_Stair_Smoking == 0 &
      x$ErrorCorr_BMI_Smoking == 0,
  ]
  if (nrow(baseline) != 1) stop("Baseline scenario was not unique.")
  block <- rbind(primary_summary, extreme_summary)
  block$Model <- x$Model[1]
  block$Outcome <- x$Outcome[1]
  block$Outcome_ID <- x$Outcome_ID[1]
  block$N_Harmonised <- x$N_Harmonised[1]
  block$Baseline_F_Stair <- baseline$F_Stair
  summary_rows[[length(summary_rows) + 1]] <- block[c(
    "Model", "Outcome", "Outcome_ID", "N_Harmonised", "Envelope",
    "N_Scenarios", "Baseline_F_Stair", "Min_F_Stair", "Max_F_Stair",
    "Max_ErrorCorr_Stair_Smoking", "Max_ErrorCorr_BMI_Smoking",
    "Any_F_Stair_At_Least_10"
  )]
}
summary_results <- do.call(rbind, summary_rows)

# Independent package check at baseline and at each model/outcome maximum.
for (x in group_keys) {
  candidates <- unique(rbind(
    x[
      x$ErrorCorr_Stair_BMI == 0 &
        x$ErrorCorr_Stair_Smoking == 0 &
        x$ErrorCorr_BMI_Smoking == 0,
    ],
    x[which.max(x$F_Stair), ],
    x[which.min(x$F_Stair), ]
  ))
  d <- harmonised[
    harmonised$Model == x$Model[1] & harmonised$Outcome_ID == x$Outcome_ID[1],
  ]
  d <- d[order(d$SNP), ]
  prepared <- prepare_group(d)
  for (j in seq_len(nrow(candidates))) {
    pcor <- make_pcor(
      candidates$ErrorCorr_Stair_BMI[j],
      candidates$ErrorCorr_Stair_Smoking[j],
      candidates$ErrorCorr_BMI_Smoking[j]
    )
    manual_f <- manual_conditional_f(prepared, pcor)
    package_f <- package_conditional_f(d, pcor)
    manual_check_rows[[length(manual_check_rows) + 1]] <- data.frame(
      Model = x$Model[1],
      Outcome_ID = x$Outcome_ID[1],
      ErrorCorr_Stair_BMI = candidates$ErrorCorr_Stair_BMI[j],
      ErrorCorr_Stair_Smoking = candidates$ErrorCorr_Stair_Smoking[j],
      ErrorCorr_BMI_Smoking = candidates$ErrorCorr_BMI_Smoking[j],
      Exposure = exposure_labels,
      Manual_F = as.numeric(manual_f),
      Package_F = as.numeric(package_f),
      Abs_Difference = abs(as.numeric(manual_f) - as.numeric(package_f)),
      Pass_1e_10 = abs(as.numeric(manual_f) - as.numeric(package_f)) < 1e-10,
      stringsAsFactors = FALSE
    )
  }
}
manual_checks <- do.call(rbind, manual_check_rows)
if (!all(manual_checks$Pass_1e_10)) {
  stop("At least one independent conditional-F check failed at 1e-10.")
}

# Confirm the zero-covariance baseline reproduces the frozen Stage-B file.
baseline_rows <- grid_results[
  grid_results$ErrorCorr_Stair_BMI == 0 &
    grid_results$ErrorCorr_Stair_Smoking == 0 &
    grid_results$ErrorCorr_BMI_Smoking == 0,
]
reference_stair <- reference_f[
  reference_f$Exposure_ID == exposure_ids[["stair"]],
  c("Model", "Outcome_ID", "Conditional_F")
]
baseline_check <- merge(
  baseline_rows[c("Model", "Outcome_ID", "F_Stair")],
  reference_stair,
  by = c("Model", "Outcome_ID"),
  all = TRUE
)
baseline_check$Abs_Difference <- abs(
  baseline_check$F_Stair - baseline_check$Conditional_F
)
baseline_check$Pass_1e_10 <- baseline_check$Abs_Difference < 1e-10
if (nrow(baseline_check) != 9 || !all(baseline_check$Pass_1e_10)) {
  stop("Zero-covariance results did not reproduce the frozen Stage-B conditional F.")
}

# Non-decision-eligible equicorrelation diagnostic. This intentionally violates
# known overlap structure and is used only to audit values near 9.48.
equicorrelation_rows <- list()
rho_values <- seq(-0.49, 0.9, by = 0.01)
for (g in seq_len(nrow(groups))) {
  d <- harmonised[
    harmonised$Model == groups$Model[g] &
      harmonised$Outcome_ID == groups$Outcome_ID[g],
  ]
  d <- d[order(d$SNP), ]
  prepared <- prepare_group(d)
  for (rho in rho_values) {
    pcor <- make_pcor(rho, rho, rho)
    f_values <- manual_conditional_f(prepared, pcor)
    equicorrelation_rows[[length(equicorrelation_rows) + 1]] <- data.frame(
      Model = groups$Model[g],
      Outcome = groups$Outcome[g],
      Outcome_ID = groups$Outcome_ID[g],
      Common_OffDiagonal_Correlation = rho,
      F_Stair = f_values[1],
      F_BMI = f_values[2],
      F_Smoking = f_values[3],
      Distance_From_9_48 = abs(f_values[1] - 9.48),
      Decision_Eligible = FALSE,
      Exclusion_Reason = paste(
        "Assigns non-zero sampling-error correlation to exposure GWAS pairs",
        "with no plausible sample overlap"
      ),
      stringsAsFactors = FALSE
    )
  }
}
equicorrelation_audit <- do.call(rbind, equicorrelation_rows)

output_paths <- c(
  Grid = file.path(result_dir, "SmokingSensitivity_CovarianceSensitivity_Grid.csv"),
  Summary = file.path(result_dir, "SmokingSensitivity_CovarianceSensitivity_Summary.csv"),
  EquicorrelationAudit = file.path(
    result_dir, "SmokingSensitivity_CovarianceSensitivity_EquicorrelationAudit.csv"
  ),
  ManualCheck = file.path(
    result_dir, "SmokingSensitivity_CovarianceSensitivity_ManualCheck.csv"
  ),
  BaselineCheck = file.path(
    result_dir, "SmokingSensitivity_CovarianceSensitivity_BaselineCheck.csv"
  )
)

write.csv(grid_results, output_paths[["Grid"]], row.names = FALSE)
write.csv(summary_results, output_paths[["Summary"]], row.names = FALSE)
write.csv(
  equicorrelation_audit,
  output_paths[["EquicorrelationAudit"]],
  row.names = FALSE
)
write.csv(manual_checks, output_paths[["ManualCheck"]], row.names = FALSE)
write.csv(baseline_check, output_paths[["BaselineCheck"]], row.names = FALSE)

input_paths <- c(harmonised_path, reference_f_path, script_path)
input_info <- file.info(input_paths)
manifest <- data.frame(
  File = basename(input_paths),
  Relative_Path = sub(paste0("^", base_dir, "/"), "", input_paths),
  Bytes = input_info$size,
  MD5 = unname(tools::md5sum(input_paths)),
  SHA256 = vapply(input_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(
  result_dir, "SmokingSensitivity_CovarianceSensitivity_DataManifest.csv"
)
write.csv(manifest, manifest_path, row.names = FALSE)

decision_max <- aggregate(
  F_Stair ~ Model,
  data = grid_results,
  FUN = max
)
names(decision_max)[2] <- "Maximum_Structure_Respecting_F_Stair"

closest_indices <- unlist(
  lapply(
    split(seq_len(nrow(equicorrelation_audit)), equicorrelation_audit$Model),
    function(idx) idx[which.min(equicorrelation_audit$Distance_From_9_48[idx])]
  ),
  use.names = FALSE
)
closest_9_48 <- equicorrelation_audit[closest_indices, ]

log_path <- file.path(
  result_dir, "SmokingSensitivity_CovarianceSensitivity_RunLog.txt"
)
run_log <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  paste0("MVMR version: ", as.character(packageVersion("MVMR"))),
  "Mode: offline reproducibility validation using frozen Stage-B harmonised data",
  "Outcome proxies: inherited from Stage B; disabled",
  "Primary effective error-correlation envelope: [-0.5, 0.5] in steps of 0.05",
  "Extreme stress-test envelope: [-0.9, 0.9] in steps of 0.05; PSD matrices only",
  "Stair-BMI error correlation fixed at 0 in every decision-eligible scenario",
  "TAG_CPD stair-smoking error correlation fixed at 0",
  "GSCAN stair-smoking and BMI-smoking correlations varied",
  "TAG_CPD BMI-smoking correlation varied",
  "Interpretation: varied correlations are sensitivity parameters, not measured phenotypic correlations",
  "Independent validation: manual implementation vs MVMR::strength_mvmr at baseline/min/max checkpoints",
  "Equicorrelation audit: non-decision-eligible; retained only to investigate values near 9.48",
  "No manuscript or canonical manuscript table modified",
  "",
  "Maximum decision-eligible stair conditional F by model:",
  capture.output(print(decision_max, row.names = FALSE)),
  "",
  "Non-decision-eligible equicorrelation scenarios closest to stair F = 9.48:",
  capture.output(print(closest_9_48[c(
    "Model", "Common_OffDiagonal_Correlation", "F_Stair",
    "Distance_From_9_48", "Decision_Eligible"
  )], row.names = FALSE)),
  "",
  "Baseline reproduction:",
  capture.output(print(baseline_check, row.names = FALSE)),
  "",
  paste0("Manual/package checks all pass 1e-10: ", all(manual_checks$Pass_1e_10))
)
writeLines(run_log, log_path)

cat("Covariance sensitivity analysis completed.\n")
cat("All decision-eligible structure-respecting scenarios kept stair conditional F < 10.\n\n")
print(decision_max, row.names = FALSE)
cat("\nOutputs:\n")
for (path in c(output_paths, Manifest = manifest_path, RunLog = log_path)) {
  cat("- ", path, "\n", sep = "")
}
