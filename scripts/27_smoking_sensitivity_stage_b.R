# ============================================================================
# Smoking sensitivity MVMR: Stage B fixed-model outcome analysis
#
# Purpose:
#   Run three prespecified, Stage-A-feasible exposure models against three
#   prespecified FinnGen COPD outcomes. Outcome proxies are disabled.
#
# Models:
#   GSCAN_CPD  = stair climbing + BMI + GSCAN cigarettes per day
#   TAG_CPD    = stair climbing + BMI + TAG cigarettes per day
#   GSCAN_INIT = stair climbing + BMI + GSCAN smoking initiation
#
# Outputs:
#   Standard no-intercept multivariable IVW estimates, an independent matrix
#   WLS check, covariance-limited conditional F-statistics, two explicitly
#   distinguished heterogeneity statistics, harmonised data, SNP audit, input
#   manifest, and run log.
#
# Safety:
#   Reads the frozen Stage-A exposure sets and writes only inside the existing
#   reanalysis workspace. It does not alter manuscripts, canonical tables, or
#   historical MVMR results.
# ============================================================================

local_lib <- file.path(getwd(), "r_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MVMR)
})

base_dir <- normalizePath(getwd(), mustWork = TRUE)
workspace_dir <- base_dir
input_dir <- file.path(workspace_dir, "inputs")
result_dir <- file.path(workspace_dir, "results")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

refresh <- "--refresh" %in% commandArgs(trailingOnly = TRUE)

models <- data.frame(
  Model = c("GSCAN_CPD", "TAG_CPD", "GSCAN_INIT"),
  Smoking_ID = c("ieu-b-25", "ieu-a-961", "ieu-b-4877"),
  Smoking_Dimension = c("Cigarettes per day", "Cigarettes per day", "Smoking initiation"),
  Overlap_Context = c(
    "Covariance-limited: GSCAN includes UK Biobank",
    "Overlap-minimised: TAG predates UK Biobank",
    "Covariance-limited: GSCAN includes UK Biobank"
  ),
  stringsAsFactors = FALSE
)

outcomes <- data.frame(
  Outcome = c("COPD hospitalization", "General COPD", "Early-onset COPD"),
  Outcome_ID = c("finn-b-COPD_HOSPITAL", "finn-b-J10_COPD", "finn-b-COPD_EARLY"),
  stringsAsFactors = FALSE
)

fixed_exposures <- c(stair = "ukb-b-298", bmi = "ieu-a-835")
exposure_labels <- c(
  "ukb-b-298" = "Stair-climbing frequency",
  "ieu-a-835" = "BMI",
  "ieu-b-25" = "Cigarettes per day (GSCAN)",
  "ieu-a-961" = "Cigarettes per day (TAG)",
  "ieu-b-4877" = "Smoking initiation (GSCAN)"
)

load_opengwas_token <- function() {
  token_files <- c(
    file.path(base_dir, "data", "opengwas_token"),
    path.expand("~/.config/opengwas_token")
  )
  token_file <- token_files[file.exists(token_files)][1]
  if (is.na(token_file)) {
    stop(
      "OpenGWAS token file not found. Expected data/opengwas_token or ",
      "~/.config/opengwas_token."
    )
  }
  token <- trimws(readLines(token_file, warn = FALSE)[1])
  if (!nzchar(token)) stop("OpenGWAS token file is empty: ", token_file)
  Sys.setenv(OPENGWAS_JWT = token)
  invisible(token_file)
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

validate_exposure_data <- function(dat, expected_ids, model) {
  required <- c(
    "SNP", "id.exposure", "effect_allele.exposure", "other_allele.exposure",
    "beta.exposure", "se.exposure", "pval.exposure"
  )
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols)) {
    stop(model, " exposure cache is missing: ", paste(missing_cols, collapse = ", "))
  }
  if (!identical(sort(unique(dat$id.exposure)), sort(expected_ids))) {
    stop(model, " exposure IDs do not match the prespecified model.")
  }
  n_snp <- length(unique(dat$SNP))
  per_id <- table(dat$id.exposure)
  if (any(per_id != n_snp) || anyDuplicated(dat[c("SNP", "id.exposure")])) {
    stop(model, " exposure cache is incomplete or duplicated.")
  }
  if (any(!is.finite(dat$beta.exposure)) || any(!is.finite(dat$se.exposure))) {
    stop(model, " exposure cache contains non-finite beta or SE values.")
  }
  invisible(TRUE)
}

validate_outcome_data <- function(dat, requested_snps, outcome_id, path) {
  required <- c(
    "SNP", "id.outcome", "effect_allele.outcome", "other_allele.outcome",
    "beta.outcome", "se.outcome", "pval.outcome"
  )
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols)) {
    stop(basename(path), " is missing: ", paste(missing_cols, collapse = ", "))
  }
  if (!all(unique(dat$id.outcome) == outcome_id)) {
    stop(basename(path), " contains an unexpected outcome ID.")
  }
  if (anyDuplicated(dat$SNP)) stop(basename(path), " contains duplicated SNP rows.")
  if (!all(dat$SNP %in% requested_snps)) {
    stop(basename(path), " contains SNPs outside the requested instrument set.")
  }
  if ("proxy.outcome" %in% names(dat)) {
    proxy_flag <- !is.na(dat$proxy.outcome) & as.character(dat$proxy.outcome) %in% c("TRUE", "T", "1")
    if (any(proxy_flag)) {
      stop(basename(path), " contains outcome proxies despite the no-proxy rule.")
    }
  }
  if (any(!is.finite(dat$beta.outcome)) || any(!is.finite(dat$se.outcome))) {
    stop(basename(path), " contains non-finite beta or SE values.")
  }
  invisible(TRUE)
}

manual_wls <- function(x, y, se_y) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  w <- 1 / as.numeric(se_y)^2
  xtwx <- crossprod(x, x * w)
  xtwy <- crossprod(x, y * w)
  beta <- as.numeric(solve(xtwx, xtwy))
  residual <- y - as.numeric(x %*% beta)
  df <- nrow(x) - ncol(x)
  sigma2 <- sum(w * residual^2) / df
  vcov <- sigma2 * solve(xtwx)
  se <- sqrt(diag(vcov))
  p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  names(beta) <- colnames(x)
  names(se) <- colnames(x)
  names(p) <- colnames(x)
  list(
    beta = beta,
    se = se,
    p = p,
    residual = residual,
    standard_q = sum(w * residual^2),
    standard_q_df = df,
    standard_q_p = pchisq(sum(w * residual^2), df = df, lower.tail = FALSE)
  )
}

extract_f_values <- function(strength_obj) {
  dat <- as.data.frame(strength_obj)
  if ("F-statistic" %in% rownames(dat)) {
    return(as.numeric(dat["F-statistic", , drop = TRUE]))
  }
  stop("Could not identify conditional F-statistics from MVMR::strength_mvmr().")
}

safe_outcome_name <- function(outcome_id) sub("^finn-b-", "", outcome_id)

estimate_rows <- list()
manual_rows <- list()
strength_rows <- list()
heterogeneity_rows <- list()
audit_rows <- list()
harmonised_rows <- list()
input_paths <- character(0)
source_log <- character(0)
token_loaded <- FALSE

for (m in seq_len(nrow(models))) {
  model <- models$Model[m]
  smoking_id <- models$Smoking_ID[m]
  exposure_ids <- c(unname(fixed_exposures), smoking_id)
  exposure_path <- file.path(
    input_dir, paste0("SmokingSensitivity_ExposureData_", model, ".csv")
  )
  if (!file.exists(exposure_path)) stop("Missing Stage-A exposure cache: ", exposure_path)
  mv_exp <- read.csv(exposure_path, stringsAsFactors = FALSE, check.names = FALSE)
  validate_exposure_data(mv_exp, exposure_ids, model)
  requested_snps <- sort(unique(mv_exp$SNP))
  input_paths <- c(input_paths, exposure_path)

  for (o in seq_len(nrow(outcomes))) {
    outcome_name <- outcomes$Outcome[o]
    outcome_id <- outcomes$Outcome_ID[o]
    outcome_short <- safe_outcome_name(outcome_id)
    outcome_path <- file.path(
      input_dir,
      paste0("SmokingSensitivity_OutcomeData_", model, "_", outcome_short, ".csv")
    )

    if (!refresh && file.exists(outcome_path)) {
      out_dat <- read.csv(outcome_path, stringsAsFactors = FALSE, check.names = FALSE)
      data_source <- "frozen direct-SNP cache"
    } else {
      if (!token_loaded) {
        load_opengwas_token()
        token_loaded <- TRUE
      }
      message("Extracting direct outcome SNPs for ", model, " / ", outcome_id, "...")
      out_dat <- tryCatch(
        TwoSampleMR::extract_outcome_data(
          snps = requested_snps,
          outcomes = outcome_id,
          proxies = FALSE
        ),
        error = function(e) {
          stop(
            "OpenGWAS direct outcome extraction failed for ", model, " / ", outcome_id,
            ": ", conditionMessage(e), call. = FALSE
          )
        }
      )
      data_source <- "OpenGWAS direct SNP extraction; proxies disabled"
      validate_outcome_data(out_dat, requested_snps, outcome_id, outcome_path)
      write.csv(out_dat, outcome_path, row.names = FALSE)
    }
    validate_outcome_data(out_dat, requested_snps, outcome_id, outcome_path)
    input_paths <- c(input_paths, outcome_path)
    source_log <- c(
      source_log,
      paste(model, outcome_id, data_source, sep = " | ")
    )

    returned_snps <- sort(unique(out_dat$SNP))
    missing_direct <- setdiff(requested_snps, returned_snps)

    mvdat <- TwoSampleMR::mv_harmonise_data(
      mv_exp, out_dat, harmonise_strictness = 2
    )
    available_ids <- colnames(mvdat$exposure_beta)
    if (!all(exposure_ids %in% available_ids)) {
      stop(model, " / ", outcome_id, " lost an exposure column during harmonisation.")
    }

    bx <- mvdat$exposure_beta[, exposure_ids, drop = FALSE]
    sebx <- mvdat$exposure_se[, exposure_ids, drop = FALSE]
    by <- mvdat$outcome_beta
    seby <- mvdat$outcome_se
    harmonised_snps <- rownames(bx)
    if (!length(harmonised_snps) || nrow(bx) <= ncol(bx)) {
      stop(model, " / ", outcome_id, " has too few harmonised instruments.")
    }

    # Standard no-intercept MV-IVW from TwoSampleMR.
    standard <- TwoSampleMR::mv_multiple(
      mvdat,
      pval_threshold = 5e-8,
      intercept = FALSE,
      instrument_specific = FALSE,
      plots = FALSE
    )$result
    standard <- standard[match(exposure_ids, standard$id.exposure), , drop = FALSE]
    if (any(is.na(standard$id.exposure))) {
      stop(model, " / ", outcome_id, " returned incomplete MV-IVW estimates.")
    }

    # Independent matrix implementation of the same weighted no-intercept model.
    manual <- manual_wls(bx, by, seby)

    # MVMR diagnostics are computed after enforcing a documented exposure order.
    mvmr_input <- MVMR::format_mvmr(
      BXGs = bx,
      BYG = by,
      seBXGs = sebx,
      seBYG = seby,
      RSID = harmonised_snps
    )
    strength <- suppressWarnings(MVMR::strength_mvmr(mvmr_input, gencov = 0))
    f_values <- extract_f_values(strength)
    package_q <- suppressWarnings(MVMR::pleiotropy_mvmr(mvmr_input, gencov = 0))

    if (length(f_values) != length(exposure_ids)) {
      stop(model, " / ", outcome_id, " returned the wrong number of conditional F-statistics.")
    }

    for (j in seq_along(exposure_ids)) {
      id <- exposure_ids[j]
      std_row <- standard[standard$id.exposure == id, , drop = FALSE]
      label <- unname(exposure_labels[id])
      estimate_rows[[length(estimate_rows) + 1]] <- data.frame(
        Model = model,
        Smoking_Dimension = models$Smoking_Dimension[m],
        Overlap_Context = models$Overlap_Context[m],
        Outcome = outcome_name,
        Outcome_ID = outcome_id,
        Exposure = label,
        Exposure_ID = id,
        N_Requested = length(requested_snps),
        N_Direct_Outcome = length(returned_snps),
        N_Harmonised = nrow(bx),
        Beta = as.numeric(std_row$b),
        SE = as.numeric(std_row$se),
        CI_Lower = as.numeric(std_row$b) - 1.96 * as.numeric(std_row$se),
        CI_Upper = as.numeric(std_row$b) + 1.96 * as.numeric(std_row$se),
        OR = exp(as.numeric(std_row$b)),
        OR_CI_Lower = exp(as.numeric(std_row$b) - 1.96 * as.numeric(std_row$se)),
        OR_CI_Upper = exp(as.numeric(std_row$b) + 1.96 * as.numeric(std_row$se)),
        P = as.numeric(std_row$pval),
        Estimator = "Standard no-intercept MV-IVW",
        Outcome_Proxy_Rule = "Disabled",
        stringsAsFactors = FALSE
      )
      manual_rows[[length(manual_rows) + 1]] <- data.frame(
        Model = model,
        Outcome_ID = outcome_id,
        Exposure_ID = id,
        Package_Beta = as.numeric(std_row$b),
        Manual_Beta = unname(manual$beta[id]),
        Abs_Beta_Difference = abs(as.numeric(std_row$b) - unname(manual$beta[id])),
        Package_SE = as.numeric(std_row$se),
        Manual_SE = unname(manual$se[id]),
        Abs_SE_Difference = abs(as.numeric(std_row$se) - unname(manual$se[id])),
        Pass_1e_10 = abs(as.numeric(std_row$b) - unname(manual$beta[id])) < 1e-10 &&
          abs(as.numeric(std_row$se) - unname(manual$se[id])) < 1e-10,
        stringsAsFactors = FALSE
      )
      strength_rows[[length(strength_rows) + 1]] <- data.frame(
        Model = model,
        Smoking_Dimension = models$Smoking_Dimension[m],
        Overlap_Context = models$Overlap_Context[m],
        Outcome = outcome_name,
        Outcome_ID = outcome_id,
        Exposure = label,
        Exposure_ID = id,
        N_Harmonised = nrow(bx),
        Conditional_F = f_values[j],
        Strength_Flag = ifelse(f_values[j] >= 10, ">=10", "<10"),
        Covariance_Assumption = "gencov = 0; reference only",
        stringsAsFactors = FALSE
      )
    }

    heterogeneity_rows[[length(heterogeneity_rows) + 1]] <- data.frame(
      Model = model,
      Smoking_Dimension = models$Smoking_Dimension[m],
      Overlap_Context = models$Overlap_Context[m],
      Outcome = outcome_name,
      Outcome_ID = outcome_id,
      N_Harmonised = nrow(bx),
      Standard_Residual_Q = manual$standard_q,
      Standard_Residual_Q_DF = manual$standard_q_df,
      Standard_Residual_Q_P = manual$standard_q_p,
      MVMR_Package_Q = as.numeric(package_q$Qstat),
      MVMR_Package_Q_DF = nrow(bx) - ncol(bx) - 1,
      MVMR_Package_Q_P = as.numeric(package_q$Qpval),
      Package_Q_Covariance_Assumption = "gencov = 0; includes exposure-SE terms",
      stringsAsFactors = FALSE
    )

    all_snps <- sort(unique(c(requested_snps, returned_snps, harmonised_snps)))
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      Model = model,
      Outcome_ID = outcome_id,
      SNP = all_snps,
      Requested = all_snps %in% requested_snps,
      Direct_Outcome_Record = all_snps %in% returned_snps,
      Harmonised_Kept = all_snps %in% harmonised_snps,
      Drop_Reason = ifelse(
        !(all_snps %in% returned_snps),
        "No direct outcome record; proxy prohibited",
        ifelse(!(all_snps %in% harmonised_snps), "Removed during harmonisation", "Kept")
      ),
      stringsAsFactors = FALSE
    )

    harmonised_rows[[length(harmonised_rows) + 1]] <- data.frame(
      Model = model,
      Smoking_Dimension = models$Smoking_Dimension[m],
      Outcome = outcome_name,
      Outcome_ID = outcome_id,
      SNP = harmonised_snps,
      Beta_Stair = bx[, fixed_exposures[["stair"]]],
      SE_Stair = sebx[, fixed_exposures[["stair"]]],
      Beta_BMI = bx[, fixed_exposures[["bmi"]]],
      SE_BMI = sebx[, fixed_exposures[["bmi"]]],
      Smoking_ID = smoking_id,
      Beta_Smoking = bx[, smoking_id],
      SE_Smoking = sebx[, smoking_id],
      Beta_Outcome = by,
      SE_Outcome = seby,
      stringsAsFactors = FALSE
    )

    message(
      model, " / ", outcome_id, ": requested ", length(requested_snps),
      ", direct outcome ", length(returned_snps), ", harmonised ", nrow(bx),
      if (length(missing_direct)) paste0(", direct-missing ", length(missing_direct)) else ""
    )
  }
}

estimates <- do.call(rbind, estimate_rows)
manual_checks <- do.call(rbind, manual_rows)
conditional_f <- do.call(rbind, strength_rows)
heterogeneity <- do.call(rbind, heterogeneity_rows)
snp_audit <- do.call(rbind, audit_rows)
harmonised <- do.call(rbind, harmonised_rows)

if (!all(manual_checks$Pass_1e_10)) {
  stop("At least one independent WLS check failed the 1e-10 tolerance.")
}
if (nrow(unique(estimates[c("Model", "Outcome_ID")])) != 9) {
  stop("Stage B did not produce all 3 models x 3 outcomes.")
}

output_paths <- c(
  Estimates = file.path(result_dir, "SmokingSensitivity_StageB_Estimates.csv"),
  ManualCheck = file.path(result_dir, "SmokingSensitivity_StageB_ManualCheck.csv"),
  ConditionalF = file.path(result_dir, "SmokingSensitivity_StageB_ConditionalF.csv"),
  Heterogeneity = file.path(result_dir, "SmokingSensitivity_StageB_Heterogeneity.csv"),
  SNPAudit = file.path(result_dir, "SmokingSensitivity_StageB_SNPAudit.csv"),
  HarmonisedData = file.path(result_dir, "SmokingSensitivity_StageB_HarmonisedData.csv")
)

write.csv(estimates, output_paths[["Estimates"]], row.names = FALSE)
write.csv(manual_checks, output_paths[["ManualCheck"]], row.names = FALSE)
write.csv(conditional_f, output_paths[["ConditionalF"]], row.names = FALSE)
write.csv(heterogeneity, output_paths[["Heterogeneity"]], row.names = FALSE)
write.csv(snp_audit, output_paths[["SNPAudit"]], row.names = FALSE)
write.csv(harmonised, output_paths[["HarmonisedData"]], row.names = FALSE)

input_paths <- unique(input_paths)
file_info <- file.info(input_paths)
manifest <- data.frame(
  File = basename(input_paths),
  Relative_Path = sub(paste0("^", base_dir, "/"), "", input_paths),
  Bytes = file_info$size,
  MD5 = unname(tools::md5sum(input_paths)),
  SHA256 = vapply(input_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(result_dir, "SmokingSensitivity_StageB_DataManifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

log_path <- file.path(result_dir, "SmokingSensitivity_StageB_RunLog.txt")
run_log <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  paste0("TwoSampleMR version: ", as.character(packageVersion("TwoSampleMR"))),
  paste0("MVMR version: ", as.character(packageVersion("MVMR"))),
  paste0("Refresh requested: ", refresh),
  "Models: GSCAN_CPD, TAG_CPD, GSCAN_INIT",
  "Outcomes: finn-b-COPD_HOSPITAL, finn-b-J10_COPD, finn-b-COPD_EARLY",
  "Outcome extraction: direct SNPs only; proxies = FALSE",
  "Estimator: standard multivariable IVW, no intercept, union instrument set",
  "Independent check: direct matrix WLS with lm-equivalent residual scaling",
  "Conditional F: MVMR::strength_mvmr(gencov = 0); covariance-limited reference",
  "Heterogeneity: standard outcome-SE residual Q plus MVMR::pleiotropy_mvmr(gencov = 0)",
  "No weak-instrument Q-minimisation run in this stage",
  "TAG_INIT excluded by the prespecified Stage-A feasibility rule",
  "",
  "Outcome data sources:",
  paste0("- ", source_log),
  "",
  "Direct-SNP and harmonisation counts:",
  capture.output(print(unique(estimates[c(
    "Model", "Outcome_ID", "N_Requested", "N_Direct_Outcome", "N_Harmonised"
  )]), row.names = FALSE)),
  "",
  "Conditional F-statistics:",
  capture.output(print(conditional_f[c(
    "Model", "Outcome_ID", "Exposure", "Conditional_F", "Strength_Flag"
  )], row.names = FALSE))
)
writeLines(run_log, log_path)

cat("Smoking sensitivity Stage B completed.\n")
cat("All 3 models x 3 outcomes were generated with outcome proxies disabled.\n")
cat("All standard MV-IVW estimates passed the independent WLS check at 1e-10.\n\n")
cat("Outputs:\n")
for (path in c(output_paths, Manifest = manifest_path, RunLog = log_path)) {
  cat("- ", path, "\n", sep = "")
}

