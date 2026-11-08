# ============================================================================
# Offline, frozen-data rerun of univariable MR analyses
#
# Purpose
#   Re-run every univariable outcome from the already-saved harmonised data,
#   using one uniform primary-analysis rule: outcome proxies are excluded.
#
# Guardrails
#   - No OpenGWAS/API calls.
#   - No re-extraction or re-harmonisation.
#   - No existing manuscript table or figure is overwritten.
#   - FinnGen outcomes must contain exactly one proxy row mapping
#     rs7851239 -> rs13293916; direct outcomes must contain no proxy rows.
# ============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(digest)
})

options(stringsAsFactors = FALSE)

source_dir <- "inputs"
result_dir <- "results"
prefix <- file.path(result_dir, "Frozen_NoProxy_UVMR")
base_seed <- 20260808L
presso_distributions <- 1000L

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  y <- toupper(trimws(as.character(x)))
  !is.na(y) & y %in% c("TRUE", "T", "1", "YES")
}

nonempty <- function(x) {
  y <- trimws(as.character(x))
  !is.na(x) & !is.na(y) & nzchar(y) & toupper(y) != "NA"
}

detect_outcome_proxy <- function(dat) {
  out <- rep(FALSE, nrow(dat))

  if ("proxy.outcome" %in% names(dat)) {
    out <- out | as_flag(dat$proxy.outcome)
  }

  if (all(c("target_snp.outcome", "proxy_snp.outcome") %in% names(dat))) {
    has_mapping <- nonempty(dat$target_snp.outcome) &
      nonempty(dat$proxy_snp.outcome) &
      dat$target_snp.outcome != dat$proxy_snp.outcome
    out <- out | has_mapping
  }

  out
}

parse_p_threshold <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x[1])) return(NA_real_)
  suppressWarnings(as.numeric(sub("^<", "", as.character(x[1]))))
}

scalar_or_na <- function(x) {
  if (length(x) == 0 || is.null(x)) return(NA)
  unname(x[1])
}

first_or_na <- function(dat, column, row_index = 1L) {
  if (is.null(dat) || !is.data.frame(dat) || !column %in% names(dat) ||
      nrow(dat) < row_index) return(NA)
  dat[[column]][row_index]
}

add_context <- function(dat, outcome, n_source, n_proxy_removed, n_primary) {
  cbind(
    data.frame(
      Outcome_Key = outcome$key,
      Outcome = outcome$label,
      Outcome_ID = outcome$id,
      Outcome_Type = outcome$type,
      Outcome_Beta_Scale = outcome$beta_scale,
      Source_Rows_MR_Keep = n_source,
      Proxies_Removed = n_proxy_removed,
      Primary_N_SNPs = n_primary,
      Analysis_Rule = "Frozen harmonised data; outcome proxies excluded",
      stringsAsFactors = FALSE
    ),
    dat
  )
}

manual_ivw <- function(dat) {
  bx <- dat$beta.exposure
  by <- dat$beta.outcome
  w <- 1 / (dat$se.outcome^2)
  denominator <- sum(w * bx^2)
  b <- sum(w * bx * by) / denominator
  q <- sum(w * (by - b * bx)^2)
  q_df <- nrow(dat) - 1L
  dispersion <- max(1, q / q_df)
  se <- sqrt(dispersion / denominator)

  data.frame(
    Manual_IVW_Beta = b,
    Manual_IVW_SE = se,
    Manual_IVW_Q = q,
    Manual_IVW_Q_df = q_df,
    Manual_Dispersion = dispersion,
    stringsAsFactors = FALSE
  )
}

run_presso <- function(dat, outcome, seed) {
  warnings_seen <- character()
  result <- tryCatch(
    withCallingHandlers(
      {
        set.seed(seed)
        mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = dat,
          NbDistribution = presso_distributions,
          SignifThreshold = 0.05
        )
      },
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    return(data.frame(
      Outcome_Key = outcome$key,
      Outcome = outcome$label,
      Outcome_ID = outcome$id,
      Primary_N_SNPs = nrow(dat),
      Seed = seed,
      NbDistribution = presso_distributions,
      Global_RSSobs = NA_real_,
      Global_P_Text = NA_character_,
      Global_P_Threshold = NA_real_,
      Raw_Beta = NA_real_, Raw_SE = NA_real_, Raw_P = NA_real_,
      Outlier_Count = NA_integer_, Outlier_SNPs = NA_character_,
      Corrected_Beta = NA_real_, Corrected_SE = NA_real_, Corrected_P = NA_real_,
      Distortion_Coefficient = NA_real_,
      Distortion_P_Text = NA_character_,
      Distortion_P_Threshold = NA_real_,
      Notes = paste("ERROR:", conditionMessage(result)),
      stringsAsFactors = FALSE
    ))
  }

  main <- result[["Main MR results"]]
  presso <- result[["MR-PRESSO results"]]
  global <- presso[["Global Test"]]
  distortion <- presso[["Distortion Test"]]

  outlier_indices <- integer()
  if (!is.null(distortion) && !is.null(distortion[["Outliers Indices"]])) {
    outlier_indices <- as.integer(distortion[["Outliers Indices"]])
    outlier_indices <- outlier_indices[
      !is.na(outlier_indices) & outlier_indices >= 1L & outlier_indices <= nrow(dat)
    ]
  }

  global_p_text <- as.character(scalar_or_na(global[["Pvalue"]]))
  distortion_p_text <- if (is.null(distortion)) NA_character_ else
    as.character(scalar_or_na(distortion[["Pvalue"]]))

  data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Primary_N_SNPs = nrow(dat),
    Seed = seed,
    NbDistribution = presso_distributions,
    Global_RSSobs = as.numeric(scalar_or_na(global[["RSSobs"]])),
    Global_P_Text = global_p_text,
    Global_P_Threshold = parse_p_threshold(global_p_text),
    Raw_Beta = as.numeric(first_or_na(main, "Causal Estimate", 1L)),
    Raw_SE = as.numeric(first_or_na(main, "Sd", 1L)),
    Raw_P = as.numeric(first_or_na(main, "P-value", 1L)),
    Outlier_Count = length(outlier_indices),
    Outlier_SNPs = if (length(outlier_indices) == 0) "" else
      paste(dat$SNP[outlier_indices], collapse = ";"),
    Corrected_Beta = as.numeric(first_or_na(main, "Causal Estimate", 2L)),
    Corrected_SE = as.numeric(first_or_na(main, "Sd", 2L)),
    Corrected_P = as.numeric(first_or_na(main, "P-value", 2L)),
    Distortion_Coefficient = if (is.null(distortion)) NA_real_ else
      as.numeric(scalar_or_na(distortion[["Distortion Coefficient"]])),
    Distortion_P_Text = distortion_p_text,
    Distortion_P_Threshold = parse_p_threshold(distortion_p_text),
    Notes = paste(unique(warnings_seen), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

outcomes <- list(
  list(key = "COPD_main", label = "COPD differential diagnosis",
       id = "ukb-d-COPD_EXCL", type = "Binary",
       beta_scale = "linear_probability", expected_proxy = 0L,
       expected_primary = 16L),
  list(key = "COPD_hospital", label = "COPD hospitalization",
       id = "finn-b-COPD_HOSPITAL", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "COPD_early", label = "Early-onset COPD",
       id = "finn-b-COPD_EARLY", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "COPD_later", label = "Late-onset COPD",
       id = "finn-b-COPD_LATER", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "COPD_general", label = "General COPD",
       id = "finn-b-J10_COPD", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "COPD_emphysema", label = "Emphysema",
       id = "finn-b-J10_EMPHYSEMA", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "COPD_insufficiency", label = "COPD respiratory insufficiency",
       id = "finn-b-COPD_INSUFFICIENCY", type = "Binary",
       beta_scale = "log_odds", expected_proxy = 1L,
       expected_primary = 15L),
  list(key = "FEV1", label = "FEV1", id = "ebi-a-GCST007432",
       type = "Continuous", beta_scale = "continuous",
       expected_proxy = 0L, expected_primary = 16L),
  list(key = "FVC", label = "FVC", id = "ebi-a-GCST007429",
       type = "Continuous", beta_scale = "continuous",
       expected_proxy = 0L, expected_primary = 16L),
  list(key = "PEF", label = "PEF", id = "ebi-a-GCST007430",
       type = "Continuous", beta_scale = "continuous",
       expected_proxy = 0L, expected_primary = 16L)
)

manifest_rows <- list()
snp_audit_rows <- list()
main_rows <- list()
manual_ivw_rows <- list()
heterogeneity_rows <- list()
pleiotropy_rows <- list()
presso_rows <- list()
loo_rows <- list()
single_rows <- list()
steiger_rows <- list()

for (i in seq_along(outcomes)) {
  outcome <- outcomes[[i]]
  source_file <- file.path(
    source_dir,
    paste0("harmonised_steiger_", outcome$key, ".csv")
  )
  stop_if_not(file.exists(source_file), paste("Missing source file:", source_file))

  dat_all <- read.csv(source_file, stringsAsFactors = FALSE, check.names = FALSE)
  stop_if_not("mr_keep" %in% names(dat_all), paste("mr_keep missing in", source_file))
  stop_if_not(all(c(
    "SNP", "beta.exposure", "se.exposure", "pval.exposure",
    "samplesize.exposure", "beta.outcome", "se.outcome", "pval.outcome",
    "samplesize.outcome", "r.exposure", "r.outcome"
  ) %in% names(dat_all)), paste("Required MR columns missing in", source_file))

  mr_keep <- as_flag(dat_all$mr_keep)
  is_proxy <- detect_outcome_proxy(dat_all)
  is_primary <- mr_keep & !is_proxy

  n_source <- sum(mr_keep)
  n_proxy_removed <- sum(mr_keep & is_proxy)
  n_primary <- sum(is_primary)

  stop_if_not(
    n_proxy_removed == outcome$expected_proxy,
    sprintf(
      "%s: expected %d proxy rows but found %d",
      outcome$key, outcome$expected_proxy, n_proxy_removed
    )
  )
  stop_if_not(
    n_primary == outcome$expected_primary,
    sprintf(
      "%s: expected %d primary SNPs but found %d",
      outcome$key, outcome$expected_primary, n_primary
    )
  )

  proxy_targets <- if ("target_snp.outcome" %in% names(dat_all))
    unique(dat_all$target_snp.outcome[mr_keep & is_proxy]) else character()
  proxy_substitutes <- if ("proxy_snp.outcome" %in% names(dat_all))
    unique(dat_all$proxy_snp.outcome[mr_keep & is_proxy]) else character()

  proxy_targets <- proxy_targets[nonempty(proxy_targets)]
  proxy_substitutes <- proxy_substitutes[nonempty(proxy_substitutes)]

  if (outcome$expected_proxy == 1L) {
    stop_if_not(
      identical(proxy_targets, "rs7851239") &&
        identical(proxy_substitutes, "rs13293916"),
      paste(
        outcome$key,
        "does not contain the expected rs7851239 -> rs13293916 proxy mapping"
      )
    )
  }

  info <- file.info(source_file)
  manifest_rows[[i]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Outcome_Type = outcome$type,
    Outcome_Beta_Scale = outcome$beta_scale,
    Source_File = source_file,
    Source_Bytes = info$size,
    Source_Modified_UTC = format(info$mtime, tz = "UTC", usetz = TRUE),
    Source_MD5 = unname(tools::md5sum(source_file)),
    Source_SHA256 = digest(source_file, algo = "sha256", file = TRUE),
    Source_Total_Rows = nrow(dat_all),
    Source_Rows_MR_Keep = n_source,
    Proxies_Removed = n_proxy_removed,
    Primary_N_SNPs = n_primary,
    Proxy_Targets = paste(proxy_targets, collapse = ";"),
    Proxy_Substitutes = paste(proxy_substitutes, collapse = ";"),
    Analysis_Rule = "Frozen harmonised data; outcome proxies excluded",
    stringsAsFactors = FALSE
  )

  snp_audit_rows[[i]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Source_Row = seq_len(nrow(dat_all)),
    SNP_Label = dat_all$SNP,
    MR_Keep = mr_keep,
    Is_Outcome_Proxy = is_proxy,
    Kept_In_Primary = is_primary,
    Exclusion_Reason = ifelse(
      !mr_keep,
      "harmonisation_mr_keep_FALSE",
      ifelse(is_proxy, "outcome_proxy_excluded_from_primary", "")
    ),
    Proxy_Target = if ("target_snp.outcome" %in% names(dat_all))
      dat_all$target_snp.outcome else NA_character_,
    Proxy_Substitute = if ("proxy_snp.outcome" %in% names(dat_all))
      dat_all$proxy_snp.outcome else NA_character_,
    stringsAsFactors = FALSE
  )

  dat <- dat_all[is_primary, , drop = FALSE]
  rownames(dat) <- NULL
  set.seed(base_seed + i)

  mr_result <- mr(dat)
  # CIs use the same reference distribution as each estimator's P value in
  # TwoSampleMR: t with nsnp - 2 df for MR-Egger, t with nsnp - 1 df for the
  # simple and weighted mode, and the normal distribution otherwise.
  z <- ifelse(
    mr_result$method == "MR Egger", qt(0.975, mr_result$nsnp - 2),
    ifelse(mr_result$method %in% c("Simple mode", "Weighted mode"),
           qt(0.975, mr_result$nsnp - 1), qnorm(0.975))
  )
  mr_result$CI_Lower <- mr_result$b - z * mr_result$se
  mr_result$CI_Upper <- mr_result$b + z * mr_result$se
  if (outcome$beta_scale == "log_odds") {
    mr_result$Effect_Measure <- "OR"
    mr_result$Effect_Estimate <- exp(mr_result$b)
    mr_result$Effect_CI_Lower <- exp(mr_result$CI_Lower)
    mr_result$Effect_CI_Upper <- exp(mr_result$CI_Upper)
  } else if (outcome$beta_scale == "linear_probability") {
    mr_result$Effect_Measure <- "Risk difference"
    mr_result$Effect_Estimate <- mr_result$b
    mr_result$Effect_CI_Lower <- mr_result$CI_Lower
    mr_result$Effect_CI_Upper <- mr_result$CI_Upper
  } else {
    mr_result$Effect_Measure <- "Beta"
    mr_result$Effect_Estimate <- mr_result$b
    mr_result$Effect_CI_Lower <- mr_result$CI_Lower
    mr_result$Effect_CI_Upper <- mr_result$CI_Upper
  }
  main_rows[[i]] <- add_context(
    mr_result, outcome, n_source, n_proxy_removed, n_primary
  )

  manual <- manual_ivw(dat)
  ivw_row <- mr_result[mr_result$method == "Inverse variance weighted", , drop = FALSE]
  stop_if_not(nrow(ivw_row) == 1L, paste("IVW result missing for", outcome$key))
  manual$TwoSampleMR_IVW_Beta <- ivw_row$b
  manual$TwoSampleMR_IVW_SE <- ivw_row$se
  manual$Beta_Absolute_Difference <- abs(manual$Manual_IVW_Beta - ivw_row$b)
  manual$SE_Absolute_Difference <- abs(manual$Manual_IVW_SE - ivw_row$se)
  manual$Pass_1e_10 <- manual$Beta_Absolute_Difference < 1e-10 &
    manual$SE_Absolute_Difference < 1e-10
  stop_if_not(manual$Pass_1e_10, paste("Manual IVW check failed for", outcome$key))
  manual_ivw_rows[[i]] <- add_context(
    manual, outcome, n_source, n_proxy_removed, n_primary
  )

  heterogeneity_rows[[i]] <- add_context(
    mr_heterogeneity(dat), outcome, n_source, n_proxy_removed, n_primary
  )
  pleiotropy_rows[[i]] <- add_context(
    mr_pleiotropy_test(dat), outcome, n_source, n_proxy_removed, n_primary
  )
  loo_rows[[i]] <- add_context(
    mr_leaveoneout(dat), outcome, n_source, n_proxy_removed, n_primary
  )
  single_rows[[i]] <- add_context(
    mr_singlesnp(dat), outcome, n_source, n_proxy_removed, n_primary
  )

  if (outcome$beta_scale == "linear_probability") {
    steiger_r_outcome <- get_r_from_bsen(
      dat$beta.outcome,
      dat$se.outcome,
      dat$samplesize.outcome
    )
    steiger_r_outcome_method <-
      "get_r_from_bsen; linear probability outcome"
  } else {
    steiger_r_outcome <- dat$r.outcome
    steiger_r_outcome_method <- if (outcome$beta_scale == "log_odds") {
      "Frozen get_r_from_lor"
    } else {
      "Frozen get_r_from_bsen"
    }
  }

  st <- mr_steiger(
    p_exp = dat$pval.exposure,
    p_out = dat$pval.outcome,
    n_exp = dat$samplesize.exposure,
    n_out = dat$samplesize.outcome,
    r_exp = dat$r.exposure,
    r_out = steiger_r_outcome
  )
  steiger_rows[[i]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Outcome_Type = outcome$type,
    Outcome_Beta_Scale = outcome$beta_scale,
    Primary_N_SNPs = n_primary,
    R2_Exposure = st$r2_exp,
    R2_Outcome = st$r2_out,
    Correct_Causal_Direction = st$correct_causal_direction,
    Steiger_P = st$steiger_test,
    Sensitivity_Ratio = st$sensitivity_ratio,
    Steiger_R_Outcome_Method = steiger_r_outcome_method,
    Analysis_Rule = "Frozen harmonised data; outcome proxies excluded",
    stringsAsFactors = FALSE
  )

  presso_rows[[i]] <- run_presso(dat, outcome, base_seed + 1000L + i)
}

manifest <- do.call(rbind, manifest_rows)
snp_audit <- do.call(rbind, snp_audit_rows)
main_results <- do.call(rbind, main_rows)
manual_ivw_results <- do.call(rbind, manual_ivw_rows)
heterogeneity <- do.call(rbind, heterogeneity_rows)
pleiotropy <- do.call(rbind, pleiotropy_rows)
presso <- do.call(rbind, presso_rows)
leave_one_out <- do.call(rbind, loo_rows)
single_snp <- do.call(rbind, single_rows)
steiger <- do.call(rbind, steiger_rows)

write.csv(manifest, paste0(prefix, "_DataManifest.csv"), row.names = FALSE, na = "")
write.csv(snp_audit, paste0(prefix, "_SNPAudit.csv"), row.names = FALSE, na = "")
write.csv(main_results, paste0(prefix, "_Main.csv"), row.names = FALSE, na = "")
write.csv(manual_ivw_results, paste0(prefix, "_ManualIVWCheck.csv"), row.names = FALSE, na = "")
write.csv(heterogeneity, paste0(prefix, "_Heterogeneity.csv"), row.names = FALSE, na = "")
write.csv(pleiotropy, paste0(prefix, "_Pleiotropy.csv"), row.names = FALSE, na = "")
write.csv(presso, paste0(prefix, "_MRPRESSO.csv"), row.names = FALSE, na = "")
write.csv(leave_one_out, paste0(prefix, "_LeaveOneOut.csv"), row.names = FALSE, na = "")
write.csv(single_snp, paste0(prefix, "_SingleSNP.csv"), row.names = FALSE, na = "")
write.csv(steiger, paste0(prefix, "_Steiger.csv"), row.names = FALSE, na = "")

ivw_summary <- main_results[
  main_results$method == "Inverse variance weighted",
  c(
    "Outcome_Key", "Outcome", "Outcome_ID", "Outcome_Type",
    "Outcome_Beta_Scale",
    "Proxies_Removed", "Primary_N_SNPs", "b", "se", "pval",
    "Effect_Measure", "Effect_Estimate", "Effect_CI_Lower", "Effect_CI_Upper"
  )
]

run_log <- c(
  "Frozen no-proxy univariable MR rerun",
  "====================================",
  paste(
    "Script:",
    "scripts/23_rerun_frozen_no_proxy_uvmr.R"
  ),
  paste("Execution time:", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste("Offline source directory:", source_dir),
  "Network/API extraction: NONE",
  "Primary rule: exclude outcome proxies from frozen harmonised data",
  paste(
    "Outcome scales:",
    "ukb-d-COPD_EXCL=linear probability; FinnGen=log odds;",
    "lung function=continuous"
  ),
  paste(
    "UKB Steiger handling:",
    "r.outcome recomputed with get_r_from_bsen; frozen log-OR conversion ignored"
  ),
  paste("Base random seed:", base_seed),
  paste("MR-PRESSO NbDistribution:", presso_distributions),
  paste("R version:", R.version.string),
  paste("TwoSampleMR version:", as.character(packageVersion("TwoSampleMR"))),
  paste("MRPRESSO version:", as.character(packageVersion("MRPRESSO"))),
  paste("digest version:", as.character(packageVersion("digest"))),
  "",
  "Primary IVW summary:",
  capture.output(print(ivw_summary, row.names = FALSE)),
  "",
  "Manual IVW validation:",
  capture.output(print(
    manual_ivw_results[, c(
      "Outcome_Key", "Primary_N_SNPs", "Beta_Absolute_Difference",
      "SE_Absolute_Difference", "Pass_1e_10"
    )],
    row.names = FALSE
  )),
  "",
  "Session information:",
  capture.output(sessionInfo())
)
writeLines(run_log, paste0(prefix, "_RunLog.txt"), useBytes = TRUE)

cat("Offline rerun completed successfully.\n")
cat(sprintf("Primary outcomes analysed: %d\n", length(outcomes)))
cat(sprintf("Proxy rows excluded: %d\n", sum(manifest$Proxies_Removed)))
cat(sprintf("All manual IVW checks passed: %s\n", all(manual_ivw_results$Pass_1e_10)))
print(ivw_summary, row.names = FALSE)
