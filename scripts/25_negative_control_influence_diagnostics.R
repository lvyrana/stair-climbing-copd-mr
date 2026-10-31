# ============================================================================
# Post hoc influence diagnostics for negative-control MR analyses
#
# Trigger
#   This diagnostic was initiated after the pre-specified full-direct IVW
#   analyses showed significant Cochran Q heterogeneity for both negative
#   controls. It was not pre-specified as a new primary hypothesis test.
#
# Scope
#   - Reads only the validated, frozen outputs of script 24.
#   - Uses all 16 directly matched instruments for each outcome.
#   - Does not delete instruments, call an API, or modify manuscript files.
#   - Influence thresholds are descriptive screening rules, not automatic
#     exclusion criteria.
# ============================================================================

suppressPackageStartupMessages({
  library(digest)
})

options(stringsAsFactors = FALSE)

result_dir <- "results"
script_path <- "scripts/25_negative_control_influence_diagnostics.R"
harmonised_path <- file.path(result_dir, "NegativeControl_Harmonised.csv")
main_path <- file.path(result_dir, "NegativeControl_Main.csv")
heterogeneity_path <- file.path(result_dir, "NegativeControl_Heterogeneity.csv")
prefix <- file.path(result_dir, "NegativeControlInfluence")
alpha_primary <- 0.025
tolerance <- 1e-10

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  y <- toupper(trimws(as.character(x)))
  !is.na(y) & y %in% c("TRUE", "T", "1", "YES")
}

file_manifest_row <- function(path, role, source_note) {
  info <- file.info(path)
  data.frame(
    File = path,
    Role = role,
    Source_Note = source_note,
    Bytes = unname(info$size),
    Modified_UTC = format(info$mtime, tz = "UTC", usetz = TRUE),
    MD5 = unname(tools::md5sum(path)),
    SHA256 = digest(file = path, algo = "sha256"),
    stringsAsFactors = FALSE
  )
}

manual_ivw <- function(dat) {
  bx <- dat$beta.exposure
  by <- dat$beta.outcome
  w <- 1 / dat$se.outcome^2
  denominator <- sum(w * bx^2)
  stop_if_not(is.finite(denominator) && denominator > 0,
              "IVW denominator is non-positive or non-finite.")
  beta <- sum(w * bx * by) / denominator
  residual <- by - beta * bx
  q_contribution <- w * residual^2
  q <- sum(q_contribution)
  q_df <- nrow(dat) - 1L
  stop_if_not(q_df > 0L, "IVW requires at least two variants.")
  dispersion_raw <- q / q_df
  dispersion_used <- max(1, dispersion_raw)
  se <- sqrt(dispersion_used / denominator)
  p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  data.frame(
    Beta = beta,
    SE = se,
    P = p,
    CI_Lower = beta - 1.96 * se,
    CI_Upper = beta + 1.96 * se,
    Q = q,
    Q_df = q_df,
    Q_P = pchisq(q, df = q_df, lower.tail = FALSE),
    Dispersion_Raw = dispersion_raw,
    Dispersion_Used = dispersion_used,
    stringsAsFactors = FALSE
  )
}

catalogued_confounders <- c(
  "rs12140153", "rs76076331", "rs6449503", "rs7651039", "rs10835211"
)

uncatalogued <- c(
  "rs3898317", "rs72820477", "rs58133129", "rs7577014",
  "rs11689201", "rs1362487", "rs4855841", "rs1829967", "rs7851239"
)

outcome_specs <- data.frame(
  Outcome_Key = c("tanning", "blonde_hair"),
  Outcome_ID = c("ukb-b-533", "ukb-d-1747_1"),
  Outcome_Label = c("Ease of skin tanning", "Natural hair colour: blonde"),
  Beta_Scale = c("SD", "linear_probability"),
  Effect_Measure = c("Beta", "Risk difference"),
  stringsAsFactors = FALSE
)

required_files <- c(harmonised_path, main_path, heterogeneity_path, script_path)
stop_if_not(all(file.exists(required_files)),
            paste("Missing required file(s):",
                  paste(required_files[!file.exists(required_files)], collapse = "; ")))

harmonised <- read.csv(harmonised_path, check.names = FALSE)
main <- read.csv(main_path, check.names = FALSE)
heterogeneity <- read.csv(heterogeneity_path, check.names = FALSE)

required_columns <- c(
  "SNP", "beta.exposure", "beta.outcome", "se.outcome",
  "Outcome_Key", "Primary_Keep"
)
stop_if_not(all(required_columns %in% names(harmonised)),
            "NegativeControl_Harmonised.csv is missing required columns.")

snp_metrics_list <- list()
loo_list <- list()
summary_list <- list()
manual_check_list <- list()

for (i in seq_len(nrow(outcome_specs))) {
  spec <- outcome_specs[i, ]
  dat <- harmonised[
    harmonised$Outcome_Key == spec$Outcome_Key &
      as_flag(harmonised$Primary_Keep),
    , drop = FALSE
  ]
  dat <- dat[order(dat$SNP), , drop = FALSE]

  stop_if_not(nrow(dat) == 16L,
              paste(spec$Outcome_Key, "does not contain exactly 16 retained SNPs."))
  stop_if_not(!anyDuplicated(dat$SNP),
              paste(spec$Outcome_Key, "contains duplicated SNPs."))
  stop_if_not(all(is.finite(dat$beta.exposure)) &&
                all(is.finite(dat$beta.outcome)) &&
                all(is.finite(dat$se.outcome)) &&
                all(dat$se.outcome > 0),
              paste(spec$Outcome_Key, "contains invalid numeric inputs."))

  full <- manual_ivw(dat)
  w <- 1 / dat$se.outcome^2
  residual <- dat$beta.outcome - full$Beta * dat$beta.exposure
  q_contribution <- w * residual^2
  denominator_piece <- w * dat$beta.exposure^2
  leverage_manual <- denominator_piece / sum(denominator_piece)
  mse <- full$Q / full$Q_df

  fit <- lm(beta.outcome ~ 0 + beta.exposure, weights = w, data = dat)
  leverage_lm <- as.numeric(hatvalues(fit))
  cook_lm <- as.numeric(cooks.distance(fit))
  standardised_residual <- as.numeric(rstandard(fit))
  cook_manual <- (q_contribution / mse) *
    leverage_manual / (1 - leverage_manual)^2

  stop_if_not(max(abs(leverage_manual - leverage_lm)) < tolerance,
              paste("Manual and lm leverage disagree for", spec$Outcome_Key))
  stop_if_not(max(abs(cook_manual - cook_lm)) < tolerance,
              paste("Manual and lm Cook distance disagree for", spec$Outcome_Key))
  stop_if_not(abs(sum(q_contribution) - full$Q) < tolerance,
              paste("SNP Q contributions do not sum to Q for", spec$Outcome_Key))

  n <- nrow(dat)
  leverage_cutoff <- 2 / n
  cook_cutoff <- 4 / n
  residual_cutoff <- 2

  snp_metrics <- data.frame(
    Outcome_Key = spec$Outcome_Key,
    Outcome_ID = spec$Outcome_ID,
    Outcome_Label = spec$Outcome_Label,
    Beta_Scale = spec$Beta_Scale,
    Effect_Measure = spec$Effect_Measure,
    SNP = dat$SNP,
    Beta_Exposure = dat$beta.exposure,
    Beta_Outcome = dat$beta.outcome,
    SE_Outcome = dat$se.outcome,
    Wald_Ratio = dat$beta.outcome / dat$beta.exposure,
    Residual = residual,
    Q_Contribution = q_contribution,
    Q_Share = q_contribution / full$Q,
    Leverage = leverage_lm,
    Cook_Distance = cook_lm,
    Standardised_Residual = standardised_residual,
    Leverage_Cutoff_2_over_n = leverage_cutoff,
    Cook_Cutoff_4_over_n = cook_cutoff,
    Absolute_Residual_Cutoff = residual_cutoff,
    High_Leverage_Flag = leverage_lm > leverage_cutoff,
    High_Cook_Flag = cook_lm > cook_cutoff,
    Large_Residual_Flag = abs(standardised_residual) > residual_cutoff,
    Any_Influence_Flag = leverage_lm > leverage_cutoff |
      cook_lm > cook_cutoff |
      abs(standardised_residual) > residual_cutoff,
    Catalogued_BMI_or_Education = dat$SNP %in% catalogued_confounders,
    GWAS_Catalog_Uncatalogued = dat$SNP %in% uncatalogued,
    Audit_Group = ifelse(
      dat$SNP %in% catalogued_confounders,
      "catalogued_BMI_or_education",
      ifelse(dat$SNP %in% uncatalogued,
             "uncatalogued_in_targeted_audit", "other_catalogued")
    ),
    stringsAsFactors = FALSE
  )
  snp_metrics$Q_Contribution_Rank <- rank(
    -snp_metrics$Q_Contribution, ties.method = "min"
  )
  snp_metrics$Cook_Distance_Rank <- rank(
    -snp_metrics$Cook_Distance, ties.method = "min"
  )
  snp_metrics$Leverage_Rank <- rank(
    -snp_metrics$Leverage, ties.method = "min"
  )
  snp_metrics <- snp_metrics[
    order(snp_metrics$Q_Contribution_Rank, snp_metrics$SNP), , drop = FALSE
  ]
  snp_metrics_list[[spec$Outcome_Key]] <- snp_metrics

  loo_rows <- lapply(seq_len(nrow(dat)), function(j) {
    omitted <- dat$SNP[j]
    reduced <- manual_ivw(dat[-j, , drop = FALSE])
    data.frame(
      Outcome_Key = spec$Outcome_Key,
      Outcome_ID = spec$Outcome_ID,
      Outcome_Label = spec$Outcome_Label,
      Beta_Scale = spec$Beta_Scale,
      Effect_Measure = spec$Effect_Measure,
      Omitted_SNP = omitted,
      Omitted_Audit_Group = ifelse(
        omitted %in% catalogued_confounders,
        "catalogued_BMI_or_education",
        ifelse(omitted %in% uncatalogued,
               "uncatalogued_in_targeted_audit", "other_catalogued")
      ),
      N_SNP = nrow(dat) - 1L,
      Full_Beta = full$Beta,
      Full_SE = full$SE,
      Full_P = full$P,
      LOO_Beta = reduced$Beta,
      LOO_SE = reduced$SE,
      LOO_P = reduced$P,
      LOO_CI_Lower = reduced$CI_Lower,
      LOO_CI_Upper = reduced$CI_Upper,
      Beta_Change = reduced$Beta - full$Beta,
      Absolute_Beta_Change = abs(reduced$Beta - full$Beta),
      Standardised_Beta_Change = abs(reduced$Beta - full$Beta) / full$SE,
      LOO_Q = reduced$Q,
      LOO_Q_df = reduced$Q_df,
      LOO_Q_P = reduced$Q_P,
      Passes_Primary_Alpha_0_025 = reduced$P < alpha_primary,
      stringsAsFactors = FALSE
    )
  })
  loo <- do.call(rbind, loo_rows)
  loo$Absolute_Beta_Change_Rank <- rank(
    -loo$Absolute_Beta_Change, ties.method = "min"
  )
  loo <- loo[order(loo$Absolute_Beta_Change_Rank, loo$Omitted_SNP), , drop = FALSE]
  loo_list[[spec$Outcome_Key]] <- loo

  top_q <- snp_metrics[order(snp_metrics$Q_Contribution_Rank), , drop = FALSE]
  top_q <- top_q[seq_len(min(3L, nrow(top_q))), , drop = FALSE]
  top_cook <- snp_metrics[order(snp_metrics$Cook_Distance_Rank), , drop = FALSE]
  top_cook <- top_cook[seq_len(min(3L, nrow(top_cook))), , drop = FALSE]
  largest_change <- loo[order(loo$Absolute_Beta_Change_Rank), , drop = FALSE]
  largest_change <- largest_change[1, , drop = FALSE]

  summary_list[[spec$Outcome_Key]] <- data.frame(
    Outcome_Key = spec$Outcome_Key,
    Outcome_ID = spec$Outcome_ID,
    Outcome_Label = spec$Outcome_Label,
    Effect_Measure = spec$Effect_Measure,
    N_SNP = n,
    Full_Beta = full$Beta,
    Full_SE = full$SE,
    Full_P = full$P,
    Full_CI_Lower = full$CI_Lower,
    Full_CI_Upper = full$CI_Upper,
    Full_Q = full$Q,
    Full_Q_df = full$Q_df,
    Full_Q_P = full$Q_P,
    LOO_Beta_Min = min(loo$LOO_Beta),
    LOO_Beta_Max = max(loo$LOO_Beta),
    LOO_P_Min = min(loo$LOO_P),
    LOO_P_Max = max(loo$LOO_P),
    N_LOO_Passes_Alpha_0_025 = sum(loo$Passes_Primary_Alpha_0_025),
    Any_LOO_Passes_Alpha_0_025 = any(loo$Passes_Primary_Alpha_0_025),
    Largest_Absolute_Beta_Change_SNP = largest_change$Omitted_SNP,
    Largest_Absolute_Beta_Change = largest_change$Absolute_Beta_Change,
    Largest_Standardised_Beta_Change = largest_change$Standardised_Beta_Change,
    Top_3_Q_Contribution_SNPs = paste(top_q$SNP, collapse = ";"),
    Top_3_Q_Shares = paste(format(top_q$Q_Share, digits = 6), collapse = ";"),
    Top_3_Cook_SNPs = paste(top_cook$SNP, collapse = ";"),
    Top_3_Cook_Distances = paste(format(top_cook$Cook_Distance, digits = 6), collapse = ";"),
    N_High_Leverage = sum(snp_metrics$High_Leverage_Flag),
    N_High_Cook = sum(snp_metrics$High_Cook_Flag),
    N_Large_Residual = sum(snp_metrics$Large_Residual_Flag),
    N_Any_Influence_Flag = sum(snp_metrics$Any_Influence_Flag),
    Post_Hoc_Trigger = "Significant full-direct IVW Cochran Q heterogeneity",
    Interpretation_Boundary = paste(
      "Descriptive influence diagnostic only; flags do not justify automatic",
      "SNP exclusion and do not create a new primary result"
    ),
    stringsAsFactors = FALSE
  )

  main_row <- main[
    main$Outcome_Key == spec$Outcome_Key &
      main$Instrument_Set == "full_direct" &
      main$method == "Inverse variance weighted",
    , drop = FALSE
  ]
  q_row <- heterogeneity[
    heterogeneity$Outcome_Key == spec$Outcome_Key &
      heterogeneity$Instrument_Set == "full_direct" &
      heterogeneity$method == "Inverse variance weighted",
    , drop = FALSE
  ]
  stop_if_not(nrow(main_row) == 1L,
              paste("Missing unique full-direct IVW row for", spec$Outcome_Key))
  stop_if_not(nrow(q_row) == 1L,
              paste("Missing unique full-direct heterogeneity row for", spec$Outcome_Key))

  manual_check_list[[spec$Outcome_Key]] <- data.frame(
    Outcome_Key = spec$Outcome_Key,
    Recalculated_Beta = full$Beta,
    Stored_TwoSampleMR_Beta = main_row$b,
    Beta_Absolute_Difference = abs(full$Beta - main_row$b),
    Recalculated_SE = full$SE,
    Stored_TwoSampleMR_SE = main_row$se,
    SE_Absolute_Difference = abs(full$SE - main_row$se),
    Recalculated_Q = full$Q,
    Stored_TwoSampleMR_Q = q_row$Q,
    Q_Absolute_Difference = abs(full$Q - q_row$Q),
    Sum_SNP_Q_Contributions = sum(q_contribution),
    Q_Sum_Absolute_Difference = abs(sum(q_contribution) - full$Q),
    Max_Leverage_Formula_Difference = max(abs(leverage_manual - leverage_lm)),
    Max_Cook_Formula_Difference = max(abs(cook_manual - cook_lm)),
    stringsAsFactors = FALSE
  )
}

snp_metrics_all <- do.call(rbind, snp_metrics_list)
loo_all <- do.call(rbind, loo_list)
summary_all <- do.call(rbind, summary_list)
manual_check_all <- do.call(rbind, manual_check_list)

stop_if_not(max(manual_check_all$Beta_Absolute_Difference) < tolerance,
            "Recalculated and stored IVW beta differ by >= 1e-10.")
stop_if_not(max(manual_check_all$SE_Absolute_Difference) < tolerance,
            "Recalculated and stored IVW SE differ by >= 1e-10.")
stop_if_not(max(manual_check_all$Q_Absolute_Difference) < tolerance,
            "Recalculated and stored IVW Q differ by >= 1e-10.")
stop_if_not(max(manual_check_all$Q_Sum_Absolute_Difference) < tolerance,
            "SNP Q contributions do not reconstruct total Q within tolerance.")

write.csv(snp_metrics_all, paste0(prefix, "_SNPMetrics.csv"),
          row.names = FALSE, na = "")
write.csv(loo_all, paste0(prefix, "_LeaveOneOut.csv"),
          row.names = FALSE, na = "")
write.csv(summary_all, paste0(prefix, "_Summary.csv"),
          row.names = FALSE, na = "")
write.csv(manual_check_all, paste0(prefix, "_ManualCheck.csv"),
          row.names = FALSE, na = "")

manifest <- rbind(
  file_manifest_row(
    harmonised_path,
    "Validated harmonised input",
    "Full-direct negative-control data generated by script 24"
  ),
  file_manifest_row(
    main_path,
    "Stored IVW cross-check source",
    "Validated negative-control main results generated by script 24"
  ),
  file_manifest_row(
    heterogeneity_path,
    "Stored heterogeneity cross-check source",
    "Validated negative-control heterogeneity generated by script 24"
  ),
  file_manifest_row(
    script_path,
    "Post hoc analysis script",
    "Triggered by significant Cochran Q heterogeneity; offline only"
  )
)
write.csv(manifest, paste0(prefix, "_DataManifest.csv"),
          row.names = FALSE, na = "")

log_lines <- c(
  "Post hoc negative-control influence diagnostics",
  paste("Run time UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R version:", R.version.string),
  paste("digest version:", as.character(packageVersion("digest"))),
  "Network/API calls: none",
  "Input: validated full-direct harmonised output from script 24",
  "Variants per outcome: 16 direct SNPs",
  "Trigger: significant full-direct IVW Cochran Q heterogeneity",
  paste("Primary alpha retained for descriptive threshold crossing check:",
        alpha_primary),
  "Influence screening: leverage > 2/n; Cook distance > 4/n; |standardised residual| > 2",
  paste(
    "Interpretation boundary:",
    "post hoc diagnostic only; no SNP is automatically excluded and no new",
    "primary hypothesis is created"
  ),
  "",
  capture.output(print(summary_all, row.names = FALSE))
)
writeLines(log_lines, paste0(prefix, "_RunLog.txt"))

cat(paste(log_lines, collapse = "\n"), "\n")
