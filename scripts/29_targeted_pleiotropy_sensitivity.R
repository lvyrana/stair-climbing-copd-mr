# ============================================================================
# Targeted pleiotropy sensitivity analysis from frozen harmonised inputs
#
# Purpose
#   Reproduce, without network access, the sensitivity analyses based on the
#   2026-08-08 targeted SNP annotation audit:
#     1) full direct-SNP IVW;
#     2) IVW after excluding five SNPs with catalogued BMI/education links;
#     3) IVW in the uncatalogued-SNP subset;
#     4) exact same-size subset permutation analysis;
#     5) formal comparison of the disjoint known-link and uncatalogued subsets;
#     6) per-SNP IVW influence diagnostics.
#
# Evidence boundary
#   The downstream sensitivity analysis is reproducible from frozen local data.
#   The external GWAS Catalog/Ensembl/OpenGWAS annotation query is NOT rerun by
#   this script. Its classifications are treated as a frozen 2026-08-08 input.
#   The raw API responses and original query code were not preserved, so the
#   external annotation itself is auditable but not independently rerunnable.
#
# Guardrails
#   - No API or network calls.
#   - No re-extraction, clumping, proxy search, or harmonisation.
#   - Outcome proxies are excluded from every primary and sensitivity analysis.
#   - No manuscript, table, or figure is overwritten.
# ============================================================================

suppressPackageStartupMessages({
  library(digest)
})

options(stringsAsFactors = FALSE)

source_dir <- "inputs"
result_dir <- "results"
annotation_file <- "inputs/instrument_annotation.csv"
analysis_date <- "2026-08-09"
annotation_date <- "2026-08-08"

prefix <- file.path(result_dir, "TargetedPleiotropy")

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  y <- toupper(trimws(as.character(x)))
  !is.na(y) & y %in% c("TRUE", "T", "1", "YES", "是")
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

manual_ivw <- function(dat) {
  stop_if_not(nrow(dat) >= 2L, "IVW requires at least two SNPs")
  bx <- as.numeric(dat$beta.exposure)
  by <- as.numeric(dat$beta.outcome)
  sey <- as.numeric(dat$se.outcome)
  stop_if_not(all(is.finite(bx) & is.finite(by) & is.finite(sey) & sey > 0),
              "Non-finite or invalid frozen MR inputs")

  w <- 1 / sey^2
  denominator <- sum(w * bx^2)
  stop_if_not(is.finite(denominator) && denominator > 0,
              "Invalid IVW denominator")
  b <- sum(w * bx * by) / denominator
  q <- sum(w * (by - b * bx)^2)
  q_df <- nrow(dat) - 1L
  dispersion <- max(1, q / q_df)
  se <- sqrt(dispersion / denominator)
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)

  list(
    b = b,
    se = se,
    ci_lo = b - 1.96 * se,
    ci_hi = b + 1.96 * se,
    p = p,
    q = q,
    q_df = q_df,
    q_p = pchisq(q, df = q_df, lower.tail = FALSE),
    dispersion = dispersion,
    denominator = denominator
  )
}

estimate_row <- function(dat, outcome, set_name, set_definition) {
  fit <- manual_ivw(dat)
  is_log_odds <- outcome$beta_scale == "log_odds"

  data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Outcome_Type = outcome$type,
    Beta_Scale = outcome$beta_scale,
    Set = set_name,
    Set_Definition = set_definition,
    N_SNPs = nrow(dat),
    SNPs = paste(sort(dat$SNP), collapse = ";"),
    Beta = fit$b,
    SE = fit$se,
    Beta_CI_Lower = fit$ci_lo,
    Beta_CI_Upper = fit$ci_hi,
    OR = if (is_log_odds) exp(fit$b) else NA_real_,
    OR_CI_Lower = if (is_log_odds) exp(fit$ci_lo) else NA_real_,
    OR_CI_Upper = if (is_log_odds) exp(fit$ci_hi) else NA_real_,
    P = fit$p,
    Q = fit$q,
    Q_df = fit$q_df,
    Q_P = fit$q_p,
    Dispersion = fit$dispersion,
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

annotation <- read.csv(
  annotation_file,
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_annotation_columns <- c(
  "rsid", "in_paper_16", "n_catalog_assoc", "n_adiposity", "n_ses_edu"
)
stop_if_not(all(required_annotation_columns %in% names(annotation)),
            "Frozen SNP annotation file is missing required columns")

paper_annotation <- annotation[as_flag(annotation$in_paper_16), , drop = FALSE]
stop_if_not(nrow(paper_annotation) == 16L,
            "Expected exactly 16 paper instruments in the frozen annotation")
stop_if_not(length(unique(paper_annotation$rsid)) == 16L,
            "Paper instrument rsIDs are not unique in the frozen annotation")

num_or_zero <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.na(y)] <- 0
  y
}

paper_annotation$Frozen_Class <- ifelse(
  num_or_zero(paper_annotation$n_adiposity) > 0 |
    num_or_zero(paper_annotation$n_ses_edu) > 0,
  "catalogued_BMI_or_education",
  ifelse(
    num_or_zero(paper_annotation$n_catalog_assoc) == 0,
    "no_catalogue_record",
    "catalogued_other"
  )
)

known_link_snps <- sort(
  paper_annotation$rsid[
    paper_annotation$Frozen_Class == "catalogued_BMI_or_education"
  ]
)
uncatalogued_snps <- sort(
  paper_annotation$rsid[
    paper_annotation$Frozen_Class == "no_catalogue_record"
  ]
)

expected_known_link_snps <- sort(c(
  "rs12140153", "rs76076331", "rs6449503", "rs7651039", "rs10835211"
))
stop_if_not(identical(known_link_snps, expected_known_link_snps),
            "Frozen BMI/education-linked SNP set differs from the audited five")
stop_if_not(length(uncatalogued_snps) == 9L,
            "Expected exactly nine uncatalogued paper instruments")

paper_annotation$Annotation_Date <- annotation_date
paper_annotation$Annotation_Status <-
  "Frozen summary is auditable; raw API responses/query code were not preserved"

manifest_rows <- list()
snp_audit_rows <- list()
subset_rows <- list()
permutation_rows <- list()
permutation_summary_rows <- list()
subset_difference_rows <- list()
influence_rows <- list()

for (outcome in outcomes) {
  input_file <- file.path(
    source_dir,
    paste0("harmonised_steiger_", outcome$key, ".csv")
  )
  stop_if_not(file.exists(input_file), paste("Missing frozen input:", input_file))

  dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  stop_if_not(all(c("SNP", "beta.exposure", "beta.outcome", "se.outcome",
                    "mr_keep") %in% names(dat)),
              paste("Missing required MR columns in", input_file))

  keep <- as_flag(dat$mr_keep)
  proxy <- detect_outcome_proxy(dat)
  n_source <- sum(keep)
  n_proxy <- sum(keep & proxy)
  primary <- dat[keep & !proxy, , drop = FALSE]

  stop_if_not(n_proxy == outcome$expected_proxy,
              paste(outcome$key, "unexpected proxy count:", n_proxy))
  stop_if_not(nrow(primary) == outcome$expected_primary,
              paste(outcome$key, "unexpected primary SNP count:", nrow(primary)))
  stop_if_not(length(unique(primary$SNP)) == nrow(primary),
              paste(outcome$key, "duplicate primary SNPs"))
  stop_if_not(all(primary$SNP %in% paper_annotation$rsid),
              paste(outcome$key, "contains an SNP absent from the frozen audit"))

  known_present <- intersect(known_link_snps, primary$SNP)
  uncatalogued_present <- intersect(uncatalogued_snps, primary$SNP)
  stop_if_not(length(known_present) == 5L,
              paste(outcome$key, "does not contain all five audited SNPs"))
  expected_uncatalogued <- if (outcome$expected_proxy == 1L) 8L else 9L
  stop_if_not(length(uncatalogued_present) == expected_uncatalogued,
              paste(outcome$key, "unexpected uncatalogued-SNP count"))

  known_only <- primary[primary$SNP %in% known_present, , drop = FALSE]
  exclude_known <- primary[!primary$SNP %in% known_present, , drop = FALSE]
  uncatalogued_only <- primary[
    primary$SNP %in% uncatalogued_present, , drop = FALSE
  ]

  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Outcome_ID = outcome$id,
    Input_File = input_file,
    Input_SHA256 = digest(input_file, algo = "sha256", file = TRUE),
    Source_MR_Keep_Rows = n_source,
    Proxies_Removed = n_proxy,
    Primary_Direct_SNPs = nrow(primary),
    Known_Link_SNPs_Present = length(known_present),
    Uncatalogued_SNPs_Present = length(uncatalogued_present),
    stringsAsFactors = FALSE
  )

  primary_annotation <- paper_annotation[
    match(primary$SNP, paper_annotation$rsid), , drop = FALSE
  ]
  stop_if_not(all(primary_annotation$rsid == primary$SNP),
              paste(outcome$key, "annotation join failed"))
  snp_audit_rows[[length(snp_audit_rows) + 1L]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    SNP = primary$SNP,
    Frozen_Class = primary_annotation$Frozen_Class,
    N_Catalog_Associations = num_or_zero(primary_annotation$n_catalog_assoc),
    N_Adiposity_Associations = num_or_zero(primary_annotation$n_adiposity),
    N_Education_SES_Associations = num_or_zero(primary_annotation$n_ses_edu),
    In_Exclude_Five_Analysis = primary$SNP %in% known_present,
    In_Uncatalogued_Only_Analysis = primary$SNP %in% uncatalogued_present,
    stringsAsFactors = FALSE
  )

  subset_rows[[length(subset_rows) + 1L]] <- estimate_row(
    primary, outcome, "full_direct_instruments",
    "All mr_keep instruments after outcome-proxy exclusion"
  )
  subset_rows[[length(subset_rows) + 1L]] <- estimate_row(
    exclude_known, outcome, "exclude_5_BMI_or_education",
    paste("Excluded:", paste(known_present, collapse = ";"))
  )
  subset_rows[[length(subset_rows) + 1L]] <- estimate_row(
    uncatalogued_only, outcome, "uncatalogued_only",
    "Only instruments with zero associations in the frozen GWAS Catalog summary"
  )
  subset_rows[[length(subset_rows) + 1L]] <- estimate_row(
    known_only, outcome, "known_link_only",
    "Only the five instruments with catalogued BMI or education/SES links"
  )

  full_fit <- manual_ivw(primary)
  observed_fit <- manual_ivw(exclude_known)
  combinations <- combn(primary$SNP, length(known_present), simplify = FALSE)
  permutation_beta <- numeric(length(combinations))

  for (i in seq_along(combinations)) {
    dropped <- combinations[[i]]
    fit_i <- manual_ivw(primary[!primary$SNP %in% dropped, , drop = FALSE])
    permutation_beta[i] <- fit_i$b
    permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
      Outcome_Key = outcome$key,
      Outcome = outcome$label,
      Combination_Index = i,
      N_Dropped = length(dropped),
      Dropped_SNPs = paste(sort(dropped), collapse = ";"),
      Kept_N_SNPs = nrow(primary) - length(dropped),
      IVW_Beta = fit_i$b,
      stringsAsFactors = FALSE
    )
  }

  observed_pct_change <-
    100 * (abs(observed_fit$b) - abs(full_fit$b)) / abs(full_fit$b)
  q_perm <- unname(quantile(permutation_beta, c(0.025, 0.975), type = 7))
  permutation_summary_rows[[length(permutation_summary_rows) + 1L]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    N_Dropped = length(known_present),
    N_Combinations = length(combinations),
    Full_Beta = full_fit$b,
    Observed_Exclude5_Beta = observed_fit$b,
    Observed_AbsoluteEffect_Change_Percent = observed_pct_change,
    P_AtLeast_As_Attenuated =
      mean(abs(permutation_beta) <= abs(observed_fit$b)),
    Permutation_Beta_2.5pct = q_perm[1],
    Permutation_Beta_97.5pct = q_perm[2],
    stringsAsFactors = FALSE
  )

  known_fit <- manual_ivw(known_only)
  uncatalogued_fit <- manual_ivw(uncatalogued_only)
  z_difference <- (uncatalogued_fit$b - known_fit$b) /
    sqrt(uncatalogued_fit$se^2 + known_fit$se^2)
  subset_difference_rows[[length(subset_difference_rows) + 1L]] <- data.frame(
    Outcome_Key = outcome$key,
    Outcome = outcome$label,
    Known_Link_N_SNPs = nrow(known_only),
    Known_Link_Beta = known_fit$b,
    Known_Link_SE = known_fit$se,
    Uncatalogued_N_SNPs = nrow(uncatalogued_only),
    Uncatalogued_Beta = uncatalogued_fit$b,
    Uncatalogued_SE = uncatalogued_fit$se,
    Z_Difference = z_difference,
    P_Difference = 2 * pnorm(abs(z_difference), lower.tail = FALSE),
    Comparison_Note = "Exploratory comparison of disjoint SNP subsets",
    stringsAsFactors = FALSE
  )

  bx <- as.numeric(primary$beta.exposure)
  by <- as.numeric(primary$beta.outcome)
  w <- 1 / as.numeric(primary$se.outcome)^2
  numerator_components <- w * bx * by
  denominator_components <- w * bx^2
  total_numerator <- sum(numerator_components)
  total_denominator <- sum(denominator_components)

  for (i in seq_len(nrow(primary))) {
    loo <- primary[-i, , drop = FALSE]
    loo_fit <- manual_ivw(loo)
    influence_rows[[length(influence_rows) + 1L]] <- data.frame(
      Outcome_Key = outcome$key,
      Outcome = outcome$label,
      SNP = primary$SNP[i],
      Frozen_Class = primary_annotation$Frozen_Class[i],
      Numerator_Component = numerator_components[i],
      Numerator_Share = if (total_numerator == 0) NA_real_ else
        numerator_components[i] / total_numerator,
      Denominator_Component = denominator_components[i],
      Precision_Share = denominator_components[i] / total_denominator,
      Full_Beta = full_fit$b,
      Leave_One_Out_Beta = loo_fit$b,
      Leave_One_Out_SE = loo_fit$se,
      Delta_Beta = loo_fit$b - full_fit$b,
      AbsoluteEffect_Change_Percent =
        100 * (abs(loo_fit$b) - abs(full_fit$b)) / abs(full_fit$b),
      stringsAsFactors = FALSE
    )
  }
}

manifest <- do.call(rbind, manifest_rows)
snp_audit <- do.call(rbind, snp_audit_rows)
subset_estimates <- do.call(rbind, subset_rows)
permutation_distribution <- do.call(rbind, permutation_rows)
permutation_summary <- do.call(rbind, permutation_summary_rows)
subset_difference <- do.call(rbind, subset_difference_rows)
snp_influence <- do.call(rbind, influence_rows)

annotation_manifest <- data.frame(
  File = annotation_file,
  SHA256 = digest(annotation_file, algo = "sha256", file = TRUE),
  Annotation_Date = annotation_date,
  Reported_Sources = "IEU OpenGWAS; Ensembl REST; GWAS Catalog",
  Reproducibility_Status = "Annotation summary fixed on the annotation date; downstream analysis is rerunnable",
  Paper_Instrument_Count = nrow(paper_annotation),
  BMI_or_Education_Link_Count = length(known_link_snps),
  No_Catalogue_Record_Count = length(uncatalogued_snps),
  stringsAsFactors = FALSE
)

output_files <- c(
  CatalogAudit = paste0(prefix, "_CatalogAudit.csv"),
  AnnotationManifest = paste0(prefix, "_AnnotationManifest.csv"),
  InputManifest = paste0(prefix, "_InputManifest.csv"),
  SNPAudit = paste0(prefix, "_SNPAudit.csv"),
  SubsetEstimates = paste0(prefix, "_SubsetEstimates.csv"),
  PermutationDistribution = paste0(prefix, "_PermutationDistribution.csv"),
  PermutationSummary = paste0(prefix, "_PermutationSummary.csv"),
  SubsetDifference = paste0(prefix, "_SubsetDifference.csv"),
  SNPInfluence = paste0(prefix, "_SNPInfluence.csv")
)

write.csv(paper_annotation, output_files[["CatalogAudit"]], row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(annotation_manifest, output_files[["AnnotationManifest"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(manifest, output_files[["InputManifest"]], row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(snp_audit, output_files[["SNPAudit"]], row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(subset_estimates, output_files[["SubsetEstimates"]], row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(permutation_distribution,
          output_files[["PermutationDistribution"]], row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(permutation_summary, output_files[["PermutationSummary"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(subset_difference, output_files[["SubsetDifference"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(snp_influence, output_files[["SNPInfluence"]], row.names = FALSE,
          fileEncoding = "UTF-8")

output_hashes <- data.frame(
  Output = names(output_files),
  File = unname(output_files),
  SHA256 = vapply(
    unname(output_files), digest, character(1), algo = "sha256", file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(output_hashes, paste0(prefix, "_OutputHashes.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

run_log <- c(
  "Targeted pleiotropy sensitivity analysis",
  paste("Analysis date:", analysis_date),
  paste("Frozen annotation date:", annotation_date),
  "Network/API calls: none",
  paste("Annotation file:", annotation_file),
  paste("Annotation SHA256:", annotation_manifest$SHA256),
  paste("Paper instruments:", nrow(paper_annotation)),
  paste("BMI/education-linked instruments:", paste(known_link_snps, collapse = ";")),
  paste("Uncatalogued instruments:", paste(uncatalogued_snps, collapse = ";")),
  "Proxy rule: all detected outcome proxy rows excluded",
  "Estimator: manual multiplicative random-effects IVW; dispersion bounded below at 1",
  "Permutation: exact enumeration of every same-size five-SNP exclusion set",
  "P_AtLeast_As_Attenuated: proportion with abs(permuted beta) <= abs(observed exclude-five beta)",
  "External annotation: fixed summary file inputs/instrument_annotation.csv",
  paste("Outcomes completed:", length(outcomes)),
  paste("Exact permutation rows:", nrow(permutation_distribution))
)
writeLines(run_log, paste0(prefix, "_RunLog.txt"), useBytes = TRUE)

message("Targeted pleiotropy sensitivity analysis completed.")
message("Outputs written with prefix: ", prefix)
