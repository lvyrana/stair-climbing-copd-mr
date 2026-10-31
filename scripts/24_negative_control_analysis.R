# ============================================================================
# Negative-control outcome analysis for stair-climbing instruments
#
# Primary diagnostic outcomes
#   1. ukb-b-533: Ease of skin tanning (ordered trait; SD scale)
#   2. ukb-d-1747_1: Natural hair colour, blonde (binary phenotype analysed
#      by the Neale Lab with linear regression; linear-probability scale)
#
# Pre-specified rules
#   - Exposure instruments are the 16 direct SNPs frozen in the existing
#     ukb-d-COPD_EXCL harmonised input.
#   - Outcome proxies are disabled. Missing direct SNPs remain missing.
#   - The full directly matched instrument set is primary.
#   - Bonferroni threshold: 0.05 / 2 = 0.025.
#   - Excluding five catalogued BMI/education variants and the uncatalogued
#     subset are secondary mechanism diagnostics, not additional primary tests.
#   - Ordinary IVW is diagnostic only because exposure and both outcomes are
#     UK Biobank-derived and individual-level overlap is unknown.
#   - Leave-one-out and SNP contribution diagnostics are emitted only if a
#     full-set IVW result has P < 0.025.
#
# Reproducibility
#   - On first run, outcome associations and metadata are downloaded once and
#     frozen under inputs/. The JWT is never written to outputs.
#   - Later runs use the frozen files and make no API calls.
#   - Use --refresh only after explicit approval to replace the frozen download.
# ============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(digest)
})

options(stringsAsFactors = FALSE)

source_dir <- "inputs"
result_dir <- "results"
script_path <- "scripts/24_negative_control_analysis.R"
exposure_path <- file.path(source_dir, "harmonised_steiger_COPD_main.csv")
outcome_cache <- file.path(source_dir, "NegativeControl_OutcomeData.csv")
metadata_cache <- file.path(source_dir, "NegativeControl_Metadata.csv")
prefix <- file.path(result_dir, "NegativeControl")
refresh <- "--refresh" %in% commandArgs(trailingOnly = TRUE)
alpha_primary <- 0.025
base_seed <- 20260809L

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  y <- toupper(trimws(as.character(x)))
  !is.na(y) & y %in% c("TRUE", "T", "1", "YES")
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA else x[1]
}

read_jwt <- function() {
  candidates <- c(
    "data/opengwas_token",
    path.expand("~/.config/opengwas_token")
  )
  existing <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  stop_if_not(length(existing) > 0L,
              "No non-empty OpenGWAS token file was found.")
  jwt <- trimws(paste(readLines(existing[1], warn = FALSE), collapse = ""))
  stop_if_not(nzchar(jwt), "OpenGWAS token file is empty after trimming.")
  jwt
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
  w <- 1 / (dat$se.outcome^2)
  denominator <- sum(w * bx^2)
  b <- sum(w * bx * by) / denominator
  q <- sum(w * (by - b * bx)^2)
  q_df <- nrow(dat) - 1L
  dispersion <- max(1, q / q_df)
  se <- sqrt(dispersion / denominator)
  data.frame(
    Manual_Beta = b,
    Manual_SE = se,
    Manual_Q = q,
    Manual_Q_df = q_df,
    Manual_Dispersion = dispersion,
    stringsAsFactors = FALSE
  )
}

format_effect <- function(dat, beta_scale) {
  dat$Beta_Scale <- beta_scale
  dat$Effect_Measure <- if (beta_scale == "linear_probability") {
    "Risk difference"
  } else {
    "Beta"
  }
  dat$Estimate <- dat$b
  dat$CI_Lower <- dat$b - 1.96 * dat$se
  dat$CI_Upper <- dat$b + 1.96 * dat$se
  dat
}

outcomes <- data.frame(
  Outcome_Key = c("tanning", "blonde_hair"),
  Outcome = c("Ease of skin tanning", "Natural hair colour: blonde"),
  Outcome_ID = c("ukb-b-533", "ukb-d-1747_1"),
  Outcome_Type = c("Categorical ordered", "Binary"),
  Beta_Scale = c("SD", "linear_probability"),
  Direction = c(
    "Higher values mean less ability to tan (towards never tan, only burn)",
    "Positive beta means a higher probability of blonde hair"
  ),
  Data_Source = c("MRC-IEU UK Biobank PHESANT", "Neale Lab UK Biobank Round 2"),
  stringsAsFactors = FALSE
)

catalogued_confounders <- c(
  "rs12140153", "rs76076331", "rs6449503", "rs7651039", "rs10835211"
)

uncatalogued <- c(
  "rs3898317", "rs72820477", "rs58133129", "rs7577014",
  "rs11689201", "rs1362487", "rs4855841", "rs1829967", "rs7851239"
)

stop_if_not(file.exists(exposure_path),
            paste("Missing frozen exposure source:", exposure_path))

frozen <- read.csv(exposure_path, check.names = FALSE)
required_exposure <- c(
  "SNP", "effect_allele.exposure", "other_allele.exposure",
  "beta.exposure", "se.exposure", "eaf.exposure", "pval.exposure",
  "samplesize.exposure", "id.exposure", "exposure"
)
stop_if_not(all(required_exposure %in% names(frozen)),
            "Frozen exposure source is missing required columns.")

exposure_dat <- unique(frozen[, required_exposure])
stop_if_not(nrow(exposure_dat) == 16L,
            paste("Expected 16 frozen exposure SNPs; found", nrow(exposure_dat)))
stop_if_not(length(unique(exposure_dat$id.exposure)) == 1L &&
              unique(exposure_dat$id.exposure) == "ukb-b-298",
            "Frozen exposure ID is not uniquely ukb-b-298.")
exposure_dat$mr_keep.exposure <- TRUE
exposure_dat$data_source.exposure <- "Frozen harmonised COPD_main input"
exposure_dat$pval_origin.exposure <- "reported"

cache_was_missing <- !file.exists(outcome_cache) || !file.exists(metadata_cache)
cache_mode <- if (refresh) {
  "refreshed_from_OpenGWAS"
} else if (cache_was_missing) {
  "created_from_OpenGWAS"
} else {
  "reused_frozen_cache"
}

if (refresh || cache_was_missing) {
  jwt <- read_jwt()
  metadata <- ieugwasr::gwasinfo(outcomes$Outcome_ID, opengwas_jwt = jwt)
  outcome_dat <- TwoSampleMR::extract_outcome_data(
    snps = exposure_dat$SNP,
    outcomes = outcomes$Outcome_ID,
    proxies = FALSE,
    opengwas_jwt = jwt
  )

  stop_if_not(is.data.frame(metadata) && nrow(metadata) == nrow(outcomes),
              "OpenGWAS metadata query did not return both outcomes.")
  stop_if_not(is.data.frame(outcome_dat) && nrow(outcome_dat) > 0L,
              "OpenGWAS outcome extraction returned no rows.")
  if ("proxy.outcome" %in% names(outcome_dat)) {
    stop_if_not(!any(as_flag(outcome_dat$proxy.outcome)),
                "Outcome extraction unexpectedly returned proxy rows.")
  }

  write.csv(metadata, metadata_cache, row.names = FALSE, na = "")
  write.csv(outcome_dat, outcome_cache, row.names = FALSE, na = "")
} else {
  metadata <- read.csv(metadata_cache, check.names = FALSE)
  outcome_dat <- read.csv(outcome_cache, check.names = FALSE)
}

stop_if_not(all(outcomes$Outcome_ID %in% unique(outcome_dat$id.outcome)),
            "Frozen outcome cache does not contain both pre-specified outcomes.")
if ("proxy.outcome" %in% names(outcome_dat)) {
  stop_if_not(!any(as_flag(outcome_dat$proxy.outcome)),
              "Frozen outcome cache contains proxy rows; refusing primary analysis.")
}

harmonised_list <- list()
audit_list <- list()
main_list <- list()
manual_list <- list()
heterogeneity_list <- list()
pleiotropy_list <- list()
contribution_list <- list()
loo_list <- list()
primary_flags <- list()

for (i in seq_len(nrow(outcomes))) {
  spec <- outcomes[i, ]
  raw_i <- outcome_dat[outcome_dat$id.outcome == spec$Outcome_ID, , drop = FALSE]
  direct_retrieved <- unique(raw_i$SNP)

  harmonised <- suppressMessages(TwoSampleMR::harmonise_data(
    exposure_dat,
    raw_i,
    action = 2
  ))
  remove_flag <- if ("remove" %in% names(harmonised)) {
    as_flag(harmonised$remove)
  } else {
    rep(FALSE, nrow(harmonised))
  }
  keep_flag <- if ("mr_keep" %in% names(harmonised)) {
    as_flag(harmonised$mr_keep)
  } else {
    !remove_flag
  }
  harmonised$Primary_Keep <- keep_flag & !remove_flag
  harmonised$Outcome_Key <- spec$Outcome_Key
  harmonised$Beta_Scale <- spec$Beta_Scale
  harmonised$Outcome_Direction <- spec$Direction
  harmonised_list[[spec$Outcome_Key]] <- harmonised

  audit_list[[spec$Outcome_Key]] <- data.frame(
    Outcome_Key = spec$Outcome_Key,
    Outcome = spec$Outcome,
    Outcome_ID = spec$Outcome_ID,
    Requested_SNP = exposure_dat$SNP,
    Direct_Retrieved = exposure_dat$SNP %in% direct_retrieved,
    Harmonised = exposure_dat$SNP %in% harmonised$SNP,
    Primary_Keep = exposure_dat$SNP %in%
      harmonised$SNP[harmonised$Primary_Keep],
    Catalogued_BMI_or_Education = exposure_dat$SNP %in% catalogued_confounders,
    GWAS_Catalog_Uncatalogued = exposure_dat$SNP %in% uncatalogued,
    Proxy_Allowed = FALSE,
    stringsAsFactors = FALSE
  )

  dat_primary <- harmonised[harmonised$Primary_Keep, , drop = FALSE]
  stop_if_not(nrow(dat_primary) >= 10L,
              paste(spec$Outcome_ID, "retained fewer than 10 direct SNPs."))

  instrument_sets <- list(
    full_direct = dat_primary,
    exclude_catalogued_confounders = dat_primary[
      !dat_primary$SNP %in% catalogued_confounders, , drop = FALSE
    ],
    uncatalogued_only = dat_primary[
      dat_primary$SNP %in% uncatalogued, , drop = FALSE
    ]
  )

  for (set_name in names(instrument_sets)) {
    dat_set <- instrument_sets[[set_name]]
    stop_if_not(nrow(dat_set) >= 3L,
                paste(spec$Outcome_ID, set_name, "has fewer than 3 SNPs."))

    set.seed(base_seed + i * 100L + match(set_name, names(instrument_sets)))
    mr_result <- suppressMessages(TwoSampleMR::mr(
      dat_set,
      method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
    ))
    mr_result <- format_effect(mr_result, spec$Beta_Scale)
    mr_result$Outcome_Key <- spec$Outcome_Key
    mr_result$Outcome_Label <- spec$Outcome
    mr_result$Instrument_Set <- set_name
    mr_result$Is_Primary_Test <- set_name == "full_direct" &
      mr_result$method == "Inverse variance weighted"
    mr_result$Primary_Alpha <- alpha_primary
    mr_result$Passes_Primary_Threshold <- mr_result$Is_Primary_Test &
      !is.na(mr_result$pval) & mr_result$pval < alpha_primary
    main_list[[paste(spec$Outcome_Key, set_name, sep = "__")]] <- mr_result

    manual <- manual_ivw(dat_set)
    ivw_row <- mr_result[mr_result$method == "Inverse variance weighted", , drop = FALSE]
    stop_if_not(nrow(ivw_row) == 1L,
                paste("IVW result missing for", spec$Outcome_ID, set_name))
    manual$Outcome_Key <- spec$Outcome_Key
    manual$Outcome_ID <- spec$Outcome_ID
    manual$Instrument_Set <- set_name
    manual$N_SNP <- nrow(dat_set)
    manual$TwoSampleMR_Beta <- ivw_row$b
    manual$TwoSampleMR_SE <- ivw_row$se
    manual$Beta_Absolute_Difference <- abs(manual$Manual_Beta - ivw_row$b)
    manual$SE_Absolute_Difference <- abs(manual$Manual_SE - ivw_row$se)
    manual_list[[paste(spec$Outcome_Key, set_name, sep = "__")]] <- manual
  }

  heterogeneity <- suppressMessages(TwoSampleMR::mr_heterogeneity(
    dat_primary,
    method_list = c("mr_ivw", "mr_egger_regression")
  ))
  heterogeneity$Outcome_Key <- spec$Outcome_Key
  heterogeneity$Instrument_Set <- "full_direct"
  heterogeneity_list[[spec$Outcome_Key]] <- heterogeneity

  pleiotropy <- suppressMessages(TwoSampleMR::mr_pleiotropy_test(dat_primary))
  pleiotropy$Outcome_Key <- spec$Outcome_Key
  pleiotropy$Instrument_Set <- "full_direct"
  pleiotropy_list[[spec$Outcome_Key]] <- pleiotropy

  full_ivw <- main_list[[paste(spec$Outcome_Key, "full_direct", sep = "__")]]
  full_ivw <- full_ivw[full_ivw$method == "Inverse variance weighted", , drop = FALSE]
  is_positive <- nrow(full_ivw) == 1L && !is.na(full_ivw$pval) &&
    full_ivw$pval < alpha_primary
  primary_flags[[spec$Outcome_Key]] <- data.frame(
    Outcome_Key = spec$Outcome_Key,
    Outcome_ID = spec$Outcome_ID,
    N_Direct_Retrieved = length(direct_retrieved),
    N_Harmonised_Keep = nrow(dat_primary),
    IVW_Beta = first_non_missing(full_ivw$b),
    IVW_SE = first_non_missing(full_ivw$se),
    IVW_P = first_non_missing(full_ivw$pval),
    Primary_Alpha = alpha_primary,
    Negative_Control_Alert = is_positive,
    Diagnostic_Only = TRUE,
    Overlap_Boundary = paste(
      "Exposure and outcome are UK Biobank-derived; individual-level overlap",
      "and exposure-outcome error covariance are unknown"
    ),
    stringsAsFactors = FALSE
  )

  if (is_positive) {
    bx <- dat_primary$beta.exposure
    by <- dat_primary$beta.outcome
    w <- 1 / (dat_primary$se.outcome^2)
    numerator_piece <- w * bx * by
    denominator_piece <- w * bx^2
    contribution <- data.frame(
      Outcome_Key = spec$Outcome_Key,
      Outcome_ID = spec$Outcome_ID,
      SNP = dat_primary$SNP,
      Wald_Ratio = by / bx,
      IVW_Numerator_Piece = numerator_piece,
      IVW_Denominator_Piece = denominator_piece,
      Absolute_Numerator_Share = abs(numerator_piece) / sum(abs(numerator_piece)),
      Denominator_Leverage = denominator_piece / sum(denominator_piece),
      Catalogued_BMI_or_Education = dat_primary$SNP %in% catalogued_confounders,
      GWAS_Catalog_Uncatalogued = dat_primary$SNP %in% uncatalogued,
      stringsAsFactors = FALSE
    )
    contribution <- contribution[
      order(-contribution$Absolute_Numerator_Share), , drop = FALSE
    ]
    contribution$Contribution_Rank <- seq_len(nrow(contribution))
    contribution_list[[spec$Outcome_Key]] <- contribution

    loo <- suppressMessages(TwoSampleMR::mr_leaveoneout(dat_primary))
    loo$Outcome_Key <- spec$Outcome_Key
    loo$Outcome_ID <- spec$Outcome_ID
    loo_list[[spec$Outcome_Key]] <- loo
  }
}

harmonised_all <- do.call(rbind, harmonised_list)
audit_all <- do.call(rbind, audit_list)
main_all <- do.call(rbind, main_list)
manual_all <- do.call(rbind, manual_list)
heterogeneity_all <- do.call(rbind, heterogeneity_list)
pleiotropy_all <- do.call(rbind, pleiotropy_list)
primary_all <- do.call(rbind, primary_flags)

stop_if_not(max(manual_all$Beta_Absolute_Difference, na.rm = TRUE) < 1e-10,
            "Manual and TwoSampleMR IVW beta differ by >= 1e-10.")
stop_if_not(max(manual_all$SE_Absolute_Difference, na.rm = TRUE) < 1e-10,
            "Manual and TwoSampleMR IVW SE differ by >= 1e-10.")

write.csv(harmonised_all, paste0(prefix, "_Harmonised.csv"), row.names = FALSE, na = "")
write.csv(audit_all, paste0(prefix, "_SNPAudit.csv"), row.names = FALSE, na = "")
write.csv(main_all, paste0(prefix, "_Main.csv"), row.names = FALSE, na = "")
write.csv(manual_all, paste0(prefix, "_ManualIVWCheck.csv"), row.names = FALSE, na = "")
write.csv(heterogeneity_all, paste0(prefix, "_Heterogeneity.csv"), row.names = FALSE, na = "")
write.csv(pleiotropy_all, paste0(prefix, "_Pleiotropy.csv"), row.names = FALSE, na = "")
write.csv(primary_all, paste0(prefix, "_PrimaryDecision.csv"), row.names = FALSE, na = "")

if (length(contribution_list) > 0L) {
  write.csv(do.call(rbind, contribution_list),
            paste0(prefix, "_SNPContributions.csv"), row.names = FALSE, na = "")
} else if (file.exists(paste0(prefix, "_SNPContributions.csv"))) {
  file.remove(paste0(prefix, "_SNPContributions.csv"))
}

if (length(loo_list) > 0L) {
  write.csv(do.call(rbind, loo_list),
            paste0(prefix, "_LeaveOneOut.csv"), row.names = FALSE, na = "")
} else if (file.exists(paste0(prefix, "_LeaveOneOut.csv"))) {
  file.remove(paste0(prefix, "_LeaveOneOut.csv"))
}

manifest <- rbind(
  file_manifest_row(exposure_path, "Exposure instrument source",
                    "Existing frozen harmonised input; exposure columns only"),
  file_manifest_row(outcome_cache, "Negative-control outcome associations",
                    "OpenGWAS direct SNP extraction; proxies disabled"),
  file_manifest_row(metadata_cache, "Negative-control outcome metadata",
                    "OpenGWAS metadata snapshot"),
  file_manifest_row(script_path, "Analysis script",
                    "Pre-specified negative-control workflow")
)
write.csv(manifest, paste0(prefix, "_DataManifest.csv"), row.names = FALSE, na = "")

log_lines <- c(
  "Negative-control outcome analysis",
  paste("Run time UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R version:", R.version.string),
  paste("TwoSampleMR version:", as.character(packageVersion("TwoSampleMR"))),
  paste("ieugwasr version:", as.character(packageVersion("ieugwasr"))),
  paste("Outcome cache mode:", cache_mode),
  "Outcome proxies: disabled",
  "Primary outcomes: ukb-b-533; ukb-d-1747_1",
  paste("Primary threshold:", alpha_primary),
  paste("Requested exposure SNPs:", nrow(exposure_dat)),
  "Blonde-hair effect scale: linear probability (Neale Lab linear regression), not odds ratio",
  paste(
    "Interpretation boundary:",
    "ordinary IVW is diagnostic only because all GWAS are UK Biobank-derived and",
    "individual-level overlap/error covariance are unknown"
  ),
  "",
  capture.output(print(primary_all, row.names = FALSE))
)
writeLines(log_lines, paste0(prefix, "_RunLog.txt"))

cat(paste(log_lines, collapse = "\n"), "\n")
