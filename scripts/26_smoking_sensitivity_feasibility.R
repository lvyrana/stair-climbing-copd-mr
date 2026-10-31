# ============================================================================
# Smoking sensitivity MVMR: Stage A feasibility and input freeze
#
# Purpose:
#   Freeze and audit four prespecified exposure models crossing smoking trait
#   (cigarettes per day vs initiation) with sample-overlap context
#   (high-power GSCAN vs pre-UK Biobank TAG).
#
# Scope:
#   This stage does not request outcome data or estimate causal effects.
#   It writes only to the existing reanalysis workspace.
# ============================================================================

local_lib <- file.path(getwd(), "r_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
})

base_dir <- normalizePath(getwd(), mustWork = TRUE)
workspace_dir <- base_dir
input_dir <- file.path(workspace_dir, "inputs")
result_dir <- file.path(workspace_dir, "results")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

refresh <- "--refresh" %in% commandArgs(trailingOnly = TRUE)
legacy_path <- file.path(
  base_dir, "results", "mvmr_european_bmi_verification", "mv_exposure_data.csv"
)

dataset_roles <- data.frame(
  ID = c("ukb-b-298", "ieu-a-835", "ieu-b-25", "ieu-a-961", "ieu-b-4877", "ieu-a-962"),
  Role = c(
    "Stair-climbing frequency",
    "BMI",
    "Cigarettes per day",
    "Cigarettes smoked per day",
    "Smoking initiation",
    "Ever vs never smoked"
  ),
  Source_Context = c(
    "UK Biobank",
    "GIANT 2015; pre-UK Biobank",
    "GSCAN 2019; includes UK Biobank",
    "TAG 2010; pre-UK Biobank",
    "GSCAN 2019; includes UK Biobank",
    "TAG 2010; pre-UK Biobank"
  ),
  stringsAsFactors = FALSE
)

models <- data.frame(
  Model = c("GSCAN_CPD", "TAG_CPD", "GSCAN_INIT", "TAG_INIT"),
  Smoking_ID = c("ieu-b-25", "ieu-a-961", "ieu-b-4877", "ieu-a-962"),
  Smoking_Dimension = c("Intensity", "Intensity", "Initiation", "Initiation"),
  Overlap_Context = c(
    "Partial exposure-exposure overlap expected: GSCAN includes UK Biobank",
    "Overlap-minimised: TAG predates UK Biobank",
    "Partial exposure-exposure overlap expected: GSCAN includes UK Biobank",
    "Overlap-minimised: TAG predates UK Biobank"
  ),
  stringsAsFactors = FALSE
)

fixed_exposures <- c(stair = "ukb-b-298", bmi = "ieu-a-835")

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
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(NA_character_)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

validate_model_data <- function(dat, expected_ids, model_name) {
  required <- c(
    "SNP", "id.exposure", "effect_allele.exposure", "other_allele.exposure",
    "beta.exposure", "se.exposure", "pval.exposure"
  )
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols)) {
    stop(model_name, " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  found_ids <- sort(unique(dat$id.exposure))
  if (!identical(found_ids, sort(expected_ids))) {
    stop(
      model_name, ": expected exposure IDs ", paste(sort(expected_ids), collapse = ", "),
      "; found ", paste(found_ids, collapse = ", ")
    )
  }

  per_id <- table(dat$id.exposure)
  n_snp <- length(unique(dat$SNP))
  if (any(per_id != n_snp)) {
    stop(
      model_name, ": exposure rows are incomplete. Unique SNPs = ", n_snp,
      "; rows by exposure = ", paste(names(per_id), per_id, collapse = "; ")
    )
  }
  if (anyDuplicated(dat[c("SNP", "id.exposure")])) {
    stop(model_name, " has duplicated SNP-exposure rows.")
  }
  if (any(!is.finite(dat$beta.exposure)) || any(!is.finite(dat$se.exposure))) {
    stop(model_name, " has non-finite beta or SE values.")
  }
  invisible(TRUE)
}

metadata_path <- file.path(input_dir, "SmokingSensitivity_GWASMetadata.csv")
source_log <- character(0)

cache_paths <- file.path(
  input_dir, paste0("SmokingSensitivity_ExposureData_", models$Model, ".csv")
)
if (refresh || !file.exists(metadata_path) || !all(file.exists(cache_paths))) {
  load_opengwas_token()
}

if (!refresh && file.exists(metadata_path)) {
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  source_log <- c(source_log, "GWAS metadata: frozen cache")
} else {
  info <- ieugwasr::gwasinfo(dataset_roles$ID)
  info <- as.data.frame(info, stringsAsFactors = FALSE)
  id_col <- intersect(c("id", "ID"), names(info))[1]
  if (is.na(id_col)) stop("OpenGWAS metadata response did not contain an ID column.")
  names(info)[names(info) == id_col] <- "ID"
  keep <- intersect(
    c(
      "ID", "trait", "population", "sample_size", "ncase", "ncontrol",
      "author", "consortium", "year", "pmid", "unit", "build", "nsnp", "note"
    ),
    names(info)
  )
  metadata <- merge(dataset_roles, info[, keep, drop = FALSE], by = "ID", all.x = TRUE)
  if (any(is.na(metadata$trait))) {
    stop("Metadata missing for: ", paste(metadata$ID[is.na(metadata$trait)], collapse = ", "))
  }
  write.csv(metadata, metadata_path, row.names = FALSE, na = "")
  source_log <- c(source_log, "GWAS metadata: downloaded from OpenGWAS and frozen")
}

expected_patterns <- c(
  "ukb-b-298" = "stair",
  "ieu-a-835" = "body mass index|bmi",
  "ieu-b-25" = "cigarettes per day",
  "ieu-a-961" = "cigarettes smoked per day",
  "ieu-b-4877" = "smoking initiation",
  "ieu-a-962" = "ever vs never smoked"
)
for (id in names(expected_patterns)) {
  trait <- metadata$trait[metadata$ID == id]
  if (length(trait) != 1 || !grepl(expected_patterns[[id]], trait, ignore.case = TRUE)) {
    stop("Dataset identity check failed for ", id, ": ", paste(trait, collapse = " | "))
  }
}

feasibility_rows <- list()
membership_rows <- list()
input_paths <- c(metadata_path)

for (i in seq_len(nrow(models))) {
  model <- models$Model[i]
  smoking_id <- models$Smoking_ID[i]
  exposure_ids <- c(unname(fixed_exposures), smoking_id)
  cache_path <- file.path(
    input_dir, paste0("SmokingSensitivity_ExposureData_", model, ".csv")
  )

  if (!refresh && file.exists(cache_path)) {
    dat <- read.csv(cache_path, stringsAsFactors = FALSE, check.names = FALSE)
    data_source <- "frozen cache"
  } else if (!refresh && model == "GSCAN_CPD" && file.exists(legacy_path)) {
    dat <- read.csv(legacy_path, stringsAsFactors = FALSE, check.names = FALSE)
    data_source <- "verified legacy 89-SNP exposure cache copied into freeze"
    validate_model_data(dat, exposure_ids, model)
    write.csv(dat, cache_path, row.names = FALSE)
  } else {
    message("Extracting multivariable instruments for ", model, "...")
    dat <- mv_extract_exposures(
      exposure_ids,
      clump_r2 = 0.001,
      clump_kb = 10000,
      harmonise_strictness = 2,
      pval_threshold = 5e-8,
      pop = "EUR"
    )
    data_source <- "downloaded from OpenGWAS and frozen"
    validate_model_data(dat, exposure_ids, model)
    write.csv(dat, cache_path, row.names = FALSE)
  }

  validate_model_data(dat, exposure_ids, model)
  input_paths <- c(input_paths, cache_path)
  source_log <- c(source_log, paste0(model, ": ", data_source))

  wide_p <- reshape(
    dat[c("SNP", "id.exposure", "pval.exposure")],
    idvar = "SNP", timevar = "id.exposure", direction = "wide"
  )
  p_cols <- paste0("pval.exposure.", exposure_ids)
  missing_p_cols <- setdiff(p_cols, names(wide_p))
  if (length(missing_p_cols)) {
    stop(model, " missing wide P-value columns: ", paste(missing_p_cols, collapse = ", "))
  }
  pmat <- as.matrix(wide_p[, p_cols, drop = FALSE])
  colnames(pmat) <- exposure_ids
  min_idx <- max.col(-pmat, ties.method = "first")
  gws_matrix <- pmat < 5e-8
  gws_list <- apply(gws_matrix, 1, function(x) {
    ids <- colnames(gws_matrix)[x]
    if (length(ids)) paste(ids, collapse = ";") else "None"
  })

  membership_rows[[length(membership_rows) + 1]] <- data.frame(
    Model = model,
    SNP = wide_p$SNP,
    Minimum_P_Exposure = exposure_ids[min_idx],
    GenomeWideSignificant_Exposures = gws_list,
    P_Stair = pmat[, "ukb-b-298"],
    P_BMI = pmat[, "ieu-a-835"],
    P_Smoking = pmat[, smoking_id],
    stringsAsFactors = FALSE
  )

  feasibility_rows[[length(feasibility_rows) + 1]] <- data.frame(
    Model = model,
    Smoking_ID = smoking_id,
    Smoking_Dimension = models$Smoking_Dimension[i],
    Overlap_Context = models$Overlap_Context[i],
    N_Unique_SNPs = nrow(wide_p),
    N_Stair_GWS = sum(pmat[, "ukb-b-298"] < 5e-8),
    N_BMI_GWS = sum(pmat[, "ieu-a-835"] < 5e-8),
    N_Smoking_GWS = sum(pmat[, smoking_id] < 5e-8),
    Min_Smoking_P = min(pmat[, smoking_id], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

feasibility <- do.call(rbind, feasibility_rows)
membership <- do.call(rbind, membership_rows)

feasibility_path <- file.path(result_dir, "SmokingSensitivity_Feasibility.csv")
membership_path <- file.path(result_dir, "SmokingSensitivity_InstrumentMembership.csv")
manifest_path <- file.path(result_dir, "SmokingSensitivity_DataManifest.csv")
log_path <- file.path(result_dir, "SmokingSensitivity_Feasibility_RunLog.txt")

write.csv(feasibility, feasibility_path, row.names = FALSE)
write.csv(membership, membership_path, row.names = FALSE)

file_info <- file.info(input_paths)
manifest <- data.frame(
  File = basename(input_paths),
  Relative_Path = sub(paste0("^", base_dir, "/"), "", input_paths),
  Bytes = file_info$size,
  Modified = format(file_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  MD5 = unname(tools::md5sum(input_paths)),
  SHA256 = vapply(input_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write.csv(manifest, manifest_path, row.names = FALSE)

run_log <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  paste0("TwoSampleMR version: ", as.character(packageVersion("TwoSampleMR"))),
  paste0("ieugwasr version: ", as.character(packageVersion("ieugwasr"))),
  paste0("Refresh requested: ", refresh),
  "Instrument threshold: P < 5e-8",
  "LD clumping: r2 < 0.001; 10000 kb; EUR reference panel",
  "Harmonisation strictness: 2",
  "No outcome data requested and no causal estimates calculated",
  "",
  "Data sources:",
  paste0("- ", source_log),
  "",
  "Feasibility summary:",
  capture.output(print(feasibility, row.names = FALSE))
)
writeLines(run_log, log_path)

cat("Smoking sensitivity feasibility freeze completed.\n")
print(feasibility, row.names = FALSE)
cat("\nOutputs:\n")
cat("- ", feasibility_path, "\n", sep = "")
cat("- ", membership_path, "\n", sep = "")
cat("- ", manifest_path, "\n", sep = "")
cat("- ", log_path, "\n", sep = "")
