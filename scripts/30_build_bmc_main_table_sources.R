# ============================================================================
# Build BMC Pulmonary Medicine main-table and Figure 1 source data
#
# Purpose
#   Create one internally consistent set of source-data files for the BMC
#   manuscript from the frozen no-proxy UVMR and verified targeted-pleiotropy
#   outputs. Historical Lung tables and figures are never overwritten.
#
# Guardrails
#   - No API/network calls.
#   - No manual transcription of estimates.
#   - COPD differential diagnosis remains on the linear-probability risk-
#     difference scale; FinnGen binary outcomes are ORs; lung function is beta.
#   - MVMR remains supplementary and is not included in the main tables/figure.
# ============================================================================

suppressPackageStartupMessages({
  library(digest)
})

options(stringsAsFactors = FALSE)

table_dir <- "tables"
figure_dir <- "tables"
result_dir <- "results"

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_csv <- function(path) {
  stop_if_not(file.exists(path), paste("Missing required input:", path))
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

old_table1_path <- file.path("inputs", "Table1_GWAS_Data_Sources.csv")
main_path <- file.path(result_dir, "Frozen_NoProxy_UVMR_Main.csv")
het_path <- file.path(result_dir, "Frozen_NoProxy_UVMR_Heterogeneity.csv")
pleio_path <- file.path(result_dir, "Frozen_NoProxy_UVMR_Pleiotropy.csv")
presso_path <- file.path(result_dir, "Frozen_NoProxy_UVMR_MRPRESSO.csv")
steiger_path <- file.path(result_dir, "Frozen_NoProxy_UVMR_Steiger.csv")
targeted_path <- file.path(result_dir, "TargetedPleiotropy_SubsetEstimates.csv")
permutation_path <- file.path(result_dir, "TargetedPleiotropy_PermutationSummary.csv")

old_table1 <- read_csv(old_table1_path)
main_all <- read_csv(main_path)
heterogeneity <- read_csv(het_path)
pleiotropy <- read_csv(pleio_path)
presso <- read_csv(presso_path)
steiger <- read_csv(steiger_path)
targeted <- read_csv(targeted_path)
permutation <- read_csv(permutation_path)

outcome_order <- c(
  "COPD_main", "COPD_hospital", "COPD_early", "COPD_later",
  "COPD_general", "COPD_emphysema", "COPD_insufficiency",
  "FEV1", "FVC", "PEF"
)
priority_order <- c("COPD_hospital", "COPD_general", "COPD_early", "COPD_main")

ivw <- main_all[main_all$method == "Inverse variance weighted", , drop = FALSE]
weighted_median <- main_all[main_all$method == "Weighted median", , drop = FALSE]
stop_if_not(nrow(ivw) == 10L, "Expected 10 frozen IVW outcome rows")
stop_if_not(all(outcome_order %in% ivw$Outcome_Key),
            "Frozen IVW outcomes differ from the expected set")

ivw$Outcome_Order <- match(ivw$Outcome_Key, outcome_order)
ivw <- ivw[order(ivw$Outcome_Order), , drop = FALSE]
weighted_median$Outcome_Order <- match(weighted_median$Outcome_Key, outcome_order)
weighted_median <- weighted_median[order(weighted_median$Outcome_Order), , drop = FALSE]

evidence_tier <- c(
  COPD_main = "Borderline and sensitivity-limited",
  COPD_hospital = "Primary supported finding",
  COPD_early = "Nominal finding",
  COPD_later = "Null",
  COPD_general = "Primary supported finding",
  COPD_emphysema = "Null; imprecise",
  COPD_insufficiency = "Null; imprecise",
  FEV1 = "Null",
  FVC = "Null",
  PEF = "Null; directional pleiotropy signal"
)

reporting_note <- c(
  COPD_main = "Report as risk difference, never as OR; IVW heterogeneity and MR-PRESSO outlier",
  COPD_hospital = "FinnGen direct-SNP estimate; one outcome proxy excluded",
  COPD_early = "Nominal only; Steiger classification TRUE but Steiger P approximately 0.994",
  COPD_later = "No statistical association",
  COPD_general = "FinnGen direct-SNP estimate; one outcome proxy excluded",
  COPD_emphysema = "Wide confidence interval; directionality classification not supportive",
  COPD_insufficiency = "Wide confidence interval; directionality classification not supportive",
  FEV1 = "No statistical association; strong heterogeneity",
  FVC = "No statistical association; strong heterogeneity",
  PEF = "No IVW association; strong heterogeneity and Egger-intercept signal"
)

# Table 1: exposure and UVMR outcomes only. Exploratory MVMR covariates remain
# in the supplementary data-source inventory.
table1_main <- old_table1[
  old_table1$Variable_Type %in% c("Exposure", "Outcome"), , drop = FALSE
]
table1_main$Primary_N_SNPs <- NA_integer_
table1_main$Proxies_Removed <- NA_integer_
table1_main$Analysis_Effect_Scale <- NA_character_
table1_main$Analysis_Position <- ifelse(
  table1_main$Variable_Type == "Exposure", "Exposure", "Main UVMR outcome"
)

for (i in seq_len(nrow(ivw))) {
  row_index <- which(table1_main$GWAS_ID == ivw$Outcome_ID[i])
  stop_if_not(length(row_index) == 1L,
              paste("Table 1 metadata match failed for", ivw$Outcome_ID[i]))
  table1_main$Primary_N_SNPs[row_index] <- ivw$Primary_N_SNPs[i]
  table1_main$Proxies_Removed[row_index] <- ivw$Proxies_Removed[i]
  table1_main$Analysis_Effect_Scale[row_index] <- ivw$Effect_Measure[i]
}
table1_main$Analysis_Effect_Scale[table1_main$Variable_Type == "Exposure"] <-
  "Per 1-SD higher genetically proxied stair-climbing frequency"

table1_supplement <- old_table1[
  old_table1$Variable_Type == "MVMR covariate", , drop = FALSE
]
table1_supplement$Analysis_Position <-
  "Exploratory supplementary MVMR covariate"

# Table 2: all ten primary no-proxy IVW estimates, with raw numeric fields for
# reproducible rendering in Word and figures.
table2 <- data.frame(
  Outcome_Order = ivw$Outcome_Order,
  Outcome_Key = ivw$Outcome_Key,
  Outcome = ivw$Outcome,
  Outcome_ID = ivw$Outcome_ID,
  Outcome_Type = ivw$Outcome_Type,
  Effect_Measure = ivw$Effect_Measure,
  Primary_N_SNPs = ivw$Primary_N_SNPs,
  Proxies_Removed = ivw$Proxies_Removed,
  Beta = ivw$b,
  SE = ivw$se,
  Beta_CI_Lower = ivw$CI_Lower,
  Beta_CI_Upper = ivw$CI_Upper,
  Effect_Estimate = ivw$Effect_Estimate,
  Effect_CI_Lower = ivw$Effect_CI_Lower,
  Effect_CI_Upper = ivw$Effect_CI_Upper,
  P = ivw$pval,
  Evidence_Tier = unname(evidence_tier[ivw$Outcome_Key]),
  Reporting_Note = unname(reporting_note[ivw$Outcome_Key]),
  stringsAsFactors = FALSE
)
stop_if_not(!any(is.na(table2$Evidence_Tier)), "Missing evidence-tier mapping")

# Full sensitivity joins used for both the compact main Table 3 and the
# supplementary all-outcome sensitivity table.
ivw_het <- heterogeneity[
  heterogeneity$method == "Inverse variance weighted",
  c("Outcome_Key", "Q", "Q_df", "Q_pval"), drop = FALSE
]
names(ivw_het)[2:4] <- c("IVW_Q", "IVW_Q_df", "IVW_Q_P")

wm <- weighted_median[
  , c("Outcome_Key", "b", "se", "pval", "Effect_Measure",
      "Effect_Estimate", "Effect_CI_Lower", "Effect_CI_Upper"), drop = FALSE
]
names(wm)[2:8] <- c(
  "Weighted_Median_Beta", "Weighted_Median_SE", "Weighted_Median_P",
  "Weighted_Median_Effect_Measure", "Weighted_Median_Effect",
  "Weighted_Median_CI_Lower", "Weighted_Median_CI_Upper"
)

egger <- pleiotropy[, c("Outcome_Key", "egger_intercept", "se", "pval"),
                    drop = FALSE]
names(egger)[2:4] <- c("Egger_Intercept", "Egger_Intercept_SE",
                       "Egger_Intercept_P")

presso_keep <- presso[
  , c("Outcome_Key", "Global_P_Text", "Outlier_Count", "Outlier_SNPs",
      "Corrected_Beta", "Corrected_SE", "Corrected_P"), drop = FALSE
]
names(presso_keep)[2:7] <- c(
  "MRPRESSO_Global_P", "MRPRESSO_Outlier_Count", "MRPRESSO_Outlier_SNPs",
  "MRPRESSO_Corrected_Beta", "MRPRESSO_Corrected_SE",
  "MRPRESSO_Corrected_P"
)

steiger_keep <- steiger[
  , c("Outcome_Key", "Correct_Causal_Direction", "Steiger_P",
      "Sensitivity_Ratio", "Steiger_R_Outcome_Method"), drop = FALSE
]

exclude5 <- targeted[
  targeted$Set == "exclude_5_BMI_or_education",
  c("Outcome_Key", "N_SNPs", "Beta", "SE", "OR", "OR_CI_Lower",
    "OR_CI_Upper", "P"), drop = FALSE
]
names(exclude5)[2:8] <- c(
  "Exclude5_N_SNPs", "Exclude5_Beta", "Exclude5_SE", "Exclude5_OR",
  "Exclude5_OR_CI_Lower", "Exclude5_OR_CI_Upper", "Exclude5_P"
)

permutation_keep <- permutation[
  , c("Outcome_Key", "Observed_AbsoluteEffect_Change_Percent",
      "P_AtLeast_As_Attenuated"), drop = FALSE
]
names(permutation_keep)[2:3] <- c(
  "Exclude5_AbsoluteEffect_Change_Percent", "Exclude5_ExactPermutation_P"
)

table3_all <- table2[
  , c("Outcome_Order", "Outcome_Key", "Outcome", "Outcome_ID",
      "Effect_Measure", "Evidence_Tier"), drop = FALSE
]
for (piece in list(wm, ivw_het, egger, presso_keep, steiger_keep,
                   exclude5, permutation_keep)) {
  table3_all <- merge(table3_all, piece, by = "Outcome_Key", all.x = TRUE,
                      sort = FALSE)
}
table3_all <- table3_all[order(table3_all$Outcome_Order), , drop = FALSE]
stop_if_not(nrow(table3_all) == 10L, "Sensitivity join did not retain 10 outcomes")

table3_main <- table3_all[
  table3_all$Outcome_Key %in% priority_order, , drop = FALSE
]
table3_main$Priority_Order <- match(table3_main$Outcome_Key, priority_order)
table3_main <- table3_main[order(table3_main$Priority_Order), , drop = FALSE]

# Figure 1 source: one source file, three scale-specific panels. Rendering must
# never put the UKB risk difference or continuous beta on an OR axis.
panel <- ifelse(
  table2$Effect_Measure == "OR", "B_FinnGen_OR",
  ifelse(table2$Effect_Measure == "Risk difference",
         "A_UKB_risk_difference", "C_lung_function_beta")
)
reference_value <- ifelse(table2$Effect_Measure == "OR", 1, 0)
figure1 <- data.frame(
  Panel = panel,
  Outcome_Order = table2$Outcome_Order,
  Outcome_Key = table2$Outcome_Key,
  Outcome = table2$Outcome,
  Outcome_ID = table2$Outcome_ID,
  Effect_Measure = table2$Effect_Measure,
  N_SNPs = table2$Primary_N_SNPs,
  Estimate = table2$Effect_Estimate,
  CI_Lower = table2$Effect_CI_Lower,
  CI_Upper = table2$Effect_CI_Upper,
  P = table2$P,
  Reference_Value = reference_value,
  Evidence_Tier = table2$Evidence_Tier,
  stringsAsFactors = FALSE
)

outputs <- c(
  BMC_Table1_GWAS_Sources_source =
    file.path(table_dir, "BMC_Table1_GWAS_Sources_source.csv"),
  BMC_Supplement_MVMR_Data_Sources_source =
    file.path(table_dir, "BMC_Supplement_MVMR_Data_Sources_source.csv"),
  BMC_Table2_Primary_UVMR_source =
    file.path(table_dir, "BMC_Table2_Primary_UVMR_source.csv"),
  BMC_Table3_Key_Sensitivity_source =
    file.path(table_dir, "BMC_Table3_Key_Sensitivity_source.csv"),
  BMC_Supplement_All_Sensitivity_source =
    file.path(table_dir, "BMC_Supplement_All_Sensitivity_source.csv"),
  BMC_Figure1_Main_MR_source =
    file.path(figure_dir, "BMC_Figure1_Main_MR_source.csv")
)

write.csv(table1_main, outputs[["BMC_Table1_GWAS_Sources_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table1_supplement,
          outputs[["BMC_Supplement_MVMR_Data_Sources_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table2, outputs[["BMC_Table2_Primary_UVMR_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table3_main, outputs[["BMC_Table3_Key_Sensitivity_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table3_all, outputs[["BMC_Supplement_All_Sensitivity_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(figure1, outputs[["BMC_Figure1_Main_MR_source"]],
          row.names = FALSE, fileEncoding = "UTF-8")

input_files <- c(
  old_table1_path, main_path, het_path, pleio_path, presso_path, steiger_path,
  targeted_path, permutation_path
)
manifest <- rbind(
  data.frame(
    Type = "Input",
    Object = basename(input_files),
    File = input_files,
    Rows = vapply(input_files, function(x) nrow(read_csv(x)), integer(1)),
    SHA256 = vapply(input_files, digest, character(1), algo = "sha256",
                    file = TRUE),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Type = "Output",
    Object = names(outputs),
    File = unname(outputs),
    Rows = c(nrow(table1_main), nrow(table1_supplement), nrow(table2),
             nrow(table3_main), nrow(table3_all), nrow(figure1)),
    SHA256 = vapply(unname(outputs), digest, character(1), algo = "sha256",
                    file = TRUE),
    stringsAsFactors = FALSE
  )
)
manifest_path <- file.path(table_dir, "BMC_Table_Figure_Source_Manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

run_log <- c(
  "BMC main table and Figure 1 source-data build",
  "Network/API calls: none",
  "Historical Lung tables/figures overwritten: no",
  "Main Table 1 rows: 11 (exposure plus 10 UVMR outcomes)",
  "Main Table 2 rows: 10 frozen no-proxy IVW outcomes",
  "Main Table 3 rows: 4 prioritized/borderline COPD outcomes",
  "Supplementary all-sensitivity rows: 10 outcomes",
  "Figure 1 source rows: 10 outcomes in three effect-scale panels",
  "MVMR position: supplementary only",
  "Effect-scale gate: UKB risk difference; FinnGen OR; lung-function beta"
)
writeLines(run_log, file.path(table_dir, "BMC_Table_Figure_Source_RunLog.txt"),
           useBytes = TRUE)

message("BMC table and Figure 1 source data completed.")
