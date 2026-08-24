###############################################################################
## s.01.2 - Error correction from the manual sheet review                    ##
##                                                                           ##
## Applies the corrections documented during E.M.'s manual revision of the   ##
## laboratory sheets to the individual records produced by s.01.1. This      ##
## step MUST run before s.01.3-dry_weight_table.R: several patches touch     ##
## fw0, which feeds the calibration of the FW0 -> DW allometric model.       ##
##                                                                           ##
## Two kinds of actions:                                                     ##
##   patch : corrected_variable != "none" -> wrong value replaced by the     ##
##           corrected value, but ONLY after verifying that the raw record   ##
##           still holds the expected wrong value (guards against silent     ##
##           drift if the raw sheets change)                                 ##
##   drop  : corrected_variable == "none"   -> uncorrectable sheet error,    ##
##           the whole acorn leaves the pipeline                             ##
##                                                                           ##
## Every action is appended to error_correction_log, which the master        ##
## exports as an audit trail (00-data/error_correction_log.csv).             ##
###############################################################################

# --- 0. Guard ----------------------------------------------------------------

if (!exists("df.germ.records")) {
  stop("df.germ.records not found. Run/source s.01.1-load_raw.R first.")
}

ID_OFFSET_PHASE2 <- 1600  # must match s.01.1

# --- 1. Correction table -----------------------------------------------------
# Manual review of records flagged by the moisture audit: hydrated acorns
# (fw_sampling > fw0), impossibly high MC, negative MC and transcription
# offsets between neighbouring rows.
#
#   orig_acorn_id      : id as written on the phase sheets
#   phase              : 1 or 2
#   corrected_variable : "fw0" | "fw_sampling" | "none" ("none" = drop acorn)
#   wrong_value        : value currently in the raw sheets (NA for drops)
#   corrected_value    : reviewed value          (NA for drops)
corrections_raw <- tribble(
  ~orig_acorn_id, ~phase, ~corrected_variable, ~wrong_value, ~corrected_value, ~correction_reason,
             441,      1,               "fw0",         2.484,            2.848, "hydrated_acorn",
             487,      1,               "none",           NA,               NA, "hydrated_acorn",
             577,      1,               "none",           NA,               NA, "negative_value",
             713,      1,               "none",           NA,               NA, "hydrated_acorn",
            1122,      1,               "fw0",         3.350,            5.350, "hydrated_acorn",
             309,      2,               "fw0",         1.404,           11.404, "impossibly_high",
             337,      2,               "none",           NA,               NA, "negative_value",
             362,      2,               "none",           NA,               NA, "negative_value",
             466,      2,               "fw0",         9.175,            4.175, "negative_value",
             566,      2,               "none",           NA,               NA, "negative_value",
             582,      2,       "fw_sampling",         2.396,            1.568, "transcription_offset",
             583,      2,       "fw_sampling",         2.690,            2.396, "transcription_offset",
             584,      2,       "fw_sampling",         1.394,            1.883, "transcription_offset",
             585,      2,       "fw_sampling",         1.568,            1.394, "transcription_offset",
            1035,      2,               "none",           NA,               NA, "hydrated_acorn",
            1310,      2,               "fw0",         1.262,           11.262, "impossibly_high",
            1350,      2,       "fw_sampling",         9.331,            4.331, "hydrated_acorn",
            1397,      2,               "none",           NA,               NA, "negative_value",
            1530,      2,               "none",           NA,               NA, "negative_value",
            1544,      2,               "none",           NA,               NA, "negative_value",
            1561,      2,               "none",           NA,               NA, "negative_value"
)

corrections <- corrections_raw |>
  mutate(acorn_id = as.character(if_else(phase == 2,
                                         orig_acorn_id + ID_OFFSET_PHASE2,
                                         orig_acorn_id)))

# --- 2. Value patches --------------------------------------------------------

patch_rows <- corrections |> filter(corrected_variable != "none")
records_patched <- df.germ.records
patch_log <- vector("list", nrow(patch_rows))

for (i in seq_len(nrow(patch_rows))) {
  r <- patch_rows[i, ]
  at_t0 <- r$corrected_variable == "fw0"

  hits <- which(records_patched$acorn_id == r$acorn_id &
                (records_patched$time == "0") == at_t0 &
                abs(records_patched$fresh_weight - r$wrong_value) < 1e-6)

  if (length(hits) != 1) {
    stop(sprintf(paste0(
      "Correction mismatch: acorn %s (%s, %s), expected fresh_weight %s, ",
      "found %d matching records. Raw data changed? Update the table."),
      r$acorn_id, r$corrected_variable, r$correction_reason,
      format(r$wrong_value), length(hits)))
  }

  old_val <- records_patched$fresh_weight[hits]
  records_patched$fresh_weight[hits] <- r$corrected_value

  patch_log[[i]] <- tibble(
    action            = "patched",
    acorn_id          = r$acorn_id,
    phase             = r$phase,
    variable          = r$corrected_variable,
    old_value         = old_val,
    new_value         = r$corrected_value,
    correction_reason = r$correction_reason
  )
}
df.germ.records <- records_patched

# --- 3. Drop uncorrectable acorns ---------------------------------------------

drop_ids <- corrections |>
  filter(corrected_variable == "none") |>
  pull(acorn_id)

n_before <- n_distinct(df.germ.records$acorn_id)
df.germ.records <- df.germ.records |> filter(!acorn_id %in% drop_ids)

drop_log <- corrections |>
  filter(corrected_variable == "none") |>
  transmute(action = "dropped_acorn", acorn_id, phase,
            variable = NA_character_,
            old_value = NA_real_, new_value = NA_real_, correction_reason)

error_correction_log <- bind_rows(patch_log, drop_log)

# --- 4. Quality-control summary -----------------------------------------------

cat("=== s.01.2 error_correction ===\n")
cat("Patches applied :", nrow(patch_rows), "\n")
cat("Acorns dropped  :", length(drop_ids),
    sprintf("(%d -> %d acorns)\n", n_before, n_distinct(df.germ.records$acorn_id)))
cat("\nAudit trail:\n")
print(as.data.frame(error_correction_log), row.names = FALSE)
