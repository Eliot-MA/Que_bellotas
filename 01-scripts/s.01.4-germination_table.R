###############################################################################
## s.01.4 - Analysis-ready germination table                                 ##
##                                                                           ##
## Joins the dry-weight results of s.01.3 with the germination outcomes and  ##
## computes the moisture content (MC) that will serve as predictor in the    ##
## sensitivity (MC50-type) analyses downstream.                              ##
##                                                                           ##
## Moisture content convention (decided with E.F.):                          ##
##   MC (%) = (FW - DW) / FW0 * 100                                          ##
##   - FW0 denominator everywhere                            ##
##   - numerator FW = fresh weight at the relevant measurement moment:       ##
##       * control acorns (t0 only)            -> FW0                        ##
##       * sampled acorns with fw_sampling     -> fw_sampling                ##
##   - DW = measured when available, otherwise predicted by the allometric   ##
##     model of s.01.2 (dry mass conserved during drying)                    ##
##                                                                           ##
## Output object: df.analysis, one row per acorn.                            ##
###############################################################################

# --- 0. Guard ---------------------------------------------------------------

if (!exists("df.acorns")) {
  stop("df.acorns not found. Run/source s.01.1-load_raw.R, ",
       "s.01.2-error_correction.R and s.01.3-dry_weight_table.R first.")
}

# --- 1. Build analysis table -------------------------------------------------

df.analysis <- df.acorns |>
  mutate(
    # Fresh weight entering the MC numerator.
    # Controls use FW0; sampled acorns need their weigh-in at removal time.
    mc_fresh_weight = if_else(is.na(sampling_time), fw0, fw_sampling),

    flag_no_fw_at_removal = !is.na(sampling_time) & is.na(fw_sampling),
    flag_no_fw0           = is.na(fw0),

    # Moisture content (% relative to initial fresh weight).
    mc = case_when(
      is.na(dw_final) | is.na(mc_fresh_weight) | mc_fresh_weight <= 0 ~ NA_real_,
      TRUE ~ (mc_fresh_weight - dw_final) / fw0 * 100
    ),

    # Suspect-moisture audit (manual review protocol, E.F.).
    # NOTE on "negative": when fw_sampling is close to the dry weight,
    # prediction error can push DW above FW, giving a small negative MC
    # near zero. Uncorrectable sheet errors of this kind are already
    # removed by s.01.2-error_correction; any remaining negative values
    # are prediction artifacts and stay flagged for the analysis stage.
    flag_weird_mc = case_when(
      is.na(mc)          ~ "mc_unavailable",
      mc > 100           ~ "impossibly_high",
      fw_sampling > fw0  ~ "hydrated_acorn",
      mc > 65            ~ "very_high",
      mc < 0             ~ "negative",
      TRUE               ~ "ok"
    )
  ) |>
  select(acorn_id, phase, species, provenance,
         sampling_time, fw0, fw_sampling,
         dw_observed, dw_model_pred, dw_final, dw_source,
         mc, mc_fresh_weight, flag_weird_mc,
         germinated, emerged, moho, observations,
         flag_no_fw_at_removal, flag_no_fw0)

# --- 2. Quality-control summary ----------------------------------------------

cat("=== s.01.4 germination_table ===\n")
cat("Acorns:", nrow(df.analysis), "\n")

cat("\nMoisture content availability:\n")
print(table(df.analysis$species,
            mc_available = !is.na(df.analysis$mc)) |> as.data.frame())

n_mc_na_no_dw <- sum(is.na(df.analysis$dw_final))
n_mc_na_no_fw <- sum(!is.na(df.analysis$dw_final) &
                     is.na(df.analysis$mc))
if (n_mc_na_no_fw > 0) {
  warning(n_mc_na_no_fw, " acorns have a dry weight but no fresh weight at ",
          "their removal time, so their moisture content cannot be computed ",
          "under the agreed convention.")
}
cat("\nMC missing because DW unavailable:", n_mc_na_no_dw, "\n")
cat("MC missing despite DW available:", n_mc_na_no_fw, "\n")

cat("\nGermination / emergence totals:\n")
print(with(df.analysis, table(germinated = factor(germinated,
                                                  levels = c(0, 1)),
                              emerged = factor(emerged,
                                               levels = c(0, 1)),
                              useNA = "ifany")))

cat("\nSuspect-moisture flags:\n")
print(table(df.analysis$species, df.analysis$flag_weird_mc))
