###############################################################################
## s.01 - Sensitivity experiment: data pipeline master                       ##
##                                                                           ##
## Sources the three child scripts in order and exports the analysis-ready   ##
## per-acorn table consumed by the MC50 / desiccation-sensitivity analyses.  ##
##                                                                           ##
##   s.01.1-load_raw.R          raw sheets -> individual germination records ##
##   s.01.2-dry_weight_table.R  one row per acorn + FW0->DW allometric model ##
##                              (also writes its own diagnostics artifacts   ##
##                               to 07-img/dw_model_diagnostics and          ##
##                               00-data/tablas_resumen)                     ##
##   s.01.3-germination_table.R moisture content + outcomes -> df.analysis   ##
##                                                                           ##
## Output:                                                                   ##
##   00-data/sensitivity_germination_long.csv  (one row per acorn)           ##
###############################################################################

source("01-scripts/s.01.1-load_raw.R")
source("01-scripts/s.01.2-dry_weight_table.R")
source("01-scripts/s.01.3-germination_table.R")

write.csv(df.analysis,
          "00-data/sensitivity_germination_long.csv",
          row.names = FALSE)

cat("\ns.01 master done:",
    nrow(df.analysis), "acorns exported to",
    "00-data/sensitivity_germination_long.csv\n")

#---

ggplot(df.analysis, aes(x = mc, y = dw_source)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red") +
  geom_vline(aes(xintercept = 60), linetype = "dashed", colour = "red")

df.check <- df.analysis |> 
  mutate(
    flag_weird_mc = case_when(
      mc > 100 ~ "impossibly_high",
      fw_sampling > fw0 ~ "hidratated_acorn",
      mc > 65 ~ "very_high", 
      mc < 0 ~ "negative", 
      TRUE ~ "dont_check"
    )
  ) |> 
  filter(flag_weird_mc %in% c("impossibly_high", "hidratated_acorn", "very_high", "negative")) |> 
  select(acorn_id, fw0, fw_sampling, dw_final, dw_source, flag_weird_mc, phase, mc)

ggplot(df.check, aes(x = fw0, y = dw_final)) +
  geom_point(aes(colour = dw_source))

# Everything okay with this variables
df.check
as.numeric(quantile(df.check$dw_final / df.check$fw0, probs = c(.05, .95)))

df.check |>
  mutate(
    orig_acorn_id = if_else(
      phase == 1, 
      acorn_id, 
      acorn_id - 1600
    )
  ) |> select(orig_acorn_id, phase, mc, fw_sampling, fw0, dw_final, flag_weird_mc)

df.checked <- tribble(
  ~orig_acorn_id, ~phase, ~variable_corrected, ~value_corrected,   ~reason_checking, 
             441,      1,               "fw0",            2.848, "hidratated_acorn",
             487,      1,              "none",               NA, "hidratated_acorn",
  
  
)
