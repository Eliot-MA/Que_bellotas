###############################################################################
## s.01 - Sensitivity experiment: data pipeline master                       ##
##                                                                           ##
## Sources the four child scripts in order and exports the analysis-ready    ##
## per-acorn table consumed by the MC50 / desiccation-sensitivity analyses.  ##
##                                                                           ##
##   s.01.1-load_raw.R          raw sheets -> individual germination records ##
##   s.01.2-error_correction.R  manual-review patches + acorn drop list      ##
##   s.01.3-dry_weight_table.R  one row per acorn + FW0->DW allometric model ##
##                              (also writes its own diagnostics artifacts   ##
##                               to 07-img/dw_model_diagnostics and          ##
##                               00-data/tablas_resumen)                     ##
##   s.01.4-germination_table.R moisture content + outcomes -> df.analysis   ##
##                                                                           ##
## Outputs:                                                                  ##
##   00-data/sensitivity_germination_long.csv  (one row per acorn)           ##
##   00-data/error_correction_log.csv          (audit trail of s.01.2)       ##
###############################################################################

source("01-scripts/s.01.1-load_raw.R")
source("01-scripts/s.01.2-error_correction.R")
source("01-scripts/s.01.3-dry_weight_table.R")
source("01-scripts/s.01.4-germination_table.R")

write.csv(df.analysis,
          "00-data/sensitivity_germination_long.csv",
          row.names = FALSE)

write.csv(error_correction_log,
          "00-data/error_correction_log.csv",
          row.names = FALSE)

cat("\ns.01 master done:",
    nrow(df.analysis), "acorns exported to",
    "00-data/sensitivity_germination_long.csv;",
    nrow(error_correction_log), "correction actions logged\n")
