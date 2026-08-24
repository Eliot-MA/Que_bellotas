###############################################################################
## s.02 - Germination GLMs: moisture content x species (+ phase modelling)   ##
##                                                                           ##
## Reads the analysis-ready table exported by the s.01 master and fits       ##
## binomial GLMs of germination success as a function of moisture content    ##
## (MC, % relative to FW0), with species interacting with MC.                ##
##                                                                           ##
## Phase modelling (E.M. protocol):                                          ##
##   Phase (drying batch) is associated with MC by design: all t0 controls    ##
##   belong to phase 1, phase 2 sampled one extra removal time, and chamber   ##
##   conditions differed between runs. Screening showed a strong adjusted     ##
##   batch effect that is heterogeneous across species (phase:species         ##
##   supported, dAICc = -18.7 vs additive phase). The phase-inclusive global ##
##   therefore carries phase*species, not just an additive phase shift.      ##
##                                                                           ##
## Globals compared (E.M. two-model plan):                                   ##
##   M2:  germ ~ mc * species                                                ##
##   M1:  germ ~ mc * species + phase * species                              ##
## Both dredged independently; M2 nests inside M1 (LRT + AICc comparison).   ##
##                                                                           ##
## Q. robur exclusion (E.M.): the RO lot failed throughout the experiment,   ##
## with near-zero germination across the entire MC gradient (suspected       ##
## lot-level viability failure unrelated to desiccation sensitivity). A      ##
## justification figure is produced BEFORE the exclusion; all downstream     ##
## models use the remaining seven species.                                   ##
##                                                                           ##
## This script stops after assumption checks (DHARMa) and predictive figures. ##
## MC50 estimation is deferred to a dedicated script; the averaged-figure     ##
## MC50 guides are graphical references read off the fitted curves, not       ##
## formal estimates.                                                          ##
##                                                                           ##
## Outputs:                                                                  ##
##   00-data/tablas_resumen/glm_dredge_with_phase.csv                        ##
##   00-data/tablas_resumen/glm_dredge_no_phase.csv                          ##
##   00-data/tablas_resumen/glm_global_comparison.csv                        ##
##   00-data/tablas_resumen/glm_selected_coefficients.csv                    ##
##   00-data/tablas_resumen/glm_dharma_tests.csv                             ##
##   07-img/germ_glm_diagnostics/*.png                                       ##
##     - germ_ro_exclusion_justification.png (all 8 species, pre-exclusion)  ##
##     - glm_pred_curves_species_phase.png (per-batch curves + raw points)   ##
##     - glm_pred_curves_species_average.png (batch-averaged curves +        ##
##       graphical MC50 guides with CI)                                      ##
###############################################################################

library(tidyverse)
library(MuMIn)
library(DHARMa)
library(emmeans)

# --- 0. Configuration --------------------------------------------------------

EXCLUDE_WEIRD_MC <- TRUE    # drop impossible predictor values (negative-MC artifact)
EXCLUDE_SPECIES  <- c("RO") # lot-level viability failure; see justification figure
SIM_SEED         <- 20260823

tablas_dir <- "00-data/tablas_resumen"
plots_dir  <- "07-img/germ_glm_diagnostics"
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Data -----------------------------------------------------------------

df.germ <- read.csv("00-data/sensitivity_germination_long.csv") |>
  mutate(
    species  = factor(species),
    phase    = recode_factor(factor(phase), `1` = "phase_1", `2` = "phase_2")
  )

cat("=== s.02 germination GLMs ===\n")
cat("Rows read:", nrow(df.germ), "\n")

df.germ <- df.germ |>
  filter(!is.na(mc), !is.na(germinated))

if (EXCLUDE_WEIRD_MC) {
  df.germ <- df.germ |> filter(flag_weird_mc %in% c("ok", "very_high"))
}

cat("Rows used for modelling:", nrow(df.germ), "\n")

cat("\nBy phase:\n")
print(
  df.germ |>
    group_by(phase) |>
    summarise(n = n(), germ_rate = round(mean(germinated), 3),
              .groups = "drop")
)

# --- 1b. Q. robur exclusion: data-based justification -------------------------

n_pre_exclusion <- nrow(df.germ)

cat("\n=== Q. robur exclusion check ===\n")
cat("Germination rate at MC >= 25% (responsive species separate clearly there):\n")
print(
  df.germ |>
    filter(mc >= 25) |>
    group_by(species) |>
    summarise(n = n(), n_germinated = sum(germinated),
              germ_rate = round(mean(germinated), 3),
              .groups = "drop")
)

p_ro_check <- ggplot(df.germ, aes(mc, germinated)) +
  geom_point(alpha = 0.12, size = 0.6,
             position = position_jitter(height = 0.03, width = 0)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, linewidth = 0.8) +
  facet_wrap(~species) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  labs(x = "Moisture content (% FW0)",
       y = "Germination (0/1)",
       title = "Species-level germination response to moisture content",
       caption = paste0(
         "Raw germination outcomes (jittered) with binomial logistic fits per ",
         "species. All species show the expected positive response to moisture ",
         "except Quercus robur, whose germination remained marginal across the ",
         "entire MC range - consistent with a lot-level viability failure ",
         "independent of desiccation sensitivity. Q. robur was therefore ",
         "excluded from all germination GLMs.")) +
  theme_bw(base_size = 12)

suppressWarnings(
  ggsave(file.path(plots_dir, "germ_ro_exclusion_justification.png"),
         p_ro_check, width = 14, height = 8, dpi = 300)
)
cat("Q. robur justification figure written\n")

df.germ <- df.germ |>
  filter(!species %in% EXCLUDE_SPECIES) |>
  mutate(species = droplevels(species))

cat("Species excluded:", paste(EXCLUDE_SPECIES, collapse = ", "),
    "| rows:", n_pre_exclusion, "->", nrow(df.germ), "\n")

# --- 2. Global models ----------------------------------------------------------

glm_global_no_phase <- glm(
  germinated ~ mc * species,
  family = binomial, data = df.germ, na.action = na.fail
)

glm_global_with_phase <- glm(
  germinated ~ mc * species + phase * species,
  family = binomial, data = df.germ, na.action = na.fail
)

# --- 3. Phase screening ---------------------------------------------------------

lrt_phase <- anova(glm_global_no_phase, glm_global_with_phase, test = "LRT")

global_comparison <- tibble(
  model        = c("germ ~ mc*species",
                   "germ ~ mc*species + phase*species"),
  df           = c(attr(logLik(glm_global_no_phase), "df"),
                   attr(logLik(glm_global_with_phase), "df")),
  logLik       = c(logLik(glm_global_no_phase), logLik(glm_global_with_phase)),
  AICc         = c(AICc(glm_global_no_phase), AICc(glm_global_with_phase)),
  lrt_deviance = c(NA, lrt_phase$Deviance[2]),
  lrt_p        = c(NA, lrt_phase[["Pr(>Chi)"]][2])
)

cat("\n=== Global-model phase comparison ===\n")
print(global_comparison)

# --- 4. Dredge -------------------------------------------------------------------

dr_with_phase <- dredge(glm_global_with_phase, rank = AICc)
dr_no_phase   <- dredge(glm_global_no_phase, rank = AICc)

cat("\n=== Dredge WITH phase (top 10 of", nrow(dr_with_phase), ") ===\n")
print(head(dr_with_phase, 10), width = 220)
cat("Models within delta-AICc < 2:", sum(dr_with_phase$delta < 2), "\n")

cat("\n=== Dredge WITHOUT phase (top 10 of", nrow(dr_no_phase), ") ===\n")
print(head(dr_no_phase, 10), width = 220)
cat("Models within delta-AICc < 2:", sum(dr_no_phase$delta < 2), "\n")

best_with_phase <- get.models(dr_with_phase, 1)[[1]]
best_no_phase   <- get.models(dr_no_phase, 1)[[1]]

cat("\nBest WITH phase formula:", deparse(formula(best_with_phase)), "\n")
cat("Best WITHOUT phase formula:", deparse(formula(best_no_phase)), "\n")

phase_in_best <- function(fit) {
  any(grepl("phase", labels(terms(formula(fit)))))
}
cat("Phase retained in best WITH-phase model?", phase_in_best(best_with_phase), "\n")

# --- 5. Selected-model coefficients ----------------------------------------------

tidy_one <- function(fit, label) {
  as.data.frame(coef(summary(fit))) |>
    rownames_to_column("term") |>
    as_tibble() |>
    rename(estimate = Estimate,
           std_error = `Std. Error`,
           z_value   = `z value`,
           p_value   = `Pr(>|z|)`) |>
    mutate(model = label, .before = 1)
}

selected_coefficients <- bind_rows(
  tidy_one(best_with_phase, "best_with_phase"),
  tidy_one(best_no_phase,   "best_no_phase")
)

cat("\n--- Coefficients: best WITH phase ---\n")
print(as.data.frame(selected_coefficients |>
                      filter(model == "best_with_phase") |>
                      select(-model)))

cat("\nPer-species phase_2 shifts (best WITH phase):\n")
print(as.data.frame(selected_coefficients |>
                      filter(model == "best_with_phase",
                             grepl("^phase", term))))

# --- 6. DHARMa assumption checks ---------------------------------------------------

run_dharma <- function(fit, label) {
  sim <- simulateResiduals(fit, n = 1000, seed = SIM_SEED)
  tu <- testUniformity(sim, plot = FALSE)$p.value
  td <- testDispersion(sim, plot = FALSE)$p.value
  to <- testOutliers(sim, plot = FALSE)$p.value
  cat(sprintf("[%s] uniformity p=%.4g | dispersion p=%.4g | outliers p=%.4g\n",
              label, tu, td, to))
  png(file.path(plots_dir, paste0("glm_dharma_", label, ".png")),
      width = 1600, height = 1600, res = 150)
  plot(sim)
  dev.off()
  tibble(model = label,
         uniformity_p = tu, dispersion_p = td, outliers_p = to)
}

cat("\n=== DHARMa simulated-residual checks ===\n")
dharma_tests <- bind_rows(
  run_dharma(best_with_phase, "best_with_phase"),
  run_dharma(best_no_phase,   "best_no_phase")
)

# --- 7. Predictive figures -------------------------------------------------------

crossing_mc <- function(x, y, target = 0.5) {
  idx <- which(y >= target)
  if (length(idx) == 0 || idx[1] == 1) return(NA_real_)
  i <- idx[1]
  x[i - 1] + (target - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
}

fmt_pct <- function(v) ifelse(is.na(v), "n.r.", sprintf("%.1f", v))

# 7a. Phase-stratified curves with raw observations ----------------------------

if (!phase_in_best(best_with_phase)) {
  cat("\nPhase dropped from the best model; skipping phase-stratified curves.\n")
} else {
  curve_grid <- df.germ |>
    group_by(species, phase) |>
    summarise(mc_min = min(mc), mc_max = max(mc), .groups = "drop") |>
    rowwise() |>
    mutate(mc = list(seq(mc_min, mc_max, length.out = 200))) |>
    ungroup() |>
    unnest(mc)

  pr <- predict(best_with_phase, newdata = curve_grid, type = "response",
                se.fit = TRUE)
  curve_grid$pred <- pr$fit
  curve_grid$lo   <- pmax(0, pr$fit - 1.96 * pr$se.fit)
  curve_grid$up   <- pmin(1, pr$fit + 1.96 * pr$se.fit)

  p_curves <- ggplot(curve_grid, aes(mc, pred, colour = phase, fill = phase)) +
    geom_point(data = df.germ,
               aes(mc, germinated, colour = phase),
               alpha = 0.12, size = 0.6,
               position = position_jitter(height = 0.03, width = 0)) +
    geom_ribbon(aes(ymin = lo, ymax = up), colour = NA, alpha = 0.15) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~species, scales = "free_x") +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
    labs(x = "Moisture content (% FW0)",
         y = "Germination probability",
         title = "Germination GLM: predicted curves per species and drying batch",
         subtitle = "Lines = model predictions with 95% CI; points = raw germination outcomes (jittered)",
         colour = "Batch", fill = "Batch") +
    theme_bw(base_size = 12)

  ggsave(file.path(plots_dir, "glm_pred_curves_species_phase.png"),
         p_curves, width = 14, height = 9, dpi = 300)
  cat("\nPhase-stratified figure written to",
      file.path(plots_dir, "glm_pred_curves_species_phase.png"), "\n")
}

# 7b. Batch-averaged sensitivity curves with graphical MC50 guides --------------
#
# Predictions marginalize the WITH-phase model over drying batches via
# emmeans (weights = "proportional": observed phase frequencies; averaging
# on the link scale, then back-transformed). The vertical guides mark where
# the 0.5 reference cuts the prediction curve and its CI ribbon: a graphical
# MC50 preview, not the formal estimate of the dedicated MC50 script.

sp_ranges <- df.germ |>
  group_by(species) |>
  summarise(mc_min = min(mc), mc_max = max(mc), .groups = "drop")

avg_grid <- bind_rows(lapply(seq_len(nrow(sp_ranges)), function(i) {
  grid_i <- seq(sp_ranges$mc_min[i], sp_ranges$mc_max[i], length.out = 200)
  emm_i <- emmeans(best_with_phase, ~ mc | species,
                   at = list(mc = grid_i),
                   weights = "proportional", type = "response")
  as.data.frame(emm_i)
})) |>
  rename(pred = prob, lo = asymp.LCL, up = asymp.UCL)

mc50_tab <- avg_grid |>
  group_by(species) |>
  summarise(
    mc50  = crossing_mc(mc, pred),
    ci_lo = crossing_mc(mc, up),
    ci_hi = crossing_mc(mc, lo),
    .groups = "drop"
  ) |>
  mutate(label = sprintf("MC50 = %s%% [%s, %s]",
                         fmt_pct(mc50), fmt_pct(ci_lo), fmt_pct(ci_hi)))

cat("\nGraphical MC50 preview (batch-marginalized, emmeans proportional weights):\n")
print(as.data.frame(mc50_tab))

p_avg <- ggplot(avg_grid, aes(mc, pred)) +
  geom_point(data = df.germ, aes(mc, germinated),
             colour = "black", alpha = 0.12, size = 0.6,
             position = position_jitter(height = 0.03, width = 0)) +
  geom_ribbon(aes(ymin = lo, ymax = up), fill = "#2c7fb8", alpha = 0.18,
              colour = NA) +
  geom_line(colour = "#2c7fb8", linewidth = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_vline(data = mc50_tab, aes(xintercept = ci_lo),
             linetype = "dashed", colour = "grey55", na.rm = TRUE) +
  geom_vline(data = mc50_tab, aes(xintercept = ci_hi),
             linetype = "dashed", colour = "grey55", na.rm = TRUE) +
  geom_vline(data = mc50_tab, aes(xintercept = mc50),
             linetype = "solid", colour = "#b2182b", linewidth = 0.7,
             na.rm = TRUE) +
  geom_text(data = mc50_tab, aes(x = -Inf, y = 0.06, label = label),
            hjust = -0.05, vjust = 0, size = 3.1, inherit.aes = FALSE) +
  facet_wrap(~species, scales = "free_x") +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  labs(x = "Moisture content (% FW0)",
       y = "Germination probability",
       title = "Moisture-sensitivity curves averaged over drying batches",
       subtitle = paste0("Predictions marginalized over drying batches ",
                         "(emmeans, proportional weights); ",
                         "red vertical = graphical MC50, dashed verticals = its ",
                         "95% CI, horizontal reference at 0.5"),
       caption = "Quercus robur excluded (lot-level viability failure)") +
  theme_bw(base_size = 12)

ggsave(file.path(plots_dir, "glm_pred_curves_species_average.png"),
       p_avg, width = 14, height = 9, dpi = 300)
cat("Batch-marginalized figure with MC50 guides written to",
    file.path(plots_dir, "glm_pred_curves_species_average.png"), "\n")

# --- 8. Exports ----------------------------------------------------------------------

write.csv(as.data.frame(dr_with_phase),
          file.path(tablas_dir, "glm_dredge_with_phase.csv"),
          row.names = FALSE)
write.csv(as.data.frame(dr_no_phase),
          file.path(tablas_dir, "glm_dredge_no_phase.csv"),
          row.names = FALSE)
write.csv(global_comparison,
          file.path(tablas_dir, "glm_global_comparison.csv"),
          row.names = FALSE)
write.csv(selected_coefficients,
          file.path(tablas_dir, "glm_selected_coefficients.csv"),
          row.names = FALSE)
write.csv(dharma_tests,
          file.path(tablas_dir, "glm_dharma_tests.csv"),
          row.names = FALSE)

cat("\ns.02 done: dredge tables, global comparison, coefficients, DHARMa checks",
    "and predictive figures written\n")
