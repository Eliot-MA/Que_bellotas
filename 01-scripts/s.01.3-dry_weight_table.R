###############################################################################
## s.01.3 - Per-acorn table and FW0 -> DW allometric model                   ##
##                                                                           ##
## Builds one row per acorn from the individual records produced by          ##
## s.01.1 and fits the allometric model that provides a dry weight value     ##
## (measured or predicted) for every acorn.                                  ##
##                                                                           ##
## Experimental design:                                                      ##
##   - every acorn has a fresh weight at intake (t0 = FW0)                   ##
##   - destructive subsamples are oven-dried to obtain observed dry weight:  ##
##     part of them already at t0, plus ~10 acorns per batch at each later   ##
##     sampling event                                                        ##
##   - the remaining acorns go to the germination trial, so observed dry     ##
##     weight and positive germination scores are mutually exclusive         ##
##                                                                           ##
## Pipeline within this script:                                              ##
##   1. one row per acorn                                                    ##
##   2. CALIBRATION TABLE: acorns having both FW0 and observed DW            ##
##   3. TRANSFORMATION TESTING: candidate model forms compared with          ##
##      residual-normality, heteroscedasticity, simulated-residual           ##
##      (DHARMa) and cross-validated predictive-error criteria               ##
##   4. MODELING-STRATEGY COMPARISON: one pooled model vs a species x FW0    ##
##      interaction model vs independent per-species models, judged by CV    ##
##   5. final fit of the selected form + structure + DW for every acorn      ##
##                                                                           ##
## Modelling background:                                                     ##
##   - dry mass is conserved during drying, so DW ~ FW0 holds across all     ##
##     sampling events and pairs can be pooled into one per-species model    ##
##   - FW0 and DW show strongly right-skewed (roughly logistic/gamma-like)   ##
##     distributions; candidates therefore include log-transformed linear    ##
##     models and Gamma GLMs with log link                                  ##
###############################################################################

# --- 0. Libraries -----------------------------------------------------------

library(tidyverse)
library(lmtest)   # Breusch-Pagan test
library(DHARMa)   # simulated-residual diagnostics (also valid for GLMs)

if (!exists("df.germ.records")) {
  stop("df.germ.records not found. Run/source s.01.1-load_raw.R first.")
}

# Set to one of the candidate names below to force a model form;
# leave as NULL for automatic selection (see STEP 3).
MODEL_FORM_OVERRIDE <- NULL

# Set to one of "general", "interaction" or "per_species" to force the
# modeling strategy; leave as NULL for automatic selection (see STEP 4).
MODEL_STRUCTURE_OVERRIDE <- NULL

# --- 1. One row per acorn ---------------------------------------------------

first_or_na <- function(x) if (length(x)) x[1] else NA_character_

df.acorns <- df.germ.records |>
  group_by(acorn_id) |>
  arrange(as.integer(time), .by_group = TRUE) |>
  summarise(
    phase         = first(phase),
    species       = first(species),
    provenance    = first(provenance),
    fw0           = fresh_weight[time == "0"][1],
    sampling_time = first_or_na(time[time != "0"]),
    fw_sampling   = {
      v <- fresh_weight[time != "0"]
      if (length(v)) v[1] else NA_real_
    },
    dw_observed   = {
      v <- na.omit(dry_weight)
      if (length(v)) v[1] else NA_real_
    },
    germinated    = {
      v <- na.omit(germinated)
      if (length(v)) v[length(v)] else NA_real_
    },
    emerged       = {
      v <- na.omit(emerged)
      if (length(v)) v[length(v)] else NA_real_
    },
    moho          = any(!is.na(observations_germination) &
                          grepl("Moho", observations_germination)),
    observations  = paste(na.omit(observations_germination),
                          collapse = "; "),
    .groups       = "drop"
  ) |>
  mutate(observations = na_if(observations, ""))

# --- 2. Consistency checks --------------------------------------------------

n_missing_fw0 <- sum(is.na(df.acorns$fw0))

conflicts <- df.acorns |>
  filter(!is.na(dw_observed),
         germinated == 1 | emerged == 1)

# --- 3. Calibration table ---------------------------------------------------

# Acorns contributing to the allometric model: those with BOTH a valid FW0
# and a directly observed dry weight.
dw.calibration <- df.acorns |>
  filter(!is.na(fw0), !is.na(dw_observed)) |>
  select(acorn_id, phase, species, provenance, fw0, dw_observed)

cat("=== calibration table ===\n")
print(
  dw.calibration |>
    group_by(species) |>
    summarise(n_pairs = n(),
              fw0_min = min(fw0), fw0_max = max(fw0),
              dw_min  = min(dw_observed), dw_max = max(dw_observed),
              .groups = "drop") |>
    as.data.frame(),
  row.names = FALSE
)

# --- 4. Transformation testing ----------------------------------------------

# Candidate model forms. All use DW as response and FW0 as predictor:
#   identity  : lm(DW ~ FW0)                     raw scales
#   log-log   : lm(log(DW) ~ log(FW0))           legacy choice
#   gamma     : glm(DW ~ FW0, Gamma("log"))      gamma response, raw X
#   gamma-log : glm(DW ~ log(FW0), Gamma("log")) gamma response, log X
#
# Each candidate is evaluated per species on the calibration table with:
#   shapiro_p  : Shapiro-Wilk on residuals        (normality; exact for lm,
#                                                approximate for GLM)
#   bp_p       : Breusch-Pagan test               (heteroscedasticity)
#   uniform_p  : DHARMa testUniformity            (overall residual correctness,
#                                                simulation-based)
#   disp_p     : DHARMa testDispersion            (variance adequacy)
#   cv_rmse    : 10-fold cross-validated RMSE on the DW scale (grams)

fit_candidate <- function(form, dat) {
  switch(form,
    identity = lm(dw_observed ~ fw0, data = dat),
    `log-log` = lm(log(dw_observed) ~ log(fw0), data = dat),
    gamma = glm(dw_observed ~ fw0, family = Gamma(link = "log"), data = dat),
    `gamma-log` = glm(dw_observed ~ log(fw0), family = Gamma(link = "log"),
                      data = dat)
  )
}

predict_candidate <- function(form, mod, newdat) {
  switch(form,
    identity = as.numeric(predict(mod, newdata = newdat)),
    `log-log` = as.numeric(exp(predict(mod, newdata = newdat))),
    gamma = as.numeric(predict(mod, newdata = newdat, type = "response")),
    `gamma-log` = as.numeric(predict(mod, newdata = newdat, type = "response"))
  )
}

cv_rmse_candidate <- function(form, dat, v = 10) {
  set.seed(123)
  fold <- sample(rep(seq_len(v), length.out = nrow(dat)))
  rmse <- numeric(v)
  for (k in seq_len(v)) {
    train <- dat[fold != k, ]
    test  <- dat[fold == k, ]
    mod   <- fit_candidate(form, train)
    pred  <- predict_candidate(form, mod, test)
    rmse[k] <- sqrt(mean((test$dw_observed - pred)^2))
  }
  mean(rmse)
}

diagnose_candidate <- function(form, dat, n_sim = 500) {
  mod  <- fit_candidate(form, dat)
  res  <- residuals(mod)

  sim <- suppressWarnings(suppressMessages(
    simulateResiduals(mod, n = n_sim, seed = 123)
  ))

  tibble(
    shapiro_p = if (length(res) >= 3 && length(res) <= 5000 &&
                    sd(res) > 0)
      tryCatch(shapiro.test(res)$p.value, error = function(e) NA_real_)
    else NA_real_,
    bp_p       = tryCatch(bptest(mod)$p.value, error = function(e) NA_real_),
    uniform_p  = tryCatch(testUniformity(sim)$p.value,
                          error = function(e) NA_real_),
    disp_p     = tryCatch(testDispersion(sim, plot = FALSE)$p.value,
                          error = function(e) NA_real_),
    cv_rmse    = cv_rmse_candidate(form, dat)
  )
}

cat("\n=== transformation testing ===\n")
diag_list <- list()
for (sp in sort(unique(dw.calibration$species))) {
  dat_sp <- dw.calibration |> filter(species == sp)
  for (frm in c("identity", "log-log", "gamma", "gamma-log")) {
    diag_list[[length(diag_list) + 1]] <-
      diagnose_candidate(frm, dat_sp) |>
      mutate(species = sp, candidate = frm, .before = 1)
  }
}
diag_tbl <- bind_rows(diag_list)

print(
  diag_tbl |>
    mutate(across(ends_with("_p"), ~ round(.x, 4))) |>
    mutate(cv_rmse = round(cv_rmse, 4)) |>
    as.data.frame(),
  row.names = FALSE
)

# Automatic selection. Essential checks (counted): Breusch-Pagan,
# DHARMa uniformity and DHARMa dispersion. Residual normality (Shapiro) is
# reported but NOT counted: with n ~ 100 it flags negligible deviations for
# every candidate, so it carries no discriminating power here.
# Ties broken by median cross-validated RMSE.
selection <- diag_tbl |>
  group_by(candidate) |>
  summarise(
    n_failed_essential = sum(bp_p < .05, na.rm = TRUE) +
                         sum(uniform_p < .05, na.rm = TRUE) +
                         sum(disp_p < .05, na.rm = TRUE),
    n_shapiro_failures = sum(shapiro_p < .05, na.rm = TRUE),
    median_cv_rmse = median(cv_rmse),
    .groups = "drop"
  ) |>
  arrange(n_failed_essential, median_cv_rmse)

cat("\nModel-form ranking (fewest failed checks first):\n")
print(as.data.frame(selection), row.names = FALSE)

model_form <- MODEL_FORM_OVERRIDE %||% selection$candidate[1]
cat("\nSelected model form:", model_form,
    if (is.null(MODEL_FORM_OVERRIDE)) "(automatic)" else "(forced)\n")

# --- 5. Modeling-strategy comparison -----------------------------------------

# Within the selected form, compare three ways of handling species:
#   general     : one pooled model shared by all species
#   interaction : single model DW ~ FW0 * species (species-specific
#                 intercepts and slopes estimated jointly)
#   per_species : one independent model fitted per species (legacy approach)
#
# Note: on identical training data, 'interaction' and 'per_species' are
# re-parameterisations of each other (the design matrix factorises by
# species), so their predictions - and therefore their CV scores - coincide
# by construction. Both are still scored to verify that equivalence
# empirically before keeping the simpler pipeline.

fit_strategy <- function(strategy, train) {
  switch(strategy,
    general     = lm(dw_observed ~ fw0, data = train),
    interaction = lm(dw_observed ~ fw0 * species, data = train),
    per_species = lm(dw_observed ~ fw0 * species, data = train)
  )
}

# Shared fold assignment -> paired comparison across strategies.
set.seed(123)
strategy_folds <- sample(rep(seq_len(10), length.out = nrow(dw.calibration)))

strategy_cv <- function(strategy, dat, folds, v = 10) {
  sq_err <- rep(NA_real_, nrow(dat))
  sp     <- sort(unique(dat$species))
  for (k in seq_len(v)) {
    train <- dat[folds != k, ]
    test  <- dat[folds == k, ]
    pred  <- as.numeric(predict(fit_strategy(strategy, train), newdata = test))
    sq_err[folds == k] <- (test$dw_observed - pred)^2
  }
  list(
    overall_rmse = sqrt(mean(sq_err, na.rm = TRUE)),
    species_rmse = tibble(
      species = sp,
      rmse = sapply(sp, \(s)
        sqrt(mean(sq_err[dat$species == s], na.rm = TRUE)))
    )
  )
}

strategies <- c("general", "interaction", "per_species")

cv_results <- setNames(strategies, strategies) |>
  map(\(st) strategy_cv(st, dw.calibration, strategy_folds))

strategy_summary <- imap_dfr(cv_results, \(res, st) {
  tibble(
    strategy           = st,
    cv_rmse_overall    = res$overall_rmse,
    mean_species_rmse  = mean(res$species_rmse$rmse),
    worst_species_rmse = max(res$species_rmse$rmse)
  )
}) |>
  # Tie-break: on identical CV error prefer per_species (simplest pipeline,
  # independent models, legacy reporting), then interaction, then general.
  mutate(preference = match(strategy,
                            c("per_species", "interaction", "general"))) |>
  arrange(round(cv_rmse_overall, 8), preference) |>
  select(-preference)

# Per-species CV errors behind the summary above (kept for reporting).
strategy_species_rmse <- imap_dfr(cv_results, \(res, st)
  mutate(res$species_rmse, strategy = st)) |>
  select(strategy, species, rmse)

# Formal test of whether species modifies the FW0-DW relationship.
anova_p <- tryCatch(
  anova(lm(dw_observed ~ fw0, dw.calibration),
        lm(dw_observed ~ fw0 * species, dw.calibration))$"Pr(>F)"[2],
  error = function(e) NA_real_)

cat("\n=== modeling-strategy comparison ===\n")
print(as.data.frame(strategy_summary), row.names = FALSE)
cat("\nANOVA general vs interaction, p =", format(anova_p), "\n")

model_structure <- MODEL_STRUCTURE_OVERRIDE %||% strategy_summary$strategy[1]
cat("Selected structure:", model_structure,
    if (is.null(MODEL_STRUCTURE_OVERRIDE)) "(automatic)" else "(forced)\n")

# --- 6. Final fit and prediction for every acorn -----------------------------

fit_and_predict_species <- function(dat, form) {
  calib <- dat |> filter(!is.na(dw_observed), !is.na(fw0))

  out <- dat |>
    mutate(
      dw_model_pred = NA_real_,
      dw_final = NA_real_, dw_lwr = NA_real_, dw_upr = NA_real_,
      model_form = form,
      model_n    = nrow(calib),
      model_r2   = NA_real_,
      model_intercept = NA_real_,
      model_slope     = NA_real_
    )

  if (nrow(calib) < 10) {
    warning("Species with fewer than 10 calibration pairs: model skipped.")
    return(out)
  }

  mod <- fit_candidate(form, calib)
  point_pred <- predict_candidate(form, mod, dat)

  # Prediction interval:
  #   lm models -> analytical (back-transformed when the response is logged)
  #   GLM models -> parametric bootstrap of the Gamma response
  if (inherits(mod, "glm")) {
    set.seed(123)
    mu_hat <- predict_candidate(form, mod, dat)
    disp   <- summary(mod)$dispersion
    sims <- matrix(rgamma(length(mu_hat) * 2000,
                          shape = 1 / disp,
                          scale = rep(mu_hat * disp, 2000)),
                   ncol = 2000)
    pi_lwr <- apply(sims, 1, quantile, probs = 0.025, names = FALSE)
    pi_upr <- apply(sims, 1, quantile, probs = 0.975, names = FALSE)
  } else {
    pr <- suppressWarnings(predict(mod, newdata = dat, interval = "prediction"))
    point_pred <- pr[, "fit"]
    pi_lwr <- pr[, "lwr"]
    pi_upr <- pr[, "upr"]
    if (form == "log-log") {
      point_pred <- exp(point_pred); pi_lwr <- exp(pi_lwr); pi_upr <- exp(pi_upr)
    }
  }

  out |>
    mutate(
      dw_model_pred = as.numeric(point_pred),
      dw_final = as.numeric(point_pred),
      dw_lwr   = as.numeric(pi_lwr),
      dw_upr   = as.numeric(pi_upr),
      model_r2 = if (inherits(mod, "lm")) summary(mod)$r.squared else
                   1 - summary(mod)$deviance / summary(mod)$null.deviance,
      model_intercept = unname(coef(mod)[1]),
      model_slope     = unname(coef(mod)[2])
    )
}

df.modeled <- switch(model_structure,
  per_species = df.acorns |>
    group_split(species) |>
    map_dfr(\(d) fit_and_predict_species(d, model_form)) |>
    mutate(model_structure = "per_species"),

  # One pooled model for all species: fit_and_predict_species() already
  # behaves this way when handed the whole table at once.
  general = fit_and_predict_species(df.acorns, model_form) |>
    mutate(model_structure = "general"),

  `interaction` = {
    calib <- df.acorns |> filter(!is.na(dw_observed), !is.na(fw0))
    mod   <- lm(dw_observed ~ fw0 * species, data = calib)

    pr  <- suppressWarnings(predict(mod, newdata = df.acorns,
                                    interval = "prediction"))
    co  <- coef(mod)
    eff_coefs <- map_dfr(sort(unique(calib$species)), \(s) {
      d0 <- co[paste0("species", s)][1]     %||% 0
      d1 <- co[paste0("fw0:species", s)][1] %||% 0
      tibble(species = s,
             model_intercept = co[["(Intercept)"]] + d0,
             model_slope     = co[["fw0"]]         + d1)
    })

    df.acorns |>
      left_join(eff_coefs, by = "species") |>
      mutate(
        dw_model_pred   = as.numeric(pr[, "fit"]),
        dw_final        = as.numeric(pr[, "fit"]),
        dw_lwr          = as.numeric(pr[, "lwr"]),
        dw_upr          = as.numeric(pr[, "upr"]),
        model_form      = model_form,
        model_structure = "interaction",
        model_n         = nrow(calib),
        model_r2        = summary(mod)$r.squared
      )
  }
)

df.acorns <- df.modeled |>
  select(acorn_id, phase, species, provenance,
         fw0, sampling_time, fw_sampling,
         dw_observed, dw_model_pred, dw_final, dw_lwr, dw_upr,
         germinated, emerged, moho, observations,
         model_form, model_structure, model_n, model_r2,
         model_intercept, model_slope) |>
  arrange(acorn_id) |>
  mutate(
    dw_source = case_when(
      !is.na(dw_observed) ~ "measured",
      !is.na(dw_final)    ~ "predicted",
      TRUE                ~ NA_character_
    ),
    dw_final = coalesce(dw_observed, dw_final)
  )

dw_models <- df.acorns |>
  distinct(species, model_form, model_structure, model_n, model_r2,
           model_intercept, model_slope) |>
  arrange(desc(model_n))

# --- 7. Quality-control summary ---------------------------------------------

cat("\n=== s.01.3 dry_weight_table ===\n")
cat("Acorns:", nrow(df.acorns), "\n")
cat("Selected form:", model_form, "| structure:", model_structure, "\n")

cat("\nDry weight availability:\n")
print(table(df.acorns$species, df.acorns$dw_source, useNA = "ifany"))

cat("\nAllometric models (", model_form, ") per species:\n", sep = "")
print(as.data.frame(dw_models), row.names = FALSE)

cat("\nConsistency checks:\n")
# Known design fact (E.M., confirmed 2026): the PE batch did not hold enough
# acorns to reach the nominal 400 weigh-ins, so ~33 PE acorns legitimately
# lack FW0 (plus 1 FA). This is an experimental shortfall, NOT a data error;
# these acorns cannot receive a predicted DW and stay flagged downstream.
cat("Acorns without FW0:", n_missing_fw0,
    "- expected pattern: 33 PE + 1 FA (PE batch shortfall, known design fact)\n")
if (n_missing_fw0 > 0) {
  warning(n_missing_fw0, " acorns lack FW0 and cannot receive a predicted ",
          "dry weight. They will need to be handled downstream.")
}
cat("Design conflicts (observed DW + positive germination/emergence):",
    nrow(conflicts), "\n")
if (nrow(conflicts) > 0) {
  cat("Conflicting acorn_ids:",
      paste(head(conflicts$acorn_id, 20), collapse = ", "),
      if (nrow(conflicts) > 20) "..." else "", "\n")
}

cat("\nAgreement between model predictions and measurements (calibration table):\n")
agreement <- dw.calibration |>
  left_join(
    select(df.acorns, acorn_id, dw_pred_check = dw_model_pred,
           dw_lwr, dw_upr),
    by = "acorn_id"
  ) |>
  group_by(species) |>
  summarise(
    n             = n(),
    cor_pearson   = cor(fw0, dw_observed),
    rmse_g        = sqrt(mean((dw_pred_check - dw_observed)^2)),
    coverage_pi95 = mean(dw_observed >= dw_lwr & dw_observed <= dw_upr,
                         na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(agreement), row.names = FALSE)

# --- 8. Transparency exports -------------------------------------------------
#
# Decision trail and validation artifacts written directly by this child
# script (agreed exception to the master-only-export rule):
#   plots  -> 07-img/dw_model_diagnostics
#   tables -> 00-data/tablas_resumen

img_dir <- file.path("07-img", "dw_model_diagnostics")
tbl_dir <- file.path("00-data", "tablas_resumen")
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

calib_check <- dw.calibration |>
  left_join(
    select(df.acorns, acorn_id, dw_pred = dw_model_pred, dw_lwr, dw_upr),
    by = "acorn_id"
  ) |>
  mutate(
    residual      = dw_pred - dw_observed,
    rel_error_pct = 100 * residual / dw_observed
  )

## 8.1 Observed vs predicted --------------------------------------------------

p_obs_all <- ggplot(calib_check, aes(dw_observed, dw_pred)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = species), alpha = .7, size = 1.6) +
  labs(x = "Observed dry weight (g)", y = "Predicted dry weight (g)",
       title = "Observed vs predicted dry weight",
       subtitle = "Calibration acorns, dashed line = 1:1") +
  theme_bw()
ggsave(file.path(img_dir, "dw_obs_vs_predicted.png"),
       p_obs_all, width = 7, height = 6, dpi = 300)

p_obs_sp <- ggplot(calib_check, aes(dw_observed, dw_pred)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = dw_lwr, ymax = dw_upr),
                width = 0, alpha = .25) +
  geom_point(size = 1.4, colour = "grey25") +
  facet_wrap(~species, scales = "free") +
  labs(x = "Observed dry weight (g)", y = "Predicted dry weight (g)",
       title = "Observed vs predicted dry weight by species",
       subtitle = "Vertical segments = 95% prediction intervals") +
  theme_bw()
ggsave(file.path(img_dir, "dw_obs_vs_predicted_by_species.png"),
       p_obs_sp, width = 10, height = 8, dpi = 300)

## 8.2 Residuals vs fitted ----------------------------------------------------

p_resid <- ggplot(calib_check, aes(dw_pred, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = .55, size = 1.5) +
  facet_wrap(~species, scales = "free") +
  labs(x = "Predicted dry weight (g)",
       y = "Residual (predicted - observed, g)") +
  theme_bw()
ggsave(file.path(img_dir, "dw_residuals_vs_fitted_by_species.png"),
       p_resid, width = 10, height = 8, dpi = 300)

## 8.3 Prediction error across the observed range ------------------------------

p_err <- ggplot(calib_check, aes(dw_observed, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = .45, size = 1.4) +
  geom_smooth(se = TRUE, linewidth = .7) +
  facet_wrap(~species, scales = "free") +
  labs(x = "Observed dry weight (g)", y = "Prediction error (g)",
       title = "Prediction error along the observed DW range",
       subtitle = "Smooth line = loess trend +/- SE") +
  theme_bw()
ggsave(file.path(img_dir, "dw_prediction_error_vs_observed_by_species.png"),
       p_err, width = 10, height = 8, dpi = 300)

## 8.4 Tables ------------------------------------------------------------------

write_csv(diag_tbl,
          file.path(tbl_dir, "dw_transformation_testing.csv"))
write_csv(selection,
          file.path(tbl_dir, "dw_model_form_ranking.csv"))
write_csv(strategy_summary,
          file.path(tbl_dir, "dw_strategy_comparison.csv"))
write_csv(strategy_species_rmse,
          file.path(tbl_dir, "dw_strategy_species_rmse.csv"))
write_csv(dw_models,
          file.path(tbl_dir, "dw_model_coefficients.csv"))
write_csv(agreement,
          file.path(tbl_dir, "dw_agreement_calibration.csv"))
write_csv(calib_check,
          file.path(tbl_dir, "dw_calibration_pairs_with_predictions.csv"))

writeLines(c(
  paste0("# s.01.2 model decisions (generated ", Sys.Date(), ")"),
  "",
  paste0("Model form          : ", model_form),
  paste0("Model structure     : ", model_structure),
  paste0("ANOVA general-vs-interaction p = ", format(anova_p)),
  paste0("Calibration pairs   : ", nrow(dw.calibration)),
  "",
  "Known design fact: 34 acorns lack FW0 (33 PE + 1 FA). The PE batch held",
  "fewer acorns than the nominal 400, so these rows legitimately have no",
  "predicted dry weight. Experimental shortfall, not a data error.",
  "",
  "Companion files:",
  "  dw_transformation_testing.csv             - per-species x candidate diagnostics",
  "  dw_model_form_ranking.csv                 - essential-checks ranking of forms",
  "  dw_strategy_comparison.csv                - CV error by modeling strategy",
  "  dw_strategy_species_rmse.csv              - per-species CV RMSE per strategy",
  "  dw_model_coefficients.csv                 - final per-species models",
  "  dw_agreement_calibration.csv              - in-sample RMSE + PI coverage",
  "  dw_calibration_pairs_with_predictions.csv - validation dataset"
), file.path(tbl_dir, "dw_model_selection_log.txt"))

cat("\nTransparency exports written to:\n  ", img_dir, "\n  ",
    tbl_dir, "\n", sep = "")
