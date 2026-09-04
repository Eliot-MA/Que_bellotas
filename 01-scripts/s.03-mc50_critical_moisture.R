###############################################################################
## s.03 - Critical moisture content (MC50) estimation                        ##
##                                                                           ##
## Estimates the moisture content at which germination probability equals    ##
## 0.5 (MC50) per species, from the batch-adjusted logistic model selected   ##
## in s.02 (germ ~ mc*species + phase; Quercus robur excluded due to lot     ##
## failure).                                                                 ##
##                                                                           ##
## Estimand (E.M. protocol): MC50 marginalized over drying batches on the    ##
## link scale with proportional batch weights - the same target as the       ##
## emmeans sensitivity curves of s.02:                                       ##
##   eta(mc | s) = b0_s + w2 * dPhase + b1_s * mc                            ##
##   MC50_s     = -(b0_s + w2 * dPhase) / b1_s                               ##
##   where w2 is the observed proportion of phase 2 in the model data.       ##
##                                                                           ##
## Inference routes (E.M. decision):                                         ##
##   1. Delta method (PRIMARY): asymptotic SE of the ratio from a            ##
##      first-order Taylor expansion using the full covariance matrix of     ##
##      every involved coefficient, Var(g) = grad' V grad (Collett 1999;     ##
##      cf. Amimi et al. 2020, adapted dose.p). With only (b0, b1) involved  ##
##      this reduces to the classic closed form.                             ##
##   2. Case-resampling bootstrap (VALIDATION): resampling stratified by     ##
##      species x phase with fixed stratum sizes (preserves the design and   ##
##      w2); percentile 95% CI. Replicates with non-convergent fits or       ##
##      biologically undefined MC50 (slope <= 0, non-finite, or              ##
##      |MC50| > 500) are discarded per species and their rates reported.    ##
##                                                                           ##
## Outputs (00-data/tablas_resumen/):                                        ##
##   mc50_delta.csv          point estimates + delta SE and CI               ##
##   mc50_bootstrap.csv      bootstrap distribution summary                  ##
##   mc50_summary.csv        side-by-side comparison                         ##
## Figure  (07-img/mc50_estimation/):                                        ##
##   mc50_forest_species.png    forest plot of MC50 with both CIs            ##
##   mc50_forest_bioclimate.png article figure: Delta MC50 coloured by       ##
##                              species bioclimate                           ##
###############################################################################

library(tidyverse)

# --- 0. Configuration ---------------------------------------------------------

EXCLUDE_SPECIES <- c("RO")   # lot-level viability failure (see s.02)
SIM_SEED        <- 20260823
BOOT_REPS       <- 2000

tablas_dir <- "00-data/tablas_resumen"

# --- 1. Data (same filters as s.02) --------------------------------------------

df.mc <- read.csv("00-data/sensitivity_germination_long.csv") |>
  mutate(
    species = factor(species),
    phase   = recode_factor(factor(phase), `1` = "phase_1", `2` = "phase_2")
  ) |>
  filter(!is.na(mc), !is.na(germinated),
         flag_weird_mc %in% c("ok", "very_high"),
         !species %in% EXCLUDE_SPECIES) |>
  mutate(species = droplevels(species))

cat("=== s.03 MC50 estimation ===\n")
cat("Rows:", nrow(df.mc),
    "| species:", paste(levels(df.mc$species), collapse = " "), "\n")

# --- 2. Model (structure frozen from the s.02 dredge winner) -----------------------

fit_mc50 <- glm(germinated ~ mc + phase + species + mc:species,
                family = binomial, data = df.mc)

# Coefficient holding the batch shift; validated once here and reused below
ph2_coef <- "phasephase_2"
stopifnot(!is.na(coef(fit_mc50)[ph2_coef]))

w2 <- mean(df.mc$phase == "phase_2")
cat(sprintf("Proportional batch weight w2 (phase_2) = %.4f\n", w2))

# --- 3. Delta method ---------------------------------------------------------------

# g = -(b0_s + w2 * dPhase) / b1_s. Gradient wrt every involved coefficient;
# Var(g) = grad' V grad reduces to the classic two-term closed form when only
# intercept and slope are involved:
#   Var(g) = Var(b0)/b1^2 + b0^2 Var(b1)/b1^4 - 2 b0 Cov(b0,b1)/b1^3
# Extended delta method: returns point estimates, SEs, CIs, AND gradient vectors
# for each species. Gradients are needed downstream to compute covariances between
# MC50 estimates of different species (required for pairwise contrasts).
#
# Reference for delta-method CI on ED50/MC50 (ratio of coefficients):
#   Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the 50 per
#   cent response dose. Statistics in Medicine, 22(12), 1977-1988.
#   doi:10.1002/sim.1368
mc50_delta <- function(fit, data, w2) {
  cf  <- coef(fit)
  V   <- vcov(fit)
  ref <- levels(data$species)[1]

  bind_rows(lapply(levels(data$species), function(sp) {
    int_part <- c("(Intercept)", if (sp != ref) paste0("species", sp))
    slp_part <- c("mc",          if (sp != ref) paste0("mc:species", sp))

    A <- sum(cf[int_part]) + w2 * cf[ph2_coef]   # marginalized intercept
    B <- sum(cf[slp_part])                       # species-specific slope
    g <- -A / B

    grad <- setNames(rep(0, length(cf)), names(cf))
    grad[int_part] <- -1 / B
    grad[ph2_coef] <- -w2 / B
    grad[slp_part] <- A / B^2

    se <- sqrt(as.numeric(t(grad) %*% V %*% grad))

    tibble(species = sp,
           mc50     = g,
           se_delta = se,
           ci_lo    = g - qnorm(0.975) * se,
           ci_hi    = g + qnorm(0.975) * se,
           grad     = list(grad))   # store gradient for pairwise covariances
  }))
}

delta_tab <- mc50_delta(fit_mc50, df.mc, w2)

cat("\n=== Delta-method MC50 (primary) ===\n")
print(as.data.frame(delta_tab))

# --- 3b. Pairwise contrasts of MC50 (multiple comparisons) --------------------------
#
# The MC50 for each species is a ratio of linear functions of the model
# coefficients: MC50_s = -A_s / B_s.  To compare MC50 between species i and j we
# need Var(MC50_i - MC50_j), which requires the covariance between the two ratio
# estimators.  This covariance is obtained from the full variance-covariance matrix
# of the model coefficients via the delta method (Faraggi et al. 2003):
#
#   Cov(MC50_i, MC50_j) = grad_i' %*% V %*% grad_j
#
# where grad_s is the gradient of MC50_s with respect to every coefficient in the
# model, and V is vcov(fit).  This is the same approach used internally by the
# drc::EDcomp function for comparing ED50 between dose-response curves.
#
# Multiple-testing protection: Sidak adjustment is appropriate for a
# pre-specified family of pairwise contrasts (analogous to a priori hypotheses;
# see → d.04.model_species.R:435-441 for rationale).

# Covariance between MC50 of two species from their gradient vectors
cov_mc50 <- function(grad_i, grad_j, V) {
  as.numeric(t(grad_i) %*% V %*% grad_j)
}

# All-pairs comparison table: for every unordered pair (i, j) compute the
# difference MC50_i - MC50_j, its SE (using the full covariance structure),
# and the unadjusted 95% CI.
pairwise_mc50 <- function(delta_tab, V) {
  sp <- delta_tab$species
  grads <- setNames(delta_tab$grad, sp)
  mc50s <- setNames(delta_tab$mc50, sp)
  n <- length(sp)

  pairs <- t(combn(sp, 2))
  bind_rows(lapply(seq_len(nrow(pairs)), function(k) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    d  <- mc50s[[i]] - mc50s[[j]]
    vi <- delta_tab$se_delta[delta_tab$species == i]^2
    vj <- delta_tab$se_delta[delta_tab$species == j]^2
    covij <- cov_mc50(grads[[i]], grads[[j]], V)
    se <- sqrt(vi + vj - 2 * covij)
    tibble(pair = paste(i, "-", j),
           difference = d, se = se,
           ci_lo = d - qnorm(0.975) * se,
           ci_hi = d + qnorm(0.975) * se)
  }))
}

V_coefs <- vcov(fit_mc50)
pw_tab  <- pairwise_mc50(delta_tab, V_coefs)

cat("\n=== Pairwise MC50 differences (unadjusted) ===\n")
print(as.data.frame(pw_tab))

# Sidak-adjusted p-values for the family of pairwise contrasts
pw_tab <- pw_tab |>
  mutate(z = abs(difference) / se,
         p_raw = 2 * pnorm(-z),
         m = nrow(pw_tab),                           # number of tests
         p_sidak = 1 - (1 - p_raw)^m,               # Sidak adjustment
         sig_005 = p_sidak < 0.05)

cat("\n=== Pairwise MC50 differences (Sidak-adjusted) ===\n")
print(as.data.frame(pw_tab |> dplyr::select(pair, difference, se, p_raw, p_sidak, sig_005)))

# --- 3c. Compact-letter display (CLD) from Sidak-adjusted pairwise CI -----------
#
# Two species share a letter when the Sidak-adjusted 95% CI of their difference
# includes zero (i.e. not significantly different).  Algorithm follows the same
# greedy approach used in d.04.model_species.R:362-388.

cld_from_sidak <- function(species_order, pw) {
  n <- length(species_order)
  sigmat <- matrix(FALSE, n, n,
                   dimnames = list(species_order, species_order))

  for (k in seq_len(nrow(pw))) {
    gr <- strsplit(pw$pair[k], " - ")[[1]]
    if (pw$sig_005[k] && gr[1] %in% species_order && gr[2] %in% species_order) {
      sigmat[gr[1], gr[2]] <- TRUE
      sigmat[gr[2], gr[1]] <- TRUE
    }
  }

  out <- setNames(rep(NA_character_, n), species_order)
  pool <- c(letters, LETTERS)
  k <- 1
  remaining <- species_order
  while (length(remaining) > 0) {
    L <- pool[k]; k <- k + 1
    group <- character(0)
    for (sp in remaining) {
      if (length(group) == 0 || all(!sigmat[sp, group]))
        group <- c(group, sp)
    }
    out[group] <- vapply(group, function(g) {
      if (is.na(out[g])) L else paste0(out[g], L)
    }, character(1))
    remaining <- setdiff(remaining, group)
  }
  out
}

# Order species by MC50 (descending) for the letter assignment
sp_order <- delta_tab$species[order(-delta_tab$mc50)]
letters_vec <- cld_from_sidak(sp_order, pw_tab)

cat("\n=== Compact letter display (Sidak-adjusted) ===\n")
print(letters_vec)

# --- 4. Stratified bootstrap (validation) ------------------------------------------

set.seed(SIM_SEED)

strata     <- interaction(df.mc$species, df.mc$phase, drop = TRUE)
idx_strata <- split(seq_len(nrow(df.mc)), strata)

boot_mat <- matrix(NA_real_, nrow = BOOT_REPS,
                   ncol = nlevels(df.mc$species),
                   dimnames = list(NULL, levels(df.mc$species)))
n_fit_fail <- 0

extract_mc50s <- function(fit_b) {
  cf  <- coef(fit_b)
  ref <- levels(df.mc$species)[1]
  vapply(levels(df.mc$species), function(sp) {
    A <- sum(cf[c("(Intercept)", if (sp != ref) paste0("species", sp))],
             na.rm = TRUE) + w2 * cf[[ph2_coef]]
    B <- sum(cf[c("mc", if (sp != ref) paste0("mc:species", sp))], na.rm = TRUE)
    g <- -A / B
    if (!is.finite(g) || B <= 0 || abs(g) > 500) return(NA_real_)
    g
  }, numeric(1))
}

for (b in seq_len(BOOT_REPS)) {
  idx_b <- unlist(lapply(idx_strata, function(ii)
    sample(ii, size = length(ii), replace = TRUE)))
  df_b <- df.mc[idx_b, ]

  fit_b <- try(suppressWarnings(
    glm(germinated ~ mc + phase + species + mc:species,
        family = binomial, data = df_b)), silent = TRUE)

  if (inherits(fit_b, "try-error") || !fit_b$converged) {
    n_fit_fail <- n_fit_fail + 1
    next
  }
  boot_mat[b, ] <- extract_mc50s(fit_b)

  if (b %% 250 == 0) cat("bootstrap progress:", b, "/", BOOT_REPS, "\n")
}

boot_tab <- bind_rows(lapply(levels(df.mc$species), function(sp) {
  vals <- boot_mat[, sp]
  ok   <- vals[is.finite(vals)]
  orig <- delta_tab$mc50[delta_tab$species == sp]
  tibble(species     = sp,
         boot_mean   = mean(ok),
         se_boot     = sd(ok),
         ci_lo       = quantile(ok, 0.025, names = FALSE),
         ci_hi       = quantile(ok, 0.975, names = FALSE),
         bias        = mean(ok) - orig,
         n_valid     = length(ok),
         pct_invalid = round(100 * (BOOT_REPS - length(ok)) / BOOT_REPS, 2))
}))

cat("\n=== Bootstrap MC50 (validation; stratified species x phase;",
    BOOT_REPS, "reps,", n_fit_fail, "non-convergent fits discarded) ===\n")
print(as.data.frame(boot_tab))

# --- 5. Side-by-side summary --------------------------------------------------------

delta_side <- delta_tab |>
  dplyr::select(species, mc50, se_delta, ci_lo_delta = ci_lo, ci_hi_delta = ci_hi)
boot_side <- boot_tab |>
  dplyr::select(species, boot_mean, se_boot,
         ci_lo_boot = ci_lo, ci_hi_boot = ci_hi,
         n_valid, pct_invalid)

summary_tab <- full_join(delta_side, boot_side, by = "species")

# Merge CLD letters (Sidak-adjusted pairwise contrasts) into the summary table
summary_tab <- summary_tab |>
  mutate(letters = unname(letters_vec[species]))

cat("\n=== MC50 comparison (%, marginalized over batches) ===\n")
print(as.data.frame(summary_tab))

# Drop the gradient (list-column) before writing to CSV; it is only needed
# in-memory for the pairwise covariance computations.
delta_tab_csv <- delta_tab |>
  dplyr::select(-grad)

write.csv(as.data.frame(delta_tab_csv),
          file.path(tablas_dir, "mc50_delta.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_tab),
          file.path(tablas_dir, "mc50_bootstrap.csv"), row.names = FALSE)
write.csv(as.data.frame(summary_tab),
          file.path(tablas_dir, "mc50_summary.csv"), row.names = FALSE)

# --- 6. Forest plot (article quality) ------------------------------------------------

img_dir <- "07-img/mc50_estimation"
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

plot_tab <- summary_tab |>
  mutate(species = factor(species, levels = levels(reorder(species, mc50))))

# Max CI upper bound per species, used to position the letters to the right
# of the error bars so they do not overlap with the graphical elements.
letter_pos <- plot_tab |>
  mutate(x_pos = ci_hi_delta + 1.0) |>          # 1 unit right of CI upper end
  dplyr::select(species, x_pos, letters)

p_mc50 <- ggplot(plot_tab, aes(x = mc50, y = species)) +
  # bootstrap percentile CI (validation layer, thin grey whiskers)
  geom_errorbar(aes(xmin = ci_lo_boot, xmax = ci_hi_boot),
                width = 0.14, colour = "grey55", linewidth = 0.35) +
  # delta CI (primary, bold black whiskers)
  geom_errorbar(aes(xmin = ci_lo_delta, xmax = ci_hi_delta),
                width = 0.28, colour = "black", linewidth = 1.1) +
  geom_point(shape = 21, size = 3.2, fill = "black", colour = "black") +
  # CLD letters (Sidak-adjusted pairwise 95% CIs)
  geom_text(data = letter_pos,
            aes(x = x_pos, y = species, label = letters),
            hjust = 0, size = 4, fontface = "bold") +
  scale_x_continuous(limits = c(24, 55), breaks = seq(25, 55, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    x = expression(MC[50] ~ "(%)"),
    y = "Species",
    caption = paste(
      "Bold whiskers: Delta-method 95% CI (primary).",
      "Thin grey whiskers: stratified-bootstrap percentile 95% CI (validation).",
      "Letters: groups from Sidak-adjusted pairwise 95% CIs on MC50",
      "(shared letter = no significant difference).",
      "MC50 marginalized over drying batches on the link scale.",
      sep = "\n")
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    plot.caption = element_text(size = 8.5, colour = "grey30",
                                hjust = 0, margin = margin(t = 6)),
    legend.position = "none"
  )

fig_path <- file.path(img_dir, "mc50_forest_species.png")
ggsave(fig_path, p_mc50, width = 140, height = 95, units = "mm", dpi = 600)
cat("Forest-plot figure written to", fig_path, "\n")

# --- 7. Article figure: MC50 by species, coloured by bioclimate --------------------

# Full names (mapping used across the repo, cf. d.05 scripts)
species_names <- c(CO = "Q. coccifera",
                   IL = "Q. ilex",
                   SU = "Q. suber",
                   FA = "Q. faginea",
                   PU = "Q. pubescens",
                   PY = "Q. pyrenaica",
                   PE = "Q. petraea")

bioclimate <- tribble(
  ~code, ~bioclimate,
  "CO",  "Mediterranean",
  "IL",  "Mediterranean",
  "SU",  "Mediterranean",
  "FA",  "Sub-Mediterranean",
  "PU",  "Sub-Mediterranean",
  "PY",  "Sub-Mediterranean",
  "PE",  "Temperate"
)

plot_bio <- summary_tab |>
  left_join(bioclimate, by = c("species" = "code")) |>
  mutate(name = unname(species_names[species])) |>
  # y order: Mediterranean (bottom) -> Sub-Mediterranean -> Temperate (top),
  # ascending MC50 within each bioclimate
  arrange(factor(bioclimate, levels = c("Mediterranean",
                                        "Sub-Mediterranean",
                                        "Temperate")), mc50) |>
  mutate(name = factor(name, levels = name))

# CLD letters positioned to the right of the CI upper bound
letter_pos_bio <- plot_bio |>
  mutate(x_pos = ci_hi_delta + 1.0) |>
  dplyr::select(name, x_pos, letters, bioclimate)

p_bio <- ggplot(plot_bio, aes(x = mc50, y = name, colour = bioclimate)) +
  geom_errorbar(aes(xmin = ci_lo_delta, xmax = ci_hi_delta),
                width = 0.28, linewidth = 1.1) +
  geom_point(size = 3.4) +
  # CLD letters (colour matches the bioclimate of each species)
  # show.legend = FALSE: keeps the letters out of the "Bioclimate" legend so
  # that only the colour key entries (Mediterranean/Sub-Mediterranean/Temperate)
  # appear there; the letters stay drawn on the plot itself.
  geom_text(data = letter_pos_bio,
            aes(x = x_pos, y = name, label = letters, colour = bioclimate),
            hjust = 0, size = 4, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(
    values = c("Mediterranean"      = "#D55E00",   # Okabe-Ito vermillion
               "Sub-Mediterranean"  = "#E69F00",   # Okabe-Ito orange
               "Temperate"          = "#0072B2"),  # Okabe-Ito blue
    name = "Bioclimate"
  ) +
  scale_x_continuous(limits = c(24, 55), breaks = seq(25, 55, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    x = expression(MC[50] ~ "(%)"),
    y = NULL,
    caption = paste0(
      "MC50: water content at 50% germination, ",
      "marginalized over drying batches on the link scale.\n",
      "Points and whiskers: Delta-method estimates with asymptotic 95% CIs.\n",
      "Letters: groups from Sidak-adjusted pairwise 95% CIs on MC50 ",
      "(shared letter = no significant difference).\n",
      "Colours denote the characteristic bioclimate of each species."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    axis.text.y = element_text(face = "italic", size = 11),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    legend.position = c(0.97, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey80",
                                     linewidth = 0.2),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 8.5, colour = "grey30",
                                hjust = 0, margin = margin(t = 6))
  )

fig_path <- file.path(img_dir, "mc50_forest_bioclimate.png")
ggsave(fig_path, p_bio, width = 140, height = 100, units = "mm", dpi = 600)
cat("Bioclimate forest-plot figure written to", fig_path, "\n")

cat("\ns.03 done: MC50 tables written to", tablas_dir,
    "| figure written to", img_dir, "\n")
