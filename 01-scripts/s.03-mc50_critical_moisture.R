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
           ci_hi    = g + qnorm(0.975) * se)
  }))
}

delta_tab <- mc50_delta(fit_mc50, df.mc, w2)

cat("\n=== Delta-method MC50 (primary) ===\n")
print(as.data.frame(delta_tab))

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
  select(species, mc50, se_delta, ci_lo_delta = ci_lo, ci_hi_delta = ci_hi)
boot_side <- boot_tab |>
  select(species, boot_mean, se_boot,
         ci_lo_boot = ci_lo, ci_hi_boot = ci_hi,
         n_valid, pct_invalid)

summary_tab <- full_join(delta_side, boot_side, by = "species")

cat("\n=== MC50 comparison (%, marginalized over batches) ===\n")
print(as.data.frame(summary_tab))

write.csv(as.data.frame(delta_tab),
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

p_mc50 <- ggplot(plot_tab, aes(x = mc50, y = species)) +
  # bootstrap percentile CI (validation layer, thin grey whiskers)
  geom_errorbar(aes(xmin = ci_lo_boot, xmax = ci_hi_boot),
                width = 0.14, colour = "grey55", linewidth = 0.35) +
  # delta CI (primary, bold black whiskers)
  geom_errorbar(aes(xmin = ci_lo_delta, xmax = ci_hi_delta),
                width = 0.28, colour = "black", linewidth = 1.1) +
  geom_point(shape = 21, size = 3.2, fill = "black", colour = "black") +
  scale_x_continuous(limits = c(24, 50), breaks = seq(25, 50, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    x = expression(MC[50] ~ "(%)"),
    y = "Species",
    caption = paste("Bold whiskers: Delta-method 95% CI (primary).\n",
                    "Thin grey whiskers: stratified-bootstrap percentile 95% CI (validation).\n",
                    "MC50 marginalized over drying batches on the link scale.")
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

p_bio <- ggplot(plot_bio, aes(x = mc50, y = name, colour = bioclimate)) +
  geom_errorbar(aes(xmin = ci_lo_delta, xmax = ci_hi_delta),
                width = 0.28, linewidth = 1.1) +
  geom_point(size = 3.4) +
  scale_colour_manual(
    values = c("Mediterranean"      = "#D55E00",   # Okabe-Ito vermillion
               "Sub-Mediterranean"  = "#E69F00",   # Okabe-Ito orange
               "Temperate"          = "#0072B2"),  # Okabe-Ito blue
    name = "Bioclimate"
  ) +
  scale_x_continuous(limits = c(24, 50), breaks = seq(25, 50, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    x = expression(MC[50] ~ "(%)"),
    y = NULL,
    caption = paste0(
      "MC50: water content at 50% germination, \n 
      marginalized over drying batches on the link scale.\n",
      "Points and whiskers: Delta-method estimates with asymptotic 95% CIs.\n",
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
