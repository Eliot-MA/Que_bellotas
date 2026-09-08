# d.03.3-bioclimate_report.R ---- patron bioclimatico en las dimensiones FAMD
# Requiere: d.03.0 (bioclimate_cols, dims) y d.03.1 (ind.bio)
if (!exists("bioclimate_cols")) source("01-scripts/d.03.0-utils.R")
if (!exists("ind.bio"))         source("01-scripts/d.03.1-load_famd.R")

# Reporte del patron bioclimatico en las dimensiones FAMD ----
## 1. Modelo mixto (bellotas anidadas en especie, bioclimate fijo) ----
library(glmmTMB)
library(car)

glmm.bio <- lapply(dims, function(dd) {
  glmmTMB::glmmTMB(reformulate("bioclimate + (1 | species)", response = dd),
                   data = ind.bio)
})
names(glmm.bio) <- dims

## 2. Tabla 1: test del efecto bioclimate por dimension (Wald chi²) + R² marginal ----
bioclimate_anova_mixed <- bind_rows(lapply(dims, function(dd) {
  a <- car::Anova(glmm.bio[[dd]])
  r <- a[row.names(a) == "bioclimate", , drop = FALSE]
  data.frame(
    dimension = dd,
    Chisq     = round(r[["Chisq"]], 2),
    df        = r[["Df"]],
    p_value   = r[["Pr(>Chisq)"]],
    marg_R2   = round(performance::r2(glmm.bio[[dd]])$R2_marginal, 3)
  )
}))

write.csv2(bioclimate_anova_mixed, "00-data/bioclimate_anova_mixed.csv", row.names = FALSE)

## 3. Tabla 2: EMMs por bioclima + letras CLD ----
bioclimate_emmeans_cld <- bind_rows(lapply(dims, function(dd) {
  emm <- emmeans::emmeans(glmm.bio[[dd]], "bioclimate")
  cf <- multcomp::cld(emm, Letters = letters) |>
    as.data.frame()
  low.col <- grep("LCL|lower.CL", names(cf), value = TRUE)
  up.col  <- grep("UCL|upper.CL", names(cf), value = TRUE)
  data.frame(dimension  = dd,
             bioclimate = cf$bioclimate,
             emmean     = cf$emmean,
             SE         = cf$SE,
             lower.CL   = cf[[low.col]],
             upper.CL   = cf[[up.col]],
             .group     = str_replace_all(cf$.group, "\\s+", ""))
}))

write.csv2(bioclimate_emmeans_cld, "00-data/bioclimate_emmeans_cld.csv", row.names = FALSE)

## 4. Tabla 3: medias (± SD) por especie y dimension ----
bioclimate_species_means <- ind.bio |>
  group_by(species, bioclimate) |>
  summarise(n = n(),
            across(all_of(dims),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        sd   = ~ sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop") |>
  arrange(factor(bioclimate, levels = levels(ind.bio$bioclimate)), species)

write.csv2(bioclimate_species_means, "00-data/bioclimate_species_means.csv", row.names = FALSE)

## 5. Tabla 4 (sensibilidad): ANOVA a nivel de especie (n = 8) ----
sp.means <- ind.bio |>
  group_by(species) |>
  summarise(across(all_of(dims), mean),
            bioclimate = unique(bioclimate),
            .groups = "drop")

bioclimate_sp_anova <- bind_rows(lapply(dims, function(dd) {
  a <- anova(lm(reformulate("bioclimate", response = dd), data = sp.means))
  r <- a["bioclimate", , drop = FALSE]
  data.frame(dimension = dd,
             F       = round(r[["F value"]], 2),
             df1     = r[["Df"]],
             df2     = a["Residuals", "Df"],
             p_value = r[["Pr(>F)"]])
}))

write.csv2(bioclimate_sp_anova, "00-data/bioclimate_species_level_anova.csv", row.names = FALSE)

## 6. Figura A: boxplot + violin por bioclima (Dim1-Dim3) con letras CLD ----
dim.long <- ind.bio |>
  pivot_longer(all_of(dims), names_to = "dimension", values_to = "value") |>
  mutate(dimension_label = factor(
    dplyr::recode(dimension,
           "Dim.1" = "Dim 1 (acorn size)",
           "Dim.2" = "Dim 2 (pericarp robustness)",
           "Dim.3" = "Dim 3 (scar / rupture)"),
    levels = c("Dim 1 (acorn size)",
               "Dim 2 (pericarp robustness)",
               "Dim 3 (scar / rupture)")))

let.df <- bioclimate_emmeans_cld |>
  mutate(dimension_label = factor(
    dplyr::recode(dimension,
           "Dim.1" = "Dim 1 (acorn size)",
           "Dim.2" = "Dim 2 (pericarp robustness)",
           "Dim.3" = "Dim 3 (scar / rupture)"),
    levels = c("Dim 1 (acorn size)",
               "Dim 2 (pericarp robustness)",
               "Dim 3 (scar / rupture)")))

p.bio.boxes <- ggplot(dim.long, aes(x = bioclimate, y = value)) +
  geom_violin(aes(fill = bioclimate), alpha = .25, colour = "grey45", linewidth = .3) +
  geom_boxplot(width = .18, alpha = .55, outlier.shape = NA,
               linewidth = .5, colour = "grey20", fill = "white") +
  geom_jitter(aes(colour = bioclimate), alpha = .12, width = .12, size = .7) +
  geom_text(data = let.df, aes(x = bioclimate, y = Inf, label = .group),
            vjust = 1.8, size = 4, colour = "black", fontface = "bold") +
  facet_wrap(~ dimension_label, scales = "free_y") +
  scale_colour_manual(values = bioclimate_cols) +
  scale_fill_manual(values = bioclimate_cols) +
  theme_minimal() +
  labs(x = NULL, y = "FAMD coordinate") +
  theme(legend.position = "none")

ggsave("07-img/FAMD_outputs/07_bioclimate_boxplots.png", p.bio.boxes, width = 9, height = 4.4)

## 7. Figura B: forest plot de EMMs por bioclima y dimension ----
p.bio.forest <- ggplot(bioclimate_emmeans_cld,
                       aes(x = emmean, y = dimension, colour = bioclimate)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL),
                  position = position_dodge2(width = 0.5, reverse = TRUE),
                  size = 0.4) +
  scale_colour_manual(values = bioclimate_cols) +
  scale_y_discrete(labels = c("Dim.1" = "Dim 1 (acorn size)",
                              "Dim.2" = "Dim 2 (pericarp robustness)",
                              "Dim.3" = "Dim 3 (scar / rupture)")) +
  theme_minimal() +
  labs(x = "Estimated marginal mean (95% CI)", y = NULL, colour = "Bioclimate")

ggsave("07-img/FAMD_outputs/08_bioclimate_emmeans_forest.png",
       p.bio.forest, width = 7, height = 4.5)

## 8. Resumen en consola ----
cat("\n===== hallazgos por bioclima (modelo mixto, 1|species) =====\n")
print(bioclimate_anova_mixed)
cat("\n--- EMMs + letras CLD ---\n")
print(bioclimate_emmeans_cld |>
        dplyr::select(dimension, bioclimate, emmean, lower.CL, upper.CL, .group))
cat("\n--- Sensibilidad a nivel de especie (n = 8) ---\n")
print(bioclimate_sp_anova)