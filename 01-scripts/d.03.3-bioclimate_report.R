# d.03.3-bioclimate_report.R ---- patron bioclimatico en las dimensiones FAMD
# Requiere: d.03.0 (bioclimate_cols, dims) y d.03.1 (ind.bio)
if (!exists("bioclimate_cols")) source("01-scripts/d.03.0-utils.R")
if (!exists("ind.bio"))         source("01-scripts/d.03.1-load_famd.R")

# Reporte del patron bioclimatico en las dimensiones FAMD ----
#
# ---------------------------------------------------------------------------
# NOTA DE DECISION (2026-10-02): el factor aleatorio es la PROCEDENCIA, no la
# especie. Motivos, todos verificados en este mismo script (secciones 2 y 6):
#
# 1. La unidad de muestreo real es la bellota dentro de una procedencia, no
#    dentro de una especie. La procedencia esta anidada en la especie, de modo
#    que (1|prov_code) ya absorbe toda la variabilidad interespecifica: los dos
#    terminos no son separables y el modelo que explicita el anidamiento
#    (1|species) + (1|species:prov_code) no mejora el ajuste (seccion 6).
# 2. El reparto de varianza (seccion 6) muestra que en Dim.1 la procedencia
#    dentro de especie explica ~46 % de la varianza frente al ~10 % de la
#    especie: el modelo con especie estaba atribuyendo a "especie" lo que en
#    realidad es variacion entre procedencias (PY2 y SU2, en particular).
# 3. El AIC favorece la procedencia en las tres dimensiones (seccion 6).
# 4. (1|prov_code) elimina los valores extremos de Dim.1 (DHARMa, seccion 2):
#    son bellotas de gran tamano de PY2 y SU2 que el modelo con especie dejaba
#    como residuos extremos. Son bellotas reales, no errores, y se conservan.
#
# CONSECUENCIA PARA LA INFERENCIA: el bioclima es un atributo de la especie, de
# modo que el tamano efectivo del contraste son 8 especies (15 procedencias), no
# 403 bellotas. Los test de Wald asintoticos de glmmTMB tratan las 403 filas como
# independientes y son optimistas. Por eso la seccion 7 anade un bootstrap por
# cluster de especie y un contraste no parametrico sobre medias de especie, que
# es la escala a la que el factor esta realmente definido. El manuscrito debe
# presentar el contraste robusto junto al asintotico.
# ---------------------------------------------------------------------------

library(glmmTMB)
library(car)
library(DHARMa)
library(e1071)   # asimetria y curtosis de los Q-residuos

dir.create("07-img/FAMD_outputs",  showWarnings = FALSE)
dir.create("07-img/FAMD_outputs/dharma", showWarnings = FALSE)
dir.create("08-reports/bioclimate", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 0. Utilidades
# ---------------------------------------------------------------------------
# DHARMa >= 0.4.0 renombro $Q.residuals a $scaledResiduals (cuantil empirico en
# [0, 1]). Se accede por el accessor para no depender de ninguna de las dos APIs.
qres <- function(sim) {
  r <- sim$Q.residuals
  if (is.null(r) || length(r) == 0) r <- sim$scaledResiduals
  r <- as.numeric(r)
  r[!is.finite(r)] <- NA_real_
  r
}
# version comparable con una normal, para QQ y asimetria
qres_z <- function(sim) stats::qnorm(pmin(pmax(qres(sim), 1e-6), 1 - 1e-6))
# evaluacion segura: devuelve NA en lugar de abortar el pipeline
try2 <- function(expr) tryCatch(expr, error = function(e) NA)

# Varianzas de los efectos aleatorios. Segun la version, VarCorr(m)$cond es una
# matriz con atributos "stddev"/"correlation" o una lista de bloques; se leen las
# diagonales de los bloques de covarianza para no depender de una sola API.
re_var <- function(m) {
  vc <- glmmTMB::VarCorr(m)
  cond <- if (!is.null(vc[["cond"]])) vc[["cond"]] else vc
  blocks <- cond
  if (is.matrix(cond) || is.array(cond)) blocks <- list(cond)
  out <- unlist(lapply(blocks, function(b) {
    b <- as.matrix(b)
    d <- diag(b)
    d[is.finite(d) & d > 0]
  }))
  if (!length(out)) NA_real_ else max(out)
}

# ---------------------------------------------------------------------------
# 1. Modelo mixto (bellotas anidadas en procedencia, bioclimate fijo)
# ---------------------------------------------------------------------------
glmm.bio <- lapply(dims, function(dd) {
  glmmTMB::glmmTMB(stats::as.formula(paste(dd, "~ bioclimate + (1 | prov_code)")),
                   data = ind.bio)
})
names(glmm.bio) <- dims

# ---------------------------------------------------------------------------
# 2. Comprobacion de residuos (DHARMa)
#    Residuos cuantiles simulados, que son el test valido para un modelo mixto;
#    los residuos de Pearson no lo son porque no incorporan el efecto aleatorio.
#    Se testea dispersion (homocedasticidad), uniformidad (normalidad) y
#    valores extremos.
# ---------------------------------------------------------------------------
cat("\n===== 2. Comprobacion de residuos (DHARMa) =====\n")

bioclimate_resid_check <- bind_rows(lapply(dims, function(dd) {
  sim <- DHARMa::simulateResiduals(glmm.bio[[dd]], n = 2000, seed = 1)
  z   <- qres_z(sim)
  data.frame(
    dimension        = dd,
    p_dispersion     = try2(DHARMa::testDispersion(sim)$p.value),
    p_uniformidad    = try2(DHARMa::testUniformity(sim)$p.value),
    p_outliers       = try2(DHARMa::testOutliers(sim)$p.value),
    n_outliers       = sum(qres(sim) <= 0.005 | qres(sim) >= 0.995, na.rm = TRUE),
    asimetria        = round(mean(e1071::skewness(z)), 3),
    exceso_curtosis  = round(e1071::kurtosis(z), 3)
  )
}))

cat("\n"); print(as.data.frame(bioclimate_resid_check), row.names = FALSE)
write.csv2(bioclimate_resid_check, "08-reports/bioclimate/bioclimate_dharma_residuals.csv",
           row.names = FALSE)

# Q-residuos medios por procedimiento: donde se concentra la desviacion
cat("\n-- Q-residuos (media, sd) por procedimiento --\n")
for (dd in dims) {
  z  <- qres_z(DHARMa::simulateResiduals(glmm.bio[[dd]], n = 2000, seed = 1))
  tb <- do.call(rbind, tapply(z, ind.bio$prov_code,
                              function(v) c(n = length(v), mean = mean(v), sd = sd(v))))
  cat("\n", dd, "\n", sep = ""); print(round(tb, 3))
}

# Figura de diagnostico
bioclimate_cols_d <- bioclimate_cols[as.character(ind.bio$bioclimate)]

# El PDF puede estar abierto en un visor, que en Windows bloquea la escritura y
# abortaria todo el pipeline. Se cae a un nombre con sufijo en ese caso.
dharma_pdf <- "07-img/FAMD_outputs/dharma/dharma_famd_dims.pdf"
dharma_dev <- tryCatch(grDevices::pdf(dharma_pdf, width = 12, height = 16),
                       error = function(e) NULL)
if (is.null(dharma_dev)) {
  dharma_pdf <- paste0(tools::file_path_sans_ext(dharma_pdf), "_",
                       format(Sys.time(), "%H%M%S"), ".pdf")
  dharma_dev <- grDevices::pdf(dharma_pdf, width = 12, height = 16)
  cat("AVISO: el PDF anterior esta bloqueado; se escribe en", dharma_pdf, "\n")
}
par(mfrow = c(3, 2))
for (dd in dims) {
  sim <- DHARMa::simulateResiduals(glmm.bio[[dd]], n = 2000, seed = 1)
  z   <- qres_z(sim)
  plot(sim, main = paste("DHARMa -", dd), xlab = "Predichos (log)",
       ylab = "Residuo cuantil", col = scales::alpha(bioclimate_cols_d, 0.30))
  points(z, ppoints(seq_along(z)), pch = 20, cex = 0.6, col = "grey30")
  abline(h = c(-2, 2), lty = 2, col = "grey60")
  hist(z, breaks = seq(min(z, -3), max(z, 3), length.out = 25),
       main = paste("Q-residuos (z) -", dd), xlab = NULL, col = "grey85", border = "white")
  qqnorm(z, main = paste("QQ -", dd)); qqline(z)
  boxplot(z ~ ind.bio$bioclimate, las = 1, ylab = "Q-residuo (z)",
          main = paste("~ bioclima -", dd))
  boxplot(z ~ ind.bio$prov_code, las = 2, ylab = "Q-residuo (z)", cex.axis = 0.6,
          main = paste("~ procedencia -", dd))
}
par(mfrow = c(1, 1))
dev.off()
cat("\nCharts ->", dharma_pdf, "\n")

# ---------------------------------------------------------------------------
# 3. Tabla 1: test del efecto bioclimate por dimension (Wald chi²) + R² marginal
# ---------------------------------------------------------------------------
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
cat("\n===== 3. Test del efecto bioclimate (Wald) =====\n")
print(as.data.frame(bioclimate_anova_mixed), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Tabla 2: EMMs por bioclima + letras CLD
# ---------------------------------------------------------------------------
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
cat("\n===== 4. EMMs por bioclima =====\n")
print(as.data.frame(bioclimate_emmeans_cld |>
                      dplyr::select(dimension, bioclimate, emmean, lower.CL,
                                    upper.CL, .group) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), function(x) round(x, 3)))))

# ---------------------------------------------------------------------------
# 5. Tabla 3: medias por especie y dimension
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 5b. Tabla 3b: medias por PROCEDENCIA (misma unidad del modelo)
# ---------------------------------------------------------------------------
bioclimate_prov_means <- ind.bio |>
  group_by(prov_code, species, bioclimate) |>
  summarise(n = n(),
            across(all_of(dims), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

write.csv2(bioclimate_prov_means, "08-reports/bioclimate/bioclimate_provenance_means.csv",
           row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Estructura de varianza: justificacion de la nota de decision
#    Compara (1|species), (1|prov_code) y el anidamiento explicito.
# ---------------------------------------------------------------------------
cat("\n===== 6. Estructura de varianza (justificacion de (1|prov_code)) =====\n")

var_specs <- list(
  "(1|species)"     = "(1|species)",
  "(1|prov_code)"   = "(1|prov_code)",
  "anidado"         = "(1|species) + (1|species:prov_code)")

variance_table <- list()
for (dd in dims) {
  for (nm in names(var_specs)) {
    m <- tryCatch(glmmTMB::glmmTMB(stats::as.formula(paste(dd, "~ 1 +", var_specs[[nm]])),
                                   data = ind.bio), error = function(e) NULL)
    if (is.null(m)) next
    v_rd <- try2(max(re_var(m)))
    variance_table[[paste(dd, nm, sep = " | ")]] <- data.frame(
      dimension = dd, estructura = nm, AIC = round(try2(AIC(m)), 1),
      logLik = round(try2(as.numeric(logLik(m))), 2),
      var_random = round(v_rd, 4),
      # glmmTMB ya no exporta isSingular(): una varianza aleatoria que colapsa a
      # ~0 indica que ese termino no es identificable en el modelo.
      singular = if (is.na(v_rd)) NA else v_rd / sigma(m)^2 < 1e-4,
      stringsAsFactors = FALSE)
  }
}
variance_table <- do.call(rbind, variance_table)
cat("\n"); print(variance_table, row.names = FALSE)
write.csv2(variance_table, "08-reports/bioclimate/bioclimate_variance_structures.csv",
           row.names = FALSE)

# Reparto porcentual de la varianza con el anidamiento explicito
cat("\n-- Reparto de varianza (modelo anidado especie:procedencia) --\n")
variance_share <- bind_rows(lapply(dims, function(dd) {
  m <- glmmTMB::glmmTMB(stats::as.formula(paste(dd, "~ 1 + (1|species) + (1|species:prov_code)")),
                        data = ind.bio)
  v_sp <- as.numeric(VarCorr(m)$cond$species[1])
  v_pv <- as.numeric(VarCorr(m)$cond$`species:prov_code`[1])
  s_b  <- sigma(m)
  tot  <- v_sp + v_pv + s_b
  data.frame(dimension = dd,
             var_especie = round(v_sp, 3), var_proc_dentro = round(v_pv, 3),
             var_bellota = round(s_b, 3),
             pct_especie = round(100 * v_sp / tot, 1),
             pct_proc_dentro = round(100 * v_pv / tot, 1),
             pct_bellota = round(100 * s_b / tot, 1))
}))
cat("\n"); print(as.data.frame(variance_share), row.names = FALSE)
write.csv2(variance_share, "08-reports/bioclimate/bioclimate_variance_share.csv",
           row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Sensibilidad: el bioclima esta definido a nivel de ESPECIE
#    (a) ANOVA de una via sobre medias por especie (n = 8)
#    (b) Kruskal-Wallis sobre esas mismas medias
#    (c) bootstrap por cluster de especie (B = 2000) del contraste Wald
# ---------------------------------------------------------------------------
cat("\n===== 7. Sensibilidad a la escala de la unidad =====\n")

# Medias por especie: la escala a la que el bioclima esta realmente definido.
sp.means <- ind.bio |>
  group_by(species) |>
  summarise(across(all_of(dims), mean),
            bioclimate = factor(first(bioclimate),
                                levels = levels(ind.bio$bioclimate)),
            .groups = "drop")

bioclimate_sp_anova <- bind_rows(lapply(dims, function(dd) {
  a <- anova(lm(stats::as.formula(paste(dd, "~ bioclimate")), data = sp.means))
  r <- a["bioclimate", , drop = FALSE]
  data.frame(dimension = dd,
             F       = round(r[["F value"]], 2),
             df1     = r[["Df"]],
             df2     = a["Residuals", "Df"],
             p_value = r[["Pr(>F)"]])
}))

write.csv2(bioclimate_sp_anova, "00-data/bioclimate_species_level_anova.csv", row.names = FALSE)

bioclimate_sp_kruskal <- bind_rows(lapply(dims, function(dd) {
  k <- kruskal.test(stats::as.formula(paste(dd, "~ bioclimate")), data = sp.means)
  data.frame(dimension = dd, statistic = round(unname(k$statistic), 2),
             df = k$parameter, p_value = k$p.value)
}))
write.csv2(bioclimate_sp_kruskal, "08-reports/bioclimate/bioclimate_species_level_kruskal.csv",
           row.names = FALSE)

# Bootstrap por cluster de especie. El re-muestreo es de ESPECIES porque el
# bioclima es constante dentro de especie: es la unidad a la que el factor
# esta definido, y las bellotas de una misma especie no son independientes.
B_boot   <- 2000
set.seed(11)
sp_all   <- unique(ind.bio$species)
boot_p <- matrix(NA_real_, nrow = B_boot, ncol = length(dims),
                 dimnames = list(NULL, dims))
for (dd in dims) {
  for (b in seq_len(B_boot)) {
    d <- ind.bio[ind.bio$species %in% sample(sp_all, length(sp_all), replace = TRUE), ]
    d$species <- factor(d$species)
    mb <- tryCatch(glmmTMB::glmmTMB(stats::as.formula(paste(dd, "~ bioclimate + (1|prov_code)")),
                                    data = d), error = function(e) NULL)
    boot_p[b, dd] <- if (is.null(mb)) NA else try2(car::Anova(mb)$`Pr(>Chisq)`[1])
  }
}

bioclimate_boot_species <- bind_rows(lapply(dims, function(dd) {
  v <- boot_p[, dd]
  data.frame(dimension = dd,
             p_asintotico = round(bioclimate_anova_mixed$p_value[
               bioclimate_anova_mixed$dimension == dd], 4),
             p_boot_mediana = round(median(v, na.rm = TRUE), 4),
             p_boot_ic_low = round(unname(quantile(v, 0.025, na.rm = TRUE)), 4),
             p_boot_ic_high = round(unname(quantile(v, 0.975, na.rm = TRUE)), 4),
             prop_p_menor_005 = round(mean(v < 0.05, na.rm = TRUE), 3),
             B = B_boot)
}))
write.csv2(bioclimate_boot_species, "08-reports/bioclimate/bioclimate_bootstrap_species.csv",
           row.names = FALSE)

cat("\n-- ANOVA / Kruskal-Wallis sobre medias de especie (n = 8) --\n")
print(as.data.frame(bioclimate_sp_anova), row.names = FALSE)
print(as.data.frame(bioclimate_sp_kruskal), row.names = FALSE)
cat("\n-- Bootstrap por cluster de especie (B =", B_boot, ") --\n")
print(as.data.frame(bioclimate_boot_species), row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. Figura A: boxplot + violin por bioclima (Dim1-Dim3) con letras CLD
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 9. Figura B: forest plot de EMMs por bioclima y dimension
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 10. Resumen en consola
# ---------------------------------------------------------------------------
cat("\n===== resumen por bioclima (modelo mixto, 1|prov_code) =====\n")
print(as.data.frame(bioclimate_anova_mixed), row.names = FALSE)
cat("\n--- EMMs + letras CLD ---\n")
print(as.data.frame(bioclimate_emmeans_cld |>
                      dplyr::select(dimension, bioclimate, emmean, lower.CL,
                                    upper.CL, .group) |>
                      dplyr::mutate(dplyr::across(where(is.numeric), round, 3))))
cat("\n--- comprobacion de residuos (DHARMa) ---\n")
print(as.data.frame(bioclimate_resid_check), row.names = FALSE)
cat("\n--- sensibilidad a la escala de unidad (especie, n = 8) ---\n")
print(as.data.frame(bioclimate_sp_kruskal), row.names = FALSE)
print(as.data.frame(bioclimate_boot_species), row.names = FALSE)
