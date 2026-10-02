# ============================================================
# d.05.1.heterogeneus_effects.R
# Analisis EXPLORATORIO de la heterogeneidad interespecifica en el efecto de
# los rasgos de la bellota sobre la tasa de desecacion.
#
# NO es el modelo de referencia. El modelo de referencia es
# d.05.3.reference_model.R (glmmTMB, con pendientes aleatorias por especie y
# sin filogenia); las figuras de efecto marginal sobre la pendiente se
# generan en d.05.5.slope_effects_figures.R.
#
# Este script justifica que la heterogeneidad por especie es real y no un
# artefacto. Los modelos con interacciones triples por especie (fit_het) son
# la version glmmTMB de la exploracion bayesiana de d.05.4_phylo_check.R:
# informativos sobre la heterogeneidad, pero no discutibles en el manuscrito.
#
# Reproduce y exporta el analisis de
#   08-reports/Heterogeneus_effects_acorn_traits.qmd
# Genera tablas en 00-data/processed/ y figuras en 07-img/:
#   - distribucion de los ejes del FAMD por especie (figura)
#   - R2 de cada eje explicado por la especie (tabla)
#   - correlacion entre ejes: global y dentro de especie (tablas)
#   - correlacion entre rasgos funcionales originales dentro de especie
#   - comparacion de modelos con/sin heterogeneidad por especie (tabla AIC/pesos)
#   - comparacion modelo "ingenuo" vs especie en random (tabla: cambio de signo)
#   - forest plot de contrastes (Delta slope Dim.alto - Dim.bajo) por especie
#   - relacion Dim.2 (pericarpo) vs contenido hidrico inicial (figura)
#   - objetos de modelos en 00-data/processed/models/
#
# La descomposicion DENTRO/ENTRE especies (within-between) se elimino de los
# scripts: daba conclusiones equivalentes a otros modelos que ya contemplan la
# heterogeneidad. Se conserva solo en los informes
# (Heterogeneus_effects_acorn_traits.qmd y d.07_brms_filogenias.qmd).
#
# Requiere en el entorno (los crea d.05.model_traits.R):
#   df, df.traits, df.t1, df.t2
# Si se ejecuta en solitario, los reconstruye con la misma receta.
# ============================================================

GEN_DASH <- "N"

if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(parameters)
  library(performance)
  library(emmeans)
  library(easystats)
})
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
suppressPackageStartupMessages(library(patchwork))

# Datos: reconstruccion en solitario por si no se ejecuta via d.05
if (!exists("df") || !exists("df.traits") || !exists("df.t1") || !exists("df.t2")) {
  # Procedencias excluidas declaradas en 01-scripts/00-config_procedencias.R
  if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")
  df.bellotas <- read.csv("00-data/processed/desiccation_traits_long.csv")
  df.famd     <- read.csv("00-data/processed/famd_ind_coord.csv")
  df <- df.bellotas |>
    dplyr::select(-X) |>
    dplyr::select(id_bellota, prov_code, tiempo_acumulado_horas, Moisture_content) |>
    # famd_ind_coord.csv tambien trae `prov_code`; se descarta para no duplicar
    # el nombre al hacer el join. Se conserva el de la tabla larga, previo al
    # filtro del FAMD.
    left_join(y = df.famd |> dplyr::select(-prov_code), by = "id_bellota") |>
    dplyr::filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
    tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
    rename(time = tiempo_acumulado_horas) |>
    mutate(
      time_s     = as.vector(scale(time)),
      species    = factor(species),
      provenance = factor(provenance),
      id_bellota = factor(id_bellota)
    )
  df.traits <- df |>
    dplyr::select(id_bellota, species, provenance, Dim.1, Dim.2, Dim.3) |>
    distinct()
  t94  <- as.vector((94 - mean(df$time)) / sd(df$time))
  df.t1 <- df |> filter(time_s < t94)
  df.t2 <- df |> filter(time_s > t94)
}

# El filtro se aplica antes de la particion temporal, de modo que PRE y POST
# comparten el mismo conjunto de procedencias.
assert_sin_procedencias_excluidas(df,    "prov_code", "d.05.1 datos")
assert_sin_procedencias_excluidas(df.t1, "prov_code", "d.05.1 fase PRE")
assert_sin_procedencias_excluidas(df.t2, "prov_code", "d.05.1 fase POST")

source("01-scripts/00-export_helpers.R")
# dir.create("07-img", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Variacion interespecifica de los ejes del FAMD
# ============================================================
p1 <- ggplot(df.traits, aes(y = species, x = Dim.1)) +
  geom_boxplot() + geom_violin(alpha = .4) + geom_jitter(alpha = .2) +
  labs(x = "Dim.1 (size)")
p2 <- ggplot(df.traits, aes(y = species, x = Dim.2)) +
  geom_boxplot() + geom_violin(alpha = .4) + geom_jitter(alpha = .2) +
  labs(x = "Dim.2 (pericarp)")
p3 <- ggplot(df.traits, aes(y = species, x = Dim.3)) +
  geom_boxplot() + geom_violin(alpha = .4) + geom_jitter(alpha = .2) +
  labs(x = "Dim.3 (scar)")

ggsave("07-img/heterogeneity_axes_distribution.png",
       p1 + p2 + p3 +
         plot_annotation(title = "Distribution of FAMD axes by species"),
       width = 13, height = 6, dpi = 300)
cat("Figura guardada: 07-img/heterogeneity_axes_distribution.png\n")

# R2 de cada eje explicado por la especie
lm.d1 <- lm(Dim.1 ~ species, data = df.traits)
lm.d2 <- lm(Dim.2 ~ species, data = df.traits)
lm.d3 <- lm(Dim.3 ~ species, data = df.traits)

axis_r2 <- tibble(
  Dim = c("Dim.1", "Dim.2", "Dim.3"),
  R2  = sapply(list(lm.d1, lm.d2, lm.d3), function(m) performance::r2(m)$R2)
)
write.csv(axis_r2, "00-data/processed/heterogeneity_axis_r2_species.csv", row.names = FALSE)
cat("Guardada: 00-data/processed/heterogeneity_axis_r2_species.csv\n")

# ============================================================
# 2. Correlacion entre ejes: global y dentro de especie
# ============================================================
cor.global <- df.traits |>
  dplyr::select(Dim.1, Dim.2, Dim.3) |>
  distinct() |>
  cor()
write.csv(cor.global, "00-data/processed/heterogeneity_cor_global.csv")

cor.species <- df |>
  group_by(species) |>
  group_modify(~{
    cormat <- cor(dplyr::select(.x, Dim.1, Dim.2, Dim.3))
    tibble(
      cor12 = cormat["Dim.1", "Dim.2"],
      cor13 = cormat["Dim.1", "Dim.3"],
      cor23 = cormat["Dim.2", "Dim.3"]
    )
  }) |>
  ungroup()
write.csv(cor.species, "00-data/processed/heterogeneity_cor_species.csv", row.names = FALSE)
cat("Guardadas: 00-data/processed/heterogeneity_cor_global.csv y 00-data/processed/heterogeneity_cor_species.csv\n")

# ============================================================
# 2b. Correlacion entre rasgos funcionales originales dentro de especies
# ============================================================
# Leer datos originales con rasgos funcionales
# Se releen del CSV largo (no del FAMD) porque se usan los rasgos en bruto, pero
# se aplican las MISMAS exclusions que el resto de la linea d.*: sin este
# filtro las correlaciones por especie incluian las bellotas de IL3.
df.traits_orig <- read.csv("00-data/processed/desiccation_traits_long.csv") |>
  dplyr::select(id_bellota, especie, prov_code, peso_seco, Volumen_estimado_cm3,
                Relacion_SV, SPM_g_cm2, Seed_Coat_Ratio, 
                Ratio_A.cicatriz_A.bellota, rajas_pericarpo) |>
  filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
  rename(species = especie,
         mass = peso_seco,
         volume = Volumen_estimado_cm3,
         SVR = Relacion_SV,
         SPM = SPM_g_cm2,
         SCR = Seed_Coat_Ratio,
         SSR = Ratio_A.cicatriz_A.bellota,
         pericarp_rupture = rajas_pericarpo) |>
  mutate(pericarp_rupture = as.numeric(as.character(pericarp_rupture))) |>
  distinct()

assert_sin_procedencias_excluidas(df.traits_orig, "prov_code",
                                   "d.05.1 rasgos originales")

species_order <- c("Quercus coccifera", "Quercus ilex", "Quercus suber",
                   "Quercus faginea", "Quercus pyrenaica", "Quercus pubescens",
                   "Quercus petraea", "Quercus robur")

# Calcular correlaciones entre rasgos funcionales por especie
cor_traits_orig <- df.traits_orig |>
  group_by(species) |>
  group_modify(~ {
    d <- dplyr::select(.x, mass, volume, SVR, SPM, SCR, SSR, pericarp_rupture)
    mat <- cor(d, use = "pairwise.complete.obs")
    # Extraer pares unicos (triangular inferior)
    pairs <- combn(colnames(mat), 2, simplify = FALSE)
    results <- lapply(pairs, function(p) {
      ct <- cor.test(d[[p[1]]], d[[p[2]]])
      tibble(
        trait1 = p[1],
        trait2 = p[2],
        cor    = ct$estimate,
        p      = ct$p.value
      )
    })
    bind_rows(results)
  }) |>
  ungroup() |>
  mutate(sig = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  ))

write.csv(cor_traits_orig, "00-data/processed/heterogeneity_cor_original_traits_by_species.csv", row.names = FALSE)
cat("Guardada: 00-data/processed/heterogeneity_cor_original_traits_by_species.csv\n")

# Grafico: heatmap de correlaciones entre rasgos originales por especie
cor_plot_orig <- cor_traits_orig |>
  mutate(
    pair = paste(trait1, trait2, sep = " - "),
    species = factor(species, levels = species_order)
  )

p_cor_traits_orig <- ggplot(cor_plot_orig, aes(x = pair, y = species, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(cor, 2), sig)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "r") +
  labs(x = "Trait pair", y = NULL,
       title = "Correlation among original functional traits, by species",
       caption = "* p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"))

ggsave("07-img/heterogeneity_cor_original_traits_heatmap.png", p_cor_traits_orig, width = 12, height = 7, dpi = 300)
cat("Figura guardada: 07-img/heterogeneity_cor_original_traits_heatmap.png\n")

# ============================================================
# 3. Efecto heterogeneo de los rasgos: modelos y comparacion
# ============================================================
# 3a. Modelos base: sin especie ("ingenuo") vs especie en el random
mm.pre.0 <- glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
                      (time_s | id_bellota), data = df.t1)
mm.pre.1 <- glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
                      (0 + time_s | species) +
                      (time_s | provenance:species) +
                      (time_s | id_bellota), data = df.t1)

mm.post.0 <- glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
                       (time_s | id_bellota), data = df.t2)
mm.post.1 <- glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
                       (0 + time_s | species) +
                       (time_s | provenance:species) +
                       (time_s | id_bellota), data = df.t2)

# 3b. Modelos que permiten un efecto heterogeneo del rasgo por especie
fit_het <- function(dat, dimv = NULL) {
  if (is.null(dimv)) {
    form <- paste0(
      "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + ",
      "time_s:Dim.1:species + time_s:Dim.2:species + time_s:Dim.3:species + ",
      "(time_s | provenance:species) + (time_s | id_bellota)"
    )
  } else {
    form <- paste0(
      "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + ",
      "time_s * ", dimv, " * species + ",
      "(time_s | provenance:species) + (time_s | id_bellota)"
    )
  }
  glmmTMB(as.formula(form), data = dat)
}

mm.pre.d1.het <- fit_het(df.t1, "Dim.1")
mm.pre.d2.het <- fit_het(df.t1, "Dim.2")
mm.pre.d3.het <- fit_het(df.t1, "Dim.3")
m.pre.all.het <- fit_het(df.t1)

mm.post.d1.het <- fit_het(df.t2, "Dim.1")
mm.post.d2.het <- fit_het(df.t2, "Dim.2")
mm.post.d3.het <- fit_het(df.t2, "Dim.3")
m.post.all.het <- fit_het(df.t2)

# 3c. Comparacion por AIC y pesos de Akaike
add_weights <- function(cmp, fase) {
  cmp |>
    as.data.frame() |>
    dplyr::mutate(DeltaAIC = AIC - min(AIC, na.rm = TRUE)) |>
    dplyr::mutate(PesoAkaike = exp(-0.5 * DeltaAIC) / sum(exp(-0.5 * DeltaAIC), na.rm = TRUE),
                  fase = fase) |>
    dplyr::arrange(DeltaAIC)
}

cmp_pre  <- add_weights(performance::compare_performance(
  mm.pre.0, mm.pre.1, mm.pre.d1.het, mm.pre.d2.het, mm.pre.d3.het, m.pre.all.het),
  "PRE (t < 94 h)")
cmp_post <- add_weights(performance::compare_performance(
  mm.post.0, mm.post.1, mm.post.d1.het, mm.post.d2.het, mm.post.d3.het, m.post.all.het),
  "POST (t > 94 h)")

cmp_all <- bind_rows(cmp_pre, cmp_post)
write.csv(cmp_all, "00-data/processed/heterogeneity_model_comparison.csv", row.names = FALSE)
cat("Guardada: 00-data/processed/heterogeneity_model_comparison.csv\n")

# 3d. Comparacion del modelo "ingenuo" vs especie en random (cambio de signo)
naive_compare <- bind_rows(
  parameters::model_parameters(mm.pre.0,  ci = 0.95) |>
    as.data.frame() |>
    dplyr::mutate(modelo = "PRE - simple (sin especie)", fase = "PRE"),
  parameters::model_parameters(mm.pre.1,  ci = 0.95) |>
    as.data.frame() |>
    dplyr::mutate(modelo = "PRE - especie en random",  fase = "PRE"),
  parameters::model_parameters(mm.post.0, ci = 0.95) |>
    as.data.frame() |>
    dplyr::mutate(modelo = "POST - simple (sin especie)", fase = "POST"),
  parameters::model_parameters(mm.post.1, ci = 0.95) |>
    as.data.frame() |>
    dplyr::mutate(modelo = "POST - especie en random",  fase = "POST")
)
write.csv(naive_compare, "00-data/processed/heterogeneity_naive_vs_random.csv", row.names = FALSE)
cat("Guardada: 00-data/processed/heterogeneity_naive_vs_random.csv\n")

# 3e. Forest plot: cambio de pendiente (Dim.alto - Dim.bajo) por especie
D1.values.species <- df.traits |>
  group_by(species) |>
  summarise(low = quantile(Dim.1, .1), high = quantile(Dim.1, .9), .groups = "drop")
D2.values.species <- df.traits |>
  group_by(species) |>
  summarise(low = quantile(Dim.2, .1), high = quantile(Dim.2, .9), .groups = "drop")
D3.values.species <- df.traits |>
  group_by(species) |>
  summarise(low = quantile(Dim.3, .1), high = quantile(Dim.3, .9), .groups = "drop")

D.global <- df.traits |>
  summarise(Dim.1 = median(Dim.1, na.rm = TRUE),
            Dim.2 = median(Dim.2, na.rm = TRUE),
            Dim.3 = median(Dim.3, na.rm = TRUE))

calc_species_contrasts <- function(model, dimension, values_species, values_global) {
  stopifnot(dimension %in% c("Dim.1", "Dim.2", "Dim.3"))
  resultados <- lapply(seq_len(nrow(values_species)), function(i) {
    sp   <- values_species$species[i]
    low  <- values_species$low[i]
    high <- values_species$high[i]

    at_values <- list(
      species = sp,
      Dim.1 = values_global$Dim.1,
      Dim.2 = values_global$Dim.2,
      Dim.3 = values_global$Dim.3
    )
    at_values[[dimension]] <- c(low, high)

    tr <- emmeans::emtrends(model, specs = dimension, var = "time_s", at = at_values)
    contr <- emmeans::contrast(tr, method = "revpairwise", adjust = "none")

    as.data.frame(contr) |>
      dplyr::left_join(
        as.data.frame(confint(contr)) |>
          dplyr::select(contrast, asymp.LCL, asymp.UCL),
        by = "contrast"
      ) |>
      dplyr::mutate(species = sp, dimension = dimension, low = low, high = high)
  })
  dplyr::bind_rows(resultados)
}

models_het <- list(
  PRE_Dim1  = mm.pre.d1.het,  PRE_Dim2  = mm.pre.d2.het,  PRE_Dim3  = mm.pre.d3.het,
  POST_Dim1 = mm.post.d1.het, POST_Dim2 = mm.post.d2.het, POST_Dim3 = mm.post.d3.het
)
model_info <- tibble(
  model     = names(models_het),
  phase     = rep(c("PRE", "POST"), each = 3),
  dimension = rep(c("Dim.1", "Dim.2", "Dim.3"), 2)
)
values_species <- list(
  "Dim.1" = D1.values.species,
  "Dim.2" = D2.values.species,
  "Dim.3" = D3.values.species
)

all_contrasts <- lapply(seq_len(nrow(model_info)), function(i) {
  calc_species_contrasts(
    model            = models_het[[model_info$model[i]]],
    dimension        = model_info$dimension[i],
    values_species   = values_species[[model_info$dimension[i]]],
    values_global    = D.global
  ) |>
    dplyr::mutate(phase = model_info$phase[i], model = model_info$model[i])
}) |>
  dplyr::bind_rows()

all_contrasts <- all_contrasts |>
  dplyr::mutate(
    species   = factor(species, levels = rev(species_order)),
    dimension = factor(dimension, levels = c("Dim.1", "Dim.2", "Dim.3")),
    phase     = factor(phase, levels = c("PRE", "POST")),
    sig       = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p_forest <- ggplot(all_contrasts, aes(x = estimate, y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), orientation = "y", width = 0) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sig), nudge_y = .2, size = 4) +
  facet_grid(phase ~ dimension, scales = "fixed") +
  labs(x = expression(Delta * " slope (Dim. high - Dim. low)"),
       y = NULL,
       caption = "x > 0 protege de la desecacion | x < 0 acelera la desecacion") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"))

ggsave("07-img/heterogeneity_forest_plot.png", p_forest, width = 11, height = 8, dpi = 300)
cat("Figura guardada: 07-img/heterogeneity_forest_plot.png\n")

saveRDS(all_contrasts, "00-data/processed/heterogeneity_contrasts.rds")
write.csv(all_contrasts, "00-data/processed/heterogeneity_contrasts.csv", row.names = FALSE)
cat("Guardado: 00-data/processed/heterogeneity_contrasts.rds y .csv\n")

# ============================================================
# 4. Hipotesis mecanistica: pericarpo (Dim.2) y contenido hidrico inicial
# ============================================================
p_dim2_t0 <- df |>
  dplyr::filter(time == 0) |>
  ggplot(aes(x = Dim.2, y = Moisture_content)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~species) +
  xlab("Pericarp robustness (Dim.2)") +
  ylab("Initial moisture content (%)")

ggsave("07-img/heterogeneity_dim2_initial_moisture.png", p_dim2_t0,
       width = 11, height = 7, dpi = 300)
cat("Figura guardada: 07-img/heterogeneity_dim2_initial_moisture.png\n")

# ============================================================
# 5. Guardado de objetos de modelos
# ============================================================
save_models(list(
  "mm.pre.0"        = mm.pre.0,
  "mm.pre.1"        = mm.pre.1,
  "mm.post.0"       = mm.post.0,
  "mm.post.1"       = mm.post.1,
  "mm.pre.d1.het"   = mm.pre.d1.het,
  "mm.pre.d2.het"   = mm.pre.d2.het,
  "mm.pre.d3.het"   = mm.pre.d3.het,
  "m.pre.all.het"   = m.pre.all.het,
  "mm.post.d1.het"  = mm.post.d1.het,
  "mm.post.d2.het"  = mm.post.d2.het,
  "mm.post.d3.het"  = mm.post.d3.het,
  "m.post.all.het"  = m.post.all.het
))

if (GEN_DASH == "Y") {
  dir.create("06-html", showWarnings = FALSE)
  modelos_dash <- list(
    mm.pre.0        = mm.pre.0,
    mm.pre.1        = mm.pre.1,
    mm.post.0       = mm.post.0,
    mm.post.1       = mm.post.1,
    mm.pre.d1.het   = mm.pre.d1.het,
    mm.pre.d2.het   = mm.pre.d2.het,
    mm.pre.d3.het   = mm.pre.d3.het,
    m.pre.all.het   = m.pre.all.het,
    mm.post.d1.het  = mm.post.d1.het,
    mm.post.d2.het  = mm.post.d2.het,
    mm.post.d3.het  = mm.post.d3.het,
    m.post.all.het  = m.post.all.het
  )
  for (nm in names(modelos_dash)) {
    tryCatch(
      easystats::model_dashboard(modelos_dash[[nm]],
                                 output_dir = "06-html/",
                                 output_file = paste0("modeldashboard_", nm, ".html")),
      error = function(e) warning("model_dashboard omitido (", nm, "): ", conditionMessage(e))
    )
  }
}

cat("\n===== d.05.1.heterogeneus_effects.R completado =====\n")



