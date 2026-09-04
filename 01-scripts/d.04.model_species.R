# Load data ----
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

# Libraries ----
library(tidyverse)
library(moments)

# Creat working dataframe ----
df <- df.bellotas |> 
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |> 
  left_join(y = df.famd, by = "id_bellota") |> 
  drop_na(Moisture_content, Dim.1, Dim.2, Dim.3) |> 
  rename(time = tiempo_acumulado_horas) |> 
  mutate(
    log.t = log(time+1), 
    sqrt.t = sqrt(time+1)
  )

# Model selection ----
library(glmmTMB)

m.orig <- glmmTMB(Moisture_content ~ time * species + (time | provenance/id_bellota), data = df)

m.log <- glmmTMB(Moisture_content ~ log.t * species + 
                   (log.t | provenance/id_bellota), 
                 data = df
                 )
# Convergence problem
# I couldn't resolve it with more iterations or alternative likelihood algorihtms

m.sqr <- glmmTMB(Moisture_content ~ sqrt.t * species + (sqrt.t|provenance/id_bellota), data = df)

## Segmented model ----
### Breakpoint detection with davies.test

library(segmented)

# Límite de búsqueda del punto de inflexión (horas)
lim_horas <- 150

# Preparar tabla de resultados
ids <- unique(df$id_bellota)
res <- data.frame(
  id_bellota = ids,
  tiempo_inflexion = NA_real_,
  davies_p = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(ids)) {
  id <- ids[i]
  dat <- subset(df, id_bellota == id)
  if (nrow(dat) < 5 || length(unique(dat$time)) < 4) next
  
  modelo <- tryCatch(lm(Moisture_content ~ time, data = dat), error = function(e) e)
  if (inherits(modelo, "error")) next
  
  x <- dat$time
  q_min_val <- min(x, na.rm = TRUE)
  alpha_low  <- mean(x <= q_min_val, na.rm = TRUE)
  alpha_high <- mean(x <= lim_horas, na.rm = TRUE)
  
  n_obs <- sum(!is.na(x))
  eps <- 1 / (10 * n_obs)
  alpha_low  <- max(alpha_low, 1 / n_obs, eps)
  alpha_high <- min(max(alpha_high, alpha_low + eps), 1 - 1 / n_obs)
  
  if (alpha_high > alpha_low + 1e-12) {
    seg.modelo <- tryCatch(
      segmented(modelo, psi = list(time = median(x, na.rm = TRUE)),
                control = seg.control(alpha = c(alpha_low, alpha_high))),
      error = function(e) e
    )
    if (inherits(seg.modelo, "error")) {
      seg.modelo <- tryCatch(
        segmented(modelo, psi = list(time = median(x, na.rm = TRUE))),
        error = function(e) e
      )
    }
  } else {
    seg.modelo <- tryCatch(
      segmented(modelo, psi = list(time = median(x, na.rm = TRUE))),
      error = function(e) e
    )
  }
  
  if (!inherits(seg.modelo, "error") && !is.null(seg.modelo$psi)) {
    res$tiempo_inflexion[i] <- as.numeric(seg.modelo$psi[1, "Est."])
  }
  
  dav <- tryCatch(davies.test(modelo, seg.Z = ~ time), error = function(e) e)
  if (!inherits(dav, "error")) res$davies_p[i] <- dav$p.value
}

puntos_inflexion <- res

# Guardar tabla por bellota con especie y procedencia (para el articulo)
# Nota: se usa un objeto aparte para no anadir columnas duplicadas a df
breakpoints_out <- puntos_inflexion |>
  left_join(
    df |> distinct(id_bellota, species, codigo),
    by = "id_bellota"
  )

write.csv(breakpoints_out, "00-data/breakpoints_per_acorn.csv", row.names = FALSE)
cat("Tabla de puntos de inflexion guardada en 00-data/breakpoints_per_acorn.csv\n")

# Añadir la columna lógica al dataframe original
df <- df %>%
  left_join(puntos_inflexion, by = "id_bellota") %>%
  mutate(punto_inflexion = time == tiempo_inflexion)

# Explore
df |> 
  dplyr::select(species, tiempo_inflexion) |>
  distinct() |> 
  ggplot(aes(x = tiempo_inflexion)) + 
  geom_histogram() +
  facet_wrap(~species)

# Scatter plot: Moisture_content ~ time con puntos coloreados antes/después del inflexión
df2 <- df %>%
  group_by(id_bellota) %>%
  mutate(fase = ifelse(time <= tiempo_inflexion, "antes", "despues")) %>%
  ungroup()

p1 <- df2 |>
  ggplot(aes(x = time, y = Moisture_content)) +
  geom_line(aes(group = id_bellota), alpha = 0.15, linewidth = 0.4) +
  geom_point(aes(color = fase), size = 0.8, alpha = .2) +
  geom_vline(aes(xintercept = 100), linetype = "dashed", colour = "red", alpha = .8) +
  scale_color_manual(values = c("antes" = "red", "despues" = "grey30"), guide = "none") +
  labs(x = "Time (hours)", y = "Moisture content (%)", 
       caption = "Points before calculated breakpoints per acorn (red).", 
       title = "Analysys of the presence of breakpoints in the desiccation curves.") +
  theme_minimal()
 
sumtable_breakpoints <- df2 |> 
  group_by(species) |> 
  distinct(id_bellota, tiempo_inflexion) |> 
  summarise(
    n = n(), 
    minimum = min(tiempo_inflexion), 
    Q05 = quantile(tiempo_inflexion, probs = .05),
    Q25 = quantile(tiempo_inflexion, probs = .25), 
    Q50 = quantile(tiempo_inflexion, probs = .5), 
    mean = mean(tiempo_inflexion), 
    Q75 = quantile(tiempo_inflexion, probs = .75), 
    Q95 = quantile(tiempo_inflexion, probs = .95),
    maximum = max(tiempo_inflexion)
  )

ggsave("07-img/breakpoint_scatter.png", p1, width = 8, height = 6)
write.csv(sumtable_breakpoints, "00-data/sumtable_breakpoints.csv", row.names = FALSE)

## Fit model

# 1) Create time segmented
umbral <- 100

df <- df |> 
  mutate(
    t_pre = time, 
    t_post = pmax(0, time - umbral)
  )

# 2) Fit model
m.seg <- glmmTMB(Moisture_content ~ t_pre * species + t_post * species + (t_pre|provenance/id_bellota), data = df)

# Compare  models
performance::compare_performance(m.orig, m.log, m.sqr, m.seg)

## Table that justifies the model selection
# library(gt)
# --- 1) Lista de modelos a comparar ---
mods <- list(
  "Original time"           = m.orig,
  "Log-transformed"         = m.log,
  "Square root transformed" = m.sqr,
  "Segmented model"         = m.seg
)

# --- 2) Función auxiliar para resumir cada modelo ---
resumen_modelo <- function(mod, nombre, vif_umbral = 10) {
  
  # R2 de Nakagawa
  r2 <- tryCatch(
    performance::r2_nakagawa(mod),
    error = function(e) NULL
  )
  
  # Colinealidad / VIF
  collin <- tryCatch(
    performance::check_collinearity(mod),
    error = function(e) NULL
  )
  
  max_vif <- if (!is.null(collin) && "VIF" %in% names(collin)) {
    max(collin$VIF, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  # Convergencia
  converged <- tryCatch(
    isTRUE(mod$fit$convergence == 0),
    error = function(e) NA
  )
  
  # Singularidad
  singular <- tryCatch(
    performance::check_singularity(mod),
    error = function(e) NA
  )
  
  tibble(
    Modelo = nombre,
    AIC = AIC(mod),
    logLik = as.numeric(logLik(mod)),
    `R2 marginal` = if (!is.null(r2)) unname(r2$R2_marginal) else NA_real_,
    `R2 condicional` = if (!is.null(r2)) unname(r2$R2_conditional) else NA_real_,
    Converge = ifelse(is.na(converged), NA, ifelse(converged, "Sí", "No")),
    Singular = ifelse(is.na(singular), NA, ifelse(singular, "Sí", "No")),
    `VIF máximo` = max_vif,
    `VIF alto` = ifelse(is.na(max_vif), NA, ifelse(max_vif >= vif_umbral, "Sí", "No"))
  )
}

# --- 3) Construir tabla base ---
tabla_comp <- imap_dfr(mods, resumen_modelo)

# --- 4) Calcular ΔAIC y pesos de Akaike ---
tabla_comp <- tabla_comp %>%
  mutate(
    `ΔAIC` = AIC - min(AIC, na.rm = TRUE),
    `Peso Akaike` = exp(-0.5 * `ΔAIC`) / sum(exp(-0.5 * `ΔAIC`), na.rm = TRUE),
    `Bandera diagnóstica` = case_when(
      Converge == "No" ~ "No converge",
      Singular == "Sí" ~ "Singular",
      `VIF alto` == "Sí" ~ "Colinealidad alta",
      TRUE ~ "OK"
    )
  ) %>%
  dplyr::select(
    Modelo, AIC, `ΔAIC`, `Peso Akaike`,
    `R2 marginal`, `R2 condicional`,
    Converge, Singular, `VIF máximo`, `VIF alto`,
    `Bandera diagnóstica`, logLik
  ) %>%
  arrange(`ΔAIC`)

# --- 5) Redondeo para presentación ---
tabla_comp_final <- tabla_comp %>%
  mutate(
    across(c(AIC, `ΔAIC`, `Peso Akaike`, `R2 marginal`, `R2 condicional`, `VIF máximo`, logLik),
           ~ round(.x, 3))
  )

# Save table
write.csv(tabla_comp_final, "00-data/model_comparisons.csv")

# Species models
## data
df.t1 <- df |> 
  filter(time < 94)

df.t2 <- df |> 
  filter(time > 94)

## Model pre
mm.pre <- glmmTMB(Moisture_content ~ time * species + (time|provenance) + (time|id_bellota), data = df.t1)

## Model post
mm.post <- glmmTMB(Moisture_content ~ time * species + (time|provenance) + (time|id_bellota), data = df.t2)

# ============================================================
# Export de resultados para el articulo / material suplementario
# ============================================================
source("01-scripts/00-export_helpers.R")

# Objetos de modelo
save_models(list(m.orig = m.orig, m.log = m.log, m.sqr = m.sqr,
                 m.seg = m.seg, mm.pre = mm.pre, mm.post = mm.post))

# Tabla de coeficientes de los modelos finales (MC ~ time * species)
coef_especie <- coef_table(list("mm.pre (t<94h)" = mm.pre, "mm.post (t>94h)" = mm.post))
write.csv(coef_especie, "00-data/coef_species_models.csv", row.names = FALSE)
cat("Tabla de coeficientes guardada en 00-data/coef_species_models.csv\n")

# Componentes de varianza + barplot por nivel
varcomp_especie <- varcomp_table(list("mm.pre (t<94h)" = mm.pre, "mm.post (t>94h)" = mm.post))
write.csv(varcomp_especie, "00-data/varcomp_species_models.csv", row.names = FALSE)
cat("Componentes de varianza guardados en 00-data/varcomp_species_models.csv\n")
plot_varcomp(varcomp_especie, "07-img/varcomp_species_models.png",
             "Varianza por nivel - modelos de especie")

# Diagnosticos de los modelos finales de especie
modelos_especie <- list(mm.pre = mm.pre, mm.post = mm.post)
plot_check_model(modelos_especie, sufijo = "species")
plot_obs_fitted(modelos_especie, sufijo = "species")
plot_dharma(modelos_especie, sufijo = "species")

# Efectos marginales de time * species (pre y post)
suppressPackageStartupMessages(library(ggeffects))
for (nm in names(modelos_especie)) {
  mod <- modelos_especie[[nm]]
  pred <- tryCatch(
    ggeffects::ggpredict(mod, terms = c("time", "species")),
    error = function(e) e
  )
  if (inherits(pred, "error")) {
    warning("ggpredict fallo para ", nm, ": ", conditionMessage(pred))
    next
  }
  p <- plot(pred) + labs(title = paste("MC ~ time * species -", nm))
  ggsave(file.path("07-img", paste0("marginal_", nm, ".png")), p, width = 9, height = 6)
  cat("Efectos marginales guardados en 07-img/marginal_", nm, ".png\n", sep = "")
}

# ============================================================
# Figura principal: tasas de desecacion pre/post por especie
# ============================================================
# Pendientes de desecacion (MC%/h) por especie estimadas con emtrends a partir
# de mm.pre y mm.post. Tasa = -pendiente (valor positivo). Especies ordenadas
# de arriba a abajo por bioclima (templado, sub-mediterraneo, mediterraneo),
# colores Okabe-Ito por bioclima (mismo codigo que el grafico de MC50).
# PRE = triangulos, POST = circulos, barra = IC95 de emtrends.
# Comparaciones multiples dentro de PRE y dentro de POST: se anaden letras de
# diferencias significativas encima de cada punto, asignadas por el metodo de
# INTERVALOS DE CONFIANZA: dos especies no comparten letra si el IC95 de su
# diferencia por pares excluye el cero. La familia de comparaciones se protege
# con el test global de la interaccion time:species (procedimiento de Fisher).
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(car))   # Anova() para el test global

# --- 4.1 Pendientes de desecacion con emtrends (por fase) ---
get_desiccation_rates <- function(model, phase) {
  emmeans::emtrends(model, ~ species, var = "time") |>
    as.data.frame() |>
    mutate(
      phase     = phase,
      rate      = -time.trend,      # valor positivo: MC%/h
      rate_low  = -asymp.UCL,       # invertir intervalo junto con el signo
      rate_high = -asymp.LCL
    ) |>
    dplyr::select(species, phase, rate, rate_low, rate_high)
}

rates <- bind_rows(
  get_desiccation_rates(mm.pre, "PRE"),
  get_desiccation_rates(mm.post, "POST")
)

# --- 4.2 Post-hoc de comparaciones multiples dentro de cada fase ---
# Metodo por INTERVALOS DE CONFIANZA. Para cada fase se obtienen los IC95 de
# las diferencias por pares entre pendientes (confint(pairs(emtrends))); dos
# especies comparten letra solo si el IC de su diferencia incluye el cero
# (no significativo). El procedimiento se protege con el test global de la
# interaccion time:species (Anova tipo II) antes de declarar diferencias.
ic_cld <- function(order, ci_df, alpha = 0.05) {
  n <- length(order)
  sigmat <- matrix(FALSE, n, n, dimnames = list(order, order))
  for (i in seq_len(nrow(ci_df))) {
    gr <- strsplit(as.character(ci_df$contrast[i]), " - ")[[1]]
    if (gr[1] %in% order && gr[2] %in% order &&
        (ci_df$asymp.LCL[i] > 0 | ci_df$asymp.UCL[i] < 0)) {
      sigmat[gr[1], gr[2]] <- TRUE
      sigmat[gr[2], gr[1]] <- TRUE
    }
  }
  out <- setNames(rep(NA_character_, n), order)
  pool <- c(letters, LETTERS)
  k <- 1
  remaining <- order
  while (length(remaining) > 0) {
    L <- pool[k]; k <- k + 1
    group <- character(0)
    for (sp in remaining) {
      if (length(group) == 0 || all(!sigmat[sp, group])) group <- c(group, sp)
    }
    out[group] <- vapply(group, function(g) {
      if (is.na(out[g])) L else paste0(out[g], L)
    }, character(1))
    remaining <- setdiff(remaining, group)
  }
  out
}

# Efectos fijos (Anova Type II) de los modelos pre/post umbral: Chisq, gl y p.
# Incluye el test global de la interaccion time:species (proteccion del procedimiento).
global_tests <- bind_rows(
  Anova(mm.pre,  type = "II") |> as.data.frame() |> rownames_to_column("term") |> mutate(phase = "PRE"),
  Anova(mm.post, type = "II") |> as.data.frame() |> rownames_to_column("term") |> mutate(phase = "POST")
) |>
  mutate(
    p_valor       = `Pr(>Chisq)`,
    significancia = case_when(
      p_valor < 0.001 ~ "***",
      p_valor < 0.01  ~ "**",
      p_valor < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) |>
  dplyr::select(phase, term, Df, Chisq, p_valor, significancia)

cat("Efectos fijos (Anova Type II) por fase:\n")
print(global_tests)

# Tabla de efectos fijos (material suplementario)
write.csv(global_tests, "00-data/anova_prepost_species.csv", row.names = FALSE)
cat("Tabla de efectos fijos guardada en 00-data/anova_prepost_species.csv\n")

# Test global de la interaccion time:species (proteccion del procedimiento)
test_interaccion <- global_tests |>
  filter(term == "time:species") |>
  mutate(sig = p_valor < 0.05)

cat("Test global time:species por fase:\n")
print(test_interaccion |> dplyr::select(phase, Chisq, Df, p_valor, sig))

# IC95 de las diferencias por pares por fase (columnas LCL/UCL por pareja)
pairs_ci <- function(model) {
  et <- emmeans::emtrends(model, ~ species, var = "time")
  ci <- as.data.frame(confint(pairs(et, adjust = "none")))
  ci
}

ic_cld_per_phase <- function(model) {
  et  <- emmeans::emtrends(model, ~ species, var = "time")
  ci  <- as.data.frame(confint(pairs(et, adjust = "none")))
  ord <- as.data.frame(et)$species[order(-as.data.frame(et)$time.trend)]  # por tasa
  tibble(species = ord, letters = unname(ic_cld(ord, ci)))
}

rates_letters <- bind_rows(
  ic_cld_per_phase(mm.pre) |> mutate(phase = "PRE"),
  ic_cld_per_phase(mm.post) |> mutate(phase = "POST")
)

# Tabla con los IC95 de las diferencias por pares (material suplementario)
ic_pairs_table <- bind_rows(
  pairs_ci(mm.pre)  |> mutate(phase = "PRE"),
  pairs_ci(mm.post) |> mutate(phase = "POST")
)
write.csv(ic_pairs_table, "00-data/desiccation_pairs_ci.csv", row.names = FALSE)
cat("Tabla de IC de comparaciones por pares guardada en 00-data/desiccation_pairs_ci.csv\n")

# ---------------------------------------------------------------------------
# 4.2b Comparaciones SELECTIVAS (inter-bioclima)
# ---------------------------------------------------------------------------
# Solo se comparan especies de DISTINTO bioclima (Mediterraneo, Sub-Mediterraneo,
# Templado). No se comparan especies del mismo bioclima (p. ej. Q. robur vs
# Q. petraea), ya que la hipotesis es bioclimatica. Sobre esta familia reducida
# de contrastes se controla la multiplicidad con Sidak (adjust = "sidak"):
#   - Tukey no es aplicable: emmeans solo permite "tukey" para el conjunto
#     completo de comparaciones por pares y lo degrada a Sidak para listas
#     arbitrarias de contrastes.
#   - Sidak es el ajuste estandar para un conjunto preespecificado de
#     contrastes (hipotesis a priori), menos conservador que Bonferroni.
# Las letras se asignan por el mismo criterio de ICs: dos especies no comparten
# letra si el IC95 (Sidak) de su diferencia selectiva excluye el cero.
bioclimate_short <- c(   # codigo corto por especie para construir los pares
  "Quercus coccifera" = "Med",
  "Quercus ilex"      = "Med",
  "Quercus suber"     = "Med",
  "Quercus faginea"   = "SubMed",
  "Quercus pyrenaica" = "SubMed",
  "Quercus pubescens" = "SubMed",
  "Quercus petraea"   = "Temp",
  "Quercus robur"     = "Temp"
)

# Construye la lista de contrastes selectivos (especies de distinto bioclima)
build_selective_contrasts <- function(species_levels) {
  idx <- setNames(seq_along(species_levels), species_levels)
  clist <- list()
  for (pr in combn(species_levels, 2, simplify = FALSE)) {
    if (bioclimate_short[[pr[1]]] != bioclimate_short[[pr[2]]]) {
      v <- rep(0, length(species_levels))
      v[idx[[pr[1]]]] <- +1
      v[idx[[pr[2]]]] <- -1
      clist[[paste0(pr[1], " - ", pr[2])]] <- v
    }
  }
  clist
}

# CLD por IC a partir de un conjunto de pares declarados significativos
cld_from_significant_pairs <- function(order, sig_pairs) {
  n <- length(order)
  sigmat <- matrix(FALSE, n, n, dimnames = list(order, order))
  for (pr in sig_pairs) {
    gr <- strsplit(pr, " - ")[[1]]
    if (gr[1] %in% order && gr[2] %in% order) {
      sigmat[gr[1], gr[2]] <- TRUE
      sigmat[gr[2], gr[1]] <- TRUE
    }
  }
  out <- setNames(rep(NA_character_, n), order)
  pool <- c(letters, LETTERS)
  k <- 1
  remaining <- order
  while (length(remaining) > 0) {
    L <- pool[k]; k <- k + 1
    group <- character(0)
    for (s in remaining) {
      if (length(group) == 0 || all(!sigmat[s, group])) group <- c(group, s)
    }
    out[group] <- vapply(group, function(x) {
      if (is.na(out[x])) L else paste0(out[x], L)
    }, character(1))
    remaining <- setdiff(remaining, group)
  }
  out
}

# Letras selectivas por fase: solo comparaciones inter-bioclima con ajuste Sidak
selective_cld_per_phase <- function(model) {
  et   <- emmeans::emtrends(model, ~ species, var = "time")
  sp   <- levels(et@grid$species)
  cn   <- build_selective_contrasts(sp)
  ci   <- as.data.frame(confint(contrast(et, method = cn, adjust = "sidak")))
  sigp <- ci$contrast[ci$asymp.LCL > 0 | ci$asymp.UCL < 0]
  ord  <- sp[order(-as.data.frame(et)$time.trend)]   # por tasa
  tibble(species = ord, letters = unname(cld_from_significant_pairs(ord, sigp)))
}

rates_letters_sel <- bind_rows(
  selective_cld_per_phase(mm.pre)  |> mutate(phase = "PRE"),
  selective_cld_per_phase(mm.post) |> mutate(phase = "POST")
)

# Tabla adicional con los IC95 de las comparaciones selectivas (Sidak)
selective_pairs_from_model <- function(model, phase) {
  et <- emmeans::emtrends(model, ~ species, var = "time")
  ci <- as.data.frame(confint(contrast(
    et, method = build_selective_contrasts(levels(et@grid$species)),
    adjust = "sidak")))
  ci$phase <- phase
  ci
}
selective_pairs_table <- bind_rows(
  selective_pairs_from_model(mm.pre,  "PRE"),
  selective_pairs_from_model(mm.post, "POST")
)
write.csv(selective_pairs_table, "00-data/desiccation_selective_pairs_ci.csv",
          row.names = FALSE)
cat("Tabla de IC de comparaciones selectivas (Sidak) guardada en 00-data/desiccation_selective_pairs_ci.csv\n")

# Anadir ambas columnas de letras a rates
rates <- rates |>
  left_join(rates_letters,    by = c("species", "phase")) |>
  left_join(rates_letters_sel, by = c("species", "phase"),
            suffix = c("", "_sel"))

# --- 4.3 Bioclima y orden del eje Y ---
bioclimate <- tribble(
  ~species,                     ~bioclimate,
  "Quercus coccifera",          "Mediterranean",
  "Quercus ilex",               "Mediterranean",
  "Quercus suber",              "Mediterranean",
  "Quercus faginea",            "Sub-Mediterranean",
  "Quercus pyrenaica",          "Sub-Mediterranean",
  "Quercus pubescens",          "Sub-Mediterranean",
  "Quercus petraea",            "Temperate",
  "Quercus robur",              "Temperate"
)

short_names <- c(
  "Quercus coccifera" = "Q. coccifera",
  "Quercus ilex"      = "Q. ilex",
  "Quercus suber"     = "Q. suber",
  "Quercus faginea"   = "Q. faginea",
  "Quercus pyrenaica" = "Q. pyrenaica",
  "Quercus pubescens" = "Q. pubescens",
  "Quercus petraea"   = "Q. petraea",
  "Quercus robur"     = "Q. robur"
)

# Orden del eje y: de arriba a abajo Templado -> Sub-mediterraneo -> Mediterraneo;
# dentro de cada bioclima por la tasa de desecacion de la fase PRE (ascendente).
rate_pre <- rates |>
  filter(phase == "PRE") |>
  dplyr::select(species, rate)

plot_tab <- rates |>
  left_join(bioclimate, by = "species") |>
  left_join(rate_pre, by = "species", suffix = c("", "_pre")) |>
  mutate(name = unname(short_names[species])) |>
  arrange(factor(bioclimate, levels = c("Mediterranean",
                                        "Sub-Mediterranean",
                                        "Temperate")), rate_pre) |>
  mutate(name = factor(name, levels = unique(name)))

# --- 4.3b Comparison intervals (flechas de comparacion) ---
# Los comparison intervals de emmeans se basan en invertir las comparaciones por
# pares: dos medias no difieren significativamente <=> sus intervalos se solapan.
# Su construccion es una aproximacion (solo exacta para <= 6 medias) y produce
# intervalos ABIERTOS (NA) en las especies extremas: la mas lenta no se acota
# hacia tasas mas bajas y la mas rapida no se acota hacia tasas mas altas.
# Aqui replicamos el algoritmo interno de emmeans ('.plot.srg') para poder
# construirlos con cualquier familia de pares (completa o selectiva).
comparison_intervals <- function(est, id1, id2, diff, LCL, UCL) {
  # est: estimaciones por especie (orden = niveles); id1/id2: indices de especie
  # de cada contraste; diff/LCL/UCL: estimacion e IC95 de cada contraste.
  neach <- length(est); npairs <- length(diff)
  del     <- (UCL - LCL) / 4
  overlap <- 2 * pmin(-LCL, UCL) / (UCL - LCL)
  involved <- lapply(seq_len(neach), function(x) union(which(id2 == x), which(id1 == x)))
  mind <- sapply(involved, function(ii) min(del[ii]))
  iden <- diag(rep(1, 2 * neach))
  lmat <- rmat <- matrix(0, nrow = npairs, ncol = neach)
  y <- numeric(npairs); v1 <- 1 - overlap
  for (i in which(!is.na(v1))) {
    wgt <- 3 + 20 * max(0, 0.5 - (1 - v1[i])^2)
    if (diff[i] > 0) lmat[i, id1[i]] <- rmat[i, id2[i]] <- wgt * v1[i]
    else            rmat[i, id1[i]] <- lmat[i, id2[i]] <- wgt * v1[i]
    y[i] <- wgt * abs(diff[i])
  }
  X  <- rbind(cbind(lmat, rmat), iden)
  yy <- c(y, rep(mind, 2)); yy[is.na(yy)] <- 0
  soln <- qr.coef(qr(X), yy); soln[is.na(soln)] <- 0
  ll <- soln[seq_len(neach)]; rl <- soln[neach + seq_len(neach)]
  rng <- suppressWarnings(range(est, na.rm = TRUE)); diffr <- diff(rng)
  ii <- which(est - ll < rng[1]); ll[ii] <- est[ii] - rng[1] + 0.02 * diffr
  ii <- which(est + rl > rng[2]); rl[ii] <- rng[2] - est[ii] + 0.02 * diffr
  ll[est < rng[1] + 1e-04 * diffr] <- NA   # lado abierto (tasa minima)
  rl[est > rng[2] - 1e-04 * diffr] <- NA   # lado abierto (tasa maxima)
  list(lcmpl = est - ll, rcmpl = est + rl)
}

# Obtiene los comparison intervals por fase y familia de pares.
#  - familia == "complete": usa pairs() con adjust "tukey" (igual que plot()).
#  - familia == "selective": usa las comparaciones inter-bioclima con adjust "sidak".
comparison_intervals_per_phase <- function(model, familia = c("complete", "selective")) {
  familia <- match.arg(familia)
  et <- emmeans::emtrends(model, ~ species, var = "time")
  sp <- levels(et@grid$species)
  est <- as.data.frame(et)$time.trend
  idx <- setNames(seq_along(sp), sp)

  if (familia == "complete") {
    ps <- as.data.frame(confint(pairs(et), level = 0.95, adjust = "tukey"))
    prs <- combn(sp, 2, simplify = FALSE)
    id1 <- sapply(prs, function(p) idx[[p[1]]]); id2 <- sapply(prs, function(p) idx[[p[2]]])
  } else {
    clist <- list(); id1 <- integer(); id2 <- integer()
    for (pr in combn(sp, 2, simplify = FALSE)) if (bioclimate_short[[pr[1]]] != bioclimate_short[[pr[2]]]) {
      v <- rep(0, length(sp)); v[idx[[pr[1]]]] <- 1; v[idx[[pr[2]]]] <- -1
      clist[[paste0(pr[1], " - ", pr[2])]] <- v
      id1 <- c(id1, idx[[pr[1]]]); id2 <- c(id2, idx[[pr[2]]])
    }
    ps <- as.data.frame(confint(contrast(et, method = clist, adjust = "sidak")))
  }
  k <- ncol(ps)
  ci <- comparison_intervals(est, id1, id2, ps[[attr(ps, "estName")]], ps[[k - 1]], ps[[k]])

  # Escala de tasa POSITIVA (endreversando el signo del emtrend).
  # Un lado abierto (NA) en emmeans es un comparison interval abierto en una
  # especie extrema: se codifica como -Inf (tasas mas lentas) o +Inf (mas rapidas).
  lpos <- -ci$rcmpl; hpos <- -ci$lcmpl
  comp_low  <- ifelse(is.na(lpos), -Inf, lpos)
  comp_high <- ifelse(is.na(hpos), +Inf, hpos)
  tibble(species = sp, comp_low = comp_low, comp_high = comp_high)
}

compints_full <- bind_rows(
  comparison_intervals_per_phase(mm.pre,  "complete") |> mutate(phase = "PRE"),
  comparison_intervals_per_phase(mm.post, "complete") |> mutate(phase = "POST")
)
compints_sel <- bind_rows(
  comparison_intervals_per_phase(mm.pre,  "selective") |> mutate(phase = "PRE"),
  comparison_intervals_per_phase(mm.post, "selective") |> mutate(phase = "POST")
)
# unir a la tabla de trazado
plot_tab <- plot_tab |>
  left_join(compints_full, by = c("species", "phase"), suffix = c("", "_full")) |>
  left_join(compints_sel,  by = c("species", "phase"), suffix = c("", "_sel"))

# --- 4.4 Generar figura ---
make_prepost_plot <- function(tab, letters_col, caption_note,
                              comp_low_col = NULL, comp_high_col = NULL,
                              dodge_w = 0.45) {
  p <- ggplot(tab, aes(x = rate, y = name, colour = bioclimate)) +
    geom_errorbar(aes(xmin = rate_low, xmax = rate_high, y = name, group = phase),
                  width = 0.28, linewidth = 0.9,
                  position = position_dodge(width = dodge_w)) +
    geom_point(aes(shape = phase, y = name, group = phase), size = 3.4,
               position = position_dodge(width = dodge_w)) +
    geom_text(aes(label = .data[[letters_col]], y = name, group = phase),
              vjust = -1.4, size = 3.2, colour = "black", show.legend = FALSE,
              position = position_dodge(width = dodge_w)) +
    scale_shape_manual(values = c("PRE" = 17, "POST" = 16),
                       name = "Phase") +
    scale_colour_manual(
      values = c("Mediterranean"     = "#D55E00",   # Okabe-Ito vermilion
                 "Sub-Mediterranean" = "#E69F00",   # Okabe-Ito orange
                 "Temperate"         = "#0072B2"),  # Okabe-Ito blue
      name = "Bioclimate"
    ) +
    labs(
      x = expression("Desiccation rate (MC%/h)"),
      y = NULL,
      caption = paste0(
        "Desiccation rate estimated with emtrends from piecewise species models\n",
        "(PRE t < 94 h; POST t > 94 h). Points: mean rate; whiskers: 95% CI.\n",
        "Rates expressed as positive values (moisture loss per hour).\n",
        caption_note,
        "Colours denote the characteristic bioclimate of each species."
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      axis.text.y = element_text(face = "italic", size = 11),
      axis.title.x = element_text(size = 12, margin = margin(t = 4)),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9.5),
      plot.caption = element_text(size = 8.5, colour = "grey30",
                                  hjust = 0, margin = margin(t = 6))
    )

  if (!is.null(comp_low_col) && !is.null(comp_high_col)) {
    cl <- tab[[comp_low_col]]; ch <- tab[[comp_high_col]]
    rate <- tab$rate; nm <- tab$name
    # limites del panel a partir de todos los valores finitos
    allv <- c(rate, tab$rate_low, tab$rate_high, cl[is.finite(cl)], ch[is.finite(ch)])
    xlo <- min(allv, na.rm = TRUE); xhi <- max(allv, na.rm = TRUE)
    span <- diff(c(xlo, xhi)); xlo <- xlo - 0.04 * span; xhi <- xhi + 0.04 * span
    # flechas a los lados: lado abierto (Inf) se extiende hasta el limite del panel
    arrow_tbl <- data.frame(
      name = nm,
      phase = rep(tab$phase, 2),
      x_start = rep(rate, 2),
      x_end   = c(ifelse(is.finite(cl), cl, xlo), ifelse(is.finite(ch), ch, xhi)),
      stringsAsFactors = FALSE
    )
    p <- p +
      geom_segment(data = arrow_tbl,
                   aes(x = x_start, xend = x_end, y = name, yend = name, group = phase),
                   colour = "grey25", linewidth = 0.55,
                   position = position_dodge(width = dodge_w),
                   arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
      coord_cartesian(xlim = c(xlo, xhi))
  }
  p
}

# Figura principal: letras de comparaciones por pares con IC95 sin ajuste
p_prepost <- make_prepost_plot(
  plot_tab, "letters",
  c("Letters: groups from pairwise 95% CIs within each phase (shared letter = no significant difference).\n")
)
ggsave("07-img/paper_prepost_desiccation.png", p_prepost,
       width = 160, height = 110, units = "mm", dpi = 600)
cat("Figura de tasas pre/post guardada en 07-img/paper_prepost_desiccation.png\n")

# Figura con comparaciones SELECTIVAS (inter-bioclima, ajuste Sidak)
p_prepost_sel <- make_prepost_plot(
  plot_tab, "letters_sel",
  c("Letters: groups from selective inter-bioclimate comparisons (only species of different bioclimate),\n",
    "95% CIs adjusted by Sidak within each phase (shared letter = no significant difference).\n")
)
ggsave("07-img/paper_prepost_desiccation_comp_sel.png", p_prepost_sel,
       width = 160, height = 110, units = "mm", dpi = 600)
cat("Figura de tasas pre/post (comparaciones selectivas) guardada en 07-img/paper_prepost_desiccation_comp_sel.png\n")

# --- Version alternativa: comparison intervals (conjunto completo, tukey) ---
p_prepost_compint <- make_prepost_plot(
  plot_tab, "letters",
  c("Letters: groups from pairwise 95% CIs within each phase (shared letter = no significant difference).\n",
    "Grey arrows: comparison intervals (inverse of pairwise comparisons, Tukey); overlapping intervals = no significant difference.\n",
    "Open ends (arrows to the plot edge) indicate an extremal species with no comparison beyond it.\n"),
  comp_low_col = "comp_low", comp_high_col = "comp_high"
)
ggsave("07-img/paper_prepost_desiccation_compint.png", p_prepost_compint,
       width = 160, height = 110, units = "mm", dpi = 600)
cat("Figura con comparison intervals guardada en 07-img/paper_prepost_desiccation_compint.png\n")

# --- Version alternativa: comparison intervals selectivos (inter-bioclima, Sidak) ---
p_prepost_sel_compint <- make_prepost_plot(
  plot_tab, "letters_sel",
  c("Letters: groups from selective inter-bioclimate comparisons (only species of different bioclimate),\n",
    "95% CIs adjusted by Sidak within each phase (shared letter = no significant difference).\n",
    "Grey arrows: comparison intervals from the same selective comparisons (Sidak); overlapping intervals = no significant difference.\n"),
  comp_low_col = "comp_low_sel", comp_high_col = "comp_high_sel"
)
ggsave("07-img/paper_prepost_desiccation_comp_sel_compint.png", p_prepost_sel_compint,
       width = 160, height = 110, units = "mm", dpi = 600)
cat("Figura con comparison intervals selectivos guardada en 07-img/paper_prepost_desiccation_comp_sel_compint.png\n")

# Tabla suplementaria de tasas de desecacion por especie (+ letras post-hoc)
write.csv(plot_tab, "00-data/desiccation_rates_prepost.csv", row.names = FALSE)
cat("Tabla de tasas guardada en 00-data/desiccation_rates_prepost.csv\n")


