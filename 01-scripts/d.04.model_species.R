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
  drop_na(Dim.1, Dim.2, Dim.3) |> 
  rename(time = tiempo_acumulado_horas) |> 
  mutate(
    log.t = log(time+1), 
    sqrt.t = sqrt(time+1)
  )

# Model selection ----
library(glmmTMB)

m.orig <- glmmTMB(Moisture_content ~ time * species + (time|id_bellota), data = df)

m.log <- glmmTMB(Moisture_content ~ log.t * species + (log.t|id_bellota), data = df)

m.sqr <- glmmTMB(Moisture_content ~ sqrt.t * species + (sqrt.t|id_bellota), data = df)

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

dir.create("07-img", showWarnings = FALSE)
ggsave("07-img/breakpoint_scatter.png", p1, width = 8, height = 6)
write.csv(sumtable_breakpoints, "00-data/sumtable_breakpoints.csv", row.names = FALSE)
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
m.seg <- glmmTMB(Moisture_content ~ t_pre * species + t_post * species + (t_pre|id_bellota), data = df)

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

