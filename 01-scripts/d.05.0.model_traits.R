# ============================================================
# d.05.model_traits.R
# Efecto de los rasgos de las bellotas sobre la tasa de desecacion.
#
# PIPELINE MAESTRO: prepara los datos comunes y lanza los tres
# sub-scripts que generan el material para el articulo:
#
#   1) d.05.1.heterogeneus_effects.R
#      Material que justifica la descomposicion del efecto DENTRO/ENTRE
#      especies (within-between): distribuciones por especie, R2 por eje,
#      correlaciones dentro de especie, comparacion de modelos con/sin
#      heterogeneidad, forest plot de contrastes por especie, modelos
#      within-between y comparacion con el modelo "ingenuo".
#
#   2) d.05.2.phylo_data.R
#      Filogenia de las 8 especies (OToL + longitudes de Grafen) y matriz
#      de covarianza A, junto con el material que justifica su eleccion
#      frente a la via V.PhyloMaker2/GBOTB (escenarios S1-S3) y el uso de
#      distancias de Grafen.
#
#   3) d.05.3_fit_models.R
#      Modelos finales a ajustar con brms. SOLO SE PRESENTAN; no se
#      ejecutan hasta validar los pasos anteriores.
#
# Rutas relativas a la raiz del repositorio.
# ============================================================

library(tidyverse)

# ---- 0. Datos ----
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd     <- read.csv("00-data/famd_ind_coord.csv")

df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  tidyr::drop_na(Moisture_content, Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    time_s     = as.vector(scale(time)),
    species    = factor(species),
    provenance = factor(provenance),
    id_bellota = factor(id_bellota)
  )

# Rasgos del FAMD a nivel de bellota (unidad muestral)
df.traits <- df |>
  dplyr::select(id_bellota, species, provenance, Dim.1, Dim.2, Dim.3) |>
  distinct()

# Separacion en fases de desecacion (punto de corte en 94 h)
t94  <- as.vector((94 - mean(df$time)) / sd(df$time))
df.t1 <- df |> filter(time_s < t94)   # fase previa (t < 94 h)
df.t2 <- df |> filter(time_s > t94)   # fase posterior (t > 94 h)

# ---- 1. Justificacion del modelo within-between ----
source("01-scripts/d.05.1.heterogeneus_effects.R")

# ---- 2. Datos filogeneticos ----
source("01-scripts/d.05.2.phylo_data.R")

# ---- 3. Modelos finales (brms) ----
# Presentados en d.05.3_fit_models.R. No se ejecutan todavia: cuando la
# filogenia (d.05.2) y la justificacion (d.05.1) esten validadas, descomentar:
# source("01-scripts/d.05.3_fit_models.R")
