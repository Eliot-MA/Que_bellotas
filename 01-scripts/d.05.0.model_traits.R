# ============================================================
# d.05.0.model_traits.R
# Efecto de los rasgos de las bellotas sobre la tasa de desecacion.
#
# PIPELINE MAESTRO: prepara los datos comunes y lanza los sub-scripts que
# generan el material para el articulo:
#
#   1) d.05.1.heterogeneus_effects.R   [exploratorio]
#      Justifica que la heterogeneidad interespecifica es real y caracteriza
#      como se distribuye: distribuciones por especie, R2 por eje,
#      correlaciones dentro de especie, comparacion con/sin heterogeneidad,
#      comparacion con el modelo "ingenuo" y correlacion de los rasgos
#      funcionales originales.
#
#   2) d.05.2.phylo_data.R             [comprobacion]
#      Filogenia de las 8 especies (arbol CROWN de Hipp et al. 2020) y matriz
#      de covarianza A. Es insumo del modelo filogenetico de comprobacion, no
#      del modelo de referencia.
#
#   3) d.05.3.reference_model.R        [MODELO DE REFERENCIA]
#      Modelo glmmTMB con pendientes aleatorias por especie y por
#      procedencia, sin filogenia. Es el modelo del articulo.
#
#   4) d.05.4_phylo_check.R            [comprobacion, opcional]
#      Modelos bayesianos (brms) con correccion filogenetica, para verificar
#      que las conclusiones del modelo de referencia se mantienen al tener en
#      cuenta el parentesco evolutivo. NO es el modelo de referencia.
#
#   5) d.05.5.slope_effects_figures.R  [figuras]
#      Efecto de cada eje del FAMD sobre la pendiente de desecacion,
#      promediado sobre todas las especies, en PRE y POST.
#
# NOTA SOBRE LA DESCOMPOSICION DENTRO/ENTRE ESPECIES (within-between): se
# elimino de los scripts porque daba conclusiones equivalentes a otros modelos
# que ya contemplan la heterogeneidad. Se conserva solo en los informes
# (Heterogeneus_effects_acorn_traits.qmd y d.07_brms_filogenias.qmd).
#
# Rutas relativas a la raiz del repositorio.
# ============================================================

library(tidyverse)

# ---- 0. Datos ----
# Las procedencias excluidas se declaran una sola vez en
# 01-scripts/00-config_procedencias.R. IL3 (Q. ilex) queda fuera de todo el
# estudio: solo 2 tiempos de muestreo por debajo de 94 h (60 obs PRE frente a
# 150 del resto) y metadatos de localidad incompletos. Se elimina de AMBAS
# fases (PRE y POST) para que las comparaciones pre-post sean validas.
if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd     <- read.csv("00-data/famd_ind_coord.csv")

df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, prov_code, tiempo_acumulado_horas, Moisture_content) |>
  # famd_ind_coord.csv tambien trae `prov_code`; se descarta para no duplicar el
  # nombre al hacer el join. Se conserva el de la tabla larga, previo al filtro
  # del FAMD.
  left_join(y = df.famd |> dplyr::select(-prov_code), by = "id_bellota") |>
  dplyr::filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
  tidyr::drop_na(Moisture_content, Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    time_s     = as.vector(scale(time)),
    species    = factor(species),
    provenance = factor(provenance),
    id_bellota = factor(id_bellota)
  )

assert_sin_procedencias_excluidas(df, "prov_code", "d.05.0 datos")

# Rasgos del FAMD a nivel de bellota (unidad muestral)
df.traits <- df |>
  dplyr::select(id_bellota, species, provenance, Dim.1, Dim.2, Dim.3) |>
  distinct()

# Separacion en fases de desecacion (punto de corte en 94 h)
t94  <- as.vector((94 - mean(df$time)) / sd(df$time))
df.t1 <- df |> filter(time_s < t94)   # fase previa (t < 94 h)
df.t2 <- df |> filter(time_s > t94)   # fase posterior (t > 94 h)

# ---- 1. Justificacion de la heterogeneidad (exploratorio) ----
source("01-scripts/d.05.1.heterogeneus_effects.R")

# ---- 2. Datos filogeneticos (comprobacion) ----
source("01-scripts/d.05.2.phylo_data.R")

# ---- 3. Modelo de referencia ----
source("01-scripts/d.05.3.reference_model.R")

# ---- 4. Comprobacion filogenetica ----
# NO se ejecuta en el flujo principal: son ~4 modelos brms de 6000
# iteraciones (varias horas de muestreo) y su unico papel es verificar que las
# conclusiones del modelo de referencia se mantienen con correccion
# filogenetica. Se lanza en solitario, con PREP_ONLY <- FALSE:
#   source("01-scripts/d.05.4_phylo_check.R")
#
# ---- 5. Figuras de efecto marginal sobre la pendiente ----
source("01-scripts/d.05.5.slope_effects_figures.R")
