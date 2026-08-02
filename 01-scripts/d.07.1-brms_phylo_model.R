# d.07.1 - Modelos filogeneticos bayesianos con brms
# Estructura propuesta (siguiendo la vignette de Burkner):
#   Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) +
#       (1 + time | gr(phylo_species, cov = A)) +  # efecto filogenetico a nivel de especie
#       (1 + time | codigo) +                       # procedencia (anidada en especie)
#       (1 + time | id_bellota)                     # efecto de bellota (anidada en procedencia)
#
# M1: filogenia + procedencia + bellota
# M2: filogenia + especie "libre" (no filogenetica) + procedencia + bellota
#     (separa el efecto filogenetico del efecto especifico no filogenetico)

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(ape)
  library(glmmTMB)
})

# ---- 1. Datos (misma construccion que d.06) ----
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    log.t = log(time + 1),
    sqrt.t = sqrt(time + 1),
    id_bellota = factor(id_bellota),
    codigo = factor(codigo),
    species = factor(species),
    # Copia de la especie para el termino filogenetico:
    # brms no permite usar el mismo factor dos veces con distinta covarianza
    phylo_species = factor(species)
  )

df.t1 <- df |> filter(time < 90)

# ---- 2. Filogenia y matriz de covarianza ----
tree <- readRDS("00-data/phylo/oak_tree.rds")
A <- readRDS("00-data/phylo/oak_vcv.rds")

# Comprobar que los niveles de species coinciden con las filas de A
cat("Niveles de species:\n"); print(levels(df.t1$species))
cat("\nTips del arbol:\n"); print(tree$tip.label)
stopifnot(all(levels(df.t1$species) %in% colnames(A)))

# Ordenar A igual que los niveles del factor (brms exige que coincidan)
A <- A[levels(df.t1$species), levels(df.t1$species)]

# ---- 3. Parametros de muestreo ----
TEST <- TRUE            # FALSE para el run "real" (mas iteraciones)
if (TEST) {
  n_iter <- 1200; n_warm <- 600; n_chain <- 4
} else {
  n_iter <- 4000; n_warm <- 2000; n_chain <- 4
}

priors <- c(
  prior(normal(40, 20), class = "Intercept"),
  prior(normal(0, 10),  class = "b"),
  prior(student_t(3, 0, 10), class = "sd"),
  prior(student_t(3, 0, 10), class = "sigma"),
  prior(lkj(2), class = "cor")
)

form1 <- bf(
  Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) +
    (1 + time | gr(phylo_species, cov = A)) +
    (1 + time | codigo) +
    (1 + time | id_bellota)
)

form2 <- bf(
  Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) +
    (1 + time | gr(phylo_species, cov = A)) +
    (1 + time | species) +
    (1 + time | codigo) +
    (1 + time | id_bellota)
)

# ---- 4. Ajuste ----
cat("\n=== M1: filogenia + procedencia + bellota ===\n")
m1 <- brm(
  form1,
  data = df.t1,
  data2 = list(A = A),
  family = gaussian(),
  prior = priors,
  iter = n_iter, warmup = n_warm, chains = n_chain, cores = n_chain,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123, refresh = 0,
  file = ifelse(TEST, "00-data/phylo/m1_phylo_test", "00-data/phylo/m1_phylo")
)

cat("\n=== M2: filogenia + especie no filogenetica + procedencia + bellota ===\n")
m2 <- brm(
  form2,
  data = df.t1,
  data2 = list(A = A),
  family = gaussian(),
  prior = priors,
  iter = n_iter, warmup = n_warm, chains = n_chain, cores = n_chain,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 124, refresh = 0,
  file = ifelse(TEST, "00-data/phylo/m2_phylo_spp_test", "00-data/phylo/m2_phylo_spp")
)

# ---- 5. Resultados ----
cat("\n################ M1 ################\n")
print(summary(m1))

cat("\n################ M2 ################\n")
print(summary(m2))

# Comparacion
cat("\n################ Comparacion LOO ################\n")
loo1 <- loo(m1, moment_match = FALSE)
loo2 <- loo(m2, moment_match = FALSE)
print(loo_compare(loo1, loo2))

# Señal filogenetica (proporcion de varianza a nivel filogenetico, M1)
# Incluye el nivel de procedencia (codigo) en el denominador
cat("\n################ Señal filogenetica (M1) ################\n")
hyp_phylo <- paste(
  "sd_phylo_species__Intercept^2 /",
  "(sd_phylo_species__Intercept^2 + sd_codigo__Intercept^2 +",
  " sd_id_bellota__Intercept^2 + sigma^2)"
)
print(hypothesis(m1, hyp_phylo, class = NULL))

saveRDS(list(m1 = m1, m2 = m2), "00-data/phylo/brms_phylo_models.rds")
cat("\nModelos guardados en 00-data/phylo/brms_phylo_models.rds\n")
