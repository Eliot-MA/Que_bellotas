# Verificacion rapida: que brms acepta la formula con gr(phylo_species, cov=A) + codigo + id_bellota
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(brms)
})

df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    id_bellota = factor(id_bellota),
    codigo = factor(codigo),
    species = factor(species),
    phylo_species = factor(species)
  )

df.t1 <- df |> filter(time < 90)
A <- readRDS("00-data/phylo/oak_vcv.rds")
A <- A[levels(df.t1$species), levels(df.t1$species)]

cat("Niveles phylo_species == filas A:", all(levels(df.t1$phylo_species) == rownames(A)), "\n")
cat("Niveles codigo:", nlevels(df.t1$codigo), "\n")

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

priors <- c(
  prior(normal(40, 20), class = "Intercept"),
  prior(normal(0, 10),  class = "b"),
  prior(student_t(3, 0, 10), class = "sd"),
  prior(student_t(3, 0, 10), class = "sigma"),
  prior(lkj(2), class = "cor")
)

cat("\n--- make_stancode M1 ---\n")
sc1 <- make_stancode(form1, data = df.t1, data2 = list(A = A),
                     family = gaussian(), prior = priors)
cat("M1 OK, stancode de", length(sc1), "lineas\n")

cat("\n--- make_stancode M2 ---\n")
sc2 <- make_stancode(form2, data = df.t1, data2 = list(A = A),
                     family = gaussian(), prior = priors)
cat("M2 OK, stancode de", length(sc2), "lineas\n")

cat("\n--- parametros aleatorios en M1 (para la CCI) ---\n")
cat(grep("sd_phylo_species__Intercept|sd_codigo__Intercept|sd_id_bellota__Intercept|sigma", sc1, value = TRUE)[1:6], sep = "\n")
