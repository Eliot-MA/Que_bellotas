# ============================================================
# d.05.3_fit_models.R
# Modelos finales para describir el efecto de los rasgos de las bellotas
# sobre la tasa de desecacion con correccion filogenetica (brms).
#
# SOLO SE PRESENTAN: por defecto no se ejecuta ningun ajuste (RUN_FIT = FALSE).
# Cuando la filogenia (d.05.2) y la justificacion del modelo dentro/entre
# especies (d.05.1) esten validadas, pon RUN_FIT <- TRUE y ejecuta este
# script (desde d.05.model_traits.R o en solitario).
#
# Modelos candidatos (a ajustar en cada fase de desecacion, t<94h y t>94h):
#   M1: filogenia + procedencia + bellota
#   M2: filogenia + especie libre + procedencia + bellota
#   M3: descomposicion dentro/entre especies + filogenia  (MODELO RECOMENDADO)
#
# Salidas (si RUN_FIT = TRUE):
#   00-data/phylo/m1_phylo, m2_phylo_spp, m3_phylo_wb  (brmsfit, con cache)
#   07-img/trace_m1.png, trace_m2.png, trace_m3.png  (trazas de las cadenas)
# ============================================================

RUN_FIT <- TRUE
SMOKE_TEST <- TRUE   # TRUE = ajuste rapido (1000 iter) para validar el pipeline
N_CORES <- 1          # 1 evita el cuelgue conocido de rstan en Windows; sube a 4 si va estable
REFRESH <- 10         # progreso visible en consola (0 = silencio total)

suppressPackageStartupMessages(library(tidyverse))

# ---- 0. Datos (si se ejecuta en solitario, reconstruir) ----
if (!exists("df") || !exists("df.t1") || !exists("df.t2")) {
  df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
  df.famd     <- read.csv("00-data/famd_ind_coord.csv")
  df <- df.bellotas |>
    dplyr::select(-X) |>
    dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
    left_join(y = df.famd, by = "id_bellota") |>
    tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
    rename(time = tiempo_acumulado_horas) |>
    mutate(
      time_s     = as.vector(scale(time)),
      species    = factor(species),
      provenance = factor(provenance),
      id_bellota = factor(id_bellota)
    )
  t94  <- as.vector((94 - mean(df$time)) / sd(df$time))
  df.t1 <- df |> filter(time_s < t94)
  df.t2 <- df |> filter(time_s > t94)
}

# Descomposicion dentro/entre especies (necesaria para M3)
add_within_between <- function(dat) {
  dat |>
    group_by(species) |>
    dplyr::mutate(
      D1_between = mean(Dim.1), D1_within = Dim.1 - D1_between,
      D2_between = mean(Dim.2), D2_within = Dim.2 - D2_between,
      D3_between = mean(Dim.3), D3_within = Dim.3 - D3_between
    ) |>
    ungroup()
}
df.wb.t1 <- add_within_between(df.t1)
df.wb.t2 <- add_within_between(df.t2)

# ---- 1. Filogenia ----
tree <- readRDS("00-data/phylo/oak_tree.rds")
A <- readRDS("00-data/phylo/oak_vcv.rds")
stopifnot(all(levels(df.t1$species) %in% colnames(A)))

# brms exige que los niveles de A coincidan con los del factor
A <- A[levels(df.t1$species), levels(df.t1$species)]

# brms no permite usar el mismo factor dos veces con distinta covarianza:
# copia dedicada al termino filogenetico
df.wb.t1$phylo_species <- factor(df.wb.t1$species)

# ---- 2. Especificacion de los modelos (presentacion) ----
form1 <- "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | codigo) +
          (1 + time_s | id_bellota)"

form2 <- "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | species) +
          (1 + time_s | codigo) +
          (1 + time_s | id_bellota)"

form3 <- "Moisture_content ~ time_s * (D1_within + D2_within + D3_within +
                                       D1_between + D2_between + D3_between) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | codigo) +
          (1 + time_s | id_bellota)"

cat("\n===== MODELOS PRESENTADOS (sin ejecutar) =====\n\n")
cat("M1 (filogenia + procedencia + bellota):\n")
cat(form1, "\n\n")
cat("M2 (filogenia + especie libre + procedencia + bellota):\n")
cat(form2, "\n\n")
cat("M3 (dentro/entre especies + filogenia) - RECOMENDADO:\n")
cat(form3, "\n\n")
cat("Los tres modelos deben ajustarse en cada fase de desecacion (df.t1 y df.t2).\n")
cat("Priors propuestos: Intercept normal(40, 20); betas normal(0, 10);\n")
cat("  sd y sigma student_t(3, 0, 10); correlaciones lkj(2).\n")
cat("Senal filogenetica (CCI) como proporcion de varianza a nivel phylo.\n\n")

# ---- 3. Ajuste (deshabilitado hasta validar pasos previos) ----
if (RUN_FIT) {
  suppressPackageStartupMessages(library(brms))

  priors <- c(
    prior(normal(40, 20), class = "Intercept"),
    prior(normal(0, 10),  class = "b"),
    prior(student_t(3, 0, 10), class = "sd"),
    prior(student_t(3, 0, 10), class = "sigma"),
    prior(lkj(2), class = "cor")
  )

  fit_model <- function(f, nm, dat, seed) {
    if (SMOKE_TEST) {
      iter_n <- 1000; warmup_n <- 500; chains_n <- 2
    } else {
      iter_n <- 4000; warmup_n <- 2000; chains_n <- 4
    }
    brms::brm(
      brms::bf(as.formula(f)),
      data = dat, data2 = list(A = A), family = gaussian(),
      prior = priors,
      iter = iter_n, warmup = warmup_n, chains = chains_n,
      cores = N_CORES,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      seed = seed, refresh = REFRESH,
      file = file.path("00-data/phylo", nm)
    )
  }

  m1 <- fit_model(form1, "m1_phylo",     df.wb.t1, seed = 123)
  m2 <- fit_model(form2, "m2_phylo_spp", df.wb.t1, seed = 124)
  m3 <- fit_model(form3, "m3_phylo_wb",  df.wb.t1, seed = 125)

  cat("\n################ Diagnostico de convergencia ################\n")
  diag_tabla <- function(fit) {
    arr <- brms::as_draws_array(fit, variable = "^b_", regex = TRUE)
    tibble(
      param    = dimnames(arr)[[3]],
      Rhat     = unname(posterior::rhat(arr)),
      Bulk_ESS = unname(posterior::ess_bulk(arr)),
      Tail_ESS = unname(posterior::ess_tail(arr))
    ) |>
      dplyr::mutate(estado = ifelse(
        Rhat < 1.01 & Bulk_ESS > 400 & Tail_ESS > 400,
        "ok", "REVISAR"
      ))
  }
  for (nm in c("m1", "m2", "m3")) {
    cat("\n-- ", toupper(nm), " --\n", sep = "")
    print(diag_tabla(get(nm)), n = Inf)
  }

  dir.create("07-img", showWarnings = FALSE, recursive = TRUE)
  for (nm in c("m1", "m2", "m3")) {
    ggplot2::ggsave(file.path("07-img", paste0("trace_", nm, ".png")),
      bayesplot::mcmc_trace(get(nm), regex_pars = "b_"),
      width = 12, height = 8, dpi = 150)
  }
  cat("\nTraceplots guardados en 07-img/trace_m1.png, trace_m2.png y trace_m3.png\n")

  cat("\n################ M3 (recomendado) ################\n")
  print(summary(m3))

  cat("\n################ Comparacion LOO ################\n")
  loo_compare(loo(m1), loo(m2), loo(m3))

  cat("\n################ Senal filogenetica (M3) ################\n")
  hyp_phylo <- paste(
    "sd_phylo_species__Intercept^2 /",
    "(sd_phylo_species__Intercept^2 + sd_codigo__Intercept^2 +",
    " sd_id_bellota__Intercept^2 + sigma^2) > 0"
  )
  print(hypothesis(m3, hyp_phylo, class = NULL))

  saveRDS(list(m1 = m1, m2 = m2, m3 = m3), "00-data/phylo/brms_phylo_models.rds")
  cat("\nModelos guardados en 00-data/phylo/brms_phylo_models.rds\n")
}
