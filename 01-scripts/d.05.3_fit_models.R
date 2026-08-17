# ============================================================
# d.05.3_fit_models.R
# Modelos finales para describir el efecto de los rasgos de las bellotas
# sobre la tasa de desecacion con correccion filogenetica (brms).
#
# TIRADA COMPLETA: iter=4000, warmup=2000, chains=4, cores=4
# Configuracion: max_treedepth=14, adapt_delta=0.99
# Tiempo estimado: ~6 horas
# ============================================================

RUN_FIT <- TRUE
SMOKE_TEST <- FALSE  # FALSE = tirada completa (4000 iter, 4 cadenas)
N_CORES <- 4          # 4 cadenas en paralelo
REFRESH <- 10         # progreso visible en consola

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
# Usar niveles comunes de ambas fases
common_species <- intersect(levels(df.t1$species), levels(df.t2$species))
A <- A[common_species, common_species]

# brms no permite usar el mismo factor dos veces con distinta covarianza:
# copia dedicada al termino filogenetico
df.wb.t1$phylo_species <- factor(df.wb.t1$species, levels = common_species)
df.wb.t2$phylo_species <- factor(df.wb.t2$species, levels = common_species)

# ---- 2. Especificacion de los modelos (presentacion) ----
form1 <- "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | codigo) +
          (1 | id_bellota)"

form2 <- "Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | species) +
          (1 + time_s | codigo) +
          (1 | id_bellota)"

form3 <- "Moisture_content ~ time_s * (D1_within + D2_within + D3_within +
                                       D1_between + D2_between + D3_between) +
          (1 + time_s | gr(phylo_species, cov = A)) +
          (1 + time_s | codigo) +
          (1 | id_bellota)"

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
      control = list(adapt_delta = 0.99, max_treedepth = 14),
      seed = seed, refresh = REFRESH,
      file = file.path("00-data/phylo", nm)
    )
  }

  # Cargar PRE ya ajustado
  pre_file <- "00-data/phylo/m3_phylo_wb_pre.rds"
  if (file.exists(pre_file)) {
    m3_pre <- readRDS(pre_file)
    cat("\nModelo PRE cargado de:", pre_file, "\n")
  } else if (exists("m3")) {
    # Si existe 'm3' en el entorno (de la ejecucion anterior), usarlo
    m3_pre <- m3
    saveRDS(m3_pre, pre_file)
    cat("\nModelo PRE recuperado de la sesion y guardado en:", pre_file, "\n")
  } else {
    stop("No se encontro el modelo PRE. Ejecuta primero la fase PRE.")
  }

  # Ajustar solo POST
  m3_post <- fit_model(form3, "m3_phylo_wb_post", df.wb.t2, seed = 126)

  cat("\n################ Diagnostico de convergencia (M3 PRE) ################\n")
  diag_tabla <- function(fit) {
    draws <- brms::as_draws_df(fit, variable = "^b_", regex = TRUE)
    b_cols <- grep("^b_", colnames(draws), value = TRUE)
    n_chains <- max(draws$.chain)
    n_iter <- nrow(draws) / n_chains
    results <- lapply(b_cols, function(col) {
      chain_vals <- split(draws[[col]], draws$.chain)
      mat <- do.call(cbind, chain_vals)
      list(
        param = col,
        Rhat = posterior::rhat(mat),
        Bulk_ESS = posterior::ess_bulk(mat),
        Tail_ESS = posterior::ess_tail(mat)
      )
    })
    tibble(
      param = sapply(results, `[[`, "param"),
      Rhat = sapply(results, `[[`, "Rhat"),
      Bulk_ESS = sapply(results, `[[`, "Bulk_ESS"),
      Tail_ESS = sapply(results, `[[`, "Tail_ESS")
    ) |>
      dplyr::mutate(estado = ifelse(
        Rhat < 1.01 & Bulk_ESS > 400 & Tail_ESS > 400,
        "ok", "REVISAR"
      ))
  }
  cat("\n-- M3 PRE (t < 94h) --\n")
  print(diag_tabla(m3_pre), n = Inf)

  cat("\n################ Diagnostico de convergencia (M3 POST) ################\n")
  cat("\n-- M3 POST (t > 94h) --\n")
  print(diag_tabla(m3_post), n = Inf)

  dir.create("07-img", showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave("07-img/trace_m3_pre.png",
    bayesplot::mcmc_trace(m3_pre, regex_pars = "b_"),
    width = 12, height = 8, dpi = 150)
  ggplot2::ggsave("07-img/trace_m3_post.png",
    bayesplot::mcmc_trace(m3_post, regex_pars = "b_"),
    width = 12, height = 8, dpi = 150)
  cat("\nTraceplots guardados en 07-img/trace_m3_pre.png y trace_m3_post.png\n")

  cat("\n################ M3 PRE (recomendado) ################\n")
  print(summary(m3_pre))

  cat("\n################ M3 POST ################\n")
  print(summary(m3_post))

  # cat("\n################ Comparacion LOO ################\n")
  # loo_compare(loo(m1), loo(m2), loo(m3))

  cat("\n################ Senal filogenetica (M3 PRE) ################\n")
  hyp_phylo <- paste(
    "sd_phylo_species__Intercept^2 /",
    "(sd_phylo_species__Intercept^2 + sd_codigo__Intercept^2 +",
    " sd_id_bellota__Intercept^2 + sigma^2) > 0"
  )
  print(hypothesis(m3_pre, hyp_phylo, class = NULL))

  cat("\n################ Senal filogenetica (M3 POST) ################\n")
  print(hypothesis(m3_post, hyp_phylo, class = NULL))

  saveRDS(list(pre = m3_pre, post = m3_post), "00-data/phylo/m3_final.rds")
  cat("\nModelos M3 finales guardados en 00-data/phylo/m3_final.rds\n")
}

