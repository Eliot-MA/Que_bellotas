# ============================================================
# d.05.z.pruebas_filo.R
# Pruebas de incorporacion de filogenia en modelos de efectos
# heterogeneos de rasgos de bellota sobre tasa de desecacion.
#
# Este script es provisional (prefijo z). Cuando el codigo este
# depurado se migrara a d.05.3_fit_models.R.
#
# Fases:
#   0. Preparacion de datos y filogenia
#   1. Smoke tests (compilacion y muestreo basico)
#
# Configuracion:
#   - Cadenas secuenciales (evita problemas con parallel en Windows)
#   - Iteraciones reducidas para smoke tests
# ============================================================

RUN_SMOKE_TEST <- TRUE
N_CHAINS_SMOKE <- 2
ITER_SMOKE     <- 200
WARMUP_SMOKE   <- 100
SEED_BASE      <- 123

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(ape)
})

dir.create("00-data/phylo", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FASE 0: Preparacion
# ============================================================
cat("\n========== FASE 0: Preparacion ==========\n")

# ---- 0.1 Cargar datos ----
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
  cat("Datos cargados:", nrow(df.t1), "obs PRE,", nrow(df.t2), "obs POST\n")
} else {
  cat("Usando dataframes existentes en el entorno\n")
}

# ---- 0.2 Cargar filogenia y matriz A ----
if (!file.exists("00-data/phylo/oak_tree.rds") || !file.exists("00-data/phylo/oak_vcv.rds")) {
  stop("Faltan archivos filogeneticos. Ejecuta primero d.05.2.phylo_data.R")
}

tree <- readRDS("00-data/phylo/oak_tree.rds")
A    <- readRDS("00-data/phylo/oak_vcv.rds")
cat("Filogenia cargada:", length(tree$tip.label), "tips\n")

# ---- 0.3 Verificar concordancia entre especies y matriz A ----
species_in_data <- levels(df$species)
species_in_A    <- colnames(A)
common_species  <- intersect(species_in_data, species_in_A)

cat("Especies en datos:", paste(species_in_data, collapse = ", "), "\n")
cat("Especies en A:", paste(species_in_A, collapse = ", "), "\n")
cat("Especies comunes:", paste(common_species, collapse = ", "), "\n")

if (length(common_species) < length(species_in_data)) {
  missing <- setdiff(species_in_data, species_in_A)
  warning("Especie(s) sin filogenia: ", paste(missing, collapse = ", "))
}

# Filtrar A a especies comunes
A <- A[common_species, common_species]

# ---- 0.4 Crear variable phylo_species (copia dedicada para termino filogenetico) ----
# brms no permite usar el mismo factor dos veces con distinta covarianza
df.t1$phylo_species <- factor(df.t1$species, levels = common_species)
df.t2$phylo_species <- factor(df.t2$species, levels = common_species)

cat("Variable phylo_species creada en df.t1 y df.t2\n")
cat("Niveles phylo_species:", levels(df.t1$phylo_species), "\n")

# ---- 0.5 Verificar que la variable no tiene NA ----
na_phylo_t1 <- sum(is.na(df.t1$phylo_species))
na_phylo_t2 <- sum(is.na(df.t2$phylo_species))
cat("NA en phylo_species PRE:", na_phylo_t1, "\n")
cat("NA en phylo_species POST:", na_phylo_t2, "\n")

if (na_phylo_t1 > 0 || na_phylo_t2 > 0) {
  stop("Hay NA en phylo_species. Revisar niveles del factor.")
}

cat("\n========== FASE 0 completada ==========\n")

# ============================================================
# FASE 1: Smoke tests
# ============================================================
cat("\n========== FASE 1: Smoke tests ==========\n")

if (!RUN_SMOKE_TEST) {
  cat("RUN_SMOKE_TEST = FALSE. Saltando smoke tests.\n")
} else {

  # ---- 1.1 Especificaciones de modelos ----

  # Modelo 1: Filogenia en pendientes + especie libre (sin interacciones triples)
  form_het_phylo <- bf(
    Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
      (0 + time_s | gr(phylo_species, cov = A)) +
      (0 + time_s | species) +
      (1 + time_s | codigo) +
      (1 | id_bellota)
  )

  # Modelo 2: Filogenia + interacciones triples explícitas
  form_het_phylo_v2 <- bf(
    Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
      time_s:Dim.1:species +
      time_s:Dim.2:species +
      time_s:Dim.3:species +
      (0 + time_s | gr(phylo_species, cov = A)) +
      (1 + time_s | codigo) +
      (1 | id_bellota)
  )

  cat("\n-- Formulas de modelos --\n")
  cat("\nM_het_1 (filogenia + especie libre, sin triples):\n")
  cat(deparse(form_het_phylo), "\n")
  cat("\nM_het_2 (filogenia + interacciones triples):\n")
  cat(deparse(form_het_phylo_v2), "\n")

  # ---- 1.2 Priors ----
  priors <- c(
    prior(normal(40, 20), class = "Intercept"),
    prior(normal(0, 10),  class = "b"),
    prior(student_t(3, 0, 10), class = "sd"),
    prior(student_t(3, 0, 10), class = "sigma"),
    prior(lkj(2), class = "cor")
  )

  # ---- 1.3 Funcion de smoke test ----
  smoke_test <- function(formula, name, data, seed) {
    cat("\n", strrep("=", 60), "\n")
    cat("SMOKE TEST:", name, "\n")
    cat(strrep("=", 60), "\n")
    cat("Formula:", deparse(formula), "\n")
    cat("Datos:", nrow(data), "obs,", nlevels(data$species), "especies\n")
    cat("Iteraciones:", ITER_SMOKE, "| Warmup:", WARMUP_SMOKE,
        "| Cadenas:", N_CHAINS_SMOKE, "\n")

    t_start <- Sys.time()

    fit <- tryCatch(
      brm(
        formula,
        data = data,
        data2 = list(A = A),
        family = gaussian(),
        prior = priors,
        iter = ITER_SMOKE,
        warmup = WARMUP_SMOKE,
        chains = N_CHAINS_SMOKE,
        cores = 1,  # Secuencial para evitar problemas en Windows
        control = list(
          adapt_delta = 0.95,
          max_treedepth = 12
        ),
        seed = seed,
        refresh = 0,  # Sin output de progreso
        silent = 2,
        backend = "cmdstanr"
      ),
      error = function(e) {
        cat("\n*** ERROR en compilacion/muestreo ***\n")
        cat(conditionMessage(e), "\n")
        return(NULL)
      }
    )

    t_end <- Sys.time()
    elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

    if (is.null(fit)) {
      cat("RESULTADO: FALLO (", round(elapsed, 1), "min)\n")
      return(invisible(NULL))
    }

    cat("\nRESULTADO: COMPILADO Y MUESTREADO (", round(elapsed, 1), "min)\n")

    # Diagnostico basico
    cat("\n-- Diagnostico de convergencia --\n")

    # Rhat y ESS para efectos fijos
    draws <- as_draws_df(fit)
    b_cols <- grep("^b_", colnames(draws), value = TRUE)

    if (length(b_cols) > 0) {
      diag_summary <- tibble(
        parametro = b_cols,
        Rhat = sapply(b_cols, function(col) {
          vals <- split(draws[[col]], draws$.chain)
          mat <- do.call(cbind, vals)
          posterior::rhat(mat)
        }),
        ESS_bulk = sapply(b_cols, function(col) {
          vals <- split(draws[[col]], draws$.chain)
          mat <- do.call(cbind, vals)
          posterior::ess_bulk(mat)
        }),
        ESS_tail = sapply(b_cols, function(col) {
          vals <- split(draws[[col]], draws$.chain)
          mat <- do.call(cbind, vals)
          posterior::ess_tail(mat)
        }),
        estado = case_when(
          Rhat > 1.05 ~ "CONVERGENCIA MALA",
          ESS_bulk < 400 ~ "ESS BAJO",
          ESS_tail < 400 ~ "ESS TAIL BAJO",
          TRUE ~ "ok"
        )
      )

      cat("\nEfectos fijos:\n")
      print(diag_summary, n = Inf)

      n_problemas <- sum(diag_summary$estado != "ok")
      if (n_problemas > 0) {
        cat("\n*** ALERTA:", n_problemas, "parametros con problemas ***\n")
      } else {
        cat("\nTodos los efectos fijos convergieron correctamente.\n")
      }
    }

    # Verificar divergencias
    n_div <- sum(draws$`.divergent__` == 1, na.rm = TRUE)
    n_total <- nrow(draws)
    cat("\nDivergencias:", n_div, "/", n_total,
        "(", round(100 * n_div / n_total, 2), "%)\n")

    # Treedepth
    treedepth_max <- max(draws$`.treedepth__`, na.rm = TRUE)
    cat("Treedepth maximo:", treedepth_max, "\n")

    # Guardar modelo
    saveRDS(fit, file.path("00-data/phylo", paste0(name, ".rds")))
    cat("\nModelo guardado en: 00-data/phylo/", name, ".rds\n", sep = "")

    return(fit)
  }

  # ---- 1.4 Ejecutar smoke tests ----
  cat("\n-- Smoke test M_het_1 (PRE) --\n")
  m_het_1_pre <- smoke_test(form_het_phylo, "m_het_1_pre_smoke", df.t1, SEED_BASE)

  cat("\n-- Smoke test M_het_2 (PRE) --\n")
  m_het_2_pre <- smoke_test(form_het_phylo_v2, "m_het_2_pre_smoke", df.t1, SEED_BASE + 1)

  # ---- 1.5 Resumen comparativo ----
  cat("\n", strrep("=", 60), "\n")
  cat("RESUMEN SMOKE TESTS\n")
  cat(strrep("=", 60), "\n")

  resultados <- list(
    M_het_1 = list(
      formula = "sin interacciones triples",
      compilado = !is.null(m_het_1_pre),
      divergencias = if (!is.null(m_het_1_pre)) {
        draws <- as_draws_df(m_het_1_pre)
        sum(draws$.divergent__ == 1, na.rm = TRUE)
      } else NA_integer_,
      treedepth = if (!is.null(m_het_1_pre)) {
        max(as_draws_df(m_het_1_pre)$.treedepth__, na.rm = TRUE)
      } else NA_integer_
    ),
    M_het_2 = list(
      formula = "con interacciones triples",
      compilado = !is.null(m_het_2_pre),
      divergencias = if (!is.null(m_het_2_pre)) {
        draws <- as_draws_df(m_het_2_pre)
        sum(draws$.divergent__ == 1, na.rm = TRUE)
      } else NA_integer_,
      treedepth = if (!is.null(m_het_2_pre)) {
        max(as_draws_df(m_het_2_pre)$.treedepth__, na.rm = TRUE)
      } else NA_integer_
    )
  )

  for (nm in names(resultados)) {
    cat("\n", nm, ":\n")
    cat("  Formula:", resultados[[nm]]$formula, "\n")
    cat("  Compilado:", resultados[[nm]]$compilado, "\n")
    cat("  Divergencias:", resultados[[nm]]$divergencias, "\n")
    cat("  Treedepth max:", resultados[[nm]]$treedepth, "\n")
  }

  # ---- 1.6 Recomendaciones ----
  cat("\n-- Recomendaciones para tirada completa --\n")

  if (!is.null(m_het_1_pre)) {
    cat("M_het_1: LISTO para tirada completa (iter=4000, warmup=2000, chains=4)\n")
  } else {
    cat("M_het_1: REVISAR antes de tirada completa (error en smoke test)\n")
  }

  if (!is.null(m_het_2_pre)) {
    cat("M_het_2: LISTO para tirada completa (iter=4000, warmup=2000, chains=4)\n")
  } else {
    cat("M_het_2: REVISAR antes de tirada completa (error en smoke test)\n")
  }

  cat("\n========== FASE 1 completada ==========\n")
}
