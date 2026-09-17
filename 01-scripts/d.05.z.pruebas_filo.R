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
#   1. Tirada completa (iter=6000, warmup=3000, chains=4 en paralelo)
#
# Configuracion:
#   - Cadenas en paralelo con cmdstanr (4 nucleos; ya no son secuenciales)
#   - Backend cmdstanr (mas rapido que rstan)
#   - adapt_delta=0.999, max_treedepth=15 (para reducir divergencias)
#   - Priors mas informativos en sd (student_t(3,0,2.5)) para estabilizar
#     grupos pequenos (8 especies, 16 codigos)
# ============================================================

N_CHAINS_FULL  <- 4
ITER_FULL      <- 6000
WARMUP_FULL    <- 3000
N_CORES_PAR    <- 4          # cadenas en paralelo
ADAPT_DELTA    <- 0.999
MAX_TREEDEPTH  <- 15
SEED_BASE      <- 123

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(ape)
})

# cmdstanr es informativo: si falla, se avisa pero el script sigue y el
# error real de brm() quedara capturado en el smoke_test.
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("cmdstanr version:", as.character(packageVersion("cmdstanr")), "\n")
  cat("cmdstan path:", cmdstanr::cmdstan_path(), "\n")
  tryCatch(
    cat("cmdstan version:", cmdstanr::cmdstan_version(), "\n"),
    error = function(e) cat("  (cmdstan no localizado)\n")
  )
} else {
  warning("cmdstanr no está instalado; los ajustes con backend='cmdstanr' fallaran.")
}

dir.create("00-data/phylo", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FASE 0: Preparacion
# ============================================================
cat("\n========== FASE 0: Preparacion ==========\n")

# ---- 0.0 Procedencias excluidas ----
# IL3 (Quercus ilex) se ELIMINA del analisis: tiene solo 60 observaciones
# en fase PRE (t < 94 h) frente a 150 en el resto de procedencias, lo que
# impide la estimacion adecuada de sus parametros. Se retira TAMBIEN de la
# fase POST para que las comparaciones pre-post se hagan sobre el mismo
# conjunto de procedencias (por eso el filtro esta ANTES de la particion).
PROCEDENCIAS_EXCLUIDAS <- "IL3"

# ---- 0.1 Cargar datos (siempre frescos; no reutilizar dataframes viejos
#   del entorno, que pueden tener columnas obsoletas) ----
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd     <- read.csv("00-data/famd_ind_coord.csv")
df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  dplyr::filter(!codigo %in% PROCEDENCIAS_EXCLUIDAS) |>
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

# Verificar que la procedencia excluida no esta en ninguna fase
stopifnot(!any(df.t1$codigo %in% PROCEDENCIAS_EXCLUIDAS),
          !any(df.t2$codigo %in% PROCEDENCIAS_EXCLUIDAS))
cat("Procedencias excluidas:", paste(PROCEDENCIAS_EXCLUIDAS, collapse = ", "), "\n")

needed_cols <- c("time_s", "species", "codigo", "id_bellota", "Dim.1", "Dim.2", "Dim.3")
missing_t1  <- setdiff(needed_cols, colnames(df.t1))
missing_t2  <- setdiff(needed_cols, colnames(df.t2))
if (length(missing_t1) > 0 || length(missing_t2) > 0) {
  stop("Faltan columnas requeridas: ",
       paste(unique(c(missing_t1, missing_t2)), collapse = ", "))
}

cat("Datos cargados:", nrow(df.t1), "obs PRE,", nrow(df.t2), "obs POST\n")

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
# FASE 1: Tirada completa
# ============================================================
cat("\n========== FASE 1: Tirada completa ==========\n")

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
    prior(student_t(3, 0, 2.5), class = "sd"),
    prior(student_t(3, 0, 2.5), class = "sigma"),
    prior(lkj(2), class = "cor")
  )

  # ---- 1.3 Funcion de smoke test ----
  smoke_test <- function(formula, name, data, seed) {
    cat("\n", strrep("=", 60), "\n")
    cat("SMOKE TEST:", name, "\n")
    cat(strrep("=", 60), "\n")
    cat("Formula:", deparse(formula), "\n")
    cat("Datos:", nrow(data), "obs,", nlevels(data$species), "especies\n")
    cat("Iteraciones:", ITER_FULL, "| Warmup:", WARMUP_FULL,
        "| Cadenas:", N_CHAINS_FULL, "| Cores en paralelo:", N_CORES_PAR,
        "\n")
    cat("adapt_delta:", ADAPT_DELTA, "| max_treedepth:", MAX_TREEDEPTH, "\n")

    t_start <- Sys.time()

    fit_err <- NULL
    fit <- tryCatch(
      brm(
        formula,
        data = data,
        data2 = list(A = A),
        family = gaussian(),
        prior = priors,
        iter = ITER_FULL,
        warmup = WARMUP_FULL,
        chains = N_CHAINS_FULL,
        cores = N_CORES_PAR,  # cadenas en paralelo con cmdstanr
        control = list(
          adapt_delta = ADAPT_DELTA,
          max_treedepth = MAX_TREEDEPTH
        ),
        seed = seed,
        refresh = 100,  # Progreso cada 100 iteraciones
        backend = "cmdstanr"
      ),
      error = function(e) {
        fit_err <<- conditionMessage(e)
        cat("\n*** ERROR en compilacion/muestreo ***\n")
        cat(fit_err, "\n")
        return(NULL)
      }
    )

    t_end <- Sys.time()
    elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

    if (is.null(fit)) {
      cat("\n!!!! FALLO:", name, " -> ", fit_err, "\n", sep = "")
      cat("RESULTADO: FALLO (", round(elapsed, 1), "min)\n", sep = "")
      return(list(fit = NULL, error = fit_err, elapsed = elapsed))
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

    # Trace plots en PNG
    cat("\n-- Generando trace plots --\n")
    trace_params <- grep(
      "^(b_|sd_|cor_|sigma)",
      posterior::variables(fit),
      value = TRUE
    )
    p_trace <- tryCatch(
      bayesplot::mcmc_trace(
        fit,
        pars = trace_params,
        facet_args = list(ncol = 4, scales = "free")
      ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7)),
      error = function(e) {
        cat("Trace plot no disponible:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(p_trace)) {
      file_png <- file.path("00-data/phylo", paste0("trace_", name, ".png"))
      n_rows <- ceiling(length(trace_params) / 4)
      ggplot2::ggsave(
        filename = file_png,
        plot     = p_trace,
        width    = 14,
        height   = max(6, 2 * n_rows),
        dpi      = 150,
        limitsize = FALSE
      )
      cat("Trace plots guardados en: ", file_png, "\n", sep = "")
    }

    return(list(fit = fit, error = NULL, elapsed = elapsed))
  }

  # ---- 1.4 Ejecutar tirada completa ----
  t_total <- Sys.time()
  cat("\n-- Tirada completa M_het_1 (PRE) --\n")
  res_het_1_pre  <- smoke_test(form_het_phylo, "m_het_1_pre", df.t1, SEED_BASE)

  cat("\n-- Tirada completa M_het_2 (PRE) --\n")
  res_het_2_pre  <- smoke_test(form_het_phylo_v2, "m_het_2_pre", df.t1, SEED_BASE + 1)

  cat("\n-- Tirada completa M_het_1 (POST) --\n")
  res_het_1_post <- smoke_test(form_het_phylo, "m_het_1_post", df.t2, SEED_BASE + 2)

  cat("\n-- Tirada completa M_het_2 (POST) --\n")
  res_het_2_post <- smoke_test(form_het_phylo_v2, "m_het_2_post", df.t2, SEED_BASE + 3)

  # ---- 1.5 Resumen comparativo ----
  cat("\n", strrep("=", 60), "\n")
  cat("RESUMEN TIRADA COMPLETA\n")
  cat(strrep("=", 60), "\n")

  resultados <- list(
    M_het_1_PRE = list(
      formula = "sin interacciones triples",
      fit = res_het_1_pre$fit,
      error = res_het_1_pre$error
    ),
    M_het_2_PRE = list(
      formula = "con interacciones triples",
      fit = res_het_2_pre$fit,
      error = res_het_2_pre$error
    ),
    M_het_1_POST = list(
      formula = "sin interacciones triples",
      fit = res_het_1_post$fit,
      error = res_het_1_post$error
    ),
    M_het_2_POST = list(
      formula = "con interacciones triples",
      fit = res_het_2_post$fit,
      error = res_het_2_post$error
    )
  )

  for (nm in names(resultados)) {
    r <- resultados[[nm]]
    cat("\n", nm, ":\n", sep = "")
    cat("  Formula:", r$formula, "\n")
    if (!is.null(r$fit)) {
      draws <- as_draws_df(r$fit)
      cat("  Compilado: TRUE\n")
      cat("  Divergencias:", sum(draws$.divergent__ == 1, na.rm = TRUE), "\n")
      cat("  Treedepth max:", max(draws$.treedepth__, na.rm = TRUE), "\n")
    } else {
      cat("  Compilado: FALSE\n")
      cat("  Error: ", r$error, "\n", sep = "")
    }
  }

  t_total_elapsed <- as.numeric(difftime(Sys.time(), t_total, units = "mins"))
  cat("\nTiempo total tirada:", round(t_total_elapsed, 1), "min\n")

  cat("\n========== FASE 1 completada ==========\n")

# ============================================================
# FASE 2: Validacion, comparacion y exportacion de resultados
#
# Carga los 4 modelos ajustados en FASE 1 y produce:
#   - resumenes (txt) y tablas de efectos fijos / varianzas (CSV)
#   - senal filogenetica (CCI) por modelo (CSV)
#   - comparacion LOO + pareto-k (CSV + figura)
#   - pendientes totales por especie para M_het_2 (CSV + forest plot)
#   - posterior predictive checks (PNG)
# Salidas en: 00-data/phylo/ y 07-img/
# ============================================================
cat("\n========== FASE 2: Validacion y comparacion ==========\n")

dir.create("07-img", showWarnings = FALSE, recursive = TRUE)

# ---- 2.1 Cargar modelos ----
# Se recargan de disco: asi la FASE 2 es independiente de la sesion
# de la FASE 1 (los RDS fueron guardados por smoke_test()).
modelos <- list(
  m_het_1_pre  = readRDS("00-data/phylo/m_het_1_pre.rds"),
  m_het_1_post = readRDS("00-data/phylo/m_het_1_post.rds"),
  m_het_2_pre  = readRDS("00-data/phylo/m_het_2_pre.rds"),
  m_het_2_post = readRDS("00-data/phylo/m_het_2_post.rds")
)

# ---- 2.2 Resumen por modelo (a pantalla y a log) ----
sink("00-data/phylo/summary_modelos.txt")
for (nm in names(modelos)) {
  cat("\n===== Summary:", nm, "=====\n")
  print(summary(modelos[[nm]]))
}
sink()

# ---- 2.3 Senal filogenetica (CCI) ----
# Unica funcion: el calculo de la CCI se repite 4 veces con denominador
# distinto segun si el modelo tiene termino de especie libre (M_het_1).
calcular_cci <- function(fit, con_especie_libre = FALSE) {
  d <- as_draws_df(fit)
  phylo_var <- d$sd_phylo_species__time_s^2
  cod_var   <- d$sd_codigo__time_s^2
  bell_var  <- d$sd_id_bellota__Intercept^2
  sigma2    <- d$sigma^2
  if (con_especie_libre) {
    sp_var <- d$sd_species__time_s^2
    var_especie <- phylo_var + sp_var
    prop_phylo_sp <- phylo_var / var_especie
  } else {
    var_especie  <- phylo_var
    prop_phylo_sp <- NA_real_
  }
  total <- var_especie + cod_var + bell_var + sigma2
  cci_phylo <- phylo_var / total
  cci_especie <- var_especie / total
  tibble(
    cci_phylo_media    = median(cci_phylo),
    cci_phylo_ic_lo    = quantile(cci_phylo, 0.05)   |> unname(),
    cci_phylo_ic_hi    = quantile(cci_phylo, 0.95)   |> unname(),
    cci_especie_media  = median(cci_especie),
    cci_especie_ic_lo  = quantile(cci_especie, 0.05) |> unname(),
    cci_especie_ic_hi  = quantile(cci_especie, 0.95) |> unname(),
    prop_phylo_en_especie = median(prop_phylo_sp)
  )
}

cci_df <- bind_rows(
  m_het_1_pre  = calcular_cci(modelos$m_het_1_pre,  con_especie_libre = TRUE),
  m_het_1_post = calcular_cci(modelos$m_het_1_post, con_especie_libre = TRUE),
  m_het_2_pre  = calcular_cci(modelos$m_het_2_pre,  con_especie_libre = FALSE),
  m_het_2_post = calcular_cci(modelos$m_het_2_post, con_especie_libre = FALSE),
  .id = "modelo"
)
write.csv(cci_df, "00-data/phylo/cci_filogenetica.csv", row.names = FALSE)
cat("\n-- Senal filogenetica (CCI) --\n")
print(cci_df)

# ---- 2.4 Efectos fijos y varianzas por nivel (CSV por modelo) ----
for (nm in names(modelos)) {
  fe <- as.data.frame(fixef(modelos[[nm]]))
  fe$parametro <- rownames(fe)
  write.csv(fe, paste0("00-data/phylo/fixef_", nm, ".csv"), row.names = FALSE)

  vc <- VarCorr(modelos[[nm]])
  sd_filas <- list()
  for (lvl in names(vc)) {
    if ("sd" %in% names(vc[[lvl]])) {
      s <- as.data.frame(vc[[lvl]]$sd)
      s$nivel  <- lvl
      s$efecto <- rownames(s)
      rownames(s) <- NULL
      sd_filas[[lvl]] <- s
    }
  }
  write.csv(dplyr::bind_rows(sd_filas),
            paste0("00-data/phylo/varcom_", nm, ".csv"), row.names = FALSE)
}
cat("\nEfectos fijos y varcom exportados a 00-data/phylo/fixef_*.csv y varcom_*.csv\n")

# ---- 2.5 Comparacion LOO + pareto-k (loop explicito por fase) ----
for (fase in c("pre", "post")) {
  cat("\n===== Comparacion LOO:", fase, "=====\n")
  nm1 <- paste0("m_het_1_", fase)
  nm2 <- paste0("m_het_2_", fase)

  loo_1 <- loo(modelos[[nm1]])
  loo_2 <- loo(modelos[[nm2]])
  tabla_loo <- loo::loo_compare(loo_1, loo_2)
  print(tabla_loo)
  write.csv(as.data.frame(tabla_loo),
            paste0("00-data/phylo/loo_comp_", fase, ".csv"))

  # pareto-k por observacion
  # brms descarto las filas con NA en Moisture_content al ajustar, de modo
  # que de df.t1/df.t2 quedan las mismas observaciones que las de pareto_k.
  dat <- (if (fase == "pre") df.t1 else df.t2) |>
    dplyr::filter(!is.na(Moisture_content))
  dat$k_1 <- loo_1$diagnostics$pareto_k
  dat$k_2 <- loo_2$diagnostics$pareto_k
  dat$pareto_status <- dplyr::case_when(
    dat$k_1 > 0.7 & dat$k_2 > 0.7 ~ "Problematic in both",
    dat$k_1 > 0.7                  ~ paste0("Problematic in ", nm1),
    dat$k_2 > 0.7                  ~ paste0("Problematic in ", nm2),
    TRUE ~ "Not problematic"
  )

  p_pareto <- ggplot(dat, aes(x = time, y = Moisture_content)) +
    geom_line(aes(group = id_bellota), alpha = 0.25) +
    geom_point(aes(colour = pareto_status), alpha = 0.6) +
    facet_wrap(~ codigo) +
    labs(x = "Time", y = "Moisture content (%)",
         colour = "Pareto k") +
    theme_classic()
  ggsave(paste0("07-img/pareto_k_", fase, ".png"),
         p_pareto, width = 14, height = 8, dpi = 150)
  cat("Figura pareto guardada en 07-img/pareto_k_", fase, ".png\n", sep = "")
}

# ---- 2.6 Pendientes totales por especie (M_het_2) ----
# La pendiente time_s:Dim de cada especie = base (coccifera) + ajuste.
# El coeficiente de la interaccion triple es siempre la DESVIACION
# respecto al nivel base del factor species.
especies        <- levels(modelos$m_het_2_pre$data$species)
especie_base    <- especies[1]           # coccifera (ordinal: brms usa el 1er nivel)
especies_ajuste <- especies[-1]

draws2      <- as_draws_df(modelos$m_het_2_pre)
pend_filas  <- list()
for (dim in c("Dim.1", "Dim.2", "Dim.3")) {
  base_col <- paste0("b_time_s:", dim)
  totales  <- list()
  totales[[especie_base]] <- draws2[[base_col]]
  for (sp in especies_ajuste) {
    adj_col <- paste0("b_time_s:", dim, ":species", gsub(" ", "", sp))
    totales[[sp]] <- draws2[[base_col]] + draws2[[adj_col]]
  }
  for (sp in especies) {
    pend_filas[[length(pend_filas) + 1]] <- tibble(
      dimension = dim,
      especie   = sp,
      pendiente = median(totales[[sp]]),
      ic_lo     = quantile(totales[[sp]], 0.025) |> unname(),
      ic_hi     = quantile(totales[[sp]], 0.975) |> unname()
    )
  }
}
pend_df <- dplyr::bind_rows(pend_filas)
write.csv(pend_df, "00-data/phylo/pendientes_por_especie_m_het_2_pre.csv",
          row.names = FALSE)

p_forest <- ggplot(pend_df, aes(x = pendiente, y = especie)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = ic_lo, xmax = ic_hi)) +
  facet_wrap(~ dimension, scales = "free_x") +
  labs(x = "Pendiente time_s:Dim (tasa de desecacion)", y = NULL) +
  theme_classic()
ggsave("07-img/forest_pendientes_m_het_2_pre.png",
       p_forest, width = 10, height = 6, dpi = 150)
cat("Pendientes por especie y forest plot guardados\n")

# ---- 2.7 Posterior predictive checks (PNG) ----
for (nm in names(modelos)) {
  p_pp <- pp_check(modelos[[nm]]) + ggplot2::ggtitle(nm)
  ggsave(paste0("07-img/pp_check_", nm, ".png"),
         p_pp, width = 8, height = 6, dpi = 150)
}
cat("Posterior predictive checks guardados en 07-img/pp_check_*.png\n")

cat("\n========== FASE 2 completada ==========\n")


