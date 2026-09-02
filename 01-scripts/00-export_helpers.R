# 00-export_helpers.R
# Funciones auxiliares para exportar resultados de modelos glmmTMB
# (tablas, objetos y figuras) de cara al articulo / material suplementario.
# Usadas por d.04.model_species.R y d.05.model_traits.R.
# Todas las rutas son relativas a la raiz del repositorio.

# Guardar objetos de modelo en 00-data/models/ ----
save_models <- function(models) {
  dir.create("00-data/models", showWarnings = FALSE, recursive = TRUE)
  for (nm in names(models)) {
    saveRDS(models[[nm]], file.path("00-data/models", paste0(nm, ".rds")))
  }
  cat("Modelos guardados en 00-data/models/:\n  ", paste(names(models), collapse = ", "), "\n", sep = "")
  invisible(NULL)
}

# Tabla de coeficientes (efectos fijos) de una lista con nombre de modelos ----
coef_table <- function(models) {
  suppressPackageStartupMessages(library(parameters))
  tabs <- lapply(names(models), function(nm) {
    parameters::model_parameters(models[[nm]], effects = "fixed", ci = 0.95) |>
      as.data.frame() |>
      dplyr::mutate(modelo = nm) |>
      dplyr::select(modelo, Parameter, Coefficient, SE, CI_low, CI_high, p)
  })
  dplyr::bind_rows(tabs)
}

# Componentes de varianza de una lista con nombre de modelos ----
varcomp_table <- function(models) {
  tabs <- lapply(names(models), function(nm) {
    mod <- models[[nm]]
    vc <- glmmTMB::VarCorr(mod)$cond
    re <- lapply(names(vc), function(g) {
      data.frame(
        nivel    = g,
        varianza = sum(diag(vc[[g]])),
        terminos = paste(rownames(vc[[g]]), collapse = ", ")
      )
    })
    re_df <- do.call(rbind, re)
    res_df <- data.frame(nivel = "Residual", varianza = sigma(mod)^2, terminos = "sigma^2")
    rbind(re_df, res_df) |>
      dplyr::mutate(modelo = nm)
  })
  dplyr::bind_rows(tabs) |>
    dplyr::select(modelo, nivel, varianza, terminos)
}

# Barplot de varianza por nivel ----
plot_varcomp <- function(vc, archivo, titulo) {
  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(vc, aes(x = reorder(nivel, varianza), y = varianza)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    facet_wrap(~modelo, scales = "free_y") +
    labs(x = NULL, y = "Varianza", title = titulo) +
    theme_minimal()
  ggplot2::ggsave(archivo, p, width = 8, height = 5)
  cat("Barplot de varianza guardado en ", archivo, "\n", sep = "")
  invisible(NULL)
}

# check_model de cada modelo de una lista ----
plot_check_model <- function(models, sufijo, dir = "07-img") {
  suppressPackageStartupMessages(library(performance))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(models)) {
    res <- tryCatch(performance::check_model(models[[nm]]), error = function(e) e)
    if (inherits(res, "error")) {
      warning("check_model fallo para ", nm, ": ", conditionMessage(res))
      next
    }
    archivo <- file.path(dir, paste0("check_model_", sufijo, "_", nm, ".png"))
    png(archivo, width = 2000, height = 1600, res = 150)
    if (!is.character(attr(res, "theme"))) {
      attr(res, "theme") <- "see::theme_lucid"
    }
    p <- tryCatch(plot(res), error = function(e) e)
    if (inherits(p, "error")) {
      dev.off()
      warning("plot(check_model) fallo para ", nm, ": ", conditionMessage(p))
      next
    }
    print(p)
    dev.off()
    cat("check_model guardado en ", archivo, "\n", sep = "")
  }
}

# Observed vs fitted de cada modelo de una lista ----
plot_obs_fitted <- function(models, sufijo, dir = "07-img") {
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(models)) {
    mod <- models[[nm]]
    mf <- if (!is.null(mod$frame)) mod$frame else tryCatch(model.frame(mod), error = function(e) NULL)
    if (is.null(mf)) {
      warning("No se pudo obtener el model frame de ", nm)
      next
    }
    d <- data.frame(observed = mf[[1]], fitted = stats::fitted(mod))
    p <- ggplot(d, aes(x = fitted, y = observed)) +
      geom_point(alpha = 0.15, size = 0.8) +
      geom_abline(slope = 1, intercept = 0, colour = "firebrick", linetype = "dashed") +
      labs(x = "Fitted", y = "Observed", title = paste("Observed vs fitted -", nm)) +
      theme_minimal()
    archivo <- file.path(dir, paste0("obs_fitted_", sufijo, "_", nm, ".png"))
    ggplot2::ggsave(archivo, p, width = 6, height = 6)
    cat("Observed vs fitted guardado en ", archivo, "\n", sep = "")
  }
}

# Graficos DHARMa de cada modelo de una lista ----
plot_dharma <- function(models, sufijo, dir = "07-img", n = 250) {
  suppressPackageStartupMessages(library(DHARMa))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(models)) {
    sim <- tryCatch(
      DHARMa::simulateResiduals(fittedModel = models[[nm]], n = n),
      error = function(e) e
    )
    if (inherits(sim, "error")) {
      warning("DHARMa fallo para ", nm, ": ", conditionMessage(sim))
      next
    }
    archivo <- file.path(dir, paste0("dharma_", sufijo, "_", nm, ".png"))
    png(archivo, width = 1400, height = 700, res = 150)
    plot(sim)
    dev.off()
    cat("DHARMa guardado en ", archivo, "\n", sep = "")
  }
}
