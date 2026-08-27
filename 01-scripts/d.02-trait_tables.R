df.bellotas <- read.csv(file = "00-data/df.bellotas.csv")

library(tidyverse) 
library(emmeans)     # Calculo de medias marginales y tendencias
library(performance) # Comprobación de supuestos
library(multcomp)    # Asignación de letras de significancia estadística
library(nlme)        # Modelos gls (generalized least squares)
library(glmmTMB)
library(e1071)   # Calcular kurtosis
library(report)  # Para hacer reportes automáticos
library(sjPlot)
library(DHARMa)
library(MuMIn)

# Crea una dataset con solo una observación por bellota y las variables de interés
## -> Si quieres otras variables cambia el siguiente vector
vars <- c("id_bellota","especie","codigo", "procedencia",
          "longitud..mm.","diametro_total","peso_seco",
          "Area_cicatriz_mm2","peso_seco_pr","Area_estimada_cm2", "SPM_g_cm2", "rajas_pericarpo")
## Elimina filas repetidas en la tabla
df_unique <- df.bellotas %>%
  distinct(id_bellota, .keep_all = TRUE) %>%
  dplyr::select(all_of(vars)) %>%
  drop_na() %>%  # elimina filas con NA en cualquier columna
  mutate(
    especie = factor(especie),
    codigo = factor(codigo),
    procedencia = factor(procedencia),
    Area_cicatriz_cm2 = Area_cicatriz_mm2 * 0.01  # pasar a cm2
  ) |> 
  rename(species = especie, 
         provenance = procedencia, 
         prov_code = codigo, 
         length_mm = longitud..mm., 
         diameter_mm = diametro_total, 
         dry_weight = peso_seco, 
         scar_area_mm2 = Area_cicatriz_mm2, 
         scar_area_cm2 = Area_cicatriz_cm2,
         pericarp_dry_weight = peso_seco_pr, 
         surface_cm2 = Area_estimada_cm2, 
         pericarp_rupture = rajas_pericarpo
         )

df <- df_unique

# 1. Fit models ----
## 1. Dry weight ----
glmm.dw.gaus <- glmmTMB(dry_weight ~ species + (1|species:provenance),
                              data = df)
glmm.dw.gamm <- glmmTMB(dry_weight ~ species + (1|species:provenance),
                              data = df, 
                              family = Gamma(link = "log"))

# compare_performance(glmm.dw.gaus, glmm.dw.gamm)
# gamma chosen
tryCatch(
  easystats::model_dashboard(glmm.dw.gamm,
                             output_dir = "06-html/",
                             output_file = "modeldashboard_dryweight.html"),
  error = function(e) warning("model_dashboard omitido (", conditionMessage(e), ")")
)

## 2. Specific pericarp mass ----
glmm.spm <- glmmTMB(SPM_g_cm2 ~ species + (1|species:provenance),
                    data = df)

tryCatch(
  easystats::model_dashboard(glmm.spm,
                             output_dir = "06-html/",
                             output_file = "modeldashboard_spm.html"),
  error = function(e) warning("model_dashboard omitido (", conditionMessage(e), ")")
)

## 3. Scar area ----
glmm.scar <- glmmTMB(scar_area_mm2 ~ species*surface_cm2 + (1|species:provenance),
                     data = df, 
                     family = Gamma(link = "log"))
options(na.action = "na.fail")
dd <- dredge(glmm.scar)
glmm.scar1 <- get.models(dd, subset = 1)[[1]]

tryCatch(
  easystats::model_dashboard(glmm.scar1,
                             output_dir = "06-html/",
                             output_file = "modeldashboard_scar.html"),
  error = function(e) warning("model_dashboard omitido (", conditionMessage(e), ")")
)

## 4. Seed coat ratio ----
glmm.scr <- glmmTMB(pericarp_dry_weight ~ dry_weight*species + 
                       (1|species:provenance),
                     data = df)

dd <- dredge(glmm.scr)
glmm.scr1 <- get.models(dd, subset = 1)[[1]]

tryCatch(
  easystats::model_dashboard(glmm.scr1,
                             output_dir = "06-html/",
                             output_file = "modeldashboard_scr.html"),
  error = function(e) warning("model_dashboard omitido (", conditionMessage(e), ")")
)

## 5. Pericarp rupture ----
df.cracks <- df |> 
  mutate(pericarp_rupture = if_else(pericarp_rupture == 2, 1, pericarp_rupture))

glmm.cracks <- glmmTMB(pericarp_rupture ~ species + (1|species:provenance), 
                      data = df.cracks , 
                      family = binomial(link = "logit"))

tryCatch(
  easystats::model_dashboard(glmm.cracks,
                             output_dir = "06-html/",
                             output_file = "modeldashboard_cracks.html"),
  error = function(e) warning("model_dashboard omitido (", conditionMessage(e), ")")
)

# 2. Generate table ----
## 1. Calculate marginal means ----
modelos <- list(glmm.dw.gamm = glmm.dw.gamm, 
                glmm.spm = glmm.spm, 
                glmm.scar1 = glmm.scar1, 
                glmm.scr1 = glmm.scr1, 
                glmm.cracks = glmm.cracks)

metadata <- data.frame(
  modelo = c("glmm.dw.gamm", 
             "glmm.spm", 
             "glmm.scar1", 
             "glmm.scr1", 
             "glmm.cracks"), 
  response = c("Dry weight",
               "Specific pericarp mass",
               "Scar area",
               "Pericarp dry weitgh",
               "Pericarp rupture"), 
  units = c("g",
            "g/cm²",
            "cm²",
            "g",
            "probability"), 
  type = as.factor(c("anova", 
           "anova", 
           "anova", 
           "ancova", 
           "binomial"))
)

library(emmeans)
library(multcomp)

extraer_metodo <- function(vec) {
  sub(" for.*", "", 
      sub("P value adjustment: ", "", 
          grep("P value adjustment", attr(vec, "mesg"), value = TRUE)))
}

# Detecta si el modelo tiene una covariable numérica además de 'factor_var'
tiene_covariable <- function(modelo, factor_var = "species") {
  datos <- model.frame(modelo)
  predictores <- names(datos)[-1]
  candidatos <- setdiff(predictores, factor_var)
  any(sapply(datos[candidatos], is.numeric))
}

tiene_interaccion <- function(modelo, factor_var = "species") {
  terminos <- attr(terms(modelo), "term.labels")
  any(grepl(paste0("(^|:)", factor_var, "(:|$)"), terminos) & grepl(":", terminos))
}

# Emmeans para modelos con interacción factor × covariable:
# evalúa el factor en los percentiles 25, 50 y 75 de la covariable continua.
emmeans_covariable_auto <- function(modelo, factor_var = "species", tipo = "response") {
  mf <- model.frame(modelo)
  preds <- names(mf)[-1]  # excluye la respuesta
  # Identificar la covariable numérica (excluye el factor y agrupamientos)
  candidatos <- setdiff(preds, factor_var)
  covar_var  <- candidatos[vapply(mf[candidatos], is.numeric, logical(1))]
  if (length(covar_var) == 0) stop("No se encontró covariable numérica en el modelo")
  covar_var <- covar_var[1]

  pct <- quantile(mf[[covar_var]], probs = c(0.25, 0.50, 0.75))
  fml <- as.formula(paste("~", factor_var, "|", covar_var))
  emm <- emmeans(modelo, fml, at = setNames(list(pct), covar_var), type = tipo)
  cld(emm, Letters = letters)
}

procesar_modelo <- function(modelo, meta_row, factor_var = "species") {
  familia <- family(modelo)
  tipo <- if (familia$family %in% c("Gamma", "binomial")) "response" else "link"
  
  hay_covariable  <- tiene_covariable(modelo, factor_var)
  hay_interaccion <- hay_covariable && tiene_interaccion(modelo, factor_var)
  
  res <- if (hay_interaccion) {
    # interacción significativa -> reportar por percentiles de la covariable
    emmeans_covariable_auto(modelo, factor_var = factor_var, tipo = tipo)
  } else {
    # sin interacción (aunque haya covariable aditiva) -> media marginal simple,
    # la covariable queda fijada en su media internamente por emmeans
    cld(emmeans(modelo, as.formula(paste("~", factor_var)), type = tipo), Letters = letters)
  }
  
  attr(res, "method")   <- extraer_metodo(res)
  attr(res, "response") <- meta_row$response
  attr(res, "units")    <- meta_row$units
  attr(res, "type")     <- meta_row$type
  attr(res, "covariable_ajustada") <- hay_covariable && !hay_interaccion  # info extra
  res
}

# Aplica a todos los modelos de una sola vez, sin duplicar bucles
emm.results <- Map(procesar_modelo, modelos, split(metadata, seq_len(nrow(metadata))))
names(emm.results) <- metadata$modelo

## 2. Creat mother table ----

library(dplyr)
library(purrr)

tidy_emm <- function(res, modelo_nombre) {
  dat <- as.data.frame(res)
  
  est_name <- attr(res, "estName")
  dat <- dat %>% dplyr::rename(estimate = dplyr::all_of(est_name))
  
  by_var <- attr(res, "by.vars")
  
  if (!is.null(by_var)) {
    valores_unicos <- sort(unique(dat[[by_var]]))
    n <- length(valores_unicos)
    etiquetas <- if (n == 2) c("Q25", "Q75") 
    else if (n == 3) c("Q25", "Q50", "Q75")
    else paste0("nivel_", seq_len(n))
    
    dat$nivel_covariable  <- etiquetas[match(dat[[by_var]], valores_unicos)]
    dat$covariable        <- by_var
    dat$valor_covariable  <- dat[[by_var]]
  } else {
    dat$nivel_covariable <- "mean"
    dat$covariable       <- NA_character_
    dat$valor_covariable <- NA_real_
  }
  
  dat %>%
    dplyr::mutate(
      modelo              = modelo_nombre,
      tipo_modelo         = as.character(attr(res, "type")),
      variable_respuesta  = attr(res, "response"),
      unidades            = attr(res, "units"),
      metodo              = attr(res, "method")
    ) %>%
    dplyr::select(modelo, tipo_modelo, variable_respuesta, unidades, species,
                  covariable, valor_covariable, nivel_covariable,
                  estimate, SE, df, asymp.LCL, asymp.UCL, .group, metodo)
}

tabla_final <- purrr::imap_dfr(emm.results, tidy_emm)

## 3. Creat publication table ----

library(dplyr)
library(tidyr)

### Change SPM units
tabla_final <- 
tabla_final |> 
  mutate(estimate = if_else(condition = variable_respuesta == "Specific pericarp mass", 
                             true = estimate*1000, false = estimate), 
         SE = if_else(condition = variable_respuesta == "Specific pericarp mass", 
                      true = SE * 1000, false = SE),
         asymp.LCL = if_else(condition = variable_respuesta == "Specific pericarp mass", 
                      true = asymp.LCL * 1000, false = asymp.LCL),
         asymp.UCL = if_else(condition = variable_respuesta == "Specific pericarp mass", 
                             true = asymp.UCL * 1000, false = asymp.UCL), 
         unidades = if_else(condition = variable_respuesta == "Specific pericarp mass", 
                            true = "mg/cm²", false = unidades))

tabla_publicacion <- tabla_final %>%
  mutate(
    # etiqueta de columna: variable (unidad) [, Q25/Q75 si aplica]
    columna = case_when(
      nivel_covariable == "mean" ~ paste0(variable_respuesta, " (", unidades, ")"),
      TRUE ~ paste0(variable_respuesta, " (", unidades, ", ", nivel_covariable, ")")
    ),
    # valor formateado: media ± SE letra
    valor_formateado = paste0(
      formatC(estimate, format = "f", digits = 2), " ± ",
      formatC(SE, format = "f", digits = 2), " ",
      trimws(.group)
    )
  ) %>%
  dplyr::select(species, columna, valor_formateado) %>%
  pivot_wider(names_from = columna, values_from = valor_formateado)

tabla_publicacion

## Save ----

write.csv2(x = tabla_final, file = "00-data/emm_traits_long.csv")
write.csv2(x = tabla_publicacion, file = "00-data/paper_traits.csv")

# ===========================================================================
# TABLAS RESUMEN: table_traitmodel_summary y table_traitmodel_effects
# ===========================================================================
# Se generan a partir de los modelos ya ajustados (lista `modelos`) y del
# `metadata` definidos arriba. No se reajustan ni modifican los modelos.
#
# Paquetes necesarios (adición a los ya cargados):
#   - parameters: model_parameters() para extraer coeficientes fijos
#   - car:        Anova() para tests globales de Wald (chi² tipo II)

library(parameters)
library(car)

# --- Funciones auxiliares reutilizables ---

# Cadena legible con las varianzas de los efectos aleatorios
.format_re <- function(model) {
  re <- tryCatch(VarCorr(model)$cond, error = function(e) NULL)
  if (is.null(re) || length(re) == 0) return(NA_character_)
  parts <- vapply(names(re), function(g) {
    mat <- re[[g]]
    paste0(g, ": σ² = ", format(mat[1, 1], digits = 4))
  }, character(1))
  paste(parts, collapse = "; ")
}

# Número de grupos de la variable de agrupación (species:provenance)
.count_groups <- function(model) {
  mod_data <- tryCatch({
    if (!is.null(model$frame)) model$frame else model.frame(model)
  }, error = function(e) NULL)
  if (is.null(mod_data)) return(NA_integer_)
  if (all(c("species", "provenance") %in% names(mod_data))) {
    return(length(unique(interaction(mod_data$species, mod_data$provenance))))
  }
  NA_integer_
}

# Codificación de significancia
.sig_code <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Niveles de referencia de factores en el modelo
.get_refs <- function(model) {
  mf <- model.frame(model)
  preds <- names(mf)[-1]  # excluye la respuesta
  refs <- character(0)
  for (nm in preds) {
    if (is.factor(mf[[nm]])) {
      refs[nm] <- levels(mf[[nm]])[1]
    }
  }
  refs
}

# Tipo de estadístico del resumen del modelo (t para gaussian, z para otros)
.get_stat_type <- function(model) {
  s <- tryCatch(summary(model), error = function(e) NULL)
  if (is.null(s) || is.null(s$coefficients$cond)) return("z")
  if ("t value" %in% colnames(s$coefficients$cond)) "t" else "z"
}

# --- 1. table_traitmodel_summary ---
# Una fila por modelo con métricas de ajuste y estructura.
# Columnas:
#   model         → nombre del objeto R
#   response      → variable respuesta (nombre legible)
#   model_type    → clase del modelo (todos GLMM aquí)
#   family        → familia y enlace del modelo
#   n             → observaciones utilizadas
#   n_groups      → unidades de agrupación (species:provenance)
#   r2_marginal   → R² marginal (performance::r2)
#   r2_conditional→ R² condicional (performance::r2)
#   AIC / BIC     → criterios de información
#   random_effects→ estructura de efectos aleatorios formateada
#   singular      → TRUE/FALSE si el modelo es singular

families_text <- sapply(modelos, function(m) {
  f <- family(m)
  paste0(f$family, "(", f$link, ")")
})

table_traitmodel_summary <- purrr::imap_dfr(modelos, function(mod, nm) {
  idx <- match(nm, metadata$modelo)

  n_obs <- nobs(mod)
  n_grp <- .count_groups(mod)

  # R² marginal y condicional
  # No aplica directamente para GLMM no gaussianos; r2() lo intenta y
  # devuelve NA si no es posible de forma fiable.
  r2_vals <- tryCatch(performance::r2(mod), error = function(e) {
    warning("[table_traitmodel_summary] r2() falló para ", nm, ": ",
            conditionMessage(e))
    list(R2_marginal = NA_real_, R2_conditional = NA_real_)
  })

  aic_val <- tryCatch(AIC(mod), error = function(e) NA_real_)
  bic_val <- tryCatch(BIC(mod), error = function(e) NA_real_)

  re_str   <- .format_re(mod)
  singular <- tryCatch(performance::check_singularity(mod), error = function(e) NA)

  tibble(
    model            = nm,
    response         = metadata$response[idx],
    model_type       = "GLMM",
    family           = families_text[idx],
    n                = n_obs,
    n_groups         = n_grp,
    r2_marginal      = r2_vals$R2_marginal,
    r2_conditional   = r2_vals$R2_conditional,
    AIC              = aic_val,
    BIC              = bic_val,
    random_effects   = re_str,
    singular         = singular
  )
})

# --- 2. table_traitmodel_effects ---
# Una fila por coeficiente individual (contrastes respecto al nivel de
# referencia) y una fila por test global de cada término del modelo.
#
# Coeficientes individuales (term_type = "coefficient"):
#   Extraídos con parameters::model_parameters(effects = "fixed", ci = 0.95).
#   estimate = β, SE = error estándar, CI_low/CI_high = IC 95%,
#   statistic = t o z según el modelo, p = p-value.
#
# Tests globales (term_type = "global_test"):
#   Extraídos con car::Anova(type = "II", test.statistic = "Chisq").
#   statistic = χ² de Wald, df = grados de libertad, p = p-value.
#
# Inferencia: se conserva el método de cada modelo.
#   - Coeficientes individuales: Wald (t para gaussian, z para otros)
#   - Tests globales: Wald χ² tipo II (car::Anova)

table_traitmodel_effects <- purrr::imap_dfr(modelos, function(mod, nm) {
  idx <- match(nm, metadata$modelo)
  resp <- metadata$response[idx]

  # --- 2a. Coeficientes individuales de efectos fijos ---
  coef_df <- tryCatch({
    params <- parameters::model_parameters(mod, effects = "fixed", ci = 0.95)
    pf <- as.data.frame(params)
    refs <- .get_refs(mod)
    stype <- .get_stat_type(mod)

    # Estadístico: puede llamarse "Statistic", "t" o "z" según versión
    stat_col <- intersect(c("Statistic", "t", "z"), names(pf))
    stat_val <- if (length(stat_col) > 0) pf[[stat_col[1]]] else rep(NA_real_, nrow(pf))

    # Grados de libertad: solo significativos para t (Inf → NA)
    df_col <- intersect(c("df_error", "df"), names(pf))
    df_val <- if (length(df_col) > 0) {
      d <- pf[[df_col[1]]]
      ifelse(is.infinite(d), NA_real_, d)
    } else {
      rep(NA_real_, nrow(pf))
    }

    # Nivel de referencia solo para terms categóricos principales
    ref_vals <- vapply(pf$Parameter, function(par) {
      if (grepl("^species", par) && !grepl(":", par)) {
        paste0("ref: ", refs["species"])
      } else {
        NA_character_
      }
    }, character(1))

    tibble(
      model          = nm,
      response       = resp,
      term           = pf$Parameter,
      term_type      = "coefficient",
      estimate       = pf$Coefficient,
      SE             = pf$SE,
      CI_low         = pf$CI_low,
      CI_high        = pf$CI_high,
      statistic      = stat_val,
      statistic_type = stype,
      df             = df_val,
      p              = pf$p,
      significance   = .sig_code(pf$p),
      reference      = ref_vals
    )
  }, error = function(e) {
    warning("[table_traitmodel_effects] Extracción de coeficientes falló para ",
            nm, ": ", conditionMessage(e))
    tibble(model = nm, response = resp, term = NA_character_,
           term_type = "coefficient",
           estimate = NA_real_, SE = NA_real_, CI_low = NA_real_,
           CI_high = NA_real_, statistic = NA_real_,
           statistic_type = NA_character_, df = NA_real_,
           p = NA_real_, significance = NA_character_,
           reference = NA_character_)
  })

  # --- 2b. Tests globales por término (Wald χ² tipo II) ---
  anova_df <- tryCatch({
    anova_res <- car::Anova(mod, type = "II", test.statistic = "Chisq")
    af <- as.data.frame(anova_res)

    # Identificar columnas (pueden variar según versión de car)
    chisq_col <- grep("Chisq|chisq|LR", names(af), value = TRUE)[1]
    p_col     <- grep("Pr\\(|p\\.value", names(af), value = TRUE)[1]
    df_col    <- grep("^Df$|^df$", names(af), value = TRUE)[1]

    if (is.na(chisq_col) || is.na(p_col)) {
      stop("No se encontraron columnas de Chisq/p en la salida de car::Anova()")
    }

    tibble(
      model          = nm,
      response       = resp,
      term           = rownames(af),
      term_type      = "global_test",
      estimate       = NA_real_,
      SE             = NA_real_,
      CI_low         = NA_real_,
      CI_high        = NA_real_,
      statistic      = af[[chisq_col]],
      statistic_type = "Chisq",
      df             = if (!is.na(df_col)) af[[df_col]] else NA_real_,
      p              = af[[p_col]],
      significance   = .sig_code(af[[p_col]]),
      reference      = NA_character_
    )
  }, error = function(e) {
    warning("[table_traitmodel_effects] car::Anova falló para ", nm, ": ",
            conditionMessage(e))
    tibble(model = nm, response = resp, term = NA_character_,
           term_type = "global_test",
           estimate = NA_real_, SE = NA_real_, CI_low = NA_real_,
           CI_high = NA_real_, statistic = NA_real_,
           statistic_type = NA_character_, df = NA_real_,
           p = NA_real_, significance = NA_character_,
           reference = NA_character_)
  })

  dplyr::bind_rows(coef_df, anova_df)
})

# --- Imprimir tablas ---
cat("\n==== table_traitmodel_summary ====\n")
print(table_traitmodel_summary)
cat("\n==== table_traitmodel_effects ====\n")
print(table_traitmodel_effects)

write.csv(table_traitmodel_summary, file = "00-data/table_traitmodel_summary.csv")
write.csv(table_traitmodel_effects, file = "00-data/table_traitmodel_effects.csv")

# --- 3. paper_table_traitmodel_effects ---
# Versión condensada de table_traitmodel_effects para el artículo.
# Una fila por término del modelo (no por coeficiente individual).
# Columnas:
#   response  → variable respuesta (nombre legible)
#   formula   → fórmula del modelo en notación compacta
#   term      → nombre del término (species, cov, species:cov)
#   test      → resultado formateado: Chisq(df) = valor***

paper_table_traitmodel_effects <- purrr::imap_dfr(modelos, function(mod, nm) {
  idx <- match(nm, metadata$modelo)
  resp <- metadata$response[idx]

  anova_res <- car::Anova(mod, type = "II", test.statistic = "Chisq")
  af <- as.data.frame(anova_res)
  anova_terms <- rownames(af)

  # Fórmula compacta: × para interacciones, + para términos aditivos
  response_var <- deparse(formula(mod)[[2]])
  formula_terms <- gsub(":", " × ", anova_terms)
  formula_str <- paste(response_var, "~", paste(formula_terms, collapse = " + "))

  # Formatear resultado: Chisq(df) = valor***
  test_strs <- vapply(seq_len(nrow(af)), function(i) {
    chisq_val <- af[["Chisq"]][i]
    df_val    <- af[["Df"]][i]
    p_val     <- af[["Pr(>Chisq)"]][i]
    stars     <- .sig_code(p_val)
    paste0("Chisq(", df_val, ") = ", format(chisq_val, digits = 4), stars)
  }, character(1))

  tibble(
    response = resp,
    formula  = formula_str,
    term     = anova_terms,
    test     = test_strs
  )
})

cat("\n==== paper_table_traitmodel_effects ====\n")
print(paper_table_traitmodel_effects)
write.csv(paper_table_traitmodel_effects, file = "00-data/paper_table_traitmodel_effects.csv", row.names = FALSE)
