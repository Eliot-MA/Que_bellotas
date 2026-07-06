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
easystats::model_dashboard(glmm.dw.gamm,
                           output_dir = "06-html/",
                           output_file = "modeldashboard_dryweight.html")

## 2. Specific pericarp mass ----
glmm.spm <- glmmTMB(SPM_g_cm2 ~ species + (1|species:provenance),
                    data = df)

easystats::model_dashboard(glmm.spm,
                           output_dir = "06-html/",
                           output_file = "modeldashboard_spm.html")

## 3. Scar area ----
glmm.scar <- glmmTMB(scar_area_mm2 ~ species*surface_cm2 + (1|species:provenance),
                     data = df, 
                     family = Gamma(link = "log"))
options(na.action = "na.fail")
dd <- dredge(glmm.scar)
glmm.scar1 <- get.models(dd, subset = 1)[[1]]

easystats::model_dashboard(glmm.scar1,
                           output_dir = "06-html/",
                           output_file = "modeldashboard_scar.html")

## 4. Seed coat ratio ----
glmm.scr <- glmmTMB(pericarp_dry_weight ~ dry_weight*species + 
                       (1|species:provenance),
                     data = df)

dd <- dredge(glmm.scr)
glmm.scr1 <- get.models(dd, subset = 1)[[1]]

easystats::model_dashboard(glmm.scr1,
                           output_dir = "06-html/",
                           output_file = "modeldashboard_scr.html")

## 5. Pericarp rupture ----
df.cracks <- df |> 
  mutate(pericarp_rupture = if_else(pericarp_rupture == 2, 1, pericarp_rupture))

glmm.cracks <- glmmTMB(pericarp_rupture ~ species + (1|species:provenance), 
                      data = df.cracks , 
                      family = binomial(link = "logit"))

easystats::model_dashboard(glmm.cracks,
                           output_dir = "06-html/",
                           output_file = "modeldashboard_cracks.html")

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
