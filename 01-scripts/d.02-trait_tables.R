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
  )

df_unique <- df_unique |> 
  rename(species = especie, 
         provenance = procedencia, 
         prov_code = codigo, 
         length_mm = longitud..mm., 
         diameter_mm = diametro_total, 
         dry_weight = peso_seco, 
         scar_area_mm2 = Area_cicatriz_mm2, 
         pericarp_dry_weight = peso_seco_pr, 
         surface_cm2 = Area_estimada_cm2, 
         pericarp_rupture = rajas_pericarpo)

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
modelos <- list(glmm.dw.gamm, glmm.spm, glmm.scar1, glmm.scr1, glmm.cracks)

modelo <- modelos[[2]]
familia <- family(modelo)

if (familia$family == "Gamma") {
  emm <- cld(emmeans(modelo, ~species), Letters = letters, type = "response")
  s <- as.data.frame(emm)
  columna <- paste0(round(s$response, 2), "±", round(s$SE, 2), s$.group)
} else {
  emm <- cld(emmeans(modelo, ~species), Letters = letters)
  s <- as.data.frame(emm)
  columna <- paste0(round(s$emmean, 2), "±", round(s$SE, 2), s$.group)
}
