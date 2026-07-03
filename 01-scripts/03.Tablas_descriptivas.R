df.bellotas <- read.csv(file = "data/df.bellotas.csv")

library(tidyverse) 
library(emmeans)     # Calculo de medias marginales y tendencias
library(performance) # Comprobación de supuestos
library(multcomp)    # Asignación de letras de significancia estadística
library(nlme)        # Modelos gls (generalized least squares)
library(glmmTMB)
library(e1071)   # Calcular kurtosis
 
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

# EDA
## Peso seco

df_unique |> 
  group_by(especie) |> 
  summarise(
    n = n(), 
    minimo = min(peso_seco),
    Q25 = quantile(peso_seco, 0.25),
    mediana = median(peso_seco),
    media = mean(peso_seco), 
    DE = sd(peso_seco), 
    Q75 = quantile(peso_seco, 0.75),
    maximo = max(peso_seco),
    IQR = IQR(peso_seco),
    kurtosis = kurtosis(peso_seco)
  )

df_unique |> 
  group_by(especie, procedencia) |> 
  summarise(
    n = n(), 
    minimo = min(peso_seco),
    Q25 = quantile(peso_seco, 0.25),
    mediana = median(peso_seco),
    media = mean(peso_seco), 
    DE = sd(peso_seco), 
    Q75 = quantile(peso_seco, 0.75),
    maximo = max(peso_seco),
    IQR = IQR(peso_seco),
    kurtosis = kurtosis(peso_seco)
  )

# Efecto fijo (especie): diferencia notable entre medias. 
# Efecto aleatorio (procedencias en especies): Si miras la segunda tabla, verás que dentro de una misma especie, la procedencia cambia drásticamente el peso. Ejemplo extremo: En Quercus pyrenaica, la procedencia de "Asturias-Lena" tiene una media de 2.55, pero la de "Salamanca-Sayago" es de 4.92.
# Probable heterocedasticidad. La Desviación Estándar (DE) varía mucho: desde 0.55 en Q. pubescens hasta 1.62 en Q. pyrenaica. Tal vez sea necesario modelar la varianza. 
# Distribución platitúrtica (curtosis negativa): las distribuciones son algo más "planas" que una normal y tienen colas menos pesadas. Sin embargo, no son valores extremos (están cerca de 0), por lo que la asunción de normalidad del modelo debería funcionar bien.
# Signos de posible distribución gamma: conforme la media aumenta (especialmente al pasar de ~1.8 a ~3.8), la Desviación Estándar (DE) se duplica o triplica. Mediana siempre mayor que media. 

hist(df_unique$peso_seco) # Se parece a una gamma
ggplot(df_unique, aes(x = peso_seco)) +
  geom_histogram() +
  facet_wrap(~especie)

## SPM

df_unique |> 
  group_by(especie) |> 
  summarise(
    n = n(), 
    minimo = min(SPM_g_cm2),
    Q25 = quantile(SPM_g_cm2, 0.25),
    mediana = median(SPM_g_cm2),
    media = mean(SPM_g_cm2), 
    DE = sd(SPM_g_cm2), 
    Q75 = quantile(SPM_g_cm2, 0.75),
    maximo = max(SPM_g_cm2),
    IQR = IQR(SPM_g_cm2),
    kurtosis = kurtosis(SPM_g_cm2)
  )

df_unique |> 
  group_by(especie, procedencia) |> 
  summarise(
    n = n(), 
    minimo = min(SPM_g_cm2),
    Q25 = quantile(SPM_g_cm2, 0.25),
    mediana = median(SPM_g_cm2),
    media = mean(SPM_g_cm2), 
    DE = sd(SPM_g_cm2), 
    Q75 = quantile(SPM_g_cm2, 0.75),
    maximo = max(SPM_g_cm2),
    IQR = IQR(SPM_g_cm2),
    kurtosis = kurtosis(SPM_g_cm2)
  )

hist(df_unique$SPM_g_cm2)

## Medias y medianas muy similares -> distribución gausiana
## Kurtosis muy elevada en Q.ilex -> posible outlier
## Pasos recomendados: 
## 1. ajustar modelo gausiano
## 2. mirar residuos y homogeneidad
## 3. si existe heterocedasticidad -> probar log(spm)
## 4. si no funciona -> probar gamma

## Area de la cicatriz


## Peso seco del pericarpo

## Rajado




# Modelos sencillos (ANOVA)
library(DHARMa)
library(performance)

## Modelo peso seco bellotas
glmm.pesoseco.gaus <- glmmTMB(peso_seco ~ especie + (1|especie:procedencia),
                              data = df_unique)
glmm.pesoseco.gamm <- glmmTMB(peso_seco ~ especie + (1|especie:procedencia),
                              data = df_unique, 
                              family = Gamma(link = "log"))

compare_performance(glmm.pesoseco.gaus, glmm.pesoseco.gamm)
# gamma sale elejido por AIC

res_peso <- simulateResiduals(glmm.pesoseco.gamm, n = 1000)
plot(res_peso)

plotResiduals(res_peso, df_unique$peso_seco)   #Quantile deviations detected; Combined adjusted quantile test significant 
plotResiduals(res_peso, df_unique$procedencia) #Within-group deviations from uniformity significant; Levene Test form homogeneity of variance significant

testDispersion(res_peso)

glmm.pesoseco.gamm.hom <- glmmTMB(peso_seco ~ especie + (1|especie:procedencia),
                              data = df_unique, 
                              family = Gamma(link = "log"))
glmm.pesoseco.gamm.het <- glmmTMB(peso_seco ~ especie + (1|especie:procedencia),
                                  data = df_unique, 
                                  family = Gamma(link = "log"), 
                                  dispformula  = ~especie)

compare_performance(glmm.pesoseco.gamm.hom, glmm.pesoseco.gamm.het)
# Ambos modelos son equivalentes en AIC, elegimos homocedástico

##
# Bootstrapping paramétrico para determinar robustez de la estimación
##

boot_fun <- function(fit) {
  fixef(fit)$cond
}

bootstrap_glmmTMB <- function(model, nsim = 1000) {
  
  coefs <- matrix(NA, nrow = nsim, ncol = length(fixef(model)$cond))
  colnames(coefs) <- names(fixef(model)$cond)
  
  for (i in 1:nsim) {
    
    # Simular nueva respuesta
    sim_data <- simulate(model)[[1]]
    
    # Reemplazar respuesta en el dataset original
    new_data <- model$frame
    response_name <- all.vars(formula(model))[1]
    new_data[[response_name]] <- sim_data
    
    # Reajustar modelo
    fit_sim <- try(
      update(model, data = new_data),
      silent = TRUE
    )
    
    # Guardar coeficientes si converge
    if (!inherits(fit_sim, "try-error")) {
      coefs[i, ] <- fixef(fit_sim)$cond
    }
  }
  
  return(coefs)
}

boot_res <- bootstrap_glmmTMB(glmm.pesoseco.gamm, nsim = 1000)

apply(boot_res, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary(mod_peso_glm)

# Bootstrap confirma los resultados de la aproximación clásica
# Los coeficientes que aparecen significativos con la aprox clasica
# tambien lo son en la validación con bootrap 

## Modelos peso específico del pericarpo
glmm.spm <- glmmTMB(SPM_g_cm2 ~ especie + (1|especie:procedencia), 
                   data = df_unique)

summary(glmm.spm)

res_spm <- simulateResiduals(glmm.spm, n = 5000)
plot(res_spm) # Levene Test positivo

glmm.logspm <- glmmTMB(log(SPM_g_cm2) ~ especie + (1|especie:procedencia), 
                    data = df_unique)

res_spm <- simulateResiduals(glmm.logspm, n = 5000)
plot(res_spm) # Levene Test positivo

glmm.sqrspm <- glmmTMB(sqrt(SPM_g_cm2) ~ especie + (1|especie:procedencia), 
                       data = df_unique)

res_spm <- simulateResiduals(glmm.sqrspm, n = 5000)
plot(res_spm) # Levene Test positivo

glmm.invspm <- glmmTMB(1/SPM_g_cm2 ~ especie + (1|especie:procedencia), 
                       data = df_unique)

res_spm <- simulateResiduals(glmm.invspm, n = 5000)
plot(res_spm) # Levene Test positivo

res_sqr <- simulateResiduals(glmm.sqrspm, n = 5000)
res_log <- simulateResiduals(glmm.logspm, n = 5000)

plot(res_sqr, main = "sqrt(SPM)") 
plot(res_log, main = "log(SPM)")  # log produce outliers

# Nada soluciona el problema, pasamos a modelos gamma

glmm.spm.gamm <- glmmTMB(SPM_g_cm2 ~ especie + (1|especie:procedencia), 
                    data = df_unique, 
                    family = Gamma(link = "log"))

res_spm <- simulateResiduals(glmm.spm.gamm, n = 5000)
plot(res_spm) # Levene Test positivo
testDispersion(res_spm)

compare_performance(glmm.spm, glmm.spm.gamm)

boot_res <- bootstrap_glmmTMB(glmm.spm.gamm, nsim = 1000)

apply(boot_res, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary(glmm.spm.gamm)


# Modelos con covariables para estandarizar (ANCOVA)
## Área de la cicatriz
mod_area <- glmmTMB(Area_cicatriz_mm2 ~ especie * Area_estimada_cm2 + (1|especie:procedencia),

               data = df_unique)
car::Anova(mod_area)

res_area <- simulateResiduals(mod_area, n = 10000)
plot(res_area) # KS positivo; Outlier test positivo; Quantile deviation detected; Combined adjusted quantile test significant

plotResiduals(res_area, df_unique$procedencia) #Within-group deviations from uniformity significant; Levene Test for homogeneity of variance significant

mod_area <- glmmTMB(Area_cicatriz_mm2 ~ especie * Area_estimada_cm2 + (1|especie:procedencia),
                    data = df_unique, 
                    family = Gamma(link = "log"))

car::Anova(mod_area)

res_area <- simulateResiduals(mod_area, n = 10000)
plot(res_area) # KS significativo; Quantile deviation detected; Combined adjusted quantile test significant

plotResiduals(res_area, df_unique$procedencia) #Within-group deviations from uniformity significant

##
# Bootstrapping paramétrico para determinar robustez de la estimación
##

boot_res <- bootstrap_glmmTMB(mod_area, nsim = 1000)

apply(boot_res, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary(mod_area)

# Bootstrap confirma los resultados de la aproximación clásica
# Los coeficientes que aparecen significativos con la aprox clasica
# tambien lo son en la validación con bootrap 

## Peso del pericarpo
mod_peri <- glmmTMB(peso_seco_pr ~ especie * peso_seco + 
                      (1|especie:procedencia), 
                    data = df_unique)
car::Anova(mod_peri)

res_peri <- simulateResiduals(mod_peri, n = 10000)
plot(res_peri) # KS positivo; Quantile deviation detected; Combined adjusted quantile test significant
plotResiduals(res_peri, df_unique$procedencia) #Within-group deviations from uniformity significant; Levene Test for homogeneity of variance significant

mod_peri <- glmmTMB(peso_seco_pr ~ especie * peso_seco + 
                      (1|especie:procedencia), 
                    data = df_unique, 
                    family = Gamma(link = "log"))

res_peri <- simulateResiduals(mod_peri, n = 10000)
plot(res_peri) # Quantile deviation detected; Combined adjusted quantile test significant
plotResiduals(res_peri, df_unique$procedencia) #Within-group deviations from uniformity significant; Levene Test for homogeneity of variance significant

##
# Bootstrapping paramétrico para determinar robustez de la estimación
##

boot_res <- bootstrap_glmmTMB(mod_peri, nsim = 1000)

apply(boot_res, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary(mod_peri)

# Bootstrap confirma los resultados de la aproximación clásica
# Los coeficientes que aparecen significativos con la aprox clasica
# tambien lo son en la validación con bootrap 

# Debido a la interaccion significativa en los modelos ancova, se tomaron las siguientes decisiones para realizar un reporte representativo de las medias marginales. 
library(sjPlot)
library(ggplot2)

# MOD SCR
plot_model(mod_peri, 
           type = "pred", 
           terms = c("peso_seco", "especie"),
           title = "Predicción de peso del pericarpo por Especie",
           axis.title = c("Peso seco (g)", "Peso del pericarpo (g)")) +
  theme_minimal()

# Regiones de prediccion sin sobre/infra prediccion
cuantiles <- quantile(x = df_unique$peso_seco, probs = seq(0.1, 0.9, 0.1))
ggplot(df_unique, aes(x = peso_seco, y = especie)) + geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  geom_vline(xintercept = cuantiles[c(3, 7)], color = "red")

cuantiles

cld(
emmeans(mod_peri, ~ especie|peso_seco, 
        at = list(peso_seco = cuantiles[3:7]))
)

# CALCULAR MEDIAS MARGINALES PARA LOS MODELOS ANOVA
library(emmeans)
library(multcomp)
library(dplyr)
library(purrr)

emm_cld_table <- function(model, trait, adjust = "tukey", digits = 2) {
  emm <- emmeans(model, ~ especie)
  
  s <- as.data.frame(summary(emm, type = "response"))
  est_name <- intersect(names(s), c("response", "emmean", "estimate"))[1]
  if (is.na(est_name)) {
    est_name <- setdiff(names(s), c("especie", "SE", "df", "lower.CL", "upper.CL",
                                    "asymp.LCL", "asymp.UCL", "t.ratio", "p.value"))[1]
  }
  
  means <- s %>%
    transmute(especie, media = .data[[est_name]], SE = SE)
  
  groups <- as.data.frame(multcomp::cld(emm, adjust = adjust, Letters = letters)) %>%
    transmute(especie, letras = gsub("\\s+", "", .group))
  
  left_join(means, groups, by = "especie") %>%
    mutate(
      `media ± SE letras` = sprintf(paste0("%.", digits, "f ± %.", digits, "f %s"),
                                    media, SE, letras),
      trait = trait
    ) %>%
    dplyr::select(trait, especie, `media ± SE letras`)
}

tabla_final <- bind_rows(
  emm_cld_table(mod_peso_glm, "Peso seco (g)"), 
  emm_cld_table(mod_area, "Area de cicatriz (mm2)")
)

tabla_final

# CALCULAR MEDIAS MARGINALES PARA LAS ANCOVA
## MOD_SCR
library(emmeans)
library(multcomp)
library(dplyr)

# Definimos los puntos de los cuantiles
puntos_evaluacion <- cuantiles[c(3, 7)]

# 1. Calculamos emmeans para ambos puntos a la vez
emm_peri <- emmeans(mod_peri, ~ especie | peso_seco, 
                    at = list(peso_seco = puntos_evaluacion), 
                    type = "response")

# 2. Obtenemos las letras de significancia (Tukey) para cada punto
# Ajustamos por punto de evaluación
letras_peri <- cld(emm_peri, Letters = letters, adjust = "tukey")

# 3. Formateamos los datos para que coincidan con tu tabla
limite <- as.numeric(cuantiles[3])


tabla_ancova_peri <- as.data.frame(letras_peri) %>%
  mutate(
    trait = "Peso seco pericarpo (g)",
    media_se = paste0(round(emmean, 2), " ± ", round(SE, 2)),
    letras = trimws(.group), # Limpiar espacios de las letras
    contexto = paste0("Peso seco (g) q", ifelse(peso_seco <= limite, "30%", "70%"))
  ) %>%
  dplyr::select(trait, especie, media_se, letras, contexto)

tabla_ancova_peri

# UNIR TODAS LAS TABLAS
library(dplyr)

# 1. Preparar tabla_ancova
# Vamos a crear la columna combinada para que sea igual a la de tabla_final
tabla_ancova_prep <- tabla_ancova_peri %>%
  mutate(`media ± SE letras` = paste(media_se, letras)) %>%
  dplyr::select(trait, especie, `media ± SE letras`, contexto)

# 2. Preparar tabla_final
# Vamos a añadir la columna 'contexto' que le falta (poniendo que es la media global)
tabla_final_prep <- tabla_final %>%
  mutate(contexto = "Promedio global") %>%
  # Aseguramos el orden de las columnas para que coincida
  dplyr::select(trait, especie, `media ± SE letras`, contexto)

# 3. Unirlas verticalmente
tabla_unida <- bind_rows(tabla_final_prep, tabla_ancova_prep)

# 4. (Opcional) Ver el resultado
print(tabla_unida)

# Crear tabla con formato cientifico
# 1. Crear el encabezado final combinando rasgo y contexto
tabla_preparada <- tabla_unida %>%
  mutate(columna_header = paste0(trait, "\n(", contexto, ")")) %>%
  dplyr::select(especie, columna_header, `media ± SE letras`)

# 2. Pivotar la tabla: Las especies se quedan en filas, 
# los rasgos + contexto pasan a ser columnas.
# Aplicamos el orden y pivotamos
orden_especies <- c("Quercus coccifera", "Quercus ilex", "Quercus suber", "Quercus faginea", "Quercus pubescens", "Quercus pyrenaica", "Quercus petraea", "Quercus robur")
tabla_publicacion <- tabla_preparada %>%
  # Convertimos especie en factor con el orden que definimos
  mutate(especie = factor(especie, levels = orden_especies)) %>%
  # Ordenamos las filas según esos niveles
  arrange(especie) %>%
  # Pasamos a formato ancho (el paso que hicimos antes)
  pivot_wider(
    names_from = columna_header, 
    values_from = `media ± SE letras`
  )

# 3. Ver el resultado
print(tabla_publicacion)

##
# Modelo de probabilidad de rajado por especie
##

df.rajas <- df_unique |> filter(rajas_pericarpo != 2) 

library(glmmTMB)

glmm.rajas <- glmmTMB(rajas_pericarpo ~ especie + (1|especie:procedencia), 
                      data = df.rajas , 
        family = binomial(link = "logit"))

VarCorr(glmm.rajas)

check_model(glmm.rajas)

library(DHARMa)
res <- simulateResiduals(glmm.rajas, n = 1000)
plot(res)

testDispersion(res)
testOutliers(res)
table(df.rajas$especie, df.rajas$rajas_pericarpo)

##
# Bootstrapping paramétrico para determinar robustez de la estimación
##

boot_res <- bootstrap_glmmTMB(glmm.rajas, nsim = 1000)

apply(boot_res, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary(glmm.rajas)

# Bootstrap confirma los resultados de la aproximación clásica
# Los coeficientes que aparecen significativos con la aprox clasica
# tambien lo son en la validación con bootrap 



# Estimar probabilidades por epsecie
library(emmeans)

emm <- emmeans(glmm.rajas, ~ especie, type = "response")
emm
pairs(emm, adjust = "tukey")

library(multcomp)

cld_res <- cld(emm, adjust = "tukey", Letters = letters)
cld_res
ggplot(cld_res, aes(x = reorder(especie, prob), y = prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  geom_text(aes(label = .group, y = prob + 0.05)) +
  coord_flip()

df.rajas |> 
  group_by(especie, procedencia) |> 
  summarise(mean_rajas = mean(rajas_pericarpo), 
            sd_rajas = sd(rajas_pericarpo), 
            n = n())

library(stringr)

cld_res2 <- cld_res %>%
  mutate(
    especie = str_trim(as.character(especie)),
    grupo = str_trim(.group)
  ) %>%
  dplyr::select(especie, prob, SE, asymp.LCL, asymp.UCL, grupo)


tabla_publicacion <- tabla_publicacion %>%
  left_join(
    cld_res2 %>% 
      dplyr::select(especie, prob, SE, grupo),
    by = "especie"
  )


tabla_publicacion <- tabla_publicacion %>%
  mutate(
    prob_fmt = paste0(
      round(prob, 2), " ± ",
      round(SE, 2), " ",
      grupo
    )
  ) |>
  dplyr::select(-prob, -SE, -grupo)



# Guardar resultados
write.csv2(tabla_unida, file = "data/tabla_descriptiva_longer.csv")
write.csv2(tabla_publicacion, file = "data/tabla_descriptiva_publicacion.csv")

write.table(
  tabla_publicacion,
  file = "data/tabla_publicacion.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
