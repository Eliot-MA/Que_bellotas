###############################################################################
## Script para el análisis del MC a la que la germinación tiene una p de 0.5 ##
###############################################################################


library(MASS)
library(DHARMa)
library(ResourceSelection)
library(tidyverse)

# 1. Cargar datos 

rD.germ <- read.csv("00-data/sensibilidad_germinacion.csv", sep=";")
str(rD.germ) # Ver estructura de los datos

# 2. Pasar a formato largo

df.germ <- rD.germ |>
  dplyr::select(matches("^(id_bellota|tiempo|DW|MC|germ|emerg)_")) |>
  pivot_longer(
    cols          = everything(),
    names_to      = c(".value", "especie"),
    names_pattern = "(.+)_([A-Z]+)$"
  )

df.germ <- df.germ |> 
  mutate(especie = recode(especie, 
                          "CO" = "Quercus coccfiera",
                          "IL" = "Quercus ilex",  
                          "FA" = "Quercus faginea",
                          "SU" = "Quercus suber",
                          "PU" = "Quercus pubescens", 
                          "PY" = "Quercus pyrenaica", 
                          "PE" = "Quercus petraea", 
                          "RO" = "Quercus robur")
  )
# quitar NA
summary(df.germ)
df.germ <- df.germ |> drop_na(germ, emerg, MC)

# 3. Exploración de los datos
## 3.1. grafico de Germinacion en función de MC facetado por especie

ggplot(df.germ, aes(x = MC, y = germ)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = binomial)) +
  facet_wrap(~especie, nrow = 2)

## 3.2. tabla de germinacion promedio y MC promedio en cada combinacion tiempo-especie

resumen <- df.germ |> 
  group_by(especie, tiempo) |> 
  summarise(n = n(), 
            media.MC = mean(MC), 
            media.germ = mean(germ)
  )

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

p <- ggplot(resumen, aes(x = media.MC, y = media.germ)) +
  geom_point(color = okabe_ito[1], alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = FALSE, color = okabe_ito[1], alhpa = 0.8) +
  geom_point(data = df.germ, aes(x = MC, y = germ), color = okabe_ito[5], alpha = 0.2) +
  geom_smooth(data = df.germ, aes(x = MC, y = germ), method = "glm", method.args = list(family = binomial), se = FALSE, color = okabe_ito[5], alhpa = 0.8) +
  facet_wrap(~ especie, nrow = 2) +
  labs(
    title    = "Germination probability as a function of moisture content",
    subtitle = "Orange: GLM fit on summarised data · Blue: GLM fit on individual data"
  ) +
  theme_minimal() +
  xlab("Moisture Content") +
  ylab("Germination ratio")

ggsave(
  filename = "07-img/germination_vs_moisture.png",
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Según la exploración gráfica los modelos glm ajustados a
# a) los datos resumidos y
# b) a los datos individuales
# se parecen bastante. 
# Parece que puede haber algunas diferencias en el MC50
# pero creo que son mínimas y, en cualquier caso, 
# entendemos la aproximación con todos los datos como más precisa

# Además: 
# La línea recta en RO indica que este lote de Q. robur
# estaba en mal estado y no germinó bien
# para ajustar los modelos se elimina del análisis

df.germ <- df.germ |> filter(especie != "Quercus robur")

# Centrar el contenido hídrico
# podría estar provocando mutlicolinealidad
MC_media <- mean(df.germ$MC)
df.germ <- df.germ |>
  mutate(
    MC_c = MC - MC_media
  )

# 4. Ajuste del modelo logístico binomial

library(glmmTMB)

# 4.1. Elección de modelos

# modelo nulo
glm.0 <- glmmTMB(germ ~ 1, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.1 <- glmmTMB(germ ~ MC, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.2 <- glmmTMB(germ ~ MC + especie, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.3 <- glmmTMB(germ ~ MC * especie, 
                 data = df.germ, 
                 family = binomial(link = "logit"))



# Comparación
library(performance)

compare_performance(glm.0, glm.1, glm.2, glm.3)

# 4.2. Revisar supuestos

# a) Vista general
check_model(glm.3)
# Parece haber problemas de colinealidad

# Centrar el contenido hídrico
# podría estar provocando mutlicolinealidad
MC_media <- mean(df.germ$MC)
df.germ <- df.germ |>
  mutate(
    MC_c = MC - MC_media
  )

glm.0 <- glmmTMB(germ ~ 1, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.1 <- glmmTMB(germ ~ MC_c, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.2 <- glmmTMB(germ ~ MC_c + especie, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

glm.3 <- glmmTMB(germ ~ MC_c * especie, 
                 data = df.germ, 
                 family = binomial(link = "logit"))

compare_performance(glm.0, glm.1, glm.2, glm.3)

check_model(glm.3)


# b) Simulación de residuos
library(DHARMa)
res <- simulateResiduals(glm.3, n = 1000)
plot(res)
testUniformity(res)
testDispersion(res)
testOutliers(res)
# Ningun test postivo
# pinta bien

# c) Linealidad en la escala logit
df.germ$pred <- predict(glm.3, type = "response")

ggplot(df.germ, aes(MC_c, germ)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess") +
  geom_line(aes(y = pred, group = especie), color = "red")


# d) colinealidad
check_collinearity(glm.3)
# Problema de colinealidad

# Estimación del MC50
## Haremos modelos ajustados separadamente por especie
## Razones
## 1) La colinealidad dejará de ser un problema importante porque desaparece el coeficiente MC*especie
## 2) Más importante, como queremos precisión en la estimación, probablemente mejore la precisión haciendo las estimaciones por separado por especie

library(dplyr)
library(purrr)
library(furrr)
library(glmmTMB)
library(parallel)

# 1) Paralelización
plan(multisession, workers = parallel::detectCores() - 1)

# 2) Función bootstrap: devuelve un vector con los MC50 de cada réplica
boot_fun_species <- function(dat, R = 5000) {
  res <- numeric(R)
  
  for (i in seq_len(R)) {
    idx <- sample(seq_len(nrow(dat)), replace = TRUE)
    d <- dat[idx, ]
    
    mod <- glmmTMB(germ ~ MC, data = d, family = binomial(link = "logit"))
    b <- fixef(mod)$cond
    
    res[i] <- -b["(Intercept)"] / b["MC"]
  }
  
  res
}

# 3) Separar por especie
df_split <- split(df.germ, df.germ$especie)



# 4) Guardar TODAS las iteraciones en boot_all
boot_all <- future_map_dfr(names(df_split), function(sp) {
  dat <- df_split[[sp]]
  res <- boot_fun_species(dat, R = 5000)
  
  tibble(
    especie = sp,
    iter = seq_along(res),
    MC50 = res
  )
}, .options = furrr_options(seed = 123))

boot_all %>%
  filter(especie == "Q. petraea") %>%
  summarise(
    sd = sd(MC50),
    min = min(MC50),
    max = max(MC50),
    n_unique = n_distinct(MC50), 
    media = mean(MC50), 
    Q02.5 = quantile(MC50, 0.025), 
    Q97.5 = quantile(MC50, 0.975)
  )

ggplot(boot_all, aes(x = MC50)) +
  geom_histogram() +
  facet_wrap(~especie)

# 5) Calcular resumen e IC95% por especie
boot_results <- boot_all %>%
  group_by(especie) %>%
  summarise(
    MC50_mean = mean(MC50, na.rm = TRUE),
    low  = quantile(MC50, 0.025, na.rm = TRUE),
    high = quantile(MC50, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 6) Guardar objetos
attr(boot_all, "R") <- 5000
attr(boot_all, "seed") <- 123
attr(boot_all, "model") <- "germ ~ MC (por especie)"

saveRDS(boot_all, "data/MC50_boot_all_5000.rds")
saveRDS(boot_results, "data/MC50_boot_results_5000.rds")

# 5) Representación
ggplot(df.germ, aes(MC, germ)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE) +
  facet_wrap(~ especie)
theme_light()
boot_results$especie <- factor(boot_results$especie, levels = c("CO", "IL", "SU", "FA", "PU", "PY", "PE"))

bioclimatic_regions <- 
  tribble(
    ~especie,       ~`bioclimatic region`,
    "Q. coccfiera",       "Mediterranean",
    "Q. ilex",       "Mediterranean",
    "Q. suber",       "Mediterranean",
    "Q. faginea",   "sub-Mediterranean",
    "Q. pubescens",   "sub-Mediterranean",
    "Q. pyrenaica",   "sub-Mediterranean",
    "Q. petraea",            "Atlantic"
  )

boot_results <- merge(x = boot_results, y = bioclimatic_regions, 
                      by = "especie")

boot_results$`bioclimatic region` <- as.factor(boot_results$`bioclimatic region`)


source("funs_publication.R")

p <- ggplot(boot_results, aes(
  y = fct_reorder(especie, as.numeric(as.factor(`bioclimatic region`))), 
  x = MC50_mean, 
  colour = `bioclimatic region`
)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = low, xmax = high), linewidth = 1, width = 0.3) +
  scale_colour_manual(values = okabe_ito) +
  labs(x = "MC50 Mean", y = NULL, colour = NULL) +
  theme_publication() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.7),
      colour = NA
    ),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 8)
  )

save_figure_dual(p, "img/meanic_germinacion")