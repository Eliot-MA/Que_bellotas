# Load data ----
library(tidyverse)
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

# Creat working dataframe ----
df <- df.bellotas |> 
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |> 
  left_join(y = df.famd, by = "id_bellota") |> 
  drop_na(Dim.1, Dim.2, Dim.3) |> 
  rename(time = tiempo_acumulado_horas) |> 
  mutate(time_s = scale(time))

# Separate data in two datasets
t94 <- as.vector((94 - mean(df$time)) / sd(df$time))

df.t1 <- df |> 
  filter(time_s < t94)

df.t2 <- df |> 
  filter(time_s > t94)

# 1) Comprobar que los PCs son ortogonales globalmente pero no dentro de especie ----
  
## Correlación global

df <- df.t1 |> 
  select(id_bellota, species, provenance, Dim.1, Dim.2, Dim.3) |> 
  distinct()

library(dplyr)

df |>
  select(Dim.1, Dim.2, Dim.3) |>
  distinct() |> 
  cor()

# Debería salir (aproximadamente)
# Dim1 Dim2 Dim3
# Dim1   1    0    0
# Dim2   0    1    0
# Dim3   0    0    1

## Correlaciones por especie
  

library(purrr)

cor_species <-
  df |> 
  group_by(species) |>
  group_modify(~{
    cormat <- cor(select(.x, Dim.1, Dim.2, Dim.3))
    
    tibble(
      cor12 = cormat["Dim.1","Dim.2"],
      cor13 = cormat["Dim.1","Dim.3"],
      cor23 = cormat["Dim.2","Dim.3"]
    )
  })

cor_species


library(GGally)

GGally::ggpairs(
  df,
  columns = c("Dim.1","Dim.2","Dim.3"),
  mapping = aes(colour = species, alpha = .5)
)



split(df, 
      df$species) |>
  purrr::walk(~GGally::ggpairs(.x,
                               columns = c("Dim.1","Dim.2","Dim.3")))
  
# 2) Buscar el posible efecto de Simpson ----

## Paso 1. Calcular la pendiente por bellota

library(broom)
pendientes <-
  df.t1 |>
  group_by(id_bellota, species, provenance,
           Dim.1, Dim.2, Dim.3) |>
  do(tidy(lm(Moisture_content ~ time_s,
             data = .))) |>
  filter(term == "time_s") |>
  ungroup() |>
  rename(slope = estimate)

# Las pendientes serán negativas.
# Cuanto más negativa, más rápida la desecación.

## Pendiente frente a Dim2 global
  

library(ggplot2)

ggplot(pendientes,
       aes(Dim.3, slope))+
  geom_point(alpha=.4)+
  geom_smooth(method="lm")

## Pendiente frente a Dim2 por especie
  
ggplot(pendientes,
       aes(Dim.3, slope,
           colour=species))+
  geom_point(alpha=.5)+
  geom_smooth(method="lm", se=FALSE)

## O mejor aún
  
ggplot(pendientes,
       aes(Dim.3, slope))+
  geom_point(alpha=.4)+
  geom_smooth(method="lm")+
  facet_wrap(~species)

modelo <- lm(slope ~ Dim.2*species, data = pendientes)

summary(modelo)


# 3) Separar efectos within y between ----
# Este análisis es probablemente el más importante.


df.within <-
  df |>
  group_by(species) |>
  mutate(
    Dim2_between = mean(Dim.2),
    Dim2_within  = Dim.2 - mean(Dim.2)
  ) |>
  ungroup()


#Comprobación


df.t1.within <- df.t1 |>
  group_by(species) |>
  mutate(
    Dim2_between = mean(Dim.2),
    Dim2_within = Dim.2 - mean(Dim.2)
  ) |>
  ungroup()


# El within debe tener media 0.

# Modelo

mm.within <- glmmTMB(
  Moisture_content ~
    time_s * (Dim2_between + Dim2_within +
                Dim.1 + Dim.3) +
    (time_s | id_bellota),
  data = df.t1.within
)

parameters::model_parameters(mm.within)


# Si
# time × Dim2_between es positivo y
# time × Dim2_within es negativo
# ya tienes la explicación del cambio de signo.

  
# 4) ¿Cuánta variación de Dim2 es entre especies? ----
## Mediante ANOVA
  

summary(
  aov(Dim.2 ~ species,
      data=df)
)

anova(
  lm(Dim.2 ~ species,
     data=df)
)


summary(
  lm(Dim.2 ~ species,
     data=df)
)
# El R² te dirá qué porcentaje de la variación de Dim2 explica la especie.

  
## Mediante ICC
library(lme4)
library(performance)

m.icc <-
  lmer(
    Dim.2 ~ (1|species),
    data=df
  )

performance::icc(m.icc)

# Si sale
# ICC = 0.92
#
  
# 5) Inspeccionar los modelos ----
## Efectos fijos

mm.pre.0 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.1 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (1|species) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.2 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (1|provenance/species) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.3 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (time_s||species) + 
            (1|provenance) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.4 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (1|species) + 
            (time_s|provenance) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.5 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (time_s||species) + 
            (time_s|provenance) + 
            (time_s|id_bellota), 
          data = df.t1)

mm.pre.6 <- 
  glmmTMB(Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + 
            (time_s||species/provenance) + 
            (time_s|id_bellota), 
          data = df.t1)


fixef(mm.pre.0)

fixef(mm.pre.2)

fixef(mm.pre.5)

## Varianzas de efectos aleatorios
  
VarCorr(mm.pre.2)

VarCorr(mm.pre.5)

## Descomposición de la varianza
  
performance::variance_decomposition(mm.pre.5)

df |> 
  group_by(species) |> 
  summarise(
    n = n(), 
    minimo = min(Dim.2), 
    q25 = quantile(Dim.2, probs = .25), 
    q50 = quantile(Dim.2, probs = .5), 
    media = mean(Dim.2), 
    q75 = quantile(Dim.2, probs = .75), 
    maximo = max(Dim.2)
  )


