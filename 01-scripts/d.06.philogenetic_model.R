# Load data ----
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

# Libraries ----
library(tidyverse)
library(moments)
library(glmmTMB)

# Creat working dataframe ----
df <- df.bellotas |> 
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |> 
  left_join(y = df.famd, by = "id_bellota") |> 
  drop_na(Dim.1, Dim.2, Dim.3) |> 
  rename(time = tiempo_acumulado_horas) |> 
  mutate(
    log.t = log(time+1), 
    sqrt.t = sqrt(time+1)
  )

# Separate data in two datasets
df.t1 <- df |> 
  filter(time < 90)

df.t2 <- df |> 
  filter(time > 90)

##
# Model before breakpoint ----
##

# Fit model without phylogenies
m.t1 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (time|id_bellota:provenance), data = df.t1)

m.t2 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (time|id_bellota:provenance), data = df.t2)

summary(m.t2)



