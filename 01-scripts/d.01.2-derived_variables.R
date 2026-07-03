# source("01.CargarDatos.R")

library(dplyr)
library(purrr)

#########################
# df.bellotas | Tiempos #
#########################

df.bellotas <- df.bellotas %>%
  arrange(id_bellota, tiempo) %>%  # Ordenamos por id_bellota y tiempo
  group_by(id_bellota) %>%
  mutate(dif_peso = peso - lag(peso)) %>%  # Calculamos la diferencia respecto al tiempo anterior
  ungroup()
df.bellotas$dif_peso[df.bellotas$tiempo == "t0"] <- 0


df.bellotas <- df.bellotas %>%
  arrange(id_bellota, tiempo) %>%  # Aseguramos el orden correcto
  group_by(id_bellota) %>%
  mutate(
    dif_tiempo_horas = as.numeric(difftime(fecha_hora, lag(fecha_hora), units = "hours")),  # Diferencia en horas
    dif_tiempo_horas = replace(dif_tiempo_horas, tiempo == "t0", 0),  # Aseguramos que en t0 sea 0
    tiempo_acumulado_horas = cumsum(dif_tiempo_horas)  # Tiempo acumulado
  ) %>%
  ungroup()

###########################
# df.bellotas | peso seco #
###########################

# peso seco de la bellota en total y sus partes
# SLA
# SCR
df.bellotas <- 
  df.bellotas |>
  mutate(
    peso_seco = peso_seco_b + peso_seco_pr + peso_seco_cot,
    peso_seco_pericarpo = peso_seco_b + peso_seco_pr,
    SLA_cm2_g = 0.3^2 * pi * n_bocados / peso_seco_b,
    SPM_g_cm2 = SLA_cm2_g^(-1),
    Seed_Coat_Ratio = peso_seco_pericarpo / peso_seco
  )

#####################
# df.bellotas | %MC #
#####################
df.bellotas <- df.bellotas %>%
  group_by(codigo, id_bellota) %>%
  mutate(peso_t0 = peso[tiempo == "t0"]) %>%
  ungroup()

df.bellotas <-
  df.bellotas |>
  mutate(
    Moisture_content = ((peso - peso_seco) / peso_t0) * 100
  )

############################
# df.bellotas | morfologia #
############################
df.bellotas <- df.bellotas |>
  mutate(
    Area_cicatriz_mm2 = cicatriz_1..mm./2 * cicatriz_2..mm./2 * pi,
    Volumen_estimado_cm3 = ((4*pi/3) * longitud..mm. * diametro_1 * diametro_2) * 0.001,
    Area_estimada_cm2 = SLA_cm2_g * peso_seco_pericarpo, 
    Ratio_A.cicatriz_A.bellota = Area_cicatriz_mm2 / (Area_estimada_cm2*100), 
    Relacion_SV = Area_estimada_cm2 / Volumen_estimado_cm3, 
    diametro_total = (diametro_1 + diametro_2) / 2
  ) 

# PESO SECO ESTIMADO (ELIMINAR)
# df.bellotas <- df.bellotas %>%
#   group_by(id_bellota) %>%
#   mutate(`Peso seco estimado (g)` = first(peso[tiempo == "t0"]) * 0.6) %>%
#   ungroup()
# 
# pesos.secos <- df.bellotas$`Peso seco estimado (g)`[df.bellotas$tiempo == "t0"]
# pesos.frescos <- df.bellotas$peso[df.bellotas$tiempo == "t0"]
# df.MC <- data.frame(1:450, pesos.secos, pesos.frescos)
# df.MC <- 
# df.MC |>
#   mutate(
#     `Peso hidrico (mg)` = (pesos.frescos - pesos.secos) * 1000
#   )

# Guardar datos 
## Tabla madre
# write.csv2(x = df.bellotas, file = "data/df.bellotas.csv")
## Solo datos de morfologia
# mi.tabla <- 
#   df.bellotas |> 
#   select(codigo, especie, numero_bellota, `longitud (mm)`, procedencia, contains(c("diametro", "cicatriz", "peso_seco")), SLA_cm2_g, SPM_g_cm2, Seed_Coat_Ratio, Area_cicatriz_mm2, Area_estimada_cm2, Volumen_estimado_cm3, Ratio_A.cicatriz_A.bellota, Relacion_SV) |> 
#   unique()
# write.csv2(x = mi.tabla, file = "data/df.bellotas.morfologia.csv")


# write.csv(df.bellotas, file = "data/df.bellotas.csv")



