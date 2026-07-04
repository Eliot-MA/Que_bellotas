
library(dplyr)
library(stringr)
library(tidyr)

# Cargar excels (raw data)
## data, wide format
rD.bellotas <- read.csv2(file = "00-data/desiccation_traits_wide.csv")
## codes used in observations
rD.codobs <- read.csv2(file = "00-data/code_observations.csv")

rD.sp.procedencias <- read.csv2(file = "00-data/procedencias_updated.csv", stringsAsFactors = FALSE)
rD.sp.procedencias <- rD.sp.procedencias |> 
  filter(Localidad != "El Pozo") |> 
  dplyr::select(ID, Procedencia, Localidad)

rD.sp.procedencias <- rD.sp.procedencias |> 
  dplyr::mutate(
    Procedencia = paste(Procedencia, Localidad, sep = "-")
  ) |> 
  dplyr::select(ID, Procedencia)
  

# Crear dataframe para el valor de fecha, hora y minuto de cada ti
## Cargar datos
patrones <- c("id_bellota", "especie", "procedencia", "fecha", "hora", "minuto")
df.ti <- rD.bellotas |>
  dplyr::select(contains(patrones))

### Normalizar variables minuto
df.ti.minuto <-
  df.ti |>
  dplyr::select(contains("minuto"))

cols <- grep("^minuto\\s*t[0-9]+$", names(df.ti.minuto), value = TRUE)

df.ti.minuto[cols] <- lapply(df.ti.minuto[cols], str_pad, width = 2, pad = "0")

### Normalizar la variable hora

cols_hora <- grep("^hora t[0-9]+$", names(df.ti), value = TRUE)

df.ti[cols_hora] <- lapply(df.ti[cols_hora], function(x) {
  str_pad(as.character(x), width = 2, pad = "0")
})


df.ti <-
df.ti |>
  dplyr::select(!contains("minuto"))

df.ti <- cbind(df.ti, df.ti.minuto)


### Crear variable Fecha y Hora para cada ti

for (i in 0:9) {
  
  df.ti[[paste0("Fecha y hora t", i)]] <-
      paste(
        df.ti[[paste0("fecha.t", i)]],
        paste(df.ti[[paste0("hora.t", i)]],
              df.ti[[paste0("minuto.t", i)]],
              sep = ":")
      )
}

df.ti <- df.ti |>
  dplyr::select(
    id_bellota,
    especie,
    procedencia,
    starts_with("Fecha y hora")
  ) |>
  mutate(
    across(
      starts_with("Fecha y hora"),
      ~ as.POSIXct(.x,
                   format = "%d/%m/%Y %H:%M",
                   tz = "Europe/Madrid")
    )
  )
  

### Hacerla tidy
df.ti <- df.ti %>%
  pivot_longer(
    cols = starts_with("Fecha y hora t"),  # Seleccionar las columnas con los tiempos
    names_to = "tiempo",                  # Nombre para la columna de tiempos
    values_to = "fecha_hora"              # Nombre para la columna de fechas y horas
  ) %>%
  mutate(tiempo = gsub("Fecha y hora ", "", tiempo))  # Eliminar "Fecha y hora " de los nombres de tiempo

# Crear un data frame para la presencia de rajas, hongos o gusanos
## Cargar datos
patrones <- c("id_bellota", "especie", "procedencia", "observaciones")
df.obs.ti <- 
  rD.bellotas |>
  dplyr::select(contains(patrones))
df.obs.ti <- 
df.obs.ti |>
  pivot_longer(
    cols = starts_with("observaciones"), 
    names_to = "tiempo_observacion", 
    values_to = "observaciones"
  )

## Crear una variable para imperfecciones del pericarpo durante el periodo de secado
## 0: Sin imperfecciones
## 1: Imperfecciones en el ápice de la bellota que, en ningún caso, llegan a la mitad de la bellota
## 2: Imperfecciones que comienzan en el ápice y superan la mitad de la bellota, en algunos casos abriendo el pericarpo completamente y dejando los cotiledones al descubierto

# Todos = 0
df.obs.ti$rajas_pericarpo <- 0

# 1 = AR
sel <- grepl(pattern = "AR", x = df.obs.ti$observaciones)
df.obs.ti$rajas_pericarpo[sel] <- 1

# 2 = BR
sel <- grepl(pattern = "BR", x = df.obs.ti$observaciones)
df.obs.ti$rajas_pericarpo[sel] <- 2

## Crear una variable para imperfecciones de los cotiledones (hongos y gusanos)
## 0: Sin imperfecciones
## 1: Presencia de moho, sin observaciones de gusanos
## 2: Presencia de un gusano, sin observaciones de moho
## 3: Presencia de moho y un gusano o presencia de dos gusanos. 
## 3: También se incluyen las bellotas podridas y las que presentaban agallas
## - las agallas crecen hacia los cotiledones, haciéndolos más pequeños (ver fotos)

# Todos = 0
df.obs.ti$cotiledon_anormal <- 0

# 1 = M 
sel <- grepl(pattern = "M", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 1

# 2 = G or gusano or AJ [2 AJ y AJ (2) aparecerán como categoria 3 al sobrescribir mas adelante]
sel <- grepl(pattern = "G", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 2

sel <- grepl(pattern = "gusano", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 2

sel <- grepl(pattern = "AJ", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 2

# 3 = "M; G" or 2 or podrida or agallas
sel <- grepl(pattern = "M; G", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 3

sel <- grepl(pattern = "2", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 3

sel <- grepl(pattern = "Podrida", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 3

sel <- grepl(pattern = "Agallas", x = df.obs.ti$observaciones)
df.obs.ti$cotiledon_anormal[sel] <- 3

## Crear varias variables para imperfecciones durante la extracción de las partes de la bellota para peso seco
## PR-ST
## CST
## B1
## Falta un trozo de pericarpo

# PR-ST
df.obs.ti$pericarpo_sin_telita <- "NO"
sel <- grepl(pattern = "PR-ST", x = df.obs.ti$observaciones)
df.obs.ti$pericarpo_sin_telita[sel] <- "SI"

# CST
df.obs.ti$circulo_sin_telita <- "NO"
sel <- grepl(pattern = "CST", x = df.obs.ti$observaciones)
df.obs.ti$circulo_sin_telita[sel] <- "SI"

# n de bocados
df.obs.ti$n_bocados <- 2
sel <- grepl(pattern = "B1", x = df.obs.ti$observaciones)
df.obs.ti$n_bocados[sel] <- 1

# pericarpo completo
df.obs.ti$pericarpo_completo <- "SI"
sel <- grepl(pattern = "Falta un trozo de pericarpo", x = df.obs.ti$observaciones)
df.obs.ti$pericarpo_completo[sel] <- "NO"



df.obs.id <- df.obs.ti %>%
  group_by(id_bellota) %>%
  summarise(
    especie = first(especie),
    procedencia = first(procedencia),
    rajas_pericarpo = max(rajas_pericarpo, na.rm = TRUE),
    cotiledon_anormal = max(cotiledon_anormal, na.rm = TRUE),
    n_bocados = min(n_bocados, na.rm = TRUE),
    pericarpo_sin_telita = first(pericarpo_sin_telita[tiempo_observacion == "observaciones t0"]),
    circulo_sin_telita = first(circulo_sin_telita[tiempo_observacion == "observaciones t0"]),
    pericarpo_completo = first(pericarpo_completo[tiempo_observacion == "observaciones t0"])
  )

# Selecciona las columnas relevantes (id_bellota y pesos en este caso)
df.pesos <- rD.bellotas %>%
  dplyr::select(!starts_with(c("fecha", "hora", "minuto", "observaciones"))) |>
  pivot_longer(
    cols = starts_with("peso.t"), # Selecciona las columnas de peso
    names_to = "tiempo",         # Nueva columna para identificar el tiempo
    values_to = "peso"           # Nueva columna para los valores de peso
  ) %>%
  mutate(tiempo = gsub("peso.t", "t", tiempo)) # Ajusta los nombres de tiempo

# Crear dataframe
df.bellotas <- merge(df.pesos, df.ti[,c("id_bellota", "tiempo", "fecha_hora")], by = c("id_bellota", "tiempo"))
df.bellotas <- 
df.bellotas |>
  dplyr::select(!peso_seco_peri)

## Introducir rajas original
### Causa: cuando estuve haciendo los análisis de los tiempos en los que aparecen las rajas del pericarpo
### me di cuenta que la variable rajas_pericarpo con la que quería trabajar aplicaba el mismo valor de
### rajas pericarpo a todas las entradas de todos los tiempos, es decir: no tenía dimensión temporal, 
### aspecto esencial para hacer un análisis de la evolución temporal de las rajas.
### Con este código introduzco una nueva variable llamada "rajas_original" que introduce
### los valores originales de la tabla, lo que me permitirá hacer el análisis temporal
### ya que en esta nueva variable si puede haber diferentes valores para distintos tiempos. 
df.bellotas <- merge(x = df.bellotas, y = df.obs.id, by = c("id_bellota", "especie", "procedencia"))
names(df.obs.ti)[names(df.obs.ti)=="rajas_pericarpo"] <- "rajas_original"
df.obs.ti <- df.obs.ti |> 
  mutate(tiempo = str_extract(tiempo_observacion, "t\\d+"))

df.bellotas <- merge(x = df.bellotas, y = df.obs.ti[,c("id_bellota", "tiempo","rajas_original")], by = c("id_bellota", "tiempo"))


# Cambiar procedencias por las definitivas
df.bellotas <- merge(x = df.bellotas, rD.sp.procedencias, 
      by.x = "codigo", by.y = "ID")
df.bellotas <- df.bellotas[, !names(df.bellotas) %in% "procedencia"]
names(df.bellotas)[names(df.bellotas) == "Procedencia"] <- "procedencia"

