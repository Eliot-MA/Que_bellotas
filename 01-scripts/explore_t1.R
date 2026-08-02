# Exploracion rapida de df.t1 (estructura de los datos para filogenia)
library(tidyverse)

df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(log.t = log(time + 1), sqrt.t = sqrt(time + 1))

df.t1 <- df |> filter(time < 90)

cat("nrow df.t1:", nrow(df.t1), "\n")
cat("n obs:", nrow(df.t1), "| n filas por bellota (min-max):\n")

# Resumen por id_bellota (¿el id se repite entre especies?)
tab <- df.t1 |>
  group_by(codigo, species) |>
  summarise(
    n_acorns = n_distinct(id_bellota),
    ids = paste(sort(unique(id_bellota)), collapse = ","),
    n_obs = n(),
    .groups = "drop"
  )
print(tab, n = 50)

cat("\n-- Numero total de bellotas por especie --\n")
print(df.t1 |> distinct(codigo, id_bellota, species) |> count(species), n = 50)

cat("\n-- id_bellota unicos GLOBALES:", length(unique(df.t1$id_bellota)), "\n")
cat("-- id_bellota unicos DENTRO de cada especie (max):",
    max(df.t1 |> group_by(species) |> summarise(n = n_distinct(id_bellota)) |> pull(n)), "\n")

cat("\n-- rango de time y MC --\n")
cat("time:", range(df.t1$time), "| MC:", range(df.t1$Moisture_content), "\n")

cat("\n-- check nesting: codigo es 1:1 con species? --\n")
print(df.t1 |> distinct(codigo, species))
