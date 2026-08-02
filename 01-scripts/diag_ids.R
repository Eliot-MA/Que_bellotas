suppressPackageStartupMessages(library(dplyr))
d <- read.csv("00-data/desiccation_traits_long.csv", stringsAsFactors = FALSE)

cat("nrow(d):", nrow(d), "\n")
cat("n_distinct(id_bellota):", n_distinct(d$id_bellota), "\n")
cat("n_distinct(codigo):", n_distinct(d$codigo), "\n")
cat("n_distinct(numero_bellota):", n_distinct(d$numero_bellota), "\n")
cat("n_distinct(codigo, numero_bellota):", d |> distinct(codigo, numero_bellota) |> nrow(), "\n")
cat("n_distinct(codigo, numero_bellota, tiempo):", d |> distinct(codigo, numero_bellota, tiempo) |> nrow(), "\n")

cat("\nhead d:\n")
print(head(d[, c("id_bellota", "codigo", "numero_bellota", "tiempo")], 15))

cat("\nid_bellota unicos por codigo:\n")
print(d |> group_by(codigo) |> summarise(n_id = n_distinct(id_bellota), n_num = n_distinct(numero_bellota), nrows = n()) |> as.data.frame())

cat("\nnumero_bellota se repite dentro del mismo codigo? (a distintas filas de tiempo):\n")
dup <- d |> group_by(codigo, numero_bellota) |> summarise(n_tiempos = n_distinct(tiempo), .groups="drop")
cat("combos codigo x numero_bellota:", nrow(dup), "| con mas de un tiempo:", sum(dup$n_tiempos > 1), "\n")
print(head(dup |> filter(n_tiempos > 1), 10) |> as.data.frame())
