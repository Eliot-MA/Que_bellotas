library(FactoMineR)  # Para PCA y FAMD
library(factoextra)  # Para visualización
library(gridExtra)   # Para dividir espacio de graficación
library(ggrepel)     # para etiquetas no sobrepuestas (opcional)
library(dplyr)
library(tidyverse)

source("01.CargarDatos.R")
source("02.CrearVariablesDerivadas.R")
df <- df.bellotas

# Variables de interés
vars_interes <- c("id_bellota", "especie", "procedencia", "codigo",
                  "Moisture_content",
                  "peso_seco", "Volumen_estimado_cm3", "Area_estimada_cm2", "Relacion_SV", # Tamaño
                  "SPM_g_cm2", "Seed_Coat_Ratio", # Robustez de pericarpo
                  "Ratio_A.cicatriz_A.bellota",   # Tamaño de la cicatriz
                  "rajas_pericarpo")   

df <- df %>%
  drop_na(all_of(vars_interes)) %>%           # Los NA dan error en el PCA
  filter(cotiledon_anormal %in% c(0, 1)) |>   # Bellotas sin defectos anómalos
  dplyr::select(all_of(vars_interes))

# Eliminar pseudorreplicación (una observación por bellota)
df.unique <- df |> dplyr::distinct(id_bellota, .keep_all = TRUE)

# Asegurar rajas como factor
df.unique$rajas_pericarpo <- factor(df.unique$rajas_pericarpo)

##################
# FAMD con rajas #
#   se quitan    #
# rajas tipo 2   #
##################

# df de solo rasgos y CON variable rajas
df.rasgos <- df.unique |> 
  select(-id_bellota, -Moisture_content, -especie, -procedencia, -codigo)

# Quitar rajas tipo 2, son demasiado pocas
df.rasgos <- df.rasgos |> 
  filter(rajas_pericarpo != 2)
df.rasgos$rajas_pericarpo <- factor(df.rasgos$rajas_pericarpo) # Reprogramar variable rajas para que solo tenga dos niveles

# FAMD
famd_rasgos <- FAMD(df.rasgos, ncp = 5, graph = FALSE)

# Unir dimensiones del famd al dataset original
# 1. Extraemos las coordenadas de las primeras 3 dimensiones
coords_famd <- as.data.frame(famd_rasgos$ind$coord[, 1:3])

# 2. Necesitamos el ID para unir. 
coords_famd$id_bellota <- df.unique |> filter(rajas_pericarpo != 2) |> pull(id_bellota)

# 3. Unimos las dimensiones a la tabla grande (df)
df_famd <- df %>%
  filter(rajas_pericarpo != 2) |> 
  left_join(coords_famd, by = "id_bellota")

# Escalar variables
df_famd_scaled <- df_famd %>%
  mutate(across(
    .cols = -c(id_bellota, especie, procedencia, codigo, Moisture_content, rajas_pericarpo, starts_with("Dim.")), 
    .fns = ~ as.vector(scale(.))
  )) |> 
  mutate(rajas_pericarpo = factor(rajas_pericarpo))

###
# PREPARAR INFORMACION PARA ARTÍCULO
# TABLAS
##

# 1) Extraer y combinar variables cuantitativas y cualitativas
library(dplyr)
library(tidyr)
library(tibble)

# --- Variables cuantitativas ---
quanti_df <- famd_rasgos$quanti.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd_rasgos$quanti.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd_rasgos$quanti.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cuantitativa")

# --- Variables cualitativas ---
quali_df <- famd_rasgos$quali.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd_rasgos$quali.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd_rasgos$quali.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cualitativa")

# --- Combinar todo ---
tabla_famd <- bind_rows(quanti_df, quali_df)

tabla_famd
# 2) Limpiar y ordenar la tabla
# Construir la tabla (variables × dimensiones con contribuciones)
# formato ancho

# Paso 1: quedarte solo con contribuciones + coordenadas
tabla_contrib <- tabla_famd %>%
  select(variable, dimension, contrib, coord)

# Paso 2: codificar el signo
# Aquí tienes varias opciones. La más limpia es añadir el signo como string:
  
  tabla_contrib <- tabla_contrib %>%
  mutate(
    signo = ifelse(coord >= 0, "+", "−"),
    contrib_signo = paste0(round(contrib, 1), signo)
  )
  
  tabla_contrib
# Paso 3: pasar a formato ancho
tabla_wide <- tabla_contrib %>%
  select(variable, dimension, contrib_signo) %>%
  pivot_wider(names_from = dimension, values_from = contrib_signo)

tabla_wide

# Paso 4: guardar informacion
# write.csv2(tabla_wide, "data/tabla_contribuciones_FAMD.csv", row.names = FALSE)
# 
# library(knitr)
# library(kableExtra)
# 
# tabla_wide_tex <- 
#   tabla_wide %>%
#   kable(format = "latex", booktabs = TRUE,
#         caption = "Contribuciones (%) de las variables a cada dimensión del FAMD. El signo indica la dirección de la asociación.",
#         align = "lccccc") %>%
#   kable_styling(latex_options = c("hold_position"))
# 
# writeLines(tabla_wide_tex, con = "latex/tabla_famd.tex")

##
# FIGURAS
##

# 1) Preparar los scores individuales del FAMD

library(dplyr)
library(ggplot2)
library(patchwork)
library(plotly)

coords <- as.data.frame(famd_rasgos$ind$coord[, 1:3])
colnames(coords) <- c("Dim1", "Dim2", "Dim3")

ides <- df.unique |> filter(rajas_pericarpo != 2) |> pull(id_bellota)
especies <- df.unique |> filter(rajas_pericarpo != 2) |> pull(especie)

coords <- cbind(ides, especies, coords)

coords <- coords %>%
  mutate(especie = factor(especies, levels = c("Quercus coccifera", "Quercus ilex", "Quercus suber", "Quercus faginea", "Quercus pubescens", "Quercus pyrenaica", "Quercus petraea", "Quercus robur")))

# 2) Función para gráficos 2D por pares de dimensiones
# Esta te permite reutilizar el mismo estilo para los tres paneles.

plot_famd_pair <- function(data, x, y) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = especie)) +
    geom_point(alpha = 0.75, size = 2) +
    stat_ellipse(level = 0.68, linewidth = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey70") +
    labs(x = x, y = y, color = "Especie") +
    theme_minimal(base_size = 12)
}
# 3) Los tres gráficos 2D
p12 <- plot_famd_pair(coords, "Dim1", "Dim2")
p13 <- plot_famd_pair(coords, "Dim1", "Dim3")
p23 <- plot_famd_pair(coords, "Dim2", "Dim3")

p12 + p13 + p23 + plot_layout(ncol = 3, guides = "collect")

# 4) Añadir centroides de especie
# Esto ayuda mucho a ver la separación global.

centroides <- coords %>%
  group_by(especie) %>%
  summarise(
    Dim1 = mean(Dim1, na.rm = TRUE),
    Dim2 = mean(Dim2, na.rm = TRUE),
    Dim3 = mean(Dim3, na.rm = TRUE),
    .groups = "drop"
  )

# Y luego, por ejemplo para Dim1–Dim2:
  
p12 <- ggplot(coords, aes(Dim1, Dim2, color = especie)) +
  geom_point(alpha = 0.1, size = 2) +
  stat_ellipse(level = 0.68, linewidth = 0.7) +
  geom_point(data = centroides, aes(Dim1, Dim2), inherit.aes = FALSE,
             shape = 4, size = 4, stroke = 1.2) +
  geom_text(data = centroides, aes(Dim1, Dim2, label = especie),
            inherit.aes = FALSE, vjust = -0.8, show.legend = FALSE) +
  coord_equal() +
  theme_minimal(base_size = 12)

p23 <- ggplot(coords, aes(Dim2, Dim3, color = especie)) +
  geom_point(alpha = 0.1, size = 2) +
  stat_ellipse(level = 0.68, linewidth = 0.7) +
  geom_point(data = centroides, aes(Dim2, Dim3), inherit.aes = FALSE,
             shape = 4, size = 4, stroke = 1.2) +
  geom_text(data = centroides, aes(Dim2, Dim3, label = especie),
            inherit.aes = FALSE, vjust = -0.8, show.legend = FALSE) +
  coord_equal() +
  theme_minimal(base_size = 12)

p13 <- ggplot(coords, aes(Dim1, Dim3, color = especie)) +
  geom_point(alpha = 0.1, size = 2) +
  stat_ellipse(level = 0.68, linewidth = 0.7) +
  geom_point(data = centroides, aes(Dim1, Dim3), inherit.aes = FALSE,
             shape = 4, size = 4, stroke = 1.2) +
  geom_text(data = centroides, aes(Dim1, Dim3, label = especie),
            inherit.aes = FALSE, vjust = -0.8, show.legend = FALSE) +
  coord_equal() +
  theme_minimal(base_size = 12)

(p12 + p13 + p23) + plot_layout(guides = "collect")



# 5) 3D interactivo con plotly

# Esto es lo más útil para explorar.

plot_ly(
  coords,
  x = ~Dim1, y = ~Dim2, z = ~Dim3,
  color = ~especie,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4, opacity = 0.75)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Dim.1"),
      yaxis = list(title = "Dim.2"),
      zaxis = list(title = "Dim.3")
    )
  )

# 6) Tabla resumen
# Tambien para explorar

resumen_coords <- coords |>
  group_by(especies) |> 
  summarise(
    n = n(), 
    minimo_Dim1 = min(Dim1),
    Q25_Dim1 = quantile(Dim1, probs = 0.25), 
    mediana_Dim1 = quantile(Dim1, probs = 0.50), 
    media_Dim1 = mean(Dim1), 
    Q75_Dim1 = quantile(Dim1, probs = 0.75), 
    maximo_Dim1 = max(Dim1),
    minimo_Dim2 = min(Dim2),
    Q25_Dim2 = quantile(Dim2, probs = 0.25), 
    mediana_Dim2 = quantile(Dim2, probs = 0.50), 
    media_Dim2 = mean(Dim2), 
    Q75_Dim2 = quantile(Dim2, probs = 0.75), 
    maximo_Dim2 = max(Dim2),
    minimo_Dim3 = min(Dim3),
    Q25_Dim3 = quantile(Dim3, probs = 0.25), 
    mediana_Dim3 = quantile(Dim3, probs = 0.50), 
    media_Dim3 = mean(Dim3), 
    Q75_Dim3 = quantile(Dim3, probs = 0.75), 
    maximo_Dim3 = max(Dim3)
  )

# 7) Explorando estrategias de agrupamiento
## 7.1. Agrupamiento por centroides
dist_mat <- dist(centroides[, -1])  # euclídea en espacio FAMD

hc <- hclust(dist_mat, method = "ward.D2")
plot(hc)

centroides$grupo <- cutree(hc, k = 3)  # prueba 2–4 grupos

## 7.2. Agrupamiento por individuos
## mas potente
## al hacerlo por individuos, el unico criterio
## son los rasgos capturados en las dimensiones
## del famd, no hay criterio de especie

kmeans_res <- kmeans(coords[,3:5], centers = 3)
coords$cluster <- kmeans_res$cluster
coords$cluster <- factor(coords$cluster)
table(coords$cluster, coords$especie)
coords$especie <- sub("^Quercus ", "Q. ", coords$especie)
centroides$especie <- sub("^Quercus ", "Q. ", centroides$especie)

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

## 7.3. Visualizar
p12 <- ggplot(coords, aes(Dim1, Dim2)) +
  geom_point(aes(colour = especie), alpha = 0.3, size = 2) +
  stat_ellipse(aes(group = cluster), level = 0.65) +
  geom_point(data = centroides, aes(Dim1, Dim2, colour = especie), size = 4, stroke = 1.2) +
  geom_text_repel(
    data = centroides,
    aes(Dim1, Dim2, label = especie),
    inherit.aes = FALSE,
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.alpha = 0.5,
    show.legend = FALSE
  ) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = okabe_ito) +
  theme(legend.position = "none")

p23 <- ggplot(coords, aes(Dim3, Dim2)) +
  geom_point(aes(colour = especie), alpha = 0.3, size = 2) +
  stat_ellipse(aes(group = cluster), level = 0.65) +
  geom_point(data = centroides, aes(Dim3, Dim2, colour = especie), size = 4, stroke = 1.2) +
  geom_text_repel(
    data = centroides,
    aes(Dim3, Dim2, label = especie),
    inherit.aes = FALSE,
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.alpha = 0.5,
    show.legend = FALSE
  ) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = okabe_ito) +
  theme(legend.position = "none")

p13 <- ggplot(coords, aes(Dim1, Dim3)) +
  geom_point(aes(colour = especie), alpha = 0.3, size = 2) +
  stat_ellipse(aes(group = cluster), level = 0.65) +
  geom_point(data = centroides, aes(Dim1, Dim3, colour = especie), size = 4, stroke = 1.2) +
  geom_text_repel(
    data = centroides,
    aes(Dim1, Dim3, label = especie),
    inherit.aes = FALSE,
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.alpha = 0.5,
    show.legend = FALSE
  ) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = okabe_ito) +
  theme(legend.position = "none")



(p12 + p23 + p13 + plot_spacer()) +
  plot_layout(ncol = 2)

plot_ly(
  coords,
  x = ~Dim1, y = ~Dim2, z = ~Dim3,
  color = ~cluster,
  type = NULL,
  mode = "markers",
  marker = list(size = 4, opacity = 0.75)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Dim.1"),
      yaxis = list(title = "Dim.2"),
      zaxis = list(title = "Dim.3")
    )
  )
