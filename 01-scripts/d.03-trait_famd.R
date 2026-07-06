library(tidyverse)

source("01-scripts/d.01-load_desiccation_exp.R")

rm(list = setdiff(ls(), "df.bellotas"))

# 1. Calculate FAMD ----
## 1.1. Create famd dataframe ----
df <- df.bellotas |>
  filter(cotiledon_anormal %in% c(0, 1)) |>
  filter(rajas_pericarpo %in% c(0, 1)) |> 
  dplyr::select(id_bellota, especie, procedencia, codigo, 
                peso_seco, Volumen_estimado_cm3, Relacion_SV, 
                SPM_g_cm2, Seed_Coat_Ratio, 
                Ratio_A.cicatriz_A.bellota, 
                rajas_pericarpo) |> 
  rename(species = especie, 
         provenance = procedencia, 
         prov_code = codigo, 
         dry_weight = peso_seco, 
         volume_cm3 = Volumen_estimado_cm3,
         SVR = Relacion_SV,
         SCR = Seed_Coat_Ratio,
         SSR = Ratio_A.cicatriz_A.bellota,
         pericarp_rupture = rajas_pericarpo
  ) |> 
  unique() |> 
  drop_na() |> 
  mutate(pericarp_rupture = as.factor(pericarp_rupture))

df |> 
  group_by(species, provenance) |> 
  summarise(
    n = n(), 
    prop = n/30
  )

## 1.2. Calculate and save famd ----

library(FactoMineR)  # Para FAMD
library(factoextra)


# Exclude id and species variables
df.famd <- df |> select(-id_bellota, -species, -provenance, -prov_code)

# Calculate famd
famd.traits <- FAMD(df.famd, graph = FALSE)

### Save famd info ----
# If not created, create famd folder
# dir.create("07-img/FAMD_outputs", showWarnings = FALSE)

library(patchwork)

# ---------------------------
# 1. Eigenvalues / Scree plot
# ---------------------------
p1 <- fviz_screeplot(famd.traits, addlabels = TRUE, ylim = c(0, 50))

ggsave("07-img/FAMD_outputs/01_scree_plot.png", p1, width = 7, height = 5)


# ---------------------------
# 2. Individuos (PC1 vs PC2)
# ---------------------------
p2.1 <- fviz_famd_ind(
  famd.traits,
  axes = c(1, 2),
  label = "none",
  col.ind = "cos2",
  gradient.cols = c("grey", "blue", "red"),
  title = "Individuals (colored by cos2)"
)

p2.2 <- fviz_famd_ind(
  famd.traits,
  axes = c(1, 3),
  label = "none",
  col.ind = "cos2",
  gradient.cols = c("grey", "blue", "red"),
  title = "Individuals (colored by cos2)"
) + labs(title = "")

p2 <- p2.1 + p2.2 + plot_layout(guides = "collect")

ggsave("07-img/FAMD_outputs/02_individuals_cos2.png", p2, width = 7, height = 6)

# ---------------------------
# 3. Variables cuantitativas y cualitativas
# ---------------------------
p3.1 <- fviz_famd_var(
  axes = c(1, 2),
  famd.traits,
  choice = "var",
  col.var = "contrib",
  gradient.cols = c("grey", "blue", "red"),
  repel = TRUE
)

p3.2 <- fviz_famd_var(
  axes = c(1, 3),
  famd.traits,
  choice = "var",
  col.var = "contrib",
  gradient.cols = c("grey", "blue", "red"),
  repel = TRUE
)

p3 <- p3.1 + p3.2

ggsave("07-img/FAMD_outputs/03_variables.png", p3, width = 7, height = 6)

# 2. Create FAMD table ----

# --- Variables cuantitativas ---
quanti_df <- famd.traits$quanti.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd.traits$quanti.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd.traits$quanti.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cuantitativa")

# --- Variables cualitativas ---
quali_df <- famd.traits$quali.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd.traits$quali.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd.traits$quali.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cualitativa")

# --- Combinar todo ---
tabla_famd <- bind_rows(quanti_df, quali_df)

tabla_famd

write.csv2(x = tabla_famd, file = "00-data/famd_long.csv")

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

# Paso 3: pasar a formato ancho
tabla_wide <- tabla_contrib %>%
  select(variable, dimension, contrib_signo) %>%
  pivot_wider(names_from = dimension, values_from = contrib_signo)

write.csv2(x = tabla_wide, file = "00-data/paper_famd.csv")

# 3. Create FAMD figure ----

## 3.1 Extraer datos del FAMD ----
library(FactoMineR)
library(ggplot2)
library(dplyr)

# Individuos
ind <- as.data.frame(famd.traits$ind$coord)

# Variables (cuantitativas + cualitativas juntas en FAMD)
# var <- as.data.frame(famd.traits$var$coord)
var <- famd.traits$quanti.var$coord
var <- as.data.frame(var)

# 3.2 Añadir Dim3 como color en individuos ----
# ind$Dim3 <- ind$Dim3

# (Ojo: en FactoMineR suele llamarse Dim.3)

# Corrigiendo forma robusta:
ind$Dim3 <- ind[, 3]

ind2 <- 
  cbind(ind, df |> dplyr::select(id_bellota, pericarp_rupture)) |> 
  mutate(dim3_bin = cut(Dim3, breaks = 2),
         dim3_bin = case_when(
           dim3_bin == "(-1.98,1.11]" ~ "low",
           dim3_bin == "(1.11,4.19]" ~ "high",
           TRUE ~ NA_character_ # Opcional: para manejar valores fuera de rango
         ), 
         dim3_bin = factor(dim3_bin, levels = c("low", "high")))

# 3.3 Escalado de variables (CLAVE para que las flechas tengan sentido) ----
# Las variables están en otra escala → hay que reescalar:
  
  mult <- min(
    (max(ind2$Dim.1) - min(ind2$Dim.1)) / (max(var$Dim.1) - min(var$Dim.1)),
    (max(ind2$Dim.2) - min(ind2$Dim.2)) / (max(var$Dim.2) - min(var$Dim.2))
  )

var_scaled <- var * mult * 0.7

# 3.4) Biplot ggplot (lo que quieres exactamente)----
p <- 
ggplot() +
  
# -------------------
# Dim3 area
# -------------------

geom_hex(
  data = ind2,
  aes(x = Dim.1, y = Dim.2, fill = dim3_bin), 
  alpha = .5
  ) +
  
# -------------------
# individuos
# -------------------
geom_point(
  data = ind2,
  aes(x = Dim.1, y = Dim.2),
  alpha = 0.2,
  size = 1.5
) +
  
  # -------------------
# flechas variables
# -------------------
geom_segment(
  data = var_scaled,
  aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
  arrow = arrow(length = unit(0.2, "cm")),
  color = "black",
  alpha = 0.6
) +
  
ggrepel::geom_text_repel(
  data = var_scaled,
  aes(x = Dim.1, y = Dim.2, label = rownames(var)),
  size = 3,
  color = "black",
  box.padding = 0.5,        # Espacio extra alrededor de la etiqueta
  point.padding = 0.5,      # Espacio extra alrededor del punto
  segment.color = "grey50"  # Dibuja una línea hacia el punto original
) +
  
# -------------------
# estética
# -------------------
scale_color_viridis_c() +
  theme_minimal() +
  
  labs(
    x = "Dim 1",
    y = "Dim 2",
    color = "Dim 3"
  )

ggsave("07-img/FAMD_outputs/paper_famd.png", p, width = 7, height = 6)

