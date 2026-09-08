# d.03.2-biplot_panel.R ---- figuras 1a (biplot FAMD), 1b (individuos por bioclima) y compuesta
if (!exists("bioclimate_cols")) source("01-scripts/d.03.0-utils.R")
if (!exists("ind.bio"))         source("01-scripts/d.03.1-load_famd.R")
library(patchwork)

# 3. Figure 1a: biplot FAMD ----

## 3.1 Extraer datos del FAMD ----
library(FactoMineR)
library(ggplot2)
library(dplyr)

# Okabe-Ito colores (dos que contrastan: azul #0072B2 y naranja #E69F00)
okabe_blue  <- "#0072B2"
okabe_orange <- "#E69F00"

# Individuos
ind <- as.data.frame(famd.traits$ind$coord)

# Variables (cuantitativas + cualitativas juntas en FAMD)
# var <- as.data.frame(famd.traits$var$coord)
var <- famd.traits$quanti.var$coord
var <- as.data.frame(var)

# Estandarizar nombres de variables
rownames(var)[rownames(var) == "SPM_g_cm2"] <- "SPM"
rownames(var)[rownames(var) == "volume_cm3"] <- "volume"
rownames(var)[rownames(var) == "dry_weight"] <- "mass"

# 3.2 A├▒adir Dim3 como color en individuos ----
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
# Las variables est├ín en otra escala ÔåÆ hay que reescalar:
  
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
# colores okabe-ito para los hex├ígonos
# -------------------
scale_fill_manual(
  values = c(low = okabe_blue, high = okabe_orange)
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
  segment.color = "grey50"  # Dibuja una l├¡nea hacia el punto original
) +
  
# -------------------
# est├®tica
# -------------------
theme_minimal() +
  
  labs(
    x = "Dim 1",
    y = "Dim 2",
    fill = "Dim 3 bin"
  )

ggsave("07-img/FAMD_outputs/paper_famd_1a.png", p, width = 7, height = 6)


## 3.2 Panel 1b ----
# Centroide de cada especie en el plano del biplot (media de Dim.1 y Dim.2)
spe.cent <- ind.bio |>
  group_by(species) |>
  summarise(Dim.1 = mean(Dim.1),
            Dim.2 = mean(Dim.2),
            bioclimate = unique(bioclimate),
            n_spots = n(),
            .groups = "drop") |>
  mutate(species_short = gsub("Quercus ", "Q. ", species))

p1b <- ggplot() +
  # individuos coloreados por bioclima
  geom_point(data = ind.bio,
             aes(x = Dim.1, y = Dim.2, colour = bioclimate),
             alpha = 0.35, size = 1.2) +
  # centroides por especie, coloreados por bioclima
  geom_point(data = spe.cent,
             aes(x = Dim.1, y = Dim.2, colour = bioclimate),
             size = 3.6) +
  # etiquetas din├ímicas de especie (ggrepel reparte autom├íticamente para evitar solapes)
  ggrepel::geom_text_repel(data = spe.cent,
                           aes(x = Dim.1, y = Dim.2, label = species_short),
                           size = 3, box.padding = 0.5, point.padding = 0.3,
                           min.segment.length = 0, segment.color = "grey50",
                           max.overlaps = Inf) +
  # flechas variables (mismo escalado que 1a)
  geom_segment(data = var_scaled,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 0.6) +
  ggrepel::geom_text_repel(data = var_scaled,
                           aes(x = Dim.1, y = Dim.2, label = rownames(var)),
                           size = 3, color = "black",
                           box.padding = 0.5, point.padding = 0.5,
                           segment.color = "grey50") +
  scale_colour_manual(values = bioclimate_cols) +
  theme_minimal() +
  labs(x = "Dim 1", y = "Dim 2", colour = "Bioclimate")

## 3.3 Fig 1a + 1b combinada ----
p_final <- p / p1b +
  plot_layout(widths = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("07-img/FAMD_outputs/paper_famd_1a_1b.png", p_final, width = 7, height = 10)

# Panel 1b por separado (misma est├®tica y tama├▒o que 1a)
ggsave("07-img/FAMD_outputs/paper_famd_1b.png", p1b, width = 7, height = 6)

