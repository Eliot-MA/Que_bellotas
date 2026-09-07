library(tidyverse)

# Si no tienes desiccation_traits_long.csv
# ejecuta primero VVV
# source("01-scripts/d.01-load_desiccation_exp.R")
# rm(list = ls())

df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")

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

## 1.1b. Correlograma de variables incluidas en el FAMD ----
library(corrplot)

df.corr <- df |>
  dplyr::select(dry_weight, volume_cm3, SVR, SPM_g_cm2, 
                SCR, SSR, pericarp_rupture) |>
  mutate(pericarp_rupture = as.numeric(as.character(pericarp_rupture)))

mat.corr <- cor(df.corr, method = "pearson", use = "pairwise.complete.obs")

rownames(mat.corr) <- colnames(mat.corr) <- c("mass", "volume", "SVR", 
                                                "SPM", "SCR", "SSR", 
                                                "pericarp rupture")

write.csv2(mat.corr, "00-data/correlation_matrix_traits.csv")

dir.create("07-img/FAMD_outputs", showWarnings = FALSE)

png("07-img/FAMD_outputs/00_correlogram_traits.png", width = 800, height = 700)
corrplot(mat.corr, 
         method = "color", 
         type = "upper", 
         addCoef.col = "black",
         tl.col = "black", 
         tl.srt = 45,
         number.cex = 0.8,
         mar = c(0, 0, 2, 0))
dev.off()
cat("Correlograma guardado en 07-img/FAMD_outputs/00_correlogram_traits.png\n")

## 1.2. Calculate and save famd ----

library(FactoMineR)  # Para FAMD
library(factoextra)


# Exclude id and species variables
df.famd <- df |> dplyr::select(-id_bellota, -species, -provenance, -prov_code)

# Calculate famd
famd.traits <- FAMD(df.famd, graph = FALSE)

df <- cbind(df, famd.traits$ind$coord[,1:5])

write.csv(x = df, "00-data/famd_ind_coord.csv")


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
  dplyr::select(variable, dimension, contrib, coord)

# Paso 2: codificar el signo
# Aquí tienes varias opciones. La más limpia es añadir el signo como string:

tabla_contrib <- tabla_contrib %>%
  mutate(
    signo = ifelse(coord >= 0, "+", "−"),
    contrib_signo = paste0(round(contrib, 1), signo)
  )

# Paso 3: pasar a formato ancho
tabla_wide <- tabla_contrib %>%
  dplyr::select(variable, dimension, contrib_signo) %>%
  pivot_wider(names_from = dimension, values_from = contrib_signo)

write.csv2(x = tabla_wide, file = "00-data/paper_famd.csv")

# 3. Create FAMD figure ----

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
# colores okabe-ito para los hexágonos
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
  segment.color = "grey50"  # Dibuja una línea hacia el punto original
) +
  
# -------------------
# estética
# -------------------
theme_minimal() +
  
  labs(
    x = "Dim 1",
    y = "Dim 2",
    fill = "Dim 3 bin"
  )

ggsave("07-img/FAMD_outputs/paper_famd_1a.png", p, width = 7, height = 6)

# 4. Figure 1b - Phase 1: exploración del agrupamiento por bioclima ----

## 4.1 Mapear el bioclimate de cada especie ----
bioclimate <- tribble(
  ~species,            ~bioclimate,
  "Quercus coccifera", "Mediterranean",
  "Quercus ilex",      "Mediterranean",
  "Quercus suber",     "Mediterranean",
  "Quercus faginea",   "Sub-Mediterranean",
  "Quercus pyrenaica", "Sub-Mediterranean",
  "Quercus pubescens", "Sub-Mediterranean",
  "Quercus petraea",   "Temperate",
  "Quercus robur",     "Temperate"
)

bioclimate_cols <- c(
  "Mediterranean"     = "#D55E00",  # Okabe-Ito vermilion
  "Sub-Mediterranean" = "#E69F00",  # Okabe-Ito orange
  "Temperate"         = "#0072B2"   # Okabe-Ito blue
)

ind.bio <- df |>
  dplyr::select(id_bellota, species, Dim.1, Dim.2, Dim.3) |>
  left_join(bioclimate, by = "species") |>
  mutate(bioclimate = factor(bioclimate,
                             levels = c("Mediterranean",
                                        "Sub-Mediterranean",
                                        "Temperate")))

## 4.2 Soporte cuantitativo: silhouette con bioclimate como grupos ----
# Valor medio >= ~0.25 => separación razonable en Dim1-Dim3; por debajo,
# solape fuerte. Es solo apoyo: la decisión final es la inspección visual (4.3).
d.bio <- dist(ind.bio |> dplyr::select(Dim.1, Dim.2, Dim.3))
sil.bio <- cluster::silhouette(as.numeric(ind.bio$bioclimate), d.bio)
sil.bio.sum <- summary(sil.bio)

sil.bio.tab <- data.frame(
  n_bellotas           = nrow(ind.bio),
  mean_silhouette      = sil.bio.sum$avg.width,
  sil_Mediterranean    = mean(sil.bio[ind.bio$bioclimate == "Mediterranean", 3]),
  sil_SubMediterranean = mean(sil.bio[ind.bio$bioclimate == "Sub-Mediterranean", 3]),
  sil_Temperate        = mean(sil.bio[ind.bio$bioclimate == "Temperate", 3])
)

# PERMANOVA opcional (si vegan está instalado)
if (requireNamespace("vegan", quietly = TRUE)) {
  set.seed(123)
  permanova <- vegan::adonis2(
    ind.bio |> dplyr::select(Dim.1, Dim.2, Dim.3) ~ bioclimate,
    data = ind.bio, permutations = 999, method = "euclidean"
  )
  sil.bio.tab$permanova_R2 <- permanova$R2[1]
  sil.bio.tab$permanova_p  <- permanova$`Pr(>F)`[1]
}

write.csv(sil.bio.tab, "00-data/bioclimate_exploration_summary.csv", row.names = FALSE)

cat("\n--- Figure 1b, Phase 1 ---\n")
cat(sprintf("Mean silhouette (bioclimate as groups, Dim1-Dim3): %.3f\n",
            sil.bio.sum$avg.width))
cat("  -> >= ~0.25 = separation razonable; por debajo = solape fuerte\n")

## 4.3 Gráficos exploratorios (decisión visual del usuario) ----
# Envolventes convexas por bioclimate para visualizar la ocupación espacial
hull.12 <- ind.bio |>
  group_by(bioclimate) |>
  slice(chull(Dim.1, Dim.2))

hull.13 <- ind.bio |>
  group_by(bioclimate) |>
  slice(chull(Dim.1, Dim.3))

p.bio.12 <- ggplot() +
  geom_polygon(data = hull.12,
               aes(x = Dim.1, y = Dim.2, fill = bioclimate, group = bioclimate),
               alpha = 0.15, colour = NA) +
  geom_point(data = ind.bio,
             aes(x = Dim.1, y = Dim.2, colour = bioclimate),
             alpha = 0.5, size = 1.2) +
  scale_colour_manual(values = bioclimate_cols) +
  scale_fill_manual(values = bioclimate_cols) +
  guides(fill = "none") +
  theme_minimal() +
  labs(x = "Dim 1", y = "Dim 2", colour = "Bioclimate")

p.bio.13 <- ggplot() +
  geom_polygon(data = hull.13,
               aes(x = Dim.1, y = Dim.3, fill = bioclimate, group = bioclimate),
               alpha = 0.15, colour = NA) +
  geom_point(data = ind.bio,
             aes(x = Dim.1, y = Dim.3, colour = bioclimate),
             alpha = 0.5, size = 1.2) +
  scale_colour_manual(values = bioclimate_cols) +
  scale_fill_manual(values = bioclimate_cols) +
  guides(fill = "none") +
  theme_minimal() +
  labs(x = "Dim 1", y = "Dim 3", colour = "Bioclimate")

p.bio <- p.bio.12 + p.bio.13 + plot_layout(guides = "collect")

ggsave("07-img/FAMD_outputs/04_bioclimate_exploration.png", p.bio, width = 12, height = 5.5)

# 5. Figure 1b - Phase 2: k-means sobre las tres primeras dimensiones ----

## 5.1 Preparar datos ----
# Se estandarizan las 3 dimensiones para que contribuyan por igual al clustering.
X.kmeans <- ind.bio |>
  dplyr::select(Dim.1, Dim.2, Dim.3) |>
  scale() |>
  as.matrix()
rownames(X.kmeans) <- ind.bio$id_bellota

## 5.2 Agrupabilidad (Hopkins) ----
set.seed(123)
H.hop <- clustertend::hopkins(X.kmeans)$H
cat(sprintf("\n[Fase 2 - k-means] Hopkins H = %.3f ", H.hop))
cat("(en clustertend: ~0.5 aleatorio, ->0 clusterable, ->1 uniforme)\n")

## 5.3 Decidir k: codo, silhouette, gap, NbClust ----
set.seed(123)
wss <- sapply(1:10, function(k) kmeans(X.kmeans, centers = k, nstart = 25)$tot.withinss)

# Elbow por "kneedle": máxima distancia perpendicular a la recta extremo-extremo
n <- length(wss)
A <- wss[n] - wss[1]; B <- -(n - 1); C <- (n - 1) * wss[1] - (wss[n] - wss[1])
perp <- abs(A * (1:n) + B * wss + C) / sqrt(A^2 + B^2)
perp[c(1, n)] <- -Inf
elbow_k <- which.max(perp)

# Silhouette medio k=2..10
avg.sil <- sapply(2:10, function(k) {
  km <- kmeans(X.kmeans, centers = k, nstart = 25)
  mean(cluster::silhouette(km$cluster, dist(X.kmeans))[, 3])
})
sil_k <- which.max(avg.sil) + 1

# Gap statistic (B=500)
set.seed(123)
gap.fit <- cluster::clusGap(X.kmeans, FUN = kmeans, K.max = 10, B = 500, nstart = 25)
gap_k <- cluster::maxSE(gap.fit$Tab[, 3], gap.fit$Tab[, 4], method = "firstSEmax")

# NbClust: voto mayoritario de ~30 índices
set.seed(123)
nb.fit <- NbClust::NbClust(X.kmeans, distance = "euclidean", min.nc = 2, max.nc = 10,
                            method = "kmeans", index = "all")
nb.votes <- table(nb.fit$Best.nc[1, ])
nb_k <- as.numeric(names(nb.votes)[which.max(nb.votes)])

k.methods <- c(elbow = elbow_k, silhouette = sil_k, gap = gap_k, NbClust = nb_k)
cat("[Fase 2] k propuesto por método:\n")
print(k.methods)

# Consenso: moda entre métodos; si empate, el que tenga mejor silhouette media
votes <- table(k.methods)
cand <- as.numeric(names(votes)[votes == max(votes)])
if (length(cand) > 1) k.final <- cand[which.max(avg.sil[cand - 1])] else k.final <- cand
cat(sprintf("[Fase 2] k final elegido: %d\n", k.final))

## 5.4 K-means final ----
set.seed(123)
km.final <- kmeans(X.kmeans, centers = k.final, nstart = 50, iter.max = 100)

ind.bio$cluster <- as.factor(km.final$cluster)

# Silhouette final
sil.final <- cluster::silhouette(km.final$cluster, dist(X.kmeans))
cat(sprintf("[Fase 2] Silhouette medio del k-means final: %.3f\n",
            summary(sil.final)$avg.width))

## 5.5 Guardar resultados ----
write.csv(ind.bio, "00-data/famd_clusters.csv", row.names = FALSE)

comp.especie <- ind.bio |>
  count(cluster, species) |>
  pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
write.csv(comp.especie, "00-data/cluster_species_composition.csv", row.names = FALSE)

comp.bioclima <- ind.bio |>
  count(cluster, bioclimate) |>
  pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
write.csv(comp.bioclima, "00-data/cluster_bioclimate_composition.csv", row.names = FALSE)

k.table <- data.frame(
  method = names(k.methods),
  k_recommended = as.integer(k.methods),
  k_final = rep(k.final, length(k.methods)),
  mean_silhouette_final = rep(summary(sil.final)$avg.width, length(k.methods)),
  hopkins_H = rep(H.hop, length(k.methods))
)
write.csv(k.table, "00-data/kmeans_k_selection.csv", row.names = FALSE)

## 5.6 Figuras de decisión (codo / silhouette / gap) ----
p.wss <- fviz_nbclust(X.kmeans, kmeans, method = "wss", k.max = 10) +
  theme_minimal() + labs(subtitle = "Elbow (WSS)")
p.sil <- fviz_nbclust(X.kmeans, kmeans, method = "silhouette", k.max = 10) +
  theme_minimal() + labs(subtitle = "Silhouette medio")
p.gap <- fviz_nbclust(X.kmeans, kmeans, method = "gap_stat", k.max = 10, nboot = 500) +
  theme_minimal() + labs(subtitle = "Gap statistic")

p.k <- (p.wss | p.sil | p.gap)
ggsave("07-img/FAMD_outputs/05_optimal_k.png", p.k, width = 15, height = 5)

# 6. Figure 1b - comparativa visual de agrupaciones k = 2, 3 y 4 ----
# Filas: k = 2, 3, 4 | Columna izquierda: Dim.1 vs Dim.2 | Columna derecha: Dim.1 vs Dim.3
grid.cols <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#E69F00", "4" = "#009E73")

p.grid <- list()
for (kk in 2:4) {
  set.seed(123)
  km.k <- kmeans(X.kmeans, centers = kk, nstart = 50, iter.max = 100)

  tmp <- ind.bio |>
    mutate(cluster = as.factor(km.k$cluster))

  cen <- tmp |>
    group_by(cluster) |>
    summarise(Dim.1 = mean(Dim.1),
              Dim.2 = mean(Dim.2),
              Dim.3 = mean(Dim.3),
              .groups = "drop")

  p12 <- ggplot(tmp, aes(x = Dim.1, y = Dim.2, colour = cluster)) +
    geom_point(alpha = 0.45, size = 1.3) +
    stat_ellipse(aes(fill = cluster, colour = cluster), geom = "polygon",
                 level = 0.65, alpha = 0.10, linewidth = 0.3) +
    geom_point(data = cen, aes(x = Dim.1, y = Dim.2, colour = cluster),
               shape = 4, size = 3, stroke = 1.1) +
    scale_colour_manual(values = grid.cols[1:kk]) +
    scale_fill_manual(values = grid.cols[1:kk]) +
    guides(fill = "none") +
    theme_minimal() +
    labs(x = "Dim 1", y = "Dim 2", colour = paste0("Cluster (k = ", kk, ")"))

  p13 <- ggplot(tmp, aes(x = Dim.1, y = Dim.3, colour = cluster)) +
    geom_point(alpha = 0.45, size = 1.3) +
    stat_ellipse(aes(fill = cluster, colour = cluster), geom = "polygon",
                 level = 0.65, alpha = 0.10, linewidth = 0.3) +
    geom_point(data = cen, aes(x = Dim.1, y = Dim.3, colour = cluster),
               shape = 4, size = 3, stroke = 1.1) +
    scale_colour_manual(values = grid.cols[1:kk]) +
    scale_fill_manual(values = grid.cols[1:kk]) +
    guides(fill = "none") +
    theme_minimal() +
    labs(x = "Dim 1", y = "Dim 3", colour = paste0("Cluster (k = ", kk, ")"))

  p.grid[[kk]] <- p12 + p13 + plot_layout(guides = "collect")
}

p.compare <- wrap_plots(p.grid[2:4], ncol = 1) +
  plot_annotation(caption = "Elipses al 65% central; X = centroide del cluster")

ggsave("07-img/FAMD_outputs/06_k2_k3_k4_grid.png", p.compare, width = 11, height = 13)

# 7. Figure 1b (definitiva) - partición k = 3 y panel 1a + 1b ----
## 7.1 Partición final k = 3 ----
# k = 3 elegido por inspección visual del usuario (compromiso Elbow/Silhouette=4 vs Gap/NbClust=1-2)
set.seed(123)
km.final3 <- kmeans(X.kmeans, centers = 3, nstart = 50, iter.max = 100)

final.cols <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#E69F00")

sil.final3 <- cluster::silhouette(km.final3$cluster, dist(X.kmeans))
cat(sprintf("[Panel 1b] Silhouette medio k=3: %.3f\n", summary(sil.final3)$avg.width))

ind.final <- ind.bio |>
  mutate(cluster = factor(km.final3$cluster, levels = 1:3))

# Guardar la partición definitiva (sobrescribe la provisional de k=4 de la sección 5)
write.csv(ind.final, "00-data/famd_clusters.csv", row.names = FALSE)

comp.especie3 <- ind.final |>
  count(cluster, species) |>
  pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
write.csv(comp.especie3, "00-data/cluster_species_composition.csv", row.names = FALSE)
cat("\n[Panel 1b] Composición por especie (k=3):\n"); print(comp.especie3)

comp.bioclima3 <- ind.final |>
  count(cluster, bioclimate) |>
  pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
write.csv(comp.bioclima3, "00-data/cluster_bioclimate_composition.csv", row.names = FALSE)

## 7.2 Elipse 65% central por Mahalanobis (2D) ----
# Recorta las distancias de Mahalanobis al centroide en el cuantil 65%:
# representa el 65% más central de cada cluster, reduciendo el solape visual.
ellipse_65 <- function(x, y, level = 0.65) {
  pts <- cbind(x, y)
  mu <- colMeans(pts)
  S <- cov(pts)
  d2 <- mahalanobis(pts, mu, S)
  q <- quantile(d2, probs = level, names = FALSE)
  ev <- eigen(S)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values * q, 0)))
  theta <- seq(0, 2 * pi, length.out = 150)
  circ <- cbind(cos(theta), sin(theta))
  ell <- t(mu + A %*% t(circ))
  data.frame(x = ell[, 1], y = ell[, 2])
}

ell.12 <- ind.final |>
  group_by(cluster) |>
  group_modify(~ ellipse_65(.x$Dim.1, .x$Dim.2))

## 7.3 Panel 1b ----
# Centroide de cada especie en el plano del biplot (media de Dim.1 y Dim.2)
spe.cent <- ind.final |>
  group_by(species) |>
  summarise(Dim.1 = mean(Dim.1),
            Dim.2 = mean(Dim.2),
            bioclimate = unique(bioclimate),
            n_spots = n(),
            .groups = "drop") |>
  mutate(species_short = gsub("Quercus ", "Q. ", species))

p1b <- ggplot() +
  # individuos coloreados por bioclima
  geom_point(data = ind.final,
             aes(x = Dim.1, y = Dim.2, colour = bioclimate),
             alpha = 0.35, size = 1.2) +
  # centroides por especie, coloreados por bioclima
  geom_point(data = spe.cent,
             aes(x = Dim.1, y = Dim.2, colour = bioclimate),
             size = 3.6) +
  # etiquetas dinámicas de especie (ggrepel reparte automáticamente para evitar solapes)
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

## 7.4 Fig 1a + 1b combinada ----
p_final <- p + p1b +
  plot_layout(widths = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("07-img/FAMD_outputs/paper_famd_1a_1b.png", p_final, width = 15, height = 7)

# Panel 1b por separado (misma estética y tamaño que 1a)
ggsave("07-img/FAMD_outputs/paper_famd_1b.png", p1b, width = 7, height = 6)

