# _legacy_clustering.R ---- ANALISIS DE AGRUPAMIENTO (ARCHIVADO, COMENTADO)
# -----------------------------------------------------------------------------
# Bloques de las antiguas secciones 4-7 de d.03 (exploracion y clustering):
# silhouette/PERMANOVA por bioclima, envolventes convexas, Hopkins, seleccion de k
# (Elbow/Silhouette/Gap/NbClust), k-means, comparativa k = 2/3/4 y elipse Mahalanobis.
# No se usan en el manuscrito final; se conservan comentados por si un revisor o
# un analisis futuro los requiere. Nada de este archivo se ejecuta.

# ## 4.2 Soporte cuantitativo: silhouette con bioclimate como grupos ----
# # Valor medio >= ~0.25 => separaci├│n razonable en Dim1-Dim3; por debajo,
# # solape fuerte. Es solo apoyo: la decisi├│n final es la inspecci├│n visual (4.3).
# d.bio <- dist(ind.bio |> dplyr::select(Dim.1, Dim.2, Dim.3))
# sil.bio <- cluster::silhouette(as.numeric(ind.bio$bioclimate), d.bio)
# sil.bio.sum <- summary(sil.bio)
# 
# sil.bio.tab <- data.frame(
#   n_bellotas           = nrow(ind.bio),
#   mean_silhouette      = sil.bio.sum$avg.width,
#   sil_Mediterranean    = mean(sil.bio[ind.bio$bioclimate == "Mediterranean", 3]),
#   sil_SubMediterranean = mean(sil.bio[ind.bio$bioclimate == "Sub-Mediterranean", 3]),
#   sil_Temperate        = mean(sil.bio[ind.bio$bioclimate == "Temperate", 3])
# )
# 
# # PERMANOVA opcional (si vegan est├í instalado)
# if (requireNamespace("vegan", quietly = TRUE)) {
#   set.seed(123)
#   permanova <- vegan::adonis2(
#     ind.bio |> dplyr::select(Dim.1, Dim.2, Dim.3) ~ bioclimate,
#     data = ind.bio, permutations = 999, method = "euclidean"
#   )
#   sil.bio.tab$permanova_R2 <- permanova$R2[1]
#   sil.bio.tab$permanova_p  <- permanova$`Pr(>F)`[1]
# }
# 
# write.csv(sil.bio.tab, "00-data/bioclimate_exploration_summary.csv", row.names = FALSE)
# 
# cat("\n--- Figure 1b, Phase 1 ---\n")
# cat(sprintf("Mean silhouette (bioclimate as groups, Dim1-Dim3): %.3f\n",
#             sil.bio.sum$avg.width))
# cat("  -> >= ~0.25 = separation razonable; por debajo = solape fuerte\n")
# 
# ## 4.3 Gr├íficos exploratorios (decisi├│n visual del usuario) ----
# # Envolventes convexas por bioclimate para visualizar la ocupaci├│n espacial
# hull.12 <- ind.bio |>
#   group_by(bioclimate) |>
#   slice(chull(Dim.1, Dim.2))
# 
# hull.13 <- ind.bio |>
#   group_by(bioclimate) |>
#   slice(chull(Dim.1, Dim.3))
# 
# p.bio.12 <- ggplot() +
#   geom_polygon(data = hull.12,
#                aes(x = Dim.1, y = Dim.2, fill = bioclimate, group = bioclimate),
#                alpha = 0.15, colour = NA) +
#   geom_point(data = ind.bio,
#              aes(x = Dim.1, y = Dim.2, colour = bioclimate),
#              alpha = 0.5, size = 1.2) +
#   scale_colour_manual(values = bioclimate_cols) +
#   scale_fill_manual(values = bioclimate_cols) +
#   guides(fill = "none") +
#   theme_minimal() +
#   labs(x = "Dim 1", y = "Dim 2", colour = "Bioclimate")
# 
# p.bio.13 <- ggplot() +
#   geom_polygon(data = hull.13,
#                aes(x = Dim.1, y = Dim.3, fill = bioclimate, group = bioclimate),
#                alpha = 0.15, colour = NA) +
#   geom_point(data = ind.bio,
#              aes(x = Dim.1, y = Dim.3, colour = bioclimate),
#              alpha = 0.5, size = 1.2) +
#   scale_colour_manual(values = bioclimate_cols) +
#   scale_fill_manual(values = bioclimate_cols) +
#   guides(fill = "none") +
#   theme_minimal() +
#   labs(x = "Dim 1", y = "Dim 3", colour = "Bioclimate")
# 
# p.bio <- p.bio.12 + p.bio.13 + plot_layout(guides = "collect")
# 
# ggsave("07-img/FAMD_outputs/04_bioclimate_exploration.png", p.bio, width = 12, height = 5.5)
# 
# # 5. Figure 1b - Phase 2: k-means sobre las tres primeras dimensiones ----
# 
# ## 5.1 Preparar datos ----
# # Se estandarizan las 3 dimensiones para que contribuyan por igual al clustering.
# X.kmeans <- ind.bio |>
#   dplyr::select(Dim.1, Dim.2, Dim.3) |>
#   scale() |>
#   as.matrix()
# rownames(X.kmeans) <- ind.bio$id_bellota
# 
# ## 5.2 Agrupabilidad (Hopkins) ----
# set.seed(123)
# H.hop <- clustertend::hopkins(X.kmeans)$H
# cat(sprintf("\n[Fase 2 - k-means] Hopkins H = %.3f ", H.hop))
# cat("(en clustertend: ~0.5 aleatorio, ->0 clusterable, ->1 uniforme)\n")
# 
# ## 5.3 Decidir k: codo, silhouette, gap, NbClust ----
# set.seed(123)
# wss <- sapply(1:10, function(k) kmeans(X.kmeans, centers = k, nstart = 25)$tot.withinss)
# 
# # Elbow por "kneedle": m├íxima distancia perpendicular a la recta extremo-extremo
# n <- length(wss)
# A <- wss[n] - wss[1]; B <- -(n - 1); C <- (n - 1) * wss[1] - (wss[n] - wss[1])
# perp <- abs(A * (1:n) + B * wss + C) / sqrt(A^2 + B^2)
# perp[c(1, n)] <- -Inf
# elbow_k <- which.max(perp)
# 
# # Silhouette medio k=2..10
# avg.sil <- sapply(2:10, function(k) {
#   km <- kmeans(X.kmeans, centers = k, nstart = 25)
#   mean(cluster::silhouette(km$cluster, dist(X.kmeans))[, 3])
# })
# sil_k <- which.max(avg.sil) + 1
# 
# # Gap statistic (B=500)
# set.seed(123)
# gap.fit <- cluster::clusGap(X.kmeans, FUN = kmeans, K.max = 10, B = 500, nstart = 25)
# gap_k <- cluster::maxSE(gap.fit$Tab[, 3], gap.fit$Tab[, 4], method = "firstSEmax")
# 
# # NbClust: voto mayoritario de ~30 ├¡ndices
# set.seed(123)
# nb.fit <- NbClust::NbClust(X.kmeans, distance = "euclidean", min.nc = 2, max.nc = 10,
#                             method = "kmeans", index = "all")
# nb.votes <- table(nb.fit$Best.nc[1, ])
# nb_k <- as.numeric(names(nb.votes)[which.max(nb.votes)])
# 
# k.methods <- c(elbow = elbow_k, silhouette = sil_k, gap = gap_k, NbClust = nb_k)
# cat("[Fase 2] k propuesto por m├®todo:\n")
# print(k.methods)
# 
# # Consenso: moda entre m├®todos; si empate, el que tenga mejor silhouette media
# votes <- table(k.methods)
# cand <- as.numeric(names(votes)[votes == max(votes)])
# if (length(cand) > 1) k.final <- cand[which.max(avg.sil[cand - 1])] else k.final <- cand
# cat(sprintf("[Fase 2] k final elegido: %d\n", k.final))
# 
# ## 5.4 K-means final ----
# set.seed(123)
# km.final <- kmeans(X.kmeans, centers = k.final, nstart = 50, iter.max = 100)
# 
# ind.bio$cluster <- as.factor(km.final$cluster)
# 
# # Silhouette final
# sil.final <- cluster::silhouette(km.final$cluster, dist(X.kmeans))
# cat(sprintf("[Fase 2] Silhouette medio del k-means final: %.3f\n",
#             summary(sil.final)$avg.width))
# 
# ## 5.5 Guardar resultados ----
# write.csv(ind.bio, "00-data/famd_clusters.csv", row.names = FALSE)
# 
# comp.especie <- ind.bio |>
#   count(cluster, species) |>
#   pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
# write.csv(comp.especie, "00-data/cluster_species_composition.csv", row.names = FALSE)
# 
# comp.bioclima <- ind.bio |>
#   count(cluster, bioclimate) |>
#   pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
# write.csv(comp.bioclima, "00-data/cluster_bioclimate_composition.csv", row.names = FALSE)
# 
# k.table <- data.frame(
#   method = names(k.methods),
#   k_recommended = as.integer(k.methods),
#   k_final = rep(k.final, length(k.methods)),
#   mean_silhouette_final = rep(summary(sil.final)$avg.width, length(k.methods)),
#   hopkins_H = rep(H.hop, length(k.methods))
# )
# write.csv(k.table, "00-data/kmeans_k_selection.csv", row.names = FALSE)
# 
# ## 5.6 Figuras de decisi├│n (codo / silhouette / gap) ----
# p.wss <- fviz_nbclust(X.kmeans, kmeans, method = "wss", k.max = 10) +
#   theme_minimal() + labs(subtitle = "Elbow (WSS)")
# p.sil <- fviz_nbclust(X.kmeans, kmeans, method = "silhouette", k.max = 10) +
#   theme_minimal() + labs(subtitle = "Silhouette medio")
# p.gap <- fviz_nbclust(X.kmeans, kmeans, method = "gap_stat", k.max = 10, nboot = 500) +
#   theme_minimal() + labs(subtitle = "Gap statistic")
# 
# p.k <- (p.wss | p.sil | p.gap)
# ggsave("07-img/FAMD_outputs/05_optimal_k.png", p.k, width = 15, height = 5)
# 
# # 6. Figure 1b - comparativa visual de agrupaciones k = 2, 3 y 4 ----
# # Filas: k = 2, 3, 4 | Columna izquierda: Dim.1 vs Dim.2 | Columna derecha: Dim.1 vs Dim.3
# grid.cols <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#E69F00", "4" = "#009E73")
# 
# p.grid <- list()
# for (kk in 2:4) {
#   set.seed(123)
#   km.k <- kmeans(X.kmeans, centers = kk, nstart = 50, iter.max = 100)
# 
#   tmp <- ind.bio |>
#     mutate(cluster = as.factor(km.k$cluster))
# 
#   cen <- tmp |>
#     group_by(cluster) |>
#     summarise(Dim.1 = mean(Dim.1),
#               Dim.2 = mean(Dim.2),
#               Dim.3 = mean(Dim.3),
#               .groups = "drop")
# 
#   p12 <- ggplot(tmp, aes(x = Dim.1, y = Dim.2, colour = cluster)) +
#     geom_point(alpha = 0.45, size = 1.3) +
#     stat_ellipse(aes(fill = cluster, colour = cluster), geom = "polygon",
#                  level = 0.65, alpha = 0.10, linewidth = 0.3) +
#     geom_point(data = cen, aes(x = Dim.1, y = Dim.2, colour = cluster),
#                shape = 4, size = 3, stroke = 1.1) +
#     scale_colour_manual(values = grid.cols[1:kk]) +
#     scale_fill_manual(values = grid.cols[1:kk]) +
#     guides(fill = "none") +
#     theme_minimal() +
#     labs(x = "Dim 1", y = "Dim 2", colour = paste0("Cluster (k = ", kk, ")"))
# 
#   p13 <- ggplot(tmp, aes(x = Dim.1, y = Dim.3, colour = cluster)) +
#     geom_point(alpha = 0.45, size = 1.3) +
#     stat_ellipse(aes(fill = cluster, colour = cluster), geom = "polygon",
#                  level = 0.65, alpha = 0.10, linewidth = 0.3) +
#     geom_point(data = cen, aes(x = Dim.1, y = Dim.3, colour = cluster),
#                shape = 4, size = 3, stroke = 1.1) +
#     scale_colour_manual(values = grid.cols[1:kk]) +
#     scale_fill_manual(values = grid.cols[1:kk]) +
#     guides(fill = "none") +
#     theme_minimal() +
#     labs(x = "Dim 1", y = "Dim 3", colour = paste0("Cluster (k = ", kk, ")"))
# 
#   p.grid[[kk]] <- p12 + p13 + plot_layout(guides = "collect")
# }
# 
# p.compare <- wrap_plots(p.grid[2:4], ncol = 1) +
#   plot_annotation(caption = "Elipses al 65% central; X = centroide del cluster")
# 
# ggsave("07-img/FAMD_outputs/06_k2_k3_k4_grid.png", p.compare, width = 11, height = 13)
# 
# # 7. Figure 1b (definitiva) - partici├│n k = 3 y panel 1a + 1b ----
# ## 7.1 Partici├│n final k = 3 ----
# # k = 3 elegido por inspecci├│n visual del usuario (compromiso Elbow/Silhouette=4 vs Gap/NbClust=1-2)
# set.seed(123)
# km.final3 <- kmeans(X.kmeans, centers = 3, nstart = 50, iter.max = 100)
# 
# final.cols <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#E69F00")
# 
# sil.final3 <- cluster::silhouette(km.final3$cluster, dist(X.kmeans))
# cat(sprintf("[Panel 1b] Silhouette medio k=3: %.3f\n", summary(sil.final3)$avg.width))
# 
# ind.final <- ind.bio |>
#   mutate(cluster = factor(km.final3$cluster, levels = 1:3))
# 
# # Guardar la partici├│n definitiva (sobrescribe la provisional de k=4 de la secci├│n 5)
# write.csv(ind.final, "00-data/famd_clusters.csv", row.names = FALSE)
# 
# comp.especie3 <- ind.final |>
#   count(cluster, species) |>
#   pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
# write.csv(comp.especie3, "00-data/cluster_species_composition.csv", row.names = FALSE)
# cat("\n[Panel 1b] Composici├│n por especie (k=3):\n"); print(comp.especie3)
# 
# comp.bioclima3 <- ind.final |>
#   count(cluster, bioclimate) |>
#   pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
# write.csv(comp.bioclima3, "00-data/cluster_bioclimate_composition.csv", row.names = FALSE)
# 
# ## 7.2 Elipse 65% central por Mahalanobis (2D) ----
# # Recorta las distancias de Mahalanobis al centroide en el cuantil 65%:
# # representa el 65% m├ís central de cada cluster, reduciendo el solape visual.
# ellipse_65 <- function(x, y, level = 0.65) {
#   pts <- cbind(x, y)
#   mu <- colMeans(pts)
#   S <- cov(pts)
#   d2 <- mahalanobis(pts, mu, S)
#   q <- quantile(d2, probs = level, names = FALSE)
#   ev <- eigen(S)
#   A <- ev$vectors %*% diag(sqrt(pmax(ev$values * q, 0)))
#   theta <- seq(0, 2 * pi, length.out = 150)
#   circ <- cbind(cos(theta), sin(theta))
#   ell <- t(mu + A %*% t(circ))
#   data.frame(x = ell[, 1], y = ell[, 2])
# }
# 
# ell.12 <- ind.final |>
#   group_by(cluster) |>
#   group_modify(~ ellipse_65(.x$Dim.1, .x$Dim.2))
# 
