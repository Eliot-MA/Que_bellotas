# d.05.5.predicted_slopes_figures.R
# Curvas de desecacion predichas por especie y pendientes estimadas,
# para valores bajos (q05) y altos (q95) de cada dimension del FAMD,
# en las fases PRE (t<94h) y POST (t>94h).
#
# Fuente de resultados: modelos brms guardados en 00-data/phylo
#   m_het_2_pre.rds / m_het_2_post.rds  (filogenia + interacciones triples)
#   m_het_1_pre.rds / m_het_1_post.rds  (filogenia + especie libre)
# Sin refit: se usan los draws almacenados.
#
# Salidas:
#   07-img/phylo_bayes_curves_M2_{pre,post}.png   curvas predichas por fase
#   07-img/phylo_bayes_curves_M2.png              grid pre+post (3 dims)
#   07-img/phylo_bayes_slopes_M2.png              pendientes estimadas por especie
#   00-data/tablas_resumen/phylo_slopes_by_species.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})

OUTDIR_IMG <- "07-img"
OUTDIR_CSV <- "00-data/tablas_resumen"
dir.create(OUTDIR_IMG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_CSV, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Datos (misma preparacion que d.05.z.pruebas_filo.R)
# ============================================================
# Procedencia IL3 excluida (ver d.05.0.model_traits.R): muy pocas obs en
# fase PRE (60 vs 150); se elimina de ambas fases para comparaciones
# pre-post validas.
PROCEDENCIAS_EXCLUIDAS <- "IL3"
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd     <- read.csv("00-data/famd_ind_coord.csv")
df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |>
  left_join(y = df.famd, by = "id_bellota") |>
  dplyr::filter(!codigo %in% PROCEDENCIAS_EXCLUIDAS) |>
  tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    time_s     = as.vector(scale(time)),
    species    = factor(species),
    provenance = factor(provenance),
    id_bellota = factor(id_bellota)
  )

TIME_M <- mean(df$time)
TIME_S <- sd(df$time)
t94    <- as.vector((94 - TIME_M) / TIME_S)
df.pre <- df |> filter(time_s < t94)
df.post<- df |> filter(time_s > t94)

species  <- levels(df$species)
spp_dots <- gsub("\\.", "", species)     # etiquetas de coeficientes fijos (sin puntos)
spp_lab  <- gsub("^Quercus\\.?", "Q. ", species)
names(spp_lab) <- species

dims     <- c("Dim.1", "Dim.2", "Dim.3")
.q      <- function(d, p) unname(quantile(df[[d]], p))
qlo      <- sapply(dims, .q, p = 0.05)
qhi      <- sapply(dims, .q, p = 0.95)
cat("Quantiles (p05/p95):\n")
print(data.frame(dim = dims, p05 = qlo, p95 = qhi))

# ============================================================
# 1. Predictor lineal y gradiente a partir de draws
# ============================================================
# lookup de columnas por especie basado en las columnas reales del modelo:
# brms sanea las etiquetas (espacio->'.', o elimina no-alfanumericos), por eso
# se compara normalizando a [a-zA-Z0-9].
sp_col <- function(nms, prefix, suffix = "") {
  hits <- nms[startsWith(nms, prefix) & endsWith(nms, suffix)]
  out <- setNames(rep(NA_character_, length(species)), species)
  if (length(hits) == 0) return(out)
  lab  <- sub(prefix, "", hits, fixed = TRUE)
  if (nzchar(suffix)) lab <- sub(suffix, "", lab, fixed = TRUE)
  keys <- gsub("[^[:alnum:]]", "", lab)
  nk   <- gsub("[^[:alnum:]]", "", species)
  m    <- match(nk, keys)
  ok   <- !is.na(m)
  out[ok] <- hits[m[ok]]
  out
}

predict_epred <- function(fit, variant, grid) {
  dd  <- posterior::as_draws_matrix(fit)
  nms <- colnames(dd)
  n   <- nrow(dd)
  get <- function(nm) if (nm %in% nms) dd[, nm] else rep(0, n)

  tD1 <- sp_col(nms, "b_time_s:Dim.1:species")
  tD2 <- sp_col(nms, "b_time_s:Dim.2:species")
  tD3 <- sp_col(nms, "b_time_s:Dim.3:species")
  rS  <- sp_col(nms, "r_species[", ",time_s]")
  rP  <- sp_col(nms, "r_phylo_species[", ",time_s]")

  b0 <- get("b_Intercept"); b_ts <- get("b_time_s")
  bD1 <- get("b_Dim.1"); bD2 <- get("b_Dim.2"); bD3 <- get("b_Dim.3")
  bt1 <- get("b_time_s:Dim.1"); bt2 <- get("b_time_s:Dim.2"); bt3 <- get("b_time_s:Dim.3")

  ts  <- grid$time_s; D1 <- grid$D1; D2 <- grid$D2; D3 <- grid$D3
  si  <- match(as.character(grid$species), species)

  M <- matrix(0.0, n, nrow(grid))
  for (j in seq_len(nrow(grid))) {
    mu <- b0 + b_ts * ts[j] +
      bD1 * D1[j] + bD2 * D2[j] + bD3 * D3[j] +
      bt1 * (ts[j] * D1[j]) + bt2 * (ts[j] * D2[j]) + bt3 * (ts[j] * D3[j])
    if (variant == "M2") {
      c1 <- tD1[si[j]]
      if (!is.na(c1)) mu <- mu + dd[, c1] * (ts[j] * D1[j])
      c2 <- tD2[si[j]]
      if (!is.na(c2)) mu <- mu + dd[, c2] * (ts[j] * D2[j])
      c3 <- tD3[si[j]]
      if (!is.na(c3)) mu <- mu + dd[, c3] * (ts[j] * D3[j])
    }
    if (variant == "M1") {
      c1 <- rS[si[j]]
      if (!is.na(c1)) mu <- mu + dd[, c1] * ts[j]
      c2 <- rP[si[j]]
      if (!is.na(c2)) mu <- mu + dd[, c2] * ts[j]
    }
    M[, j] <- mu
  }
  M
}

grad <- function(fit, variant, grid) {
  dd  <- posterior::as_draws_matrix(fit)
  nms <- colnames(dd)
  n   <- nrow(dd)
  get <- function(nm) if (nm %in% nms) dd[, nm] else rep(0, n)

  tD1 <- sp_col(nms, "b_time_s:Dim.1:species")
  tD2 <- sp_col(nms, "b_time_s:Dim.2:species")
  tD3 <- sp_col(nms, "b_time_s:Dim.3:species")
  rS  <- sp_col(nms, "r_species[", ",time_s]")
  rP  <- sp_col(nms, "r_phylo_species[", ",time_s]")

  b_ts <- get("b_time_s")
  bt1 <- get("b_time_s:Dim.1"); bt2 <- get("b_time_s:Dim.2"); bt3 <- get("b_time_s:Dim.3")

  D1 <- grid$D1; D2 <- grid$D2; D3 <- grid$D3
  si  <- match(as.character(grid$species), species)

  G <- matrix(0.0, n, nrow(grid))
  for (j in seq_len(nrow(grid))) {
    g <- b_ts + bt1 * D1[j] + bt2 * D2[j] + bt3 * D3[j]
    if (variant == "M2") {
      c1 <- tD1[si[j]]
      if (!is.na(c1)) g <- g + dd[, c1] * D1[j]
      c2 <- tD2[si[j]]
      if (!is.na(c2)) g <- g + dd[, c2] * D2[j]
      c3 <- tD3[si[j]]
      if (!is.na(c3)) g <- g + dd[, c3] * D3[j]
    }
    if (variant == "M1") {
      c1 <- rS[si[j]]
      if (!is.na(c1)) g <- g + dd[, c1]
      c2 <- rP[si[j]]
      if (!is.na(c2)) g <- g + dd[, c2]
    }
    G[, j] <- g
  }
  G
}

summ <- function(M, grid, probs = c(0.10, 0.50, 0.90)) {
  q <- t(apply(M, 2, quantile, probs = probs))
  grid |> mutate(lwr = q[, 1], med = q[, 2], upr = q[, 3])
}

# ============================================================
# 2. Rejilla de escenarios (especie x dim x nivel q05/q95)
# ============================================================
scenarios <- function(ts_range, n_ts = 30) {
  ts_seq <- seq(ts_range[1], ts_range[2], length.out = n_ts)
  rows <- bind_rows(lapply(dims, function(d) {
    bind_rows(
      tibble(dim = d, level = "q05", D1 = if (d == "Dim.1") qlo[d] else 0,
             D2 = if (d == "Dim.2") qlo[d] else 0, D3 = if (d == "Dim.3") qlo[d] else 0),
      tibble(dim = d, level = "q95", D1 = if (d == "Dim.1") qhi[d] else 0,
             D2 = if (d == "Dim.2") qhi[d] else 0, D3 = if (d == "Dim.3") qhi[d] else 0))
  }))
  expand_grid(species = species, time_s = ts_seq, rows) |>
    mutate(hr = TIME_M + time_s * TIME_S)
}

sc_pre  <- scenarios(range(df.pre$time_s))
sc_post <- scenarios(range(df.post$time_s))

# ============================================================
# 3. Curvas y pendientes (modelo M_het_2)
# ============================================================
fit_pre  <- readRDS(file.path("00-data/phylo", "m_het_2_pre.rds"))
fit_post <- readRDS(file.path("00-data/phylo", "m_het_2_post.rds"))

curves <- bind_rows(
  summ(predict_epred(fit_pre,  "M2", sc_pre),  sc_pre)  |> mutate(phase = "PRE"),
  summ(predict_epred(fit_post, "M2", sc_post), sc_post) |> mutate(phase = "POST")
)

sc_slope <- bind_rows(sc_pre, sc_post) |>
  distinct(species, D1, D2, D3, dim, level)
slopes <- bind_rows(
  summ(grad(fit_pre,  "M2", sc_slope), sc_slope) |> mutate(phase = "PRE"),
  summ(grad(fit_post, "M2", sc_slope), sc_slope) |> mutate(phase = "POST")
) |>
  mutate(slope_hr = med / TIME_S, lwr_hr = lwr / TIME_S, upr_hr = upr / TIME_S) |>
  mutate(spp = unname(spp_lab[species]))

cat("curvas:", nrow(curves), "| pendientes:", nrow(slopes), "\n")

# ============================================================
# 4. Figuras
# ============================================================
lab_dim  <- c(Dim.1 = "Dim.1 (tamano)",
              Dim.2 = "Dim.2 (pericarpo)",
              Dim.3 = "Dim.3 (cicatriz)")
cols_spp <- setNames(RColorBrewer::brewer.pal(8, "Dark2"), unname(spp_lab))

curves_plot <- curves |>
  mutate(level = ifelse(level == "q05", "p05", "p95"),
         spp = unname(spp_lab[species])) |>
  ggplot(aes(hr, med, color = spp, fill = spp, linetype = level)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.55) +
  facet_grid(phase ~ dim, scales = "free_x", labeller = labeller(dim = lab_dim)) +
  labs(x = "Tiempo (h)", y = "Humedad predicha (%)",
       color = "Especie", fill = "Especie", linetype = "Dimension =") +
  scale_color_manual(values = cols_spp) +
  scale_fill_manual(values = cols_spp) +
  theme_classic(base_size = 11) +
  theme(legend.key.width = unit(1.3, "cm"))

ggsave(file.path(OUTDIR_IMG, "phylo_bayes_curves_M2.png"), curves_plot,
       width = 12, height = 7.5, dpi = 150)

for (ph in c("PRE", "POST")) {
  p <- curves %>%
    mutate(level = ifelse(level == "q05", "p05", "p95"),
           spp = unname(spp_lab[species])) %>%
    filter(phase == ph) %>%
    ggplot(aes(hr, med, color = spp, fill = spp, linetype = level)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~ dim, ncol = 1, labeller = labeller(dim = lab_dim), scales = "free_x") +
    labs(x = "Tiempo (h)", y = "Humedad predicha (%)",
         color = "Especie", fill = "Especie", linetype = "Dimension =") +
    scale_color_manual(values = cols_spp) +
    scale_fill_manual(values = cols_spp) +
    theme_classic(base_size = 11) +
    theme(legend.key.width = unit(1.3, "cm"))
  ggsave(file.path(OUTDIR_IMG, paste0("phylo_bayes_curves_M2_", tolower(ph), ".png")),
         p, width = 6.5, height = 12, dpi = 150)
}

p_slopes <- slopes |>
  ggplot(aes(spp, slope_hr, color = level)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_linerange(aes(ymin = lwr_hr, ymax = upr_hr),
                 position = position_dodge(width = 0.55), linewidth = 0.7) +
  geom_point(aes(shape = level), position = position_dodge(width = 0.55), size = 2) +
  facet_grid(phase ~ dim, labeller = labeller(dim = lab_dim)) +
  labs(x = NULL, y = "Pendiente estimada (% h\u207b\u00b9)",
       color = "Dimension =", shape = "Dimension =") +
  scale_color_manual(values = c(q05 = "#0072B2", q95 = "#D55E00")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR_IMG, "phylo_bayes_slopes_M2.png"), p_slopes,
       width = 11, height = 7, dpi = 150)

# ============================================================
# 5. Tabla CSV
# ============================================================
slopes |>
  dplyr::select(phase, dim, species, level, slope_hr, lwr_hr, upr_hr) |>
  arrange(phase, dim, species, level) |>
  write_csv(file.path(OUTDIR_CSV, "phylo_slopes_by_species.csv"))

cat("Hecho. Figuras en", OUTDIR_IMG, "y tabla en", OUTDIR_CSV, "\n")