# ============================================================
# d.05.5.slope_effects_figures.R
# Figuras y tabla del EFECTO MARGINAL de cada eje del FAMD sobre la PENDIENTE
# de desecacion, en las fases PRE (t < 94 h) y POST (t > 94 h).
#
# MODELO: m.ref.pre / m.ref.post (glmmTMB), los de referencia ajustados en
#   01-scripts/d.05.3.reference_model.R
#
#   Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
#     (0 + time_s | species) + (1 + time_s | prov_code) + (1 | id_bellota)
#
# POR QUE NO HAY PENDIENTES POR ESPECIE
#   El modelo de referencia lleva la especie como PENDIENTE ALEATORIA, no como
#   coeficiente fijo por especie. Como las pendientes aleatorias tienen media
#   cero, la prediccion marginal que devuelve emmeans YA es la pendiente
#   promediada sobre todas las especies: no hace falta estimar ni resumir una
#   pendiente por especie. Es la lectura que conecta con la pregunta del
#   articulo (los rasgos cambian la velocidad de desecacion) y evita
#   discutir 8 efectos por rasgo sin informacion para hacerlo.
#
#   La heterogeneidad por especie NO desaparece: vive en la varianza del
#   termino aleatorio, que se exporta en 00-data/reference_model_varcomp.csv y
#   se justifica en d.05.1.heterogeneus_effects.R. Si esa varianza de
#   pendiente fuera casi nula, d.05.3 lo avisa antes de llegar aqui.
#
# PROTOCOLO DE VALORES
#   Deciles 10 y 90 de cada eje, igual que en Heterogeneus_effects_acorn_traits.qmd,
#   para que el contraste no extrapole fuera del rango observado. El resto de
#   ejes se fija en su mediana.
#
# Sin refit: se usan los modelos guardados en 00-data/models/.
#
# Salidas:
#   07-img/slope_effect_reference.png        efecto sobre la pendiente (principal)
#   07-img/slope_effect_reference_pre.png    idem, solo PRE
#   07-img/slope_effect_reference_post.png   idem, solo POST
#   07-img/curves_reference.png             curvas predichas a p10 / p90
#   00-data/tablas_resumen/reference_slope_effects.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(emmeans)
})

OUTDIR_IMG <- "07-img"
OUTDIR_CSV <- "00-data/tablas_resumen"
dir.create(OUTDIR_IMG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_CSV, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Datos
# ============================================================
if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd     <- read.csv("00-data/famd_ind_coord.csv")
df <- df.bellotas |>
  dplyr::select(-X) |>
  dplyr::select(id_bellota, prov_code, tiempo_acumulado_horas, Moisture_content) |>
  # famd_ind_coord.csv tambien trae `prov_code`; se descarta para no duplicar el
  # nombre al hacer el join.
  left_join(y = df.famd |> dplyr::select(-prov_code), by = "id_bellota") |>
  dplyr::filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
  tidyr::drop_na(Dim.1, Dim.2, Dim.3) |>
  rename(time = tiempo_acumulado_horas) |>
  mutate(
    time_s     = as.vector(scale(time)),
    species    = factor(species),
    provenance = factor(provenance),
    prov_code  = factor(prov_code),
    id_bellota = factor(id_bellota)
  )

assert_sin_procedencias_excluidas(df, "prov_code", "d.05.5 datos")

TIME_M <- mean(df$time)
TIME_S <- sd(df$time)
t94    <- as.vector((94 - TIME_M) / TIME_S)

dims <- c("Dim.1", "Dim.2", "Dim.3")
lab_dim <- c(Dim.1 = "Dim.1 (tamano)",
             Dim.2 = "Dim.2 (pericarpo)",
             Dim.3 = "Dim.3 (cicatriz)")

# Deciles 10 y 90 por eje, y mediana de cada eje para los puntos de referencia.
q10  <- sapply(dims, function(d) unname(quantile(df[[d]], 0.10)))
q90  <- sapply(dims, function(d) unname(quantile(df[[d]], 0.90)))
meds <- sapply(dims, function(d) unname(median(df[[d]], na.rm = TRUE)))

cat("Valores de referencia por eje (p10 / mediana / p90):\n")
print(data.frame(dim = dims, p10 = q10, mediana = meds, p90 = q90))

# ============================================================
# 1. Cargar los modelos de referencia
# ============================================================
rutas <- c(PRE = "00-data/models/m.ref.pre.rds",
           POST = "00-data/models/m.ref.post.rds")
faltantes <- rutas[!file.exists(rutas)]
if (length(faltantes) > 0) {
  stop("Faltan los modelos de referencia: ", paste(faltantes, collapse = ", "),
       ". Ejecuta primero 01-scripts/d.05.3.reference_model.R", call. = FALSE)
}
modelos <- setNames(lapply(rutas, readRDS), names(rutas))

# ============================================================
# 2. Efecto de cada eje sobre la pendiente, promediado sobre especies
# ============================================================
# emtrends con `specs = dimension` y `var = "time_s"` devuelve la pendiente de
# la humedad frente al tiempo en el valor del eje indicado. Al no incluir
# `species` en las especificaciones, emmeans marginaliza sobre el termino
# aleatorio de especie: el resultado es la pendiente de la poblacion, es decir
# la media sobre todas las especies.
#
# El contraste entre p10 y p90 es el cambio de pendiente atribuible al rasgo,
# en % de humedad por unidad de tiempo escalada. Se divide por TIME_S para
# pasarlo a %/h, que es la unidad con sentido biologico.
efecto_por_eje <- function(mod, fase) {
  map_dfr(dims, function(d) {
    at <- as.list(meds)
    names(at) <- dims
    at[[d]] <- c(q10[d], q90[d])

    tr <- emmeans::emtrends(mod, specs = d, var = "time_s", at = at)
    co <- emmeans::contrast(tr, method = "revpairwise", adjust = "none")

    as.data.frame(co) |>
      # emmeans anade `dim.coe` cuando la tendencia tiene mas de un
      # coeficiente; aqui solo interesa el contraste de la pendiente.
      dplyr::select(-dplyr::any_of("dim.coe")) |>
      dplyr::mutate(dim = d, phase = fase) |>
      dplyr::rename(estimacion = em.trend, se = SE,
                    ic_lo = asymp.LCL, ic_hi = asymp.UCL) |>
      dplyr::mutate(
        # de % por sd(horas) a % por hora
        estimacion_h = estimacion / TIME_S,
        ic_lo_h     = ic_lo / TIME_S,
        ic_hi_h     = ic_hi / TIME_S
      )
  })
}

efectos <- bind_rows(
  map_dfr(names(modelos), function(f) efecto_por_eje(modelos[[f]], f))
)

write_csv(efectos, file.path(OUTDIR_CSV, "reference_slope_effects.csv"))
cat("\nEfecto de cada eje sobre la pendiente (promediado sobre especies):\n")
print(efectos |> dplyr::select(phase, dim, estimacion_h, ic_lo_h, ic_hi_h))

# ============================================================
# 3. Curvas predichas a p10 y p90 de cada eje
# ============================================================
# Curvas de humedad predicha a lo largo del tiempo, fijando un eje en p10 o
# p90 y los demas en su mediana. De nuevo se marginaliza sobre especie.
curvas_por_eje <- function(mod, fase, ts_seq) {
  map_dfr(dims, function(d) {
    map_dfr(c("p10", "p90"), function(niv) {
      val <- if (niv == "p10") q10[d] else q90[d]
      at  <- as.list(meds)
      names(at) <- dims
      at[[d]] <- val

      emmeans::emmeans(mod, ~ time_s, at = at) |>
        as.data.frame() |>
        dplyr::mutate(
          dim   = d,
          level = niv,
          phase = fase,
          hr    = TIME_M + time_s * TIME_S,
          med   = emmean,
          lwr   = asyp.LCL,
          upr   = asyp.UCL
        )
    })
  })
}

rango_pre  <- range(df$time_s[df$time_s < t94])
rango_post <- range(df$time_s[df$time_s > t94])

curvas <- bind_rows(
  curvas_por_eje(modelos[["PRE"]],  "PRE",  seq(rango_pre[1],  rango_pre[2],  length.out = 30)),
  curvas_por_eje(modelos[["POST"]], "POST", seq(rango_post[1], rango_post[2], length.out = 30))
)

# ============================================================
# 4. Figuras
# ============================================================
# Figura principal: efecto de cada eje sobre la pendiente.
p_efecto <- efectos |>
  dplyr::mutate(dim = factor(dim, levels = dims)) |>
  ggplot(aes(x = estimacion_h, y = dim)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ic_lo_h, xmax = ic_hi_h),
                 height = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ phase, ncol = 2) +
  scale_y_discrete(labels = lab_dim) +
  labs(x = expression(paste("Cambio de la pendiente (",
                            Delta * " p90 - p10, % h"^-1 * ")")),
       y = NULL) +
  theme_classic(base_size = 11)

ggsave(file.path(OUTDIR_IMG, "slope_effect_reference.png"), p_efecto,
       width = 8.5, height = 3.6, dpi = 300)

for (ph in c("PRE", "POST")) {
  p <- p_efecto + facet_wrap(~ phase, ncol = 1) +
    ggtitle(paste0("Fase ", ph))
  ggsave(file.path(OUTDIR_IMG, paste0("slope_effect_reference_",
                                     tolower(ph), ".png")),
         p, width = 7, height = 3.2, dpi = 300)
}

# Curvas predichas.
p_curvas <- curvas |>
  dplyr::mutate(dim = factor(dim, levels = dims),
                level = factor(level, levels = c("p10", "p90"))) |>
  ggplot(aes(hr, med, colour = level, fill = level, linetype = level)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(phase ~ dim, scales = "free_x",
             labeller = labeller(dim = lab_dim)) +
  scale_colour_manual(values = c(p10 = "#0072B2", p90 = "#D55E00")) +
  scale_fill_manual(values = c(p10 = "#0072B2", p90 = "#D55E00")) +
  labs(x = "Tiempo (h)", y = "Humedad predicha (%)",
       colour = "Nivel del eje", fill = "Nivel del eje",
       linetype = "Nivel del eje") +
  theme_classic(base_size = 11) +
  theme(legend.key.width = unit(1.3, "cm"))

ggsave(file.path(OUTDIR_IMG, "curves_reference.png"), p_curvas,
       width = 12, height = 7.5, dpi = 300)

# ============================================================
# 5. Nube de puntos observada (contexto visual de las curvas predichas)
# ============================================================
p_obs <- df |>
  dplyr::mutate(phase = ifelse(time_s < t94, "PRE", "POST")) |>
  ggplot(aes(time, Moisture_content)) +
  geom_point(alpha = 0.06, size = 0.6, colour = "grey55") +
  facet_wrap(~ phase, ncol = 2) +
  labs(x = "Tiempo (h)", y = "Humedad (%)") +
  theme_classic(base_size = 11)

ggsave(file.path(OUTDIR_IMG, "observed_reference.png"), p_obs,
       width = 8, height = 3.4, dpi = 300)

cat("\nHecho. Figuras en", OUTDIR_IMG, "y tabla en", OUTDIR_CSV, "\n")
