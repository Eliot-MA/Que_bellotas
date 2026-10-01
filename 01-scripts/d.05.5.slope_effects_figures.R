# ============================================================
# d.05.5.slope_effects_figures.R
# Figuras y tabla del EFECTO de cada eje del FAMD sobre la PENDIENTE de
# desecacion, en las fases PRE (t < 94 h) y POST (t > 94 h).
#
# MODELO: m.ref.pre / m.ref.post (glmmTMB), los de referencia ajustados en
#   01-scripts/d.05.3.reference_model.R
#
#   Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
#     (0 + time_s | species) + (1 + time_s | prov_code) + (1 | id_bellota)
#
# POR QUE SE LEN LOS COEFICIENTES DE LA INTERACCION
#   Cuando la especie estaba en los FIJOS, con una pendiente por especie, el
#   coeficiente `time_s:Dim.x` describia solo a la especie de referencia y la
#   pendiente de cualquier otra especie salia sumando coeficientes a mano. Por
#   eso antes se usaban pendientes marginales con emmeans.
#
#   Ahora la especie es una pendiente ALEATORIA y `prov_code` esta anidada en
#   `species`, de modo que los terminos aleatorios ya capturan toda la
#   heterogeneidad entre especies. No queda especie de referencia que explicar:
#   el coeficiente de la interaccion es directamente el efecto del rasgo sobre
#   la pendiente de la poblacion. Es una lectura, no un calculo.
#
#   La heterogeneidad por especie NO desaparece, vive en la varianza del termino
#   aleatorio (00-data/reference_model_varcomp.csv, d.05.1). No hace falta
#   resumir una pendiente por especie.
#
# Sin refit: se usan los modelos guardados en 00-data/models/.
#
# Salidas:
#   07-img/slope_effect_reference.png        efecto sobre la pendiente (principal)
#   07-img/slope_effect_reference_pre.png    idem, solo PRE
#   07-img/slope_effect_reference_post.png   idem, solo POST
#   07-img/curves_reference.png             curvas predichas a p10 / p90
#   07-img/observed_reference.png           nube de puntos observada
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
# 2. Efecto de cada eje sobre la pendiente: coeficientes de la interaccion
# ============================================================
# POR QUE SE LEEN LOS COEFICIENTES Y NO PENDIENTES MARGINALES
#   En el modelo anterior la especie estaba en los FIJOS, con una pendiente por
#   especie. El coeficiente de la interaccion describia solo a la especie de
#   referencia, y obtener la pendiente de otra especie obligaba a sumar a mano
#   su coeficiente mas el de la interaccion. Era aritmetica mental, no lectura.
#
#   Aqui la especie es una pendiente ALEATORIA, y ademas `prov_code` esta
#   anidada en `species`, de modo que los terminos aleatorios ya capturan toda
#   la heterogeneidad entre especies. No queda ninguna especie de referencia que
#   explicar: el coeficiente `time_s:Dim.x` ES el efecto del rasgo sobre la
#   pendiente de la poblacion, directamente y sin correcciones.
#
#   La heterogeneidad por especie sigue existiendo, pero vive en la varianza del
#   termino aleatorio, no en los coeficientes fijos (ver d.05.1 y
#   reference_model_varcomp.csv). Por eso no hace falta resumir una pendiente
#   por especie.
#
# UNIDADES
#   `time_s` esta escalado, asi que el coeficiente esta en % de humedad por
#   desviacion estandar de tiempo. Se divide por TIME_S para obtener %/h, que es
#   la unidad con sentido biologico. El valor de Dim.x es la variacion por una
#   desviacion estandar de ese eje del FAMD.
#
# El tramo p10-p90 se calcula aparte solo como contexto de magnitud (cuanto
# cambio de pendiente cubre el recorrido del 10 al 90 percentil del eje). No es
# el efecto que se interpreta: es el coeficiente multiplicado por el rango
# observado del eje.
#
# `intervalos_de` lo usan las curvas de la seccion 4.

# Los limites del intervalo tambien cambian de nombre segun el metodo
# (asymp.LCL/asymp.UCL frente a lower.CL/upper.CL). Si no estan, se recalculan
# con el cuantil t usando los grados de libertad que da el propio emmeans, y no
# con 1.96: sobre glmmTMB los intervalos son t y usar la normal acortaria los
# limites alatorios.
intervalos_de <- function(df, col_est) {
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(df))) {
    return(list(lo = df$asymp.LCL, hi = df$asymp.UCL))
  }
  if (all(c("lower.CL", "upper.CL") %in% names(df))) {
    return(list(lo = df$lower.CL, hi = df$upper.CL))
  }
  q <- if ("df" %in% names(df)) stats::qt(0.975, df$df) else stats::qnorm(0.975)
  list(lo = df[[col_est]] - q * df$SE, hi = df[[col_est]] + q * df$SE)
}

efecto_por_eje <- function(mod, fase) {
  suppressPackageStartupMessages(library(parameters))
  tab <- parameters::model_parameters(mod, effects = "fixed", ci = 0.95) |>
    as.data.frame()

  # La interaccion de interes se llama `time_s:Dim.x`. Se busca por patron y no
  # por nombre exacto para no depender de como glmmTMB etiquete la columna.
  es_interaccion <- grepl("^time_s:Dim\\.[123]$", tab$Parameter)
  if (!any(es_interaccion)) {
    stop("No encuentro las interacciones time_s:Dim.x en ", fase,
         ". Coeficientes disponibles: ",
         paste(tab$Parameter, collapse = ", "), call. = FALSE)
  }
  tab <- tab[es_interaccion, ]

  tibble::tibble(
    phase       = fase,
    dim         = sub("^time_s:", "", tab$Parameter),
    estimacion  = tab$Coefficient,
    se          = tab$SE,
    ic_lo       = tab$CI_low,
    ic_hi       = tab$CI_high,
    p           = tab$p,
    # de % por sd(horas) a % por hora
    estimacion_h = tab$Coefficient / TIME_S,
    ic_lo_h     = tab$CI_low / TIME_S,
    ic_hi_h     = tab$CI_high / TIME_S
  ) |>
    # Contexto: que parte del efecto cubre el recorrido p10-p90 del eje.
    dplyr::mutate(
      p90_menos_p10_h = estimacion_h * (q90[dims] - q10[dims])[dim]
    )
}

efectos <- bind_rows(
  map_dfr(names(modelos), function(f) efecto_por_eje(modelos[[f]], f))
) |>
  dplyr::mutate(dim = factor(dim, levels = dims))

write_csv(efectos, file.path(OUTDIR_CSV, "reference_slope_effects.csv"))
cat("\nEfecto de cada eje sobre la pendiente (coeficiente time_s:Dim.x):\n")
print(efectos |>
        dplyr::select(phase, dim, estimacion_h, ic_lo_h, ic_hi_h, p) |>
        as.data.frame())

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

      em <- as.data.frame(
        emmeans::emmeans(mod, ~ time_s, at = c(at, list(time_s = ts_seq)))
      )
      # `emmean` y los limites tampoco tienen nombre fijo segun el metodo.
      col_est <- intersect(c("emmean", "emestimate", "estimate"), names(em))[1]
      ic <- intervalos_de(em, col_est)

      tibble::tibble(
        dim   = d,
        level = niv,
        phase = fase,
        hr    = TIME_M + em$time_s * TIME_S,
        med   = em[[col_est]],
        lwr   = ic$lo,
        upr   = ic$hi
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
# Cada punto es el coeficiente `time_s:Dim.x`, es decir el cambio de pendiente
# por una desviacion estandar de ese eje. El cero vertical es "el rasgo no
# cambia la velocidad de desecacion".
p_efecto <- efectos |>
  ggplot(aes(x = estimacion_h, y = dim)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ic_lo_h, xmax = ic_hi_h),
                 height = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ phase, ncol = 2) +
  scale_y_discrete(labels = lab_dim) +
  labs(x = expression(paste("Cambio de la pendiente (",
                            Delta * " por 1 sd del eje, % h"^-1 * ")")),
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
