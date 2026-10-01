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
#   07-img/curves_reference.png             curvas predichas a p10 / p90 con las
#                                           observaciones individuales de fondo
#   07-img/cumulative_loss_reference.png    diferencia de perdida acumulada entre
#                                           p90 y p10 (solo efectos significativos)
#   00-data/tablas_resumen/reference_slope_effects.csv
#   00-data/tablas_resumen/reference_cumulative_loss.csv
#   00-data/tablas_resumen/reference_predicted_mc_window.csv
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
TIME_CORTE <- 94
t94    <- as.vector((TIME_CORTE - TIME_M) / TIME_S)

dims <- c("Dim.1", "Dim.2", "Dim.3")
# Etiquetas de las facetas en ingles (figuras destined to the paper).
lab_dim <- c(Dim.1 = "Dim.1 (size)",
             Dim.2 = "Dim.2 (pericarp)",
             Dim.3 = "Dim.3 (scar)")

# Deciles 10 y 90 por eje, y mediana de cada eje para los puntos de referencia.
q10  <- sapply(dims, function(d) unname(quantile(df[[d]], 0.10)))
q90  <- sapply(dims, function(d) unname(quantile(df[[d]], 0.90)))
meds <- sapply(dims, function(d) unname(median(df[[d]], na.rm = TRUE)))

# Desviacion estandar de cada eje. El efecto de un rasgo se reporta por DESVIACION
# ESTANDAR del eje, no por unidad, porque solo asi se pueden comparar entre ejes:
# el coeficiente crudo esta en "% por unidad de Dim.x", y las unidades de los tres
# ejes no son intercambiables si sus dispersiones difieren.
#
# Las coordenadas del FAMD NO estan estandarizadas: su desviacion estandar es la
# raiz de su autovalor (1.674, 1.297 y 1.089 para Dim.1, Dim.2 y Dim.3), no 1. Es
# un error de interpretacion leer "1 unidad de Dim.x" como "1 desviacion estandar de
# Dim.x": el factor es 1.67 para el eje mas pesado y 1.09 para el mas ligero, de modo
# que la lectura por unidad subestima el efecto del tamaño respecto a los demas y
# exagera el de la cicatriz.
sd_ejes <- sapply(dims, function(d) sd(df[[d]]))

cat("Valores de referencia por eje (p10 / mediana / p90 / sd):\n")
print(data.frame(dim = dims, p10 = q10, mediana = meds, p90 = q90, sd = sd_ejes))
cat("\nNota: la desviacion de cada eje no es 1. Un incremento de UNA UNIDAD del eje\n",
    "no equivale a una desviacion estandar; el efecto por sd es el efecto por unidad\n",
    "multiplicado por la sd del eje.\n", sep = "")

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
#   Hay dos conversiones encadenadas y conviene no confundirlas.
#
#   1. `time_s` esta escalado, asi que el coeficiente viene en % de humedad por
#      desviacion estandar de tiempo. Dividir por TIME_S lo pasa a %/h, que es
#      la unidad con sentido biologico. Es una conversion de unidades.
#
#   2. El valor de Dim.x NO esta estandarizado. La coordenada de cada eje tiene
#      desviacion estandar igual a la raiz de su autovalor (1.674, 1.297 y 1.089
#      para Dim.1, Dim.2 y Dim.3), no 1. Por tanto el coeficiente por unidad de
#      Dim.x no es el efecto de una desviacion estandar del rasgo, sino el de una
#      unidad arbitraria de eje, y las unidades de los tres ejes no son
#      intercambiables entre si.
#
#   Por eso el efecto que se interpreta y se grafica es el POR DESVIACION
#   ESTANDAR: `estimacion_h * sd_ejes[dim]`. Es la unica lectura que permite
#   comparar entre ejes, que es lo que hace el analisis. El error de saltarse
#   este factor no es pequeno ni neutro: 1.67 para Dim.1 frente a 1.09 para
#   Dim.3, de modo que la lectura por unidad subestima al eje del tamaño y
#   exagera al de la cicatriz, comprimiendo la diferencia entre ambos.
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
    # Lectura por DESVIACION ESTANDAR del eje, que es la comparable entre ejes.
    # La de por unidad se conserva en las columnas `_h`: es la que se necesita
    # para los tramos p10-p90, y sirve de trazabilidad del calculo.
    dplyr::mutate(
      sd_eje        = sd_ejes[dim],
      estimacion_sd = estimacion_h * sd_eje,
      ic_lo_sd      = ic_lo_h     * sd_eje,
      ic_hi_sd      = ic_hi_h     * sd_eje
    ) |>
    # Contexto: que parte del efecto cubre el recorrido p10-p90 del eje.
    # Se mantiene en la escala por unidad porque el rango p10-p90 ya es una
    # cantidad observada de cada eje, comparable entre si sin factor adicional.
    dplyr::mutate(
      p90_menos_p10_h = estimacion_h * (q90[dims] - q10[dims])[dim]
    )
}

efectos <- bind_rows(
  map_dfr(names(modelos), function(f) efecto_por_eje(modelos[[f]], f))
) |>
  # `phase` como factor con orden explicito: si fuera character, facet_wrap la
  # ordenaria alfabeticamente y POST (P-O-S-T) apareceria antes que PRE.
  dplyr::mutate(dim = factor(dim, levels = dims),
                phase = factor(phase, levels = c("PRE", "POST")))

write_csv(efectos, file.path(OUTDIR_CSV, "reference_slope_effects.csv"))
cat("\nEfecto de cada eje sobre la pendiente (coeficiente time_s:Dim.x).\n",
    "Se lee por desviacion estandar del eje, que es la escala comparable\n",
    "entre ejes; las columnas por unidad quedan en el CSV solo como\n",
    "trazabilidad del calculo.\n", sep = "")
print(efectos |>
        dplyr::select(phase, dim, sd_eje, estimacion_sd, ic_lo_sd, ic_hi_sd, p) |>
        as.data.frame())

# ============================================================
# 3. Curvas predichas a p10 y p90 de cada eje
# ============================================================
# Curvas de humedad predicha a lo largo del tiempo, fijando un eje en p10 o
# p90 y los demas en su mediana. De nuevo se marginaliza sobre especie.
# Curvas predichas a p10 y p90 de cada eje.
#
# PRE y POST se dibujan en la MISMA faceta, sobre la misma escala de tiempo,
# para que se lea la desecacion como un proceso continuo y no como dos
# analisis separados. Cada fase aporta el tramo de su propio rango observado:
# la curva PRE cubre desde el inicio hasta 94 h y la curva POST desde 94 h
# hasta el final. Las curvas no se extrapolan fuera de ese tramo, porque el
# modelo de cada fase no esta ajustado para esos valores de tiempo.
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
        time_s = em$time_s,
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

# Secuencia sobre el rango temporal completo. Se anaden explicitamente el
# ultimo tiempo de PRE y el primero de POST: son los puntos por los que cada
# curva empieza y termina tras el recorte, y sin ellos geom_line los
# descartaria al no haber ningun valor adyacente dentro del rango.
secuencia_completa <- sort(unique(c(
  seq(rango_pre[1], rango_post[2], length.out = 60),
  rango_pre[2], rango_post[1], t94
)))


curvas <- bind_rows(
  curvas_por_eje(modelos[["PRE"]],  "PRE",  secuencia_completa),
  curvas_por_eje(modelos[["POST"]], "POST", secuencia_completa)
) |>
  # Cada curva se recorta a su propio rango observado. Sin esto, emmeans
  # devolveria la prediccion del modelo PRE para tiempos de la fase POST, que es
  # extrapolacion: ese modelo no esta ajustado para esos valores.
  #
  # El recorte es elemento a elemento, asi que no hace falta agrupar. La
  # condicion se evalua con `case_when` y no con `if`: `if` exige una condicion
  # de longitud 1 y aqui `phase` tiene un valor por fila, que es exactamente lo
  # que hace que aborte con "the condition has length > 1".
  dplyr::mutate(
    en_rango = dplyr::case_when(
      phase == "PRE"  ~ time_s <= rango_pre[2],
      phase == "POST" ~ time_s >= rango_post[1],
      .default = FALSE
    ),
    med = ifelse(en_rango, med, NA_real_),
    lwr = ifelse(en_rango, lwr, NA_real_),
    upr = ifelse(en_rango, upr, NA_real_)
  )

# ============================================================
# 4. Perdida acumulada de agua entre p10 y p90, dentro de cada fase
# ============================================================
# POR QUE NO SE INTEGRA HASTA EL UMBRAL DE 94 h
#   El corte en 94 h es ANALITICO: separa las fases de los modelos. No es un
#   tiempo de muestreo. Entre 45.2 h y 139.67 h no hay ninguna observacion, de
#   modo que el umbral cae dentro de un hueco de ~94 h. Acumular la perdida
#   "hasta las 94 h" obligaria a proyectar el modelo sobre un tramo que ninguna
#   de las dos fases tiene ajustado: es la misma extrapolacion que las curvas de
#   la seccion 3 evitan por construccion.
#   La ventana de integracion es por tanto la REAL de cada fase, que es la
#   unica que su propio modelo puede integrar sin salirse de sus datos.
#
# POR QUE PUNTOS DE HUMEDAD Y NO %/h
#   Un efecto de 0.0027 %/h no dice cuanta agua se pierde: depende de cuanto
#   dura la fase. La perdida acumulada en puntos de humedad es comparable con la
#   desecacion total observada y responde a la pregunta en contenido de agua, no
#   en velocidad. Ademas la fase POST dura ~5 veces mas que la PRE, asi que un
#   efecto por hora mas pequeno puede dar una perdida acumulada comparable.
#
# POR QUE LA DIFERENCIA SE SACA DEL COEFICIENTE Y NO DE UN CONTRASTE
#   El modelo es lineal en `time_s`, luego MC(d,t) = b0 + (b_t + b_int*d)*t.
#   Al restar la prediccion en p10 de la de p90 se cancelan los terminos que no
#   dependen de d, y la diferencia de perdida queda
#
#       dif_pp = (b_int / TIME_S) * (p90 - p10) * duracion_fase
#
#   es decir `p90_menos_p10_h` multiplicado por la duracion observada. Es una
#   multiplicacion por una constante positiva, asi que el intervalo del
#   coeficiente se traslada tal cual y solo hace falta el error estandar de la
#   interaccion, que ya se leyo en la seccion 2.

ventana_fase <- tibble::tibble(
  # `phase` es factor con los mismos niveles que en `efectos`: el join de la
  # seccion siguiente empareja factor contra factor, no factor contra texto.
  phase      = factor(c("PRE", "POST"), levels = c("PRE", "POST")),
  t_ini_h    = c(min(df$time[df$time_s < t94]), min(df$time[df$time_s > t94])),
  t_fin_h    = c(max(df$time[df$time_s < t94]), max(df$time[df$time_s > t94]))
) |>
  dplyr::mutate(
    duracion_h = t_fin_h - t_ini_h,
    t_ini_s    = (t_ini_h - TIME_M) / TIME_S,
    t_fin_s    = (t_fin_h - TIME_M) / TIME_S
  )

cat("\nVentana observada de cada fase (el umbral de 94 h cae en el hueco entre\n",
    "las dos; no se puede integrar hasta 94 h sin extrapolar):\n", sep = "")
print(as.data.frame(ventana_fase |> dplyr::select(phase, t_ini_h, t_fin_h, duracion_h)))

acumulado <- efectos |>
  dplyr::left_join(ventana_fase, by = "phase") |>
  dplyr::mutate(
    rango_dim = (q90[dims] - q10[dims])[dim],
    dif_pp    = estimacion_h * rango_dim * duracion_h,
    ic_lo_pp  = ic_lo_h     * rango_dim * duracion_h,
    ic_hi_pp  = ic_hi_h     * rango_dim * duracion_h,
    # p90 retiene MAS humedad si la diferencia es positiva; si es negativa, p90
    # pierde mas agua que p10.
    lectura = dplyr::case_when(
      dif_pp > 0 ~ "p90 retains more moisture",
      TRUE        ~ "p90 loses more moisture"
    )
  ) |>
  # Solo las combinaciones con efecto significativo. En las otras dos el
  # intervalo ya demuestra que cualquier efecto remanente es despreciable, y
  # multiplicar un efecto no significativo por la duracion de la fase genera un
  # numero mayor que parece informativo y no lo es.
  dplyr::filter(p < 0.05) |>
  dplyr::mutate(
    dim  = factor(dim, levels = dims),
    phase = factor(phase, levels = c("PRE", "POST"))
  )

write_csv(acumulado, file.path(OUTDIR_CSV, "reference_cumulative_loss.csv"))
cat("\nPerdida acumulada: diferencia entre p90 y p10 (puntos de humedad).\n",
    "Positivo = las bellotas de p90 terminan con mas agua que las de p10.\n", sep = "")
print(as.data.frame(acumulado |>
  dplyr::select(phase, dim, duracion_h, dif_pp, ic_lo_pp, ic_hi_pp, p, lectura)))

# Predicciones de humedad al inicio y al final de la ventana de cada fase, con
# cada eje en p10 y p90 y los otros dos en su mediana. Da el contexto que
# necesita la tabla anterior: no solo cuanto separa a los extremos del eje, sino
# cuantos puntos de humedad se pierden en total.
mc_ventana <- function(mod, fase, t_ini_s, t_fin_s) {
  bind_rows(map_dfr(dims, function(d) {
    map_dfr(c("p10", "p90"), function(niv) {
      at <- as.list(meds)
      names(at) <- dims
      at[[d]] <- if (niv == "p10") q10[d] else q90[d]

      em <- as.data.frame(emmeans::emmeans(mod, ~ time_s,
                                           at = c(at, list(time_s = c(t_ini_s, t_fin_s)))))
      col_est <- intersect(c("emmean", "emestimate", "estimate"), names(em))[1]
      ic <- intervalos_de(em, col_est)

      tibble::tibble(
        dim      = d,
        level    = niv,
        phase    = fase,
        # El tiempo se recupera de la escala, no del orden de las filas de
        # emmeans: si el orden de la retícula cambiara, el par de tiempos
        # quedaria descolgado de sus predicciones.
        tiempo_h = TIME_M + em$time_s * TIME_S,
        mc       = em[[col_est]],
        lo       = ic$lo,
        hi       = ic$hi
      )
    })
  }))
}

mc_ventana_largo <- bind_rows(
  mc_ventana(modelos[["PRE"]],  "PRE",  ventana_fase$t_ini_s[1], ventana_fase$t_fin_s[1]),
  mc_ventana(modelos[["POST"]], "POST", ventana_fase$t_ini_s[2], ventana_fase$t_fin_s[2])
)

perdida_por_nivel <- mc_ventana_largo |>
  dplyr::group_by(phase, dim, level) |>
  dplyr::arrange(tiempo_h, .by_group = TRUE) |>
  dplyr::summarise(
    mc_ini     = dplyr::first(mc),
    mc_fin     = dplyr::last(mc),
    perdida_pp = dplyr::first(mc) - dplyr::last(mc),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dim   = factor(dim, levels = dims),
    phase = factor(phase, levels = c("PRE", "POST")),
    level = factor(level, levels = c("p10", "p90"))
  )

write_csv(perdida_por_nivel,
          file.path(OUTDIR_CSV, "reference_predicted_mc_window.csv"))
cat("\nHumedad predicha al inicio y al final de cada fase, y perdida total:\n")
print(as.data.frame(perdida_por_nivel |>
  dplyr::select(phase, dim, level, mc_ini, mc_fin, perdida_pp)))

# ============================================================
# 5. Figuras
# ============================================================
# Figura principal: efecto de cada eje sobre la pendiente.
# Cada punto es el coeficiente `time_s:Dim.x` expresado por desviacion estandar
# de ese eje y por hora. El cero vertical es "el rasgo no cambia la velocidad de
# desecacion".
#
# Se grafica por DESVIACION ESTANDAR y no por unidad porque el objetivo de la
# figura es comparar los tres ejes entre si, y solo la escala por sd es comparable:
# las coordenadas del FAMD no estan estandarizadas (sd = raiz del autovalor) y sus
# dispersiones difieren, de modo que la escala por unidad premia al eje mas
# disperso sin que ello signifique mayor influencia biologica.
p_efecto <- efectos |>
  ggplot(aes(x = estimacion_sd, y = dim)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ic_lo_sd, xmax = ic_hi_sd),
                 height = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ phase, ncol = 2) +
  scale_y_discrete(labels = lab_dim) +
  labs(x = expression(paste("Change in desiccation rate (",
                            Delta * " per 1 sd of the axis, % h"^-1 * ")")),
       y = NULL) +
  theme_classic(base_size = 11)

ggsave(file.path(OUTDIR_IMG, "slope_effect_reference.png"), p_efecto,
       width = 8.5, height = 3.6, dpi = 300)

# Perdida acumulada entre los extremos del eje, en puntos de humedad. Es la
# misma magnitud que la figura anterior pero en la unidad en la que se lee la
# desecacion: cuantos puntos de humedad separa a p90 de p10 al final de la fase.
# Solo aparecen las combinaciones con p < 0.05; las demas ya se ha mostrado que
# son despreciables, y dibujarlas daria a un intervalo que cruza el cero la
# apariencia de un resultado.
p_acumulado <- acumulado |>
  ggplot(aes(x = dif_pp, y = dim)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ic_lo_pp, xmax = ic_hi_pp),
                 height = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ phase, ncol = 2) +
  scale_y_discrete(labels = lab_dim) +
  labs(x = expression(paste("Difference in water loss between p90 and p10 ",
                            "(percentage points)")),
       y = NULL,
       subtitle = paste0("Significant trait effects only (p < 0.05). ",
                         "Positive values: p90 acorns retain more water")) +
  theme_classic(base_size = 11)

ggsave(file.path(OUTDIR_IMG, "cumulative_loss_reference.png"), p_acumulado,
       width = 8.5, height = 3.6, dpi = 300)

for (ph in c("PRE", "POST")) {
  p <- p_efecto + facet_wrap(~ phase, ncol = 1) +
    ggtitle(paste0("Phase ", ph))
  ggsave(file.path(OUTDIR_IMG, paste0("slope_effect_reference_",
                                     tolower(ph), ".png")),
         p, width = 7, height = 3.2, dpi = 300)
}

# Figura secundaria: curvas predichas. PRE y POST comparten faceta para leer el proceso completo:
# la curva PRE llega hasta 94 h y la POST arranca ahi. La linea vertical marca
# el corte. El color sigue siendo el nivel del eje (p10/p90); la fase se
# distingue por el tramo de tiempo que ocupa cada curva.
#
# `group = interaction(level, phase)` es lo que impide que geom_line una PRE con
# POST. Sin el, ggplot agruparia solo por `level` (la unica variable discreta
# mapeada) y la linea cruzaria la region intermedia, uniendo el ultimo punto de
# PRE con el primero de POST. Ese tramo se deja en blanco a proposito: son horas
# con muy poca densidad de datos y con las que ninguna de las dos fases esta
# ajustada, asi que no se puede predecir con seguridad. Los valores ya vienen en
# NA desde el recorte por rango, y el grupo garantiza que aunque se solaparan en
# el tiempo seguirian siendo dos segmentos separados.
#
# Las observaciones individuales van de fondo en la misma figura. Se dibujan con
# `inherit.aes = FALSE` porque si heredaran el aes principal (`colour = level`,
# `group`) cada punto recibiria el color de un nivel que no le corresponde y se
# agruparia por `dim`+`level`, que no es lo que se quiere. Sin faceta propia: se
# reparten en las tres paneles por el tiempo que ocupan, que es lo que hace falta
# para juzgar si las curvas siguen los datos.
p_curvas <- curvas |>
  dplyr::mutate(dim = factor(dim, levels = dims),
                level = factor(level, levels = c("p10", "p90")),
                phase = factor(phase, levels = c("PRE", "POST"))) |>
  ggplot(aes(hr, med, colour = level, fill = level, linetype = level,
             group = interaction(level, phase))) +
  geom_point(data = df, aes(time, Moisture_content), inherit.aes = FALSE,
             alpha = 0.05, size = 0.5, colour = "grey60") +
  geom_vline(xintercept = TIME_CORTE, linetype = 3, colour = "grey45") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  facet_wrap(~ dim, nrow = 1, labeller = labeller(dim = lab_dim)) +
  scale_colour_manual(values = c(p10 = "#0072B2", p90 = "#D55E00")) +
  scale_fill_manual(values = c(p10 = "#0072B2", p90 = "#D55E00")) +
  scale_x_continuous(breaks = seq(0, max(df$time), by = 24)) +
  labs(x = "Time (h)", y = "Moisture content (%)",
       title = "Desiccation in two phases (points: individual observations)",
       subtitle = paste0("Dashed line: phase boundary at ", TIME_CORTE, " h"),
       colour = "Axis level", fill = "Axis level",
       linetype = "Axis level") +
  theme_classic(base_size = 11) +
  theme(legend.key.width = unit(1.3, "cm"),
        strip.background = element_blank())

ggsave(file.path(OUTDIR_IMG, "curves_reference.png"), p_curvas,
       width = 13, height = 4.4, dpi = 300)

# Las observaciones individuales ya van de fondo en la figura de curvas, asi que
# no hay figura aparte de nube de puntos: era el mismo contenido duplicado y sin
# las curvas, que es justo lo que hace falta para juzgar si las ajustan.

cat("\nHecho. Figuras en", OUTDIR_IMG, "y tabla en", OUTDIR_CSV, "\n")
