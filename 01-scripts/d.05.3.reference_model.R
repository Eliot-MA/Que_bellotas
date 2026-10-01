# ============================================================
# d.05.3.reference_model.R
# MODELO DE REFERENCIA del articulo.
#
# Efecto de los rasgos funcionales de la bellota (ejes del FAMD) sobre la tasa
# de desecacion, con pendientes aleatorias por ESPECIE y por PROCEDENCIA, sin
# correccion filogenetica. Ajuste persistente con glmmTMB.
#
#   Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) +
#     (0 + time_s | species) +     # pendiente por especie, SIN intercepto
#     (1 + time_s | prov_code) +   # pendiente por procedencia (nivel basal)
#     (1 | id_bellota)             # bellota anidada, solo intercepto
#
# Se ajusta en las dos fases de desecacion (PRE t < 94 h, POST t > 94 h).
#
# Decisiones de diseno que conviene no deshacer sin evidencia:
#
#   1. El termino de especie es de PENDIENTE ALEATORIA, no un coeficiente fijo
#      por especie. Permite que la heterogeneidad exista (la exploracion de
#      d.05.1 muestra que es real) sin estimar un efecto del rasgo por cada
#      especie, que no es discutible con la informacion disponible.
#
#   2. El termino de especie es `(0 + time_s | species)`, SIN intercepto.
#      Esto no es una preferencia estetica, es una consecuencia de la
#      anidacion: `prov_code` esta dentro de `species` (15 codigos en 8
#      especies), de modo que cualquier intercepto de especie es ya
#      representable como intercepto de procedencia y anadir ambos deja una
#      dimension redundante. La variante con intercepto `(1 + time_s | species)`
#      esta ajustada mas abajo COMO `m.ref_int`. En PRE no converge: Hessian no
#      definida positiva, logLik NULL y AIC NA. En POST si converge, con aviso de
#      ajuste singular y AIC finito (8077.673), asi que aparece en la tabla de
#      comparacion marcado como singular. Es la evidencia directa de que la
#      dimension sobrante es real. Como
#      `prov_code` aporta el nivel basal entre especies, perder el intercepto
#      de especie no pierde informacion: en `m.ref_int.pre` la varianza de
#      intercepto de especie sale en 4.16e-06, o sea cero. En cambio la
#      varianza de PENDIENTE por especie si es real (3.38 en PRE, 1.08 en POST),
#      y por eso el termino se queda. La comparacion AIC de la seccion 2 da un
#      empate (PRE DeltaAIC 1.47 a favor de `m.ref_naive`), de modo que la
#      decision es conceptual y no sale del ajuste; el comentario sobre
#      `form_ref` mas abajo desarrolla el matiz.
#
#   3. NO se incluye filogenia. d.05.4_phylo_check.R muestra que la senal
#      filogenetica es debil e indistinguible de cero (CCI ~ 0.007) y que la
#      interpretacion no cambia al incorporarla. Se conserva como
#      comprobacion de robustez, no como modelo de referencia.
#
# Salidas:
#   00-data/reference_model_comparison.csv   comparacion AIC/pesos (PRE y POST)
#   00-data/reference_model_coef.csv        coeficientes fijos
#   00-data/reference_model_varcomp.csv     componentes de varianza
#   00-data/reference_model_varcomp_detalle.csv  varianzas separadas + ICC
#   00-data/reference_model_anova.csv       Anova por tipo III
#   00-data/models/m.ref*.rds               objetos de modelo
#   06-html/modeldashboard_m.ref*.html      dashboards (si GEN_DASH == "Y")
# ============================================================

GEN_DASH <- "N"

if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(car)         # Anova() para los tests globales
  library(performance) # compare_performance()
})

# ---- 0. Datos ----
# Reutiliza el dataframe del master si existe; si no, lo reconstruye con la
# misma receta que d.05.0.model_traits.R.
if (!exists("df") || !exists("df.t1") || !exists("df.t2")) {
  df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
  df.famd     <- read.csv("00-data/famd_ind_coord.csv")
  df <- df.bellotas |>
    dplyr::select(-X) |>
    dplyr::select(id_bellota, prov_code, tiempo_acumulado_horas, Moisture_content) |>
    # famd_ind_coord.csv tambien trae `prov_code`; se descarta para no duplicar
    # el nombre al hacer el join. Se conserva el de la tabla larga, previo al
    # filtro del FAMD.
    left_join(y = df.famd |> dplyr::select(-prov_code), by = "id_bellota") |>
    filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
    tidyr::drop_na(Moisture_content, Dim.1, Dim.2, Dim.3) |>
    rename(time = tiempo_acumulado_horas) |>
    mutate(
      time_s     = as.vector(scale(time)),
      species    = factor(species),
      provenance = factor(provenance),
      prov_code  = factor(prov_code),
      id_bellota = factor(id_bellota)
    )
  t94  <- as.vector((94 - mean(df$time)) / sd(df$time))
  df.t1 <- df |> filter(time_s < t94)
  df.t2 <- df |> filter(time_s > t94)
}

# `prov_code` debe ser factor incluso si el dataframe viene del master, que solo
# factoriza species / provenance / id_bellota.
df.t1$prov_code <- factor(df.t1$prov_code)
df.t2$prov_code <- factor(df.t2$prov_code)

# El termino aleatorio (1 + time_s | prov_code) debe correr sobre el mismo
# conjunto de procedencias que el resto de la linea d.*.
assert_sin_procedencias_excluidas(df,    "prov_code", "d.05.3 datos")
assert_sin_procedencias_excluidas(df.t1, "prov_code", "d.05.3 fase PRE")
assert_sin_procedencias_excluidas(df.t2, "prov_code", "d.05.3 fase POST")

source("01-scripts/00-export_helpers.R")

# ============================================================
# 1. Especificacion de los modelos
# ============================================================
# El cuerpo comun: interaccion tiempo x cada eje del FAMD, pendiente por
# procedencia y bellota anidada. El termino de especie se inyecta como
# parametro porque es justo lo que se quiere contrastar.
form_con_especie <- function(termino_especie) {
  paste0("Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + ",
         termino_especie, " + ",
         "(1 + time_s | prov_code) + ",
         "(1 | id_bellota)")
}

# ---- m.ref: EL MODELO DE REFERENCIA ----
#
# DECISION (2026): el termino de especie va SIN intercepto, `(0 + time_s | species)`.
#
# Motivo, que no es una preferencia de escritura sino una consecuencia de la
# estructura de los datos:
#
#   `prov_code` esta ANIDADA dentro de `species` (15 codigos en 8 especies): cada
#   codigo de procedencia pertenece a una sola especie. Por tanto cualquier
#   intercepto de especie ya es representable como intercepto de procedencia, y
#   estimarlos los dos a la vez deja una dimension redundante.
#
#   Esa redundancia no es teorica: la variante con intercepto de especie
#   (`m.ref_int` mas abajo) falla en PRE con "non-positive-definite Hessian
#   matrix", donde glmmTMB deja logLik en NULL y su AIC es NA, sin comparacion
#   posible. En POST converge, pero con "singular convergence (7)" y un AIC
#   finito (8077.673), asi que se mantiene en la comparacion marcado como
#   singular. El aviso de convergencia de este script comprueba ambas fases y
#   lo dice en voz alta en lugar de renormalizar pesos sobre los modelos que si
#   funcionan.
#
#   Que el nivel basal entre especies no se pierda: lo aporta `prov_code`. Al
#   haber al menos una procedencia por especie, el intercepto aleatorio de
#   procedencia ya absorbe la varianza de intercepto entre especies. Los datos
#   lo confirman de forma inequivoca: en `m.ref_int.pre` la varianza de
#   intercepto de especie sale en 4.16e-06, es decir, exactamente cero, porque
#   `prov_code` la absorbe entera. Lo propio de la especie, y lo que el modelo
#   conserva, es la PENDIENTE: cada especie se seca con su propia velocidad.
#
#   Comprobacion por AIC, y por que NO decide aqui: `m.ref` sale por DETRAS de
#   `m.ref_naive` (PRE DeltaAIC = 1.47 a favor de naive). Es un empate, no una
#   victoria del modelo con heterogeneidad. Pero la varianza de la pendiente por
#   especie esta estimada y no es nula: 3.38 en PRE (SD ~ 1.84) y 1.08 en POST
#   (SD ~ 1.04). O sea, el termino hace su trabajo; lo que ocurre es que con
#   solo 8 especies una varianza entre grupos apenas mueve la verosimilitud
#   marginal, de modo que AIC cobra sus dos grados de libertad sin devolver nada.
#
#   DECISION: se conserva `m.ref`. Es una eleccion deliberada a favor de no
#   asumir que las especies se secan igual, no una victoria del ajuste. La
#   honestidad del articulo va en decir la estimacion con su incertidumbre
#   (con 8 grupos no es precisa), no en presentarla como un efecto demostrado.
form_ref <- form_con_especie("(0 + time_s | species)")

# ---- m.ref_int: variante NO identificable, se conserva por registro ----
#
# Se mantiene ajustada a proposito, aunque no sea la referencia y aunque no
# converja, porque es la evidencia de por que la referencia no lleva intercepto
# de especie. Si algun dia cambia el diseno (mas especies, o procedencia no
# anidada), esta es la variante que habria que reevaluar.
form_ref_int <- form_con_especie("(1 + time_s | species)")

# ---- m.ref_naive: linea base sin especie ----
# Ignora la heterogeneidad interespecifica. Sirve para medir cuanto cuesta no
# modelarla (el "modelo ingenuo" de d.05.1).
form_ref_naive <- paste0("Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3) + ",
                         "(1 + time_s | prov_code) + ",
                         "(1 | id_bellota)")

cat("\n===== MODELOS AJUSTADOS =====\n\n")
cat("m.ref  (REFERENCIA, pendiente por especie sin intercepto):\n")
cat(form_ref, "\n\n")
cat("m.ref_int  (con intercepto de especie; no identificable):\n")
cat(form_ref_int, "\n\n")
cat("m.ref_naive  (linea base, sin especie):\n")
cat(form_ref_naive, "\n\n")

ajustar <- function(formula, data) {
  glmmTMB(as.formula(formula), data = data)
}

cat("\n-- Ajuste PRE (t < 94 h) --\n")
cat("--   m.ref (referencia)\n")
m.ref.pre       <- ajustar(form_ref, df.t1)
cat("--   m.ref_int (con intercepto de especie; se espera que falle en PRE)\n")
m.ref_int.pre   <- ajustar(form_ref_int, df.t1)
cat("--   m.ref_naive (linea base)\n")
m.ref_naive.pre <- ajustar(form_ref_naive, df.t1)

cat("\n-- Ajuste POST (t > 94 h) --\n")
cat("--   m.ref (referencia)\n")
m.ref.post       <- ajustar(form_ref, df.t2)
cat("--   m.ref_int (con intercepto de especie; se espera singular en POST)\n")
m.ref_int.post   <- ajustar(form_ref_int, df.t2)
cat("--   m.ref_naive (linea base)\n")
m.ref_naive.post <- ajustar(form_ref_naive, df.t2)

# ============================================================
# 1b. Estado de convergencia
# ============================================================
# Se informa antes de comparar nada. Un "convergence problem" de glmmTMB
# significa que la Hessian no es definida positiva o el ajuste es singular. En el
# primer caso logLik queda en NULL y el modelo no admite comparacion; en el
# segundo el ajuste puede conservar una verosimilitud valida, asi que se marca
# como singular sin excluirlo. Conviene verlo de forma explicita y no perderlo en
# el log.
#
# Que `m.ref_int` de problemas NO es un fallo del script ni del optimizador: es
# la manifestacion de la redundancia entre el intercepto de especie y el de
# procedencia descrita al definir las formulas. Se espera, y se comprueba.
estado_convergencia <- function(mod, etiqueta) {
  singular <- tryCatch(glmmTMB::isSingular(mod), error = function(e) NA)
  util     <- tryCatch(is.finite(stats::AIC(mod)), error = function(e) FALSE)
  cat(sprintf("  %-18s singular=%-5s AIC=%s\n",
              etiqueta, ifelse(is.na(singular), "?", singular),
              if (util) sprintf("%.1f", stats::AIC(mod)) else "NA (no utilizable)"))
  invisible(util)
}

cat("\n===== ESTADO DE CONVERGENCIA =====\n")
estado <- list(
  ref.pre        = estado_convergencia(m.ref.pre,        "m.ref.pre"),
  ref_int.pre    = estado_convergencia(m.ref_int.pre,    "m.ref_int.pre"),
  naive.pre      = estado_convergencia(m.ref_naive.pre,  "m.ref_naive.pre"),
  ref.post       = estado_convergencia(m.ref.post,       "m.ref.post"),
  ref_int.post   = estado_convergencia(m.ref_int.post,   "m.ref_int.post"),
  naive.post     = estado_convergencia(m.ref_naive.post, "m.ref_naive.post")
)

if (!estado$ref.pre || !estado$ref.post) {
  stop(paste0(
    "\nEL MODELO DE REFERENCIA m.ref NO HA CONVERGIDO.\n",
    "  Sin verosimilitud no hay AIC, ni comparacion, ni Anova, ni figura.\n",
    "  No se sigue: hay que resolver la causa antes de continuar.\n"), call. = FALSE)
}

if (!estado$ref_int.pre || !estado$ref_int.post) {
  cat(paste0(
    "\n  m.ref_int falla en al menos una fase, como se esperaba. Es coherente\n",
    "  con la anidacion de `prov_code` en `species`: el intercepto de especie\n",
    "  es redundante con el de procedencia. La referencia sigue siendo m.ref.\n",
    "  Ver el comentario de DECISION sobre form_ref.\n"))
  if (!estado$ref_int.pre && estado$ref_int.post) {
    cat(paste0(
      "  PRE queda fuera de la comparacion (sin verosimilitud). POST si\n",
      "  entra, con AIC finito y marcado como singular: el aviso no lo\n",
      "  invalida, pero conviene tenerlo presente al leer la tabla.\n"))
  }
}

# ============================================================
# 2. Comparacion de modelos (AIC y pesos de Akaike)
# ============================================================
# Un modelo sin verosimilitud no admite comparacion: glmmTMB deja logLik en
# NULL y AIC() devuelve NA. Meterlo daria pesos de Akaike renormalizados sobre
# un subconjunto, que es exactamente el modo sutil de que un modelo fallido
# desaparezca del papel sin dejar rastro. Aqui se filtra y se dice en voz alta
# que modelo se ha dejado fuera y por que.
#
# El filtro es `AIC` finito, no "sin avisos": un ajuste singular pero con
# verosimilitud (m.ref_int.post) entra en la tabla marcado como singular.
es_utilizable <- function(mod) {
  if (is.null(mod)) return(FALSE)
  a <- tryCatch(stats::AIC(mod), error = function(e) NA_real_)
  is.finite(a)
}

comparar_fase <- function(fase, ...) {
  mods <- list(...)
  ok   <- vapply(mods, es_utilizable, logical(1))
  if (any(!ok)) {
    cat(sprintf("\n  [%s] AVISO: %s no tienen AIC y quedan fuera de la comparacion.\n",
                fase, paste(names(mods)[!ok], collapse = ", ")))
    cat(sprintf("  [%s] Los pesos de Akaike de abajo suman 1 SOLO entre los modelos listados.\n",
                fase))
  }
  if (!any(ok)) {
    cat(sprintf("\n  [%s] Ningun modelo de la fase es utilizable.\n", fase))
    return(NULL)
  }
  cmp <- performance::compare_performance(mods[ok]) |>
    as.data.frame() |>
    dplyr::mutate(
      DeltaAIC   = AIC - min(AIC, na.rm = TRUE),
      PesoAkaike = exp(-0.5 * DeltaAIC) / sum(exp(-0.5 * DeltaAIC), na.rm = TRUE),
      fase       = fase
    ) |>
    dplyr::arrange(DeltaAIC)
  cmp
}

cmp_pre  <- comparar_fase("PRE (t < 94 h)",
                          m.ref_naive.pre = m.ref_naive.pre,
                          m.ref_int.pre   = m.ref_int.pre,
                          m.ref.pre       = m.ref.pre)
cmp_post <- comparar_fase("POST (t > 94 h)",
                          m.ref_naive.post = m.ref_naive.post,
                          m.ref_int.post   = m.ref_int.post,
                          m.ref.post       = m.ref.post)

cmp_all <- bind_rows(cmp_pre, cmp_post)
write.csv(cmp_all, "00-data/reference_model_comparison.csv", row.names = FALSE)
cat("\nGuardada: 00-data/reference_model_comparison.csv\n\n")
print(cmp_all)

# ============================================================
# 3. Aviso sobre la pendiente por especie
# ============================================================
# La referencia lleva `(0 + time_s | species)`, asi que la heterogeneidad que
# este termino puede incorporar es la PENDIENTE. Merece la pena vigilarla porque
# el ajuste por si solo casi no la premia: en PRE sale AIC 1.47 por debajo de
# `m.ref_naive` y aun asi la varianza estimada es 3.38 (SD ~ 1.84). O sea, el
# termino esta haciendo su trabajo y a la vez el criterio de informacion no lo
# distingue de no tenerlo. Conviene saberlo antes de escribir el articulo.
#
# Se lee la diagonal con diag() porque si una varianza colapsa a cero glmmTMB
# devuelve una matriz 1x1 sin rownames y pedir ["Intercept", 1] aborta con
# "subindice fuera de los limites" (mismo criterio que varcomp_table()).
avisar_pendiente_especie <- function(mod, etiqueta) {
  vc <- glmmTMB::VarCorr(mod)$cond
  if (!("species" %in% names(vc))) {
    cat(sprintf("  [%s] AVISO: no hay termino aleatorio para `species`.\n", etiqueta))
    return(invisible(NULL))
  }
  v <- diag(vc[["species"]])
  # Con `(0 + time_s | species)` hay una sola componente: la varianza de la
  # pendiente. Si se mantuviera el intercepto seria la segunda de la diagonal.
  sd_pend <- sqrt(abs(v[1]))
  sd_res  <- sqrt(glmmTMB::sigma(mod)^2)
  cat(sprintf("  [%s] sd(species)[time_s] = %.4f   (sigma = %.4f, ratio = %.2f)\n",
              etiqueta, sd_pend, sd_res, sd_pend / sd_res))
  # El umbral es RELATIVO a la sigma, no absoluto: con 8 especies un valor
  # pequeno en unidades absolutas puede ser una fraccion grande del ruido, y al
  # reves. Se avisa cuando la heterogeneidad entre especies es despreciable
  # frente al error residual.
  if (is.finite(sd_pend) && sd_res > 0 && (sd_pend / sd_res) < 0.5) {
    cat(sprintf(paste0(
      "  [%s] AVISO: la pendiente por especie aporta menos de la mitad de la\n",
      "  desviacion residual. Las figuras marginales de d.05.5 resumiran mal\n",
      "  la heterogeneidad de d.05.1. Comparar con m.ref_naive antes de\n",
      "  interpretar las cifras del articulo.\n"),
      etiqueta))
  } else {
    cat(sprintf(
      "  [%s] La pendiente por especie es apreciable frente al error residual.\n",
      etiqueta))
  }
  invisible(NULL)
}

cat("\n-- Pendiente por especie (heterogeneidad) --\n")
avisar_pendiente_especie(m.ref.pre,  "PRE")
avisar_pendiente_especie(m.ref.post, "POST")

# ============================================================
# 4. Exportacion de resultados
# ============================================================
modelos_ref <- list(
  "m.ref.pre"        = m.ref.pre,
  "m.ref.post"       = m.ref.post,
  "m.ref_int.pre"    = m.ref_int.pre,
  "m.ref_int.post"   = m.ref_int.post,
  "m.ref_naive.pre"  = m.ref_naive.pre,
  "m.ref_naive.post" = m.ref_naive.post
)

coef_ref <- coef_table(modelos_ref)
write.csv(coef_ref, "00-data/reference_model_coef.csv", row.names = FALSE)
cat("Guardada: 00-data/reference_model_coef.csv\n")

varcomp_ref <- varcomp_table(modelos_ref)
write.csv(varcomp_ref, "00-data/reference_model_varcomp.csv", row.names = FALSE)
cat("Guardada: 00-data/reference_model_varcomp.csv\n")

# varcomp_table() exporta `varianza = sum(diag(vc[[nivel]]))`, o sea el trazo de
# la matriz de covarianza. Para un grupo con pendiente e intercepto eso mezcla
# dos varianzas en un numero que no se puede leer como "varianza de la
# pendiente". Se anade una tabla en formato largo, una fila por componente.
componentes_icc <- function(modelos) {
  # Nombres de las componentes de un grupo, deducidos de la formula. Hace falta
  # porque un grupo con una sola componente, como `species` con
  # `(0 + time_s | species)`, produce una matriz 1x1 cuyos rownames no son
  # fiables; sin esto su varianza se perderia en la exportacion.
  nombres_componentes <- function(mod, grupo) {
    fallback <- function(k) paste0("comp", seq_len(k))
    tr <- tryCatch(terms(formula(mod)), error = function(e) NULL)
    if (is.null(tr)) return(NULL)
    for (et in attr(tr, "term.labels")) {
      partes <- strsplit(et, "|", fixed = TRUE)
      if (length(partes) != 2) next
      lhs <- trimws(partes[[1]]); rhs <- trimws(partes[[2]])
      # el factor de agrupacion puede venir como `a/b`; se compara el ultimo
      # segmento, que es el grupo de este VarCorr
      rhs_seg <- trimws(strsplit(rhs, "/", fixed = TRUE)[[1]])
      if (length(rhs_seg) == 0 || tail(rhs_seg, 1) != grupo) next
      vars <- trimws(strsplit(lhs, "+", fixed = TRUE)[[1]])
      vars <- vars[nzchar(vars) & vars != "0"]
      if (!"1" %in% vars) vars <- c("Intercept", vars)
      return(vars)
    }
    NULL
  }

  filas <- list()
  for (nm in names(modelos)) {
    mod <- modelos[[nm]]
    if (!es_utilizable(mod)) next
    vc <- glmmTMB::VarCorr(mod)$cond
    if (is.null(vc)) next
    s2 <- glmmTMB::sigma(mod)^2

    # varianzas de todas las componentes del modelo, para el peso relativo
    todas <- unlist(lapply(vc, function(m) sum(diag(m))))
    total <- sum(todas, na.rm = TRUE) + s2

    for (g in names(vc)) {
      m   <- vc[[g]]
      d   <- diag(m)
      rn  <- rownames(m)
      comp <- if (!is.null(rn) && !anyNA(rn) && all(nzchar(rn))) {
        rn
      } else {
        nc <- nombres_componentes(mod, g)
        if (is.null(nc) || length(nc) != length(d)) paste0("comp", seq_along(d)) else nc
      }
      if (length(comp) != length(d)) comp <- paste0("comp", seq_along(d))
      for (i in seq_along(d)) {
        filas[[length(filas) + 1]] <- data.frame(
          modelo      = nm,
          nivel       = g,
          componente  = comp[i],
          varianza    = d[i],
          sd          = sqrt(abs(d[i])),
          sigma2      = s2,
          peso_total  = if (total > 0) d[i] / total else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, filas)
}

comp_ref <- componentes_icc(modelos_ref)
write.csv(comp_ref, "00-data/reference_model_varcomp_detalle.csv", row.names = FALSE)
cat("Guardada: 00-data/reference_model_varcomp_detalle.csv\n")
cat("\nComponentes por nivel (varianza, sd y peso sobre el total del modelo):\n")
print(comp_ref)

# El Anova solo se calcula si m.ref ha convergido; con una Hessian no definida
# positiva car::Anova() falla o devuelve cifras sin sentido, que es peor que no
# devolver nada. El ajuste previo ya detiene el script en ese caso, asi que
# aqui solo queda el intento por si glmmTMB deviviera un AIC finito con
# advertencia.
if (es_utilizable(m.ref.pre) && es_utilizable(m.ref.post)) {
  anova_ref <- bind_rows(
    car::Anova(m.ref.pre,  type = "III") |> as.data.frame() |>
      rownames_to_column("term") |> mutate(phase = "PRE"),
    car::Anova(m.ref.post, type = "III") |> as.data.frame() |>
      rownames_to_column("term") |> mutate(phase = "POST")
  )
  write.csv(anova_ref, "00-data/reference_model_anova.csv", row.names = FALSE)
  cat("Guardada: 00-data/reference_model_anova.csv\n")
} else {
  cat("Omitido reference_model_anova.csv: m.ref no ha convergido.\n")
}

save_models(modelos_ref)

# ============================================================
# 5. Diagnosticos
# ============================================================
# Solo sobre los modelos utilizables. Hacer diagnostico de residuos de un
# ajuste que no ha convergido produce imagenes bonitas de nada, y ademas
# varios de estos helpers fallan al pedir la verosimilitud.
diagnostico <- modelos_ref[vapply(modelos_ref, es_utilizable, logical(1))]
if (length(diagnostico) < length(modelos_ref)) {
  cat(sprintf("\n  Diagnosticos omitidos para: %s\n",
              paste(setdiff(names(modelos_ref), names(diagnostico)), collapse = ", ")))
}
dir.create("07-img", showWarnings = FALSE, recursive = TRUE)
plot_check_model(diagnostico, sufijo = "ref")
plot_obs_fitted(diagnostico,  sufijo = "ref")
plot_dharma(diagnostico,         sufijo = "ref")
plot_varcomp(varcomp_ref, "07-img/varcomp_reference_model.png",
             "Varianza por nivel - modelo de referencia")

# ============================================================
# 6. Dashboards (opcional)
# ============================================================
if (GEN_DASH == "Y") {
  dir.create("06-html", showWarnings = FALSE)
  for (nm in names(diagnostico)) {
    tryCatch(
      easystats::model_dashboard(diagnostico[[nm]],
                                 output_dir = "06-html/",
                                 output_file = paste0("modeldashboard_", nm, ".html")),
      error = function(e)
        warning("model_dashboard omitido (", nm, "): ", conditionMessage(e))
    )
  }
}

# ============================================================
# 7. Resumen a pantalla
# ============================================================
cat("\n################ Modelo de referencia: efectos fijos ################\n")
print(coef_ref |> filter(modelo %in% c("m.ref.pre", "m.ref.post")))

cat("\n################ Modelo de referencia: varianza por nivel ################\n")
print(varcomp_ref |> filter(modelo %in% c("m.ref.pre", "m.ref.post")))

cat("\n################ Todos los modelos: varianza por nivel ################\n")
print(varcomp_ref)

cat("\n===== d.05.3.reference_model.R completado =====\n")
cat("Las figuras de efecto marginal sobre la pendiente se generan en\n")
cat("01-scripts/d.05.5.slope_effects_figures.R a partir de estos modelos.\n")
