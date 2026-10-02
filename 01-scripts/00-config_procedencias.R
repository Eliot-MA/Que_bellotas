# 00-config_procedencias.R ---- declaracion unica de las procedencias excluidas
#
# Afecta SOLO a la linea d.* (desecacion). NO debe cargarse desde la linea
# s.* (germinacion): ahi se utilizo una sola procedencia por especie, de modo
# que no existe la necesidad de distinguir procedencias dentro de especie y la
# exclusion a nivel de especie (lote RO en s.02/s.03) no se gestiona aqui.
#
# Motivo de la exclusion de IL3 (Quercus ilex, lote "El Serranillo"):
#   1. Desequilibrio muestral en la fase previa. IL3 tiene solo 2 tiempos de
#      muestreo por debajo de 94 h (t0 = 0 h y t1 ~ 23,8 h) frente a 5 en el
#      resto de procedencias: 60 observaciones PRE frente a 150. La pendiente
#      PRE de cada bellota IL3 se estima con 2 puntos, sin grados de libertad.
#   2. Metadatos de procedencia incompletos. En 00-data/rawd/procedencias_updated.csv
#      la fila de IL3 reutiliza la cabecera como datos
#      ("IL3;Q. ilex;El Serranillo;Localidad;Procedencia;Latitud;..."), por lo
#      que d.01.1 (paste(Procedencia, Localidad, sep = "-")) produce
#      procedencia = "Procedencia-Localidad": sin localidad, sin coordenadas
#      y sin altitud. No puede entrar en analisis espaciales ni bioclimaticos
#      ni describirse en el manuscrito.
#
# Decision adoptada: la procedencia se excluye de TODO el estudio, incluida la
# tabla de rasgos (d.02) y el calculo del FAMD (d.03.1). Los ejes del FAMD se
# recalculan sobre el conjunto real de procedencias (15 codigos, 403 bellotas).
#
# Se aplica ANTES de la particion temporal para que PRE y POST compartan el
# mismo conjunto de procedencias y las comparaciones pre-post sean validas.

PROCEDENCIAS_EXCLUIDAS <- "IL3"

# ---- Helper de verificacion ------------------------------------------------
# Falla ruidosamente si alguna procedencia excluida sobrevive en un data.frame.
# Se llama en cada punto del pipeline donde el filtro debe haber surtido efecto,
# para que un olvido no pase inadvertido hasta la redaccion.
#   col: columna que contiene el codigo de procedencia ("prov_code")
assert_sin_procedencias_excluidas <- function(df, col = "prov_code", contexto = NULL) {
  etiqueta <- if (is.null(contexto)) deparse(substitute(df)) else contexto
  if (!col %in% names(df))
    stop("[", etiqueta, "] no existe la columna '", col, "'")
  presentes <- intersect(unique(as.character(df[[col]])), PROCEDENCIAS_EXCLUIDAS)
  if (length(presentes) > 0)
    stop("[", etiqueta, "] se han encontrado procedencias excluidas: ",
         paste(presentes, collapse = ", "))
  invisible(TRUE)
}

# Resumen de composicion por procedencia, para dejar traza en consola.
reportar_composicion_procedencias <- function(df, col = "prov_code", etiqueta = "df") {
  cod <- as.character(df[[col]])
  tb <- as.data.frame(table(prov_code = cod), stringsAsFactors = FALSE)
  names(tb) <- c(col, "n_observaciones")
  tb$n_observaciones <- as.integer(tb$n_observaciones)
  tb$n_bellotas <- vapply(tb[[col]], function(cc)
    length(unique(df[["id_bellota"]][cod == cc])), integer(1))
  cat("\n[", etiqueta, "] composicion por procedencia (",
      nrow(tb), " codigos, ", sum(tb$n_observaciones), " observaciones)\n", sep = "")
  print(tb, row.names = FALSE)
  invisible(tb)
}
