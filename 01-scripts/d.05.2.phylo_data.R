# ============================================================
# d.05.2.phylo_data.R
# Construccion de la filogenia de las 8 especies de Quercus.
#
# Fuente principal: arbol CROWN de Hipp et al. (2020) "Global oak phylogeny of
# the genus Quercus" (Nature Plants 6: 1113-1119) - repositorio
# andrew-hipp/global-oaks-2019. Es un filogenoma de Quercus (RAD-seq, ~functional
# por miles de marcadores) datado con ~8 fosiles. Comparado con OToL (que carece
# de longitudes de rama) y con GBOTB (que injerta ad hoc especies ausentes),
# el crown tree es la via mas informativa: todos los tips de nuestras 8
# especies existen y con datum.
#
# Genera y guarda:
#   - arbol crown descargado y PODADO a las 8 especies (sin injertos)
#   - matriz de covarianza filogenetica A = vcv.phylo
#   - verificacion de coincidencia de nombres (tips crown vs spp del analisis)
#   - via descartada V.PhyloMaker2/GBOTB (escenarios S1-S3) solo como cotejo
#   - metricas comparativas para justificar (a) la eleccion del crown tree
#     frente a OToL/GBOTB y (b) el uso distancias de Grafen si procediera
#
# Salidas:
#   00-data/phylo/oak_crown.tre                   arbol crown crudo (descarga)
#   00-data/phylo/tip_check_summary.csv           coincidencia de nombres
#   00-data/phylo/oak_tree.rds | .nwk            arbol crown podado (principal)
#   00-data/phylo/oak_vcv.rds                    matriz A
#   00-data/phylo/gbotb_scenarios.rds            arboles S1-S3 + estado de especies
#   00-data/phylo/tree_comparison_summary.csv    metricas de justificacion
#   07-img/oak_phylo.png
#   07-img/phylomaker_scenarios.png
#   07-img/phylo_trees_comparison.png
# ============================================================

library(ape)

spp <- c("Quercus petraea", "Quercus robur", "Quercus faginea",
         "Quercus coccifera", "Quercus ilex", "Quercus pyrenaica",
         "Quercus suber", "Quercus pubescens")

dir.create("00-data/phylo", showWarnings = FALSE, recursive = TRUE)
dir.create("07-img", showWarnings = FALSE, recursive = TRUE)

RECREATE_PHYLO <- FALSE   # TRUE para re-descargar y regenerar todo
DEMO_GBOTB     <- TRUE    # FALSE para no ejecutar V.PhyloMaker2 (tarda ~1.5 min)
RECREATE_OTOL  <- FALSE   # TRUE para generar el arbol OToL de cotejo (requiere rotl)

# ============================================================
# 1. Arbol CROWN de Hipp et al. (2020) — fuente principal
# ============================================================
# El crown tree es un filogenoma de Quercus (RAD-seq) datado con ~8 fosiles.
# Su nodo raiz (crown age) marca el inicio de la diversificacion de Quercus,
# lo cual es lo relevante para un analisis comparativo de 8 especies del
# genero. No se necesita el stem tree ni longitudes sinteticas de Grafen.

CROWN_URL <- "https://raw.githubusercontent.com/andrew-hipp/global-oaks-2019/master/ANALYSES/2019-06_globalOaks-gitUpdate/OUT/ANALYSIS.PRODUCTS/tr.singletons.correlated.1.taxaGrepCrown.tre"

if (!file.exists("00-data/phylo/oak_tree.rds") || RECREATE_PHYLO) {

  # 1a. Descargar el crown tree (si no existe o se pide re-descarga)
  crown_dest <- "00-data/phylo/oak_crown.tre"
  if (!file.exists(crown_dest) || RECREATE_PHYLO) {
    download.file(CROWN_URL, crown_dest, mode = "wb")
    cat("Crown tree descargado de Hipp et al. (2020)\n")
  }

  # 1b. Cargar con ape
  crown_full <- ape::read.tree(crown_dest)
  cat("Crown tree cargado:", length(crown_full$tip.label), "tips\n")
  cat("Branch lengths:", !is.null(crown_full$edge.length), "\n")
  cat("Ultrametric:", ape::is.ultrametric(crown_full), "\n")

  # 1c. Normalizar labels del crown tree
  # Los tips vienen como "Quercus_petraea|NA|NA|1982.0337|QUE001513"
  # (genero_especie|pais|prov|coleccion|codigo). Necesitamos "Quercus petraea"
  # Normalizar a "Genero especie" (2 palabras, primera mayuscula, resto minuscula)
  # para que coincidan EXACTAMENTE con nuestro vector spp.
  # e.g. Quercus_petraea|NA|NA|1982.0337|QUE001513 -> "Quercus petraea"
  normalize_tip <- function(x) {
    sp  <- sub("\\|.*$", "", x)        # quitar todo despues del primer |
    sp  <- gsub("_", " ", sp)          # guion bajo -> espacio
    parts <- strsplit(trimws(sp), "\\s+")[[1]]
    if (length(parts) < 2) return(sp)
    genus   <- paste0(toupper(substr(parts[1], 1, 1)),
                      tolower(substr(parts[1], 2, nchar(parts[1]))))
    epithet <- tolower(parts[-1])
    paste(c(genus, epithet), collapse = " ")
  }

  tip_sp <- sapply(crown_full$tip.label, normalize_tip, USE.NAMES = FALSE)
  names(tip_sp) <- crown_full$tip.label

  # 1d. Verificar que las 8 especies estan presentes
  cat("-- Verificacion de nombres de tips (crown tree) --\n")
  check_df <- data.frame(
    expected   = spp,
    n_tips     = sapply(spp, function(sp) sum(tip_sp == sp)),
    stringsAsFactors = FALSE
  )
  print(check_df)
  missing <- spp[check_df$n_tips == 0]
  if (length(missing) > 0) {
    stop("Especies no encontradas en crown tree: ",
         paste(missing, collapse = ", "))
  }

  # Guardar verificacion como CSV
  write.csv(check_df, "00-data/phylo/tip_check_summary.csv", row.names = FALSE)

  # 1e. Elegir 1 tip por especie: preferir la accesion con metadatos completos
  #     (label mas largo = mas pipes). Robur e ilex tienen 2 accesiones.
  keep_tips <- character(length(spp))
  names(keep_tips) <- spp
  for (sp in spp) {
    labels <- crown_full$tip.label[tip_sp == sp]
    keep_tips[sp] <- labels[which.max(nchar(labels))]
  }
  cat("\nTips seleccionados (1 por especie):\n")
  print(keep_tips)

  # 1f. Podar a las 8 especies y renombrar
  tree <- ape::keep.tip(crown_full, unname(keep_tips))
  tree$tip.label <- spp[match(tree$tip.label, unname(keep_tips))]

  cat("\nArbol crown podado a 8 especies:\n")
  cat("Tips:", paste(tree$tip.label, collapse = ", "), "\n")
  cat("Branch lengths:", !is.null(tree$edge.length), "\n")
  cat("Ultrametric:", ape::is.ultrametric(tree), "\n")

  ape::write.tree(tree, file = "00-data/phylo/oak_tree.nwk")
  saveRDS(tree, "00-data/phylo/oak_tree.rds")

} else {
  cat("Cargando arbol crown ya construido\n")
  tree <- readRDS("00-data/phylo/oak_tree.rds")
}

# Matriz de varianza-covarianza filogenetica
A <- ape::vcv.phylo(tree)
saveRDS(A, "00-data/phylo/oak_vcv.rds")

# Figura del arbol crown
png("07-img/oak_phylo.png", width = 1200, height = 800, res = 150)
plot(tree, cex = 1.1)
title("Quercus crown tree — Hipp et al. (2020)")
dev.off()
cat("Guardado: 00-data/phylo/oak_tree.rds, oak_tree.nwk, oak_vcv.rds, 07-img/oak_phylo.png\n")

# ============================================================
# 2. Cotejo: OToL + Grafen (opcional, requiere rotl)
# ============================================================
# Via secundaria de cotejo documentada en el metodo: subarbol inducido de
# Open Tree of Life con longitudes de Grafen. El crown tree es el principal;
# este bloque solo verifica que la topologia OToL es coherente con la del
# crown tree. Habilitar con RECREATE_OTOL <- TRUE (tarda ~30 s). Si rotl no
# esta instalado, se salta sin error.
if ((!file.exists("00-data/phylo/oak_otol.rds") || RECREATE_OTOL) &&
    requireNamespace("rotl", quietly = TRUE)) {

  suppressPackageStartupMessages(library(rotl))

  # 2a. Resolucion de nombres en la taxonomia de OToL
  resolved <- rotl::tnrs_match_names(spp, context_name = "Vascular plants")
  saveRDS(resolved, "00-data/phylo/otol_resolution.rds")

  cat("-- Resolucion de nombres en OToL (cotejo) --\n")
  print(resolved[, c("search_string", "unique_name", "ott_id", "is_synonym", "flags")])
  stopifnot(!any(is.na(resolved$ott_id)))   # las 8 especies deben resolverse

  # 2b. Subarbol inducido
  tree_otol <- rotl::tol_induced_subtree(ott_ids = resolved$ott_id[!is.na(resolved$ott_id)])

  # 2c. Limpiar labels y hacerlos coincidir con spp
  strip <- function(x) gsub("_ott\\d+$", "", x)
  tip_clean <- strip(tree_otol$tip.label)
  map <- setNames(resolved$unique_name, as.character(resolved$ott_id))
  ott_from_tip <- gsub(".*_ott(\\d+)$", "\\1", tree_otol$tip.label)
  tip_species <- map[ott_from_tip]
  tree_otol$tip.label <- unname(ifelse(is.na(tip_species), tip_clean, tip_species))

  missing <- setdiff(spp, tree_otol$tip.label)
  stopifnot(length(missing) == 0)

  # 2d. OToL no trae longitudes de rama: derivacion de Grafen
  if (is.null(tree_otol$edge.length)) {
    cat("OToL sin branch lengths; aplicando Grafen (compute.brlen)\n")
    tree_otol <- ape::compute.brlen(tree_otol, method = "Grafen")
  }

  ape::write.tree(tree_otol, file = "00-data/phylo/oak_otol.nwk")
  saveRDS(tree_otol, "00-data/phylo/oak_otol.rds")

} else if (file.exists("00-data/phylo/oak_otol.rds")) {
  cat("Cargando arbol OToL de cotejo ya construido\n")
  tree_otol <- readRDS("00-data/phylo/oak_otol.rds")
} else {
  cat("OToL de cotejo omitido (rotl no disponible o RECREATE_OTOL = FALSE)\n")
  tree_otol <- NULL
}

# ============================================================
# 3. Via descartada: V.PhyloMaker2 / GBOTB (escenarios S1-S3)
# ============================================================
if (!file.exists("00-data/phylo/gbotb_scenarios.rds") || DEMO_GBOTB) {
  if (!requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
    cat("V.PhyloMaker2 no instalado; omitiendo cotejo GBOTB\n")
    res_sc <- list(S1 = NULL, S2 = NULL, S3 = NULL)
    status <- data.frame(species = spp, status = NA_character_)
  } else {
  suppressPackageStartupMessages(library(V.PhyloMaker2))

  sp_list <- data.frame(
    species = spp,
    genus   = "Quercus",
    family  = "Fagaceae",
    stringsAsFactors = FALSE
  )

  run_demo <- function(sc) {
    cat("\n-- escenario", sc, "--\n")
    res <- tryCatch(
      V.PhyloMaker2::phylo.maker(sp.list = sp_list,
                                 tree = GBOTB.extended.TPL,
                                 scenarios = sc),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      cat("ERROR:", conditionMessage(res), "\n")
      return(invisible(NULL))
    }
    phy <- res[[1]]
    cat("tips:", length(phy$tip.label),
        "| branch lengths:", !is.null(phy$edge.length), "\n")
    invisible(res)
  }

  res_sc <- lapply(c("S1", "S2", "S3"), run_demo)
  names(res_sc) <- c("S1", "S2", "S3")

  if (is.null(res_sc$S3)) {
    warning("Escenario S3 fallo; no se puede guardar el estado de especies.")
    status <- data.frame(species = spp, status = NA_character_)
  } else {
    status <- res_sc$S3$species.list[, c("species", "status")]
    print(status)
  }

  saveRDS(list(
    trees  = lapply(res_sc, function(r) if (is.null(r)) NULL else r[[1]]),
    status = status
  ), "00-data/phylo/gbotb_scenarios.rds")

  # Figura: como resuelve cada escenario el arbol de Quercus
  png("07-img/phylomaker_scenarios.png", width = 2100, height = 900, res = 150)
  par(mfrow = c(1, 3), mar = c(2, 1, 3, 1))
  for (sc in names(res_sc)) {
    if (is.null(res_sc[[sc]])) next
    phy <- res_sc[[sc]][[1]]
    tip_col <- ifelse(gsub("_", " ", phy$tip.label) %in%
                        status$species[status$status == "bind"],
                      "firebrick", "black")
    plot(phy, cex = 1.1, tip.color = tip_col, label.offset = 0.1)
    title(paste("V.PhyloMaker2 - escenario", sc))
  }
  mtext("rojo = especies no presentes en GBOTB (injertadas ad hoc)",
        side = 1, line = -1, outer = TRUE)
  dev.off()
  cat("Guardado: 00-data/phylo/gbotb_scenarios.rds y 07-img/phylomaker_scenarios.png\n")
  }

} else {
  gbotb <- readRDS("00-data/phylo/gbotb_scenarios.rds")
  res_sc <- gbotb$trees
  status <- gbotb$status
}

# Normalizar a una lista de arboles (el run fresco devuelve el objeto completo
# de phylo.maker con el arbol en [[1]]; la cache guarda ya solo los arboles).
trees_gbotb <- lapply(res_sc, function(r) {
  if (is.null(r)) return(NULL)
  if (inherits(r, "phylo")) r else r[[1]]
})

# ============================================================
# 4. Comparacion de arboles y justificacion de la eleccion
# ============================================================
# El arbol crown (Hipp et al. 2020) es la referencia. Se compara contra
# GBOTB (V.PhyloMaker2, escenarios S1-S3) y, si existe, OToL+Grafen para
# justificar que el crown tree no necesita injertar ninguna especie.

tidy_labels <- function(phy) {
  phy$tip.label <- gsub("_", " ", phy$tip.label)
  phy
}
tree_crown_ref <- tidy_labels(tree)
gbotb_ok <- !is.null(trees_gbotb) && !all(sapply(trees_gbotb, is.null))
otol_ok  <- !is.null(tree_otol)

cands <- list()
if (otol_ok) cands[["OToL_Grafen"]] <- tidy_labels(tree_otol)
if (gbotb_ok) {
  for (sc in names(trees_gbotb)) {
    phy <- trees_gbotb[[sc]]
    if (!is.null(phy)) cands[[sc]] <- tidy_labels(phy)
  }
}

if (length(cands) > 0) {

  bind_sp <- if (exists("status")) status$species[status$status == "bind"] else character(0)

  rf <- sapply(cands, function(phy) {
    tryCatch(ape::dist.topo(tree_crown_ref, phy, method = "PH85"),
             error = function(e) NA_real_)
  })
  coph <- sapply(cands, function(phy) {
    tryCatch(cor(as.vector(stats::cophenetic(tree_crown_ref)),
                 as.vector(stats::cophenetic(phy))),
             error = function(e) NA_real_)
  })

  tree_summary <- tibble::tibble(
    arbol                 = names(cands),
    tips                  = sapply(cands, function(phy) length(phy$tip.label)),
    branch_lengths        = sapply(cands, function(phy) !is.null(phy$edge.length)),
    ultrametric           = sapply(cands, function(phy) isTRUE(ape::is.ultrametric(phy))),
    RF_distance_crown     = unname(rf),
    cophenetic_cor_crown  = unname(coph)
  )

  tree_summary <- tibble::add_row(
    tree_summary,
    arbol                 = "Crown_Hipp2020",
    tips                  = length(tree_crown_ref$tip.label),
    branch_lengths        = !is.null(tree_crown_ref$edge.length),
    ultrametric           = isTRUE(ape::is.ultrametric(tree_crown_ref)),
    RF_distance_crown     = 0,
    cophenetic_cor_crown  = 1
  )

  write.csv(tree_summary, "00-data/phylo/tree_comparison_summary.csv", row.names = FALSE)
  cat("Guardada: 00-data/phylo/tree_comparison_summary.csv\n")
  print(tree_summary)

  # Figura comparativa: crown + (OToL opcional) + escenarios GBOTB
  n_panels <- 1 + length(cands)
  png("07-img/phylo_trees_comparison.png",
      width = 400 * n_panels, height = 900, res = 150)
  par(mfrow = c(1, n_panels), mar = c(2, 1, 3, 1))
  plot(tree_crown_ref, cex = 1.1)
  title("Crown (Hipp 2020)")
  for (nm in names(cands)) {
    phy <- cands[[nm]]
    bind_pos <- vapply(phy$tip.label, function(t) any(t %in% bind_sp), logical(1))
    tip_col <- ifelse(bind_pos, "firebrick", "black")
    plot(phy, cex = 1.1, tip.color = tip_col, label.offset = 0.1)
    title(nm)
  }
  mtext("rojo = especies injertadas ad hoc en GBOTB", side = 1, line = -1, outer = TRUE)
  dev.off()
  cat("Figura guardada: 07-img/phylo_trees_comparison.png\n")

  cat("\nJustificacion: el crown tree (Hipp et al. 2020) resuelve las 8 especies\n")
  cat("sin injertos, con longitudes de rama reales datadas con fosiles.\n")
  if (gbotb_ok) {
    cat("GBOTB injerta ad hoc las 3 especies ausentes (robur, ilex, pubescens),\n")
    cat("que cambian de posicion segun el escenario S1/S2/S3.\n")
  }
}

cat("\n===== d.05.2.phylo_data.R completado =====\n")
