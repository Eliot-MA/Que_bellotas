# ============================================================
# d.05.2.phylo_data.R
# Construccion de la filogenia de las 8 especies de Quercus.
#
# Genera y guarda:
#   - arbol de Open Tree of Life (subarbol inducido) con longitudes de Grafen
#   - matriz de covarianza filogenetica A = vcv.phylo
#   - arboles de la via descartada V.PhyloMaker2/GBOTB (escenarios S1-S3)
#   - metricas comparativas para justificar (a) la eleccion de OToL frente a
#     GBOTB y (b) el uso de distancias de Grafen
#
# Salidas:
#   00-data/phylo/otol_resolution.rds            resolucion de nombres en OToL
#   00-data/phylo/oak_tree.rds | .nwk            arbol OToL + Grafen
#   00-data/phylo/oak_vcv.rds                    matriz A
#   00-data/phylo/gbotb_scenarios.rds            arboles S1-S3 + estado de especies
#   00-data/phylo/tree_comparison_summary.csv    metricas de justificacion
#   07-img/oak_phylo.png
#   07-img/phylomaker_scenarios.png
#   07-img/phylo_trees_comparison.png
# ============================================================

library(ape)
library(rotl)

spp <- c("Quercus petraea", "Quercus robur", "Quercus faginea",
         "Quercus coccifera", "Quercus ilex", "Quercus pyrenaica",
         "Quercus suber", "Quercus pubescens")

dir.create("00-data/phylo", showWarnings = FALSE, recursive = TRUE)
dir.create("07-img", showWarnings = FALSE, recursive = TRUE)

RECREATE_PHYLO <- FALSE   # TRUE para re-descargar de OToL
DEMO_GBOTB     <- TRUE    # FALSE para no ejecutar V.PhyloMaker2 (tarda ~1.5 min)

# ============================================================
# 1. Arbol de OToL + longitudes de Grafen
# ============================================================
if (!file.exists("00-data/phylo/oak_tree.rds") || RECREATE_PHYLO) {

  # 1a. Resolucion de nombres en la taxonomia de OToL
  resolved <- rotl::tnrs_match_names(spp, context_name = "Vascular plants")
  saveRDS(resolved, "00-data/phylo/otol_resolution.rds")

  cat("-- Resolucion de nombres en OToL --\n")
  print(resolved[, c("search_string", "unique_name", "ott_id", "is_synonym", "flags")])
  stopifnot(!any(is.na(resolved$ott_id)))   # las 8 especies deben resolverse

  # 1b. Subarbol inducido
  tree <- rotl::tol_induced_subtree(ott_ids = resolved$ott_id[!is.na(resolved$ott_id)])

  # 1c. Limpiar labels y hacerlos coincidir con los nombres originales
  strip <- function(x) gsub("_ott\\d+$", "", x)
  tip_clean <- strip(tree$tip.label)
  map <- setNames(resolved$unique_name, as.character(resolved$ott_id))
  ott_from_tip <- gsub(".*_ott(\\d+)$", "\\1", tree$tip.label)
  tip_species <- map[ott_from_tip]
  tree$tip.label <- unname(ifelse(is.na(tip_species), tip_clean, tip_species))

  missing <- setdiff(spp, tree$tip.label)
  stopifnot(length(missing) == 0)

  # 1d. JUSTIFICACION DE GRAFEN: el subarbol inducido de OToL no trae
  #     longitudes de rama; Grafen las deriva a partir de la profundidad de
  #     los nodos (ultrametrico y estandar en ausencia de longitudes).
  otol_had_lengths <- !is.null(tree$edge.length)
  cat("OToL inducido con branch lengths:", otol_had_lengths, "\n")
  if (!otol_had_lengths) {
    cat("Aplicando Grafen (compute.brlen) para derivar longitudes de rama\n")
    tree <- ape::compute.brlen(tree, method = "Grafen")
  }

  ape::write.tree(tree, file = "00-data/phylo/oak_tree.nwk")
  saveRDS(tree, "00-data/phylo/oak_tree.rds")

} else {
  cat("Cargando arbol OToL ya construido\n")
  tree <- readRDS("00-data/phylo/oak_tree.rds")
  otol_had_lengths <- FALSE
}

# Matriz de varianza-covarianza filogenetica
A <- ape::vcv.phylo(tree)
saveRDS(A, "00-data/phylo/oak_vcv.rds")

# Figura del arbol OToL
png("07-img/oak_phylo.png", width = 1200, height = 800, res = 150)
plot(tree, cex = 1.1)
title("Open Tree of Life - Quercus subtree (Grafen branch lengths)")
dev.off()
cat("Guardado: 00-data/phylo/oak_tree.rds, oak_tree.nwk, oak_vcv.rds y 07-img/oak_phylo.png\n")

# ============================================================
# 2. Via descartada: V.PhyloMaker2 / GBOTB (escenarios S1-S3)
# ============================================================
if (!file.exists("00-data/phylo/gbotb_scenarios.rds") || DEMO_GBOTB) {
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
# 3. Comparacion de arboles y justificacion de la eleccion
# ============================================================
tidy_labels <- function(phy) {
  phy$tip.label <- gsub("_", " ", phy$tip.label)
  phy
}
tree_otol <- tidy_labels(tree)
gbotb_ok <- !is.null(trees_gbotb) && !all(sapply(trees_gbotb, is.null))

if (gbotb_ok) {
  cands <- lapply(trees_gbotb, function(phy) {
    if (is.null(phy)) return(NULL)
    tidy_labels(phy)
  })

  rf <- sapply(cands, function(phy) {
    if (is.null(phy)) return(NA_real_)
    tryCatch(ape::dist.topo(tree_otol, phy, method = "PH85"),
             error = function(e) NA_real_)
  })
  coph <- sapply(cands, function(phy) {
    if (is.null(phy)) return(NA_real_)
    tryCatch(cor(as.vector(ape::cophenetic(tree_otol)),
                 as.vector(ape::cophenetic(phy))),
             error = function(e) NA_real_)
  })

  tree_summary <- tibble::tibble(
    arbol              = names(cands),
    tips               = sapply(cands, function(phy) if (is.null(phy)) NA_integer_ else length(phy$tip.label)),
    branch_lengths     = sapply(cands, function(phy) if (is.null(phy)) NA else !is.null(phy$edge.length)),
    ultrametric        = sapply(cands, function(phy) if (is.null(phy)) NA else isTRUE(ape::is.ultrametric(phy))),
    RF_distance_OToL   = unname(rf),
    cophenetic_cor_OToL = unname(coph)
  )

  tree_summary <- tibble::add_row(
    tree_summary,
    arbol = "OToL_Grafen",
    tips = length(tree_otol$tip.label),
    branch_lengths = !is.null(tree_otol$edge.length),
    ultrametric = isTRUE(ape::is.ultrametric(tree_otol)),
    RF_distance_OToL = 0,
    cophenetic_cor_OToL = 1
  )

  # Justificacion de Grafen: el subarbol inducido de OToL carecia de longitudes
  tree_summary <- tibble::add_column(
    tree_summary,
    otol_without_branch_lengths = ifelse(tree_summary$arbol == "OToL_Grafen",
                                         otol_had_lengths == FALSE, NA),
    .after = "branch_lengths"
  )

  write.csv(tree_summary, "00-data/phylo/tree_comparison_summary.csv", row.names = FALSE)
  cat("Guardada: 00-data/phylo/tree_comparison_summary.csv\n")
  print(tree_summary)

  # Figura comparativa de los 4 arboles candidatos
  png("07-img/phylo_trees_comparison.png", width = 2600, height = 900, res = 150)
  par(mfrow = c(1, 4), mar = c(2, 1, 3, 1))
  plot(tree_otol, cex = 1.1)
  title("OToL + Grafen")
  for (sc in names(cands)) {
    phy <- cands[[sc]]
    if (is.null(phy)) next
    tip_col <- ifelse(phy$tip.label %in% status$species[status$status == "bind"],
                      "firebrick", "black")
    plot(phy, cex = 1.1, tip.color = tip_col, label.offset = 0.1)
    title(paste("V.PhyloMaker2 -", sc))
  }
  mtext("rojo = especies injertadas ad hoc en GBOTB", side = 1, line = -1, outer = TRUE)
  dev.off()
  cat("Figura guardada: 07-img/phylo_trees_comparison.png\n")

  cat("\nJustificacion: la via GBOTB no lanza error, pero injerta ad hoc las\n")
  cat("3 especies ausentes (robur, ilex, pubescens), que cambian de posicion\n")
  cat("segun el escenario S1/S2/S3. OToL resuelve las 8 especies sin injertos\n")
  cat("y, al carecer de longitudes de rama, se completan con Grafen.\n")
}

cat("\n===== d.05.2.phylo_data.R completado =====\n")
