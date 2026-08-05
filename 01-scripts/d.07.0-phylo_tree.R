# d.07.0 - Construccion de la filogenia de las 8 especies de Quercus
# Estrategia: resolver nombres en Open Tree of Life y extraer el subarbol inducido
library(rotl)
library(ape)

spp <- c(
  "Quercus petraea", "Quercus robur", "Quercus faginea",
  "Quercus coccifera", "Quercus ilex", "Quercus pyrenaica",
  "Quercus suber", "Quercus pubescens"
)

# 1) Resolver nombres en la taxonomia de OToL
resolved <- tnrs_match_names(spp, context_name = "Vascular plants")
print(resolved[, c("search_string", "unique_name", "ott_id", "is_synonym", "flags")])

cat("\n-- especies NO resueltas --\n")
print(resolved |> dplyr::filter(is.na(ott_id)))

# 2) Subarbol inducido
tree <- tol_induced_subtree(ott_ids = resolved$ott_id[!is.na(resolved$ott_id)])
cat("\n-- tips del arbol inducido --\n")
print(tree$tip.label)

# 3) Limpiar labels y hacer coincidir con los nombres originales
strip <- function(x) gsub("_ott\\d+$", "", x)
tip_clean <- strip(tree$tip.label)

# mapa ott -> nombre de especie original
map <- setNames(resolved$unique_name, as.character(resolved$ott_id))
# labels con ott_id
ott_from_tip <- gsub(".*_ott(\\d+)$", "\\1", tree$tip.label)
tip_species <- map[ott_from_tip]
tree$tip.label <- unname(ifelse(is.na(tip_species), tip_clean, tip_species))

cat("\n-- tips tras mapear a nombres de especie --\n")
print(tree$tip.label)

# Comprobar que estan las 8
missing <- setdiff(spp, tree$tip.label)
cat("\nEspecies ausentes:", if (length(missing)) paste(missing, collapse = ", ") else "ninguna", "\n")

# 4) Longitudes de rama: el arbol de OToL inducido no tiene branch lengths.
# Usamos el metodo de Grafen (ultrametrico, estandar en ausencia de longitudes).
if (is.null(tree$edge.length)) {
  cat("\n-- Sin branch lengths: aplicando Grafen (compute.brlen) --\n")
  tree <- ape::compute.brlen(tree, method = "Grafen")
}

# Topologia resultante
cat("\n-- Newick --\n")
cat(ape::write.tree(tree), "\n")

# Guardar
dir.create("00-data/phylo", showWarnings = FALSE, recursive = TRUE)
ape::write.tree(tree, file = "00-data/phylo/oak_tree.nwk")
saveRDS(tree, "00-data/phylo/oak_tree.rds")

cat("\n-- Resumen del arbol --\n")
print(tree)
cat("\n-- Matriz de var-cov filogenetica (A) --\n")
A <- ape::vcv.phylo(tree)
print(A)
saveRDS(A, "00-data/phylo/oak_vcv.rds")

# Plot
png("07-img/oak_phylo.png", width = 1200, height = 800, res = 150)
plot(tree, cex = 1.1)
title("Open Tree of Life - Quercus subtree (Grafen branch lengths)")
dev.off()
cat("\nGuardado en 00-data/phylo/ y 07-img/oak_phylo.png\n")

# ============================================================
# APENDICE - Verificacion de la via descartada (mega-arbol GBOTB)
# Premisa a comprobar: "intentar extraer la filogenia de Quercus por
# la via de V.PhyloMaker2/GBOTB da error". Resultado empirico: NO da
# error; phylo.maker completa y devuelve un arbol con branch lengths.
# El problema es silencioso: las especies ausentes de GBOTB se injertan
# ad hoc ("bind") y la topologia intra-Quercus resultante es arbitraria.
# Por eso se descarto esta via en favor de OToL + Grafen.
# ============================================================
DEMO_GBOTB <- TRUE   # FALSE para omitir este bloque (tarda ~1.5 min)
if (DEMO_GBOTB) {
  cat("\n\n===== VERIFICACION: via GBOTB (V.PhyloMaker2) =====\n")
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
      V.PhyloMaker2::phylo.maker(sp.list = sp_list, tree = GBOTB.extended.TPL, scenarios = sc),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      cat("ERROR:", conditionMessage(res), "\n")
      return(invisible(NULL))
    }
    phy <- res[[1]]
    cat("tips:", length(phy$tip.label),
        "| nodos internos:", phy$Nnode,
        "| con branch lengths:", !is.null(phy$edge.length), "\n")
    cat("newick:", ape::write.tree(phy), "\n")
    invisible(res)
  }

  res_sc <- lapply(c("S1", "S2", "S3"), run_demo)
  names(res_sc) <- c("S1", "S2", "S3")

  cat("\nEstado de cada especie:\n")
  cat("  'prune' = ya estaba en GBOTB (hereda la topologia del mega-arbol)\n")
  cat("  'bind'  = no estaba en GBOTB, se injerta ad hoc\n")
  print(res_sc$S3$species.list[, c("species", "status")])
  bind_spp <- res_sc$S3$species.list$species[res_sc$S3$species.list$status == "bind"]

  # Guardar arboles de los escenarios para reproductibilidad del informe
  saveRDS(list(
    trees  = lapply(res_sc, function(r) r[[1]]),
    status = res_sc$S3$species.list[, c("species", "status")]
  ), "00-data/phylo/gbotb_scenarios.rds")

  # Figura: como resuelve cada escenario el arbol de Quercus
  dir.create("07-img", showWarnings = FALSE, recursive = TRUE)
  png("07-img/phylomaker_scenarios.png", width = 2100, height = 900, res = 150)
  par(mfrow = c(1, 3), mar = c(2, 1, 3, 1))
  for (sc in names(res_sc)) {
    phy <- res_sc[[sc]][[1]]
    tip_col <- ifelse(gsub("_", " ", phy$tip.label) %in% bind_spp, "firebrick", "black")
    plot(phy, cex = 1.1, tip.color = tip_col, label.offset = 0.1)
    title(paste("V.PhyloMaker2 - escenario", sc))
  }
  mtext("rojo = especies no presentes en GBOTB (injertadas ad hoc)", side = 1, line = -1, outer = TRUE)
  dev.off()
  cat("\nFigura guardada en 07-img/phylomaker_scenarios.png\n")

  cat("\nConclusion de la verificacion: la via GBOTB NO lanza error, pero devuelve\n")
  cat("una filogenia intra-Quercus arbitraria: las 3 especies 'bind' (robur, ilex,\n")
  cat("pubescens) cambian de posicion segun el escenario S1/S2/S3 (en S1/S3 quedan\n")
  cat("en politomia en el nodo de Quercus y en S2 se emparejan con un congenere\n")
  cat("aleatorio) y la topologia resultante puede contradecir filogenias de\n")
  cat("referencia (p.ej. Q. suber queda entre los robles blancos). Por ello se\n")
  cat("descarto y se uso OToL + Grafen.\n")
}
