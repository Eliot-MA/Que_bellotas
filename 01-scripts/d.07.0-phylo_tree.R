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
