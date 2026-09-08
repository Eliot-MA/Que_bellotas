# d.03-trait_famd.R ---- orquestador del analisis FAMD
# Reproduce todo el flujo en orden:
#   d.03.0 utils (constantes) -> d.03.1 carga/FAMD -> d.03.2 figuras 1a/1b -> d.03.3 reporte bioclima
# Los analisis de agrupamiento (k-means, silhouette, k-grid, etc.) se retiraron
# del manuscrito; quedan archivados y comentados en 01-scripts/_legacy_clustering.R.

library(tidyverse)

source("01-scripts/d.03.0-utils.R")
source("01-scripts/d.03.1-load_famd.R")
source("01-scripts/d.03.2-biplot_panel.R")
source("01-scripts/d.03.3-bioclimate_report.R")