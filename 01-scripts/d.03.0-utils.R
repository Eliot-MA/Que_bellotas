# d.03.0-utils.R ---- constantes y utilidades compartidas por los hijos de d.03
# Cargado desde d.03.1, d.03.2 y d.03.3 (guardas idempotentes).

# Dimensiones FAMD que se interpretan en el manuscrito
dims <- c("Dim.1", "Dim.2", "Dim.3")

# Clasificacion bioclimatica de cada especie
bioclimate <- tibble::tribble(
  ~species,            ~bioclimate,
  "Quercus coccifera", "Mediterranean",
  "Quercus ilex",      "Mediterranean",
  "Quercus suber",     "Mediterranean",
  "Quercus faginea",   "Sub-Mediterranean",
  "Quercus pyrenaica", "Sub-Mediterranean",
  "Quercus pubescens", "Sub-Mediterranean",
  "Quercus petraea",   "Temperate",
  "Quercus robur",     "Temperate"
)

# Paleta Okabe-Ito: vermilion/orange/blue
bioclimate_cols <- c(
  "Mediterranean"     = "#D55E00",  # Okabe-Ito vermilion
  "Sub-Mediterranean" = "#E69F00",  # Okabe-Ito orange
  "Temperate"         = "#0072B2"   # Okabe-Ito blue
)

okabe_blue   <- "#0072B2"
okabe_orange <- "#E69F00"