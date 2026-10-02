if (!exists("bioclimate_cols")) source("01-scripts/d.03.0-utils.R")
if (!exists("PROCEDENCIAS_EXCLUIDAS")) source("01-scripts/00-config_procedencias.R")

library(tidyverse)

# Si no tienes desiccation_traits_long.csv
# ejecuta primero VVV
# source("01-scripts/d.01-load_desiccation_exp.R")
# rm(list = ls())

df.bellotas <- read.csv("00-data/processed/desiccation_traits_long.csv")

# 1. Calculate FAMD ----
## 1.1. Create famd dataframe ----
df <- df.bellotas |>
  filter(cotiledon_anormal %in% c(0, 1)) |>
  filter(rajas_pericarpo %in% c(0, 1)) |> 
  dplyr::select(id_bellota, especie, procedencia, prov_code,
                peso_seco, Volumen_estimado_cm3, Relacion_SV, 
                SPM_g_cm2, Seed_Coat_Ratio, 
                Ratio_A.cicatriz_A.bellota, 
                rajas_pericarpo) |> 
  rename(species = especie, 
         provenance = procedencia, 
         dry_weight = peso_seco, 
         volume_cm3 = Volumen_estimado_cm3,
         SVR = Relacion_SV,
         SCR = Seed_Coat_Ratio,
         SSR = Ratio_A.cicatriz_A.bellota,
         pericarp_rupture = rajas_pericarpo
  ) |> 
  unique() |> 
  # Procedencias excluidas ANTES de calcular el FAMD: los ejes deben reflejar
  # el conjunto real de procedencias del estudio (15 codigos), no uno que
  # incluye una procedencia que luego se descarta de los modelos.
  filter(!prov_code %in% PROCEDENCIAS_EXCLUIDAS) |>
  drop_na() |> 
  mutate(pericarp_rupture = as.factor(pericarp_rupture))

assert_sin_procedencias_excluidas(df, "prov_code", "FAMD (matriz de entrada)")
reportar_composicion_procedencias(df, "prov_code", "FAMD (matriz de entrada)")

df |> 
  group_by(species, provenance) |> 
  summarise(
    n = n(), 
    prop = n/30
  )

## 1.1b. Correlograma de variables incluidas en el FAMD ----
library(corrplot)

df.corr <- df |>
  dplyr::select(dry_weight, volume_cm3, SVR, SPM_g_cm2, 
                SCR, SSR, pericarp_rupture) |>
  mutate(pericarp_rupture = as.numeric(as.character(pericarp_rupture)))

mat.corr <- cor(df.corr, method = "pearson", use = "pairwise.complete.obs")

rownames(mat.corr) <- colnames(mat.corr) <- c("mass", "volume", "SVR", 
                                                "SPM", "SCR", "SSR", 
                                                "pericarp rupture")

write.csv2(mat.corr, "00-data/processed/correlation_matrix_traits.csv")

dir.create("07-img/FAMD_outputs", showWarnings = FALSE)

png("07-img/FAMD_outputs/00_correlogram_traits.png", width = 800, height = 700)
corrplot(mat.corr, 
         method = "color", 
         type = "upper", 
         addCoef.col = "black",
         tl.col = "black", 
         tl.srt = 45,
         number.cex = 0.8,
         mar = c(0, 0, 2, 0))
dev.off()
cat("Correlograma guardado en 07-img/FAMD_outputs/00_correlogram_traits.png\n")

## 1.2. Calculate and save famd ----

library(FactoMineR)  # Para FAMD
library(factoextra)


# Exclude id and species variables
df.famd <- df |> dplyr::select(-id_bellota, -species, -provenance, -prov_code)

# Calculate famd
# Los ejes se calculan sobre `df`, que ya viene filtrado por
# PROCEDENCIAS_EXCLUIDAS: el espacio de rasgos refleja el conjunto real de
# procedencias del estudio.
famd.traits <- FAMD(df.famd, graph = FALSE)

# Coordenadas individuales: se seleccionan por nombre y se renombran Dim.1..Dim.5
# para depender solo del contenido y no del formato de $ind de FactoMineR.
coord_ind <- as.data.frame(famd.traits$ind$coord)
coord_ind <- coord_ind[, grep("^Dim\\.", colnames(coord_ind)), drop = FALSE]
coord_ind <- coord_ind[, 1:5, drop = FALSE]
names(coord_ind) <- paste0("Dim.", 1:5)

df <- cbind(df, coord_ind)

assert_sin_procedencias_excluidas(df, "prov_code", "famd_ind_coord.csv")
# row.names = FALSE: por defecto write.csv escribe los nombres de fila como
# primera columna sin cabecera, que al releer el CSV aparece como "H1".
write.csv(x = df, "00-data/processed/famd_ind_coord.csv", row.names = FALSE)
cat("Coordenadas FAMD guardadas en 00-data/processed/famd_ind_coord.csv (",
    nrow(df), " bellotas, ", length(PROCEDENCIAS_EXCLUIDAS) ,
    " procedencia(s) excluida(s): ",
    paste(PROCEDENCIAS_EXCLUIDAS, collapse = ", "), ")\n", sep = "")


### Save famd info ----
# If not created, create famd folder
# dir.create("07-img/FAMD_outputs", showWarnings = FALSE)

library(patchwork)

# ---------------------------
# 1. Eigenvalues / Scree plot
# ---------------------------
p1 <- fviz_screeplot(famd.traits, addlabels = TRUE, ylim = c(0, 50))

ggsave("07-img/FAMD_outputs/01_scree_plot.png", p1, width = 7, height = 5)


# ---------------------------
# 2. Individuos (PC1 vs PC2)
# ---------------------------
p2.1 <- fviz_famd_ind(
  famd.traits,
  axes = c(1, 2),
  label = "none",
  col.ind = "cos2",
  gradient.cols = c("grey", "blue", "red"),
  title = "Individuals (colored by cos2)"
)

p2.2 <- fviz_famd_ind(
  famd.traits,
  axes = c(1, 3),
  label = "none",
  col.ind = "cos2",
  gradient.cols = c("grey", "blue", "red"),
  title = "Individuals (colored by cos2)"
) + labs(title = "")

p2 <- p2.1 + p2.2 + plot_layout(guides = "collect")

ggsave("07-img/FAMD_outputs/02_individuals_cos2.png", p2, width = 7, height = 6)

# ---------------------------
# 3. Variables cuantitativas y cualitativas
# ---------------------------
p3.1 <- fviz_famd_var(
  axes = c(1, 2),
  famd.traits,
  choice = "var",
  col.var = "contrib",
  gradient.cols = c("grey", "blue", "red"),
  repel = TRUE
)

p3.2 <- fviz_famd_var(
  axes = c(1, 3),
  famd.traits,
  choice = "var",
  col.var = "contrib",
  gradient.cols = c("grey", "blue", "red"),
  repel = TRUE
)

p3 <- p3.1 + p3.2

ggsave("07-img/FAMD_outputs/03_variables.png", p3, width = 7, height = 6)

# 2. Create FAMD table ----

# --- Variables cuantitativas ---
quanti_df <- famd.traits$quanti.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd.traits$quanti.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd.traits$quanti.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cuantitativa")

# --- Variables cualitativas ---
quali_df <- famd.traits$quali.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to = "dimension", values_to = "coord") %>%
  left_join(
    famd.traits$quali.var$contrib %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "contrib"),
    by = c("variable", "dimension")
  ) %>%
  left_join(
    famd.traits$quali.var$cos2 %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "dimension", values_to = "cos2"),
    by = c("variable", "dimension")
  ) %>%
  mutate(tipo = "cualitativa")

# --- Combinar todo ---
tabla_famd <- bind_rows(quanti_df, quali_df)

tabla_famd

write.csv2(x = tabla_famd, file = "00-data/processed/famd_long.csv")

# 2) Limpiar y ordenar la tabla
# Construir la tabla (variables ├ù dimensiones con contribuciones)
# formato ancho

# Paso 1: quedarte solo con contribuciones + coordenadas
tabla_contrib <- tabla_famd %>%
  dplyr::select(variable, dimension, contrib, coord)

# Paso 2: codificar el signo
# Aqu├¡ tienes varias opciones. La m├ís limpia es a├▒adir el signo como string:

tabla_contrib <- tabla_contrib %>%
  mutate(
    # Signo en ASCII plano: el guion unicode estaba corrupto en el codigo
    # fuente y write.csv2 lo exportaba como mojibake en los CSV de salida.
    signo = ifelse(coord >= 0, "+", "-"),
    contrib_signo = paste0(round(contrib, 1), signo)
  )

# Paso 3: pasar a formato ancho
tabla_wide <- tabla_contrib %>%
  dplyr::select(variable, dimension, contrib_signo) %>%
  pivot_wider(names_from = dimension, values_from = contrib_signo)

write.csv2(x = tabla_wide, file = "00-data/processed/paper_famd.csv")


## 3. Tabla de individuos con bioclimate (compartida con d.03.2 y d.03.3) ----
# prov_code viaja aqui porque d.03.3 modela la procedencia como unidad aleatoria
# (ver la nota de decision al inicio de d.03.3).
ind.bio <- df |>
  dplyr::select(id_bellota, species, prov_code, Dim.1, Dim.2, Dim.3) |>
  left_join(bioclimate, by = "species") |>
  mutate(bioclimate = factor(bioclimate,
                             levels = c("Mediterranean",
                                        "Sub-Mediterranean",
                                        "Temperate")))

assert_sin_procedencias_excluidas(ind.bio, "prov_code", "ind.bio (tabla de individuos)")

