library(FactoMineR)  # Para PCA y FAMD
library(factoextra)  # Para visualización
library(gridExtra)   # Para dividir espacio de graficación
library(ggrepel)     # para etiquetas no sobrepuestas (opcional)
library(dplyr)

df.orig <- df

# Variables de interés
vars_interes <- c("id_bellota", "especie", "procedencia", "codigo",
                  "Moisture_content", "tiempo_acumulado_horas",
                  "peso_seco", "Volumen_estimado_cm3", "Area_estimada_cm2", "Relacion_SV", # Tamaño
                  "SPM_g_cm2", "Seed_Coat_Ratio", # Robustez de pericarpo
                  "Ratio_A.cicatriz_A.bellota",   # Tamaño de la cicatriz
                  "rajas_pericarpo")   

df <- df %>%
  drop_na(all_of(vars_interes)) %>%           # Los NA dan error en el PCA
  filter(cotiledon_anormal %in% c(0, 1)) |>   # Bellotas sin defectos anómalos
  dplyr::select(all_of(vars_interes))

# Eliminar pseudorreplicación (una observación por bellota)
df.unique <- df |> dplyr::distinct(id_bellota, .keep_all = TRUE)

# Asegurar rajas como factor
df.unique$rajas_pericarpo <- factor(df.unique$rajas_pericarpo)

##################
# FAMD con rajas #
#   se quitan    #
# rajas tipo 2   #
##################

# df de solo rasgos y CON variable rajas
df.rasgos <- df.unique |> 
  select(-id_bellota, -Moisture_content, -especie, -procedencia, -codigo)

# Quitar rajas tipo 2, son demasiado pocas
df.rasgos <- df.rasgos |> 
  filter(rajas_pericarpo != 2)
df.rasgos$rajas_pericarpo <- factor(df.rasgos$rajas_pericarpo) # Reprogramar variable rajas para que solo tenga dos niveles

# FAMD
famd_rasgos <- FAMD(df.rasgos, ncp = 5, graph = FALSE)

# Unir dimensiones del famd al dataset original
# 1. Extraemos las coordenadas de las primeras 3 dimensiones
coords_famd <- as.data.frame(famd_rasgos$ind$coord[, 1:3])

# 2. Necesitamos el ID para unir. 
coords_famd$id_bellota <- df.unique |> filter(rajas_pericarpo != 2) |> pull(id_bellota)

# 3. Unimos las dimensiones a la tabla grande (df)
df_famd <- df %>%
  filter(rajas_pericarpo != 2) |> 
  left_join(coords_famd, by = "id_bellota")

# Escalar variables
df_famd_scaled <- df_famd %>%
  mutate(across(
    .cols = -c(id_bellota, especie, procedencia, codigo, Moisture_content, rajas_pericarpo, starts_with("Dim.")), 
    .fns = ~ as.vector(scale(.))
  )) |> 
  mutate(rajas_pericarpo = factor(rajas_pericarpo))


