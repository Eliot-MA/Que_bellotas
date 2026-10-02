source("01-scripts/d.01.1-load_data.R")
source("01-scripts/d.01.2-derived_variables.R")

# Frontera de nombres del dato derivado: el crudo llama a esta columna "codigo"
# y "procedencia" es la localidad textual. Aqui se renombra a "prov_code" para
# que todo el pipeline d.* (y las aserciones de 00-config_procedencias.R) use
# un identificador en ingles. "procedencia" se conserva sin tocar.
df.bellotas <- df.bellotas |> rename(prov_code = codigo)

write.csv(x = df.bellotas, "00-data/processed/desiccation_traits_long.csv")
