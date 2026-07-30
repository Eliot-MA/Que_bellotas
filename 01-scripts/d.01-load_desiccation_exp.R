source("01-scripts/d.01.1-load_data.R")
source("01-scripts/d.01.2-derived_variables.R")
write.csv(x = df.bellotas, "00-data/desiccation_traits_long.csv")
