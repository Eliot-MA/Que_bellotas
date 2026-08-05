# Load data ----
library(tidyverse)
df.bellotas <- read.csv("00-data/desiccation_traits_long.csv")
df.famd <- read.csv("00-data/famd_ind_coord.csv")

# Creat working dataframe ----
df <- df.bellotas |> 
  dplyr::select(-X) |>
  dplyr::select(id_bellota, codigo, tiempo_acumulado_horas, Moisture_content) |> 
  left_join(y = df.famd, by = "id_bellota") |> 
  drop_na(Dim.1, Dim.2, Dim.3) |> 
  rename(time = tiempo_acumulado_horas) |> 
  mutate(time = scale(time))

# Separate data in two datasets
df.t1 <- df |> 
  filter(time < 94)

df.t2 <- df |> 
  filter(time > 94)

# Fit models
## Pre 94 hours
library(glmmTMB)
library(parameters)
library(performance)

## Quiero testar structuras de efectos aleatorios: 

mm.pre.1 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (1|species) + (time|provenance) + (time|id_bellota), data = df.t1)

mm.pre.2 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (time||species) + (time|provenance) + (time|id_bellota), data = df.t1)

compare_parameters(
  mm.pre.1, mm.pre.2,
  ci = 0.95,
  select = "{estimate}{stars} [{ci}]", standardize = TRUE
)

compare_performance(mm.pre.1, mm.pre.2)

## Post 94 hours
mm.post.1 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (1|species) + (time|provenance) + (time|id_bellota), data = df.t2)

mm.post.2 <- glmmTMB(Moisture_content ~ time * (Dim.1 + Dim.2 + Dim.3) + (time||species) + (time|provenance) + (time|id_bellota), data = df.t2)

compare_parameters(
  mm.post.1, mm.post.2,
  ci = 0.95,
  select = "{estimate}{stars} [{ci}]", standardize = TRUE
)

check_model(mm.post.2)

# ============================================================
# Export de resultados para el articulo / material suplementario
# ============================================================
source("01-scripts/00-export_helpers.R")

# Objetos de modelo (todas las variantes comparadas)
save_models(list("mm.pre.1" = mm.pre.1, "mm.pre.2" = mm.pre.2,
                 "mm.post.1" = mm.post.1, "mm.post.2" = mm.post.2))

# Comparacion de estructuras de efectos aleatorios (compare_performance)
cmp_rasgos <- dplyr::bind_rows(
  performance::compare_performance(mm.pre.1, mm.pre.2) |>
    as.data.frame() |>
    dplyr::mutate(dataset = "pre (t<94h)"),
  performance::compare_performance(mm.post.1, mm.post.2) |>
    as.data.frame() |>
    dplyr::mutate(dataset = "post (t>94h)")
)
write.csv(cmp_rasgos, "00-data/model_comparisons_traits.csv", row.names = FALSE)
cat("Comparacion de modelos de rasgos guardada en 00-data/model_comparisons_traits.csv\n")

# Tabla de coeficientes de los modelos finales de rasgos (.2)
coef_rasgos <- coef_table(list("mm.pre.2 (t<94h)" = mm.pre.2, "mm.post.2 (t>94h)" = mm.post.2))
write.csv(coef_rasgos, "00-data/coef_traits_models.csv", row.names = FALSE)
cat("Tabla de coeficientes guardada en 00-data/coef_traits_models.csv\n")

# Componentes de varianza + barplot por nivel
varcomp_rasgos <- varcomp_table(list("mm.pre.2 (t<94h)" = mm.pre.2, "mm.post.2 (t>94h)" = mm.post.2))
write.csv(varcomp_rasgos, "00-data/varcomp_traits_models.csv", row.names = FALSE)
cat("Componentes de varianza guardados en 00-data/varcomp_traits_models.csv\n")
plot_varcomp(varcomp_rasgos, "07-img/varcomp_traits_models.png",
             "Varianza por nivel - modelos de rasgos")

# Diagnosticos de los modelos finales de rasgos (.2)
modelos_rasgos <- list("mm.pre.2" = mm.pre.2, "mm.post.2" = mm.post.2)
plot_check_model(modelos_rasgos, sufijo = "traits")
plot_obs_fitted(modelos_rasgos, sufijo = "traits")
plot_dharma(modelos_rasgos, sufijo = "traits")

# Efectos marginales de time * (Dim.1 + Dim.2 + Dim.3) (pre y post)
suppressPackageStartupMessages(library(ggeffects))
suppressPackageStartupMessages(library(patchwork))
for (nm in names(modelos_rasgos)) {
  mod <- modelos_rasgos[[nm]]
  pl <- lapply(c("Dim.1", "Dim.2", "Dim.3"), function(dimv) {
    pred <- tryCatch(
      ggeffects::ggpredict(mod, terms = c("time", paste0(dimv, " [-1, 0, 1]"))),
      error = function(e) e
    )
    if (inherits(pred, "error")) {
      warning("ggpredict fallo para ", dimv, " en ", nm, ": ", conditionMessage(pred))
      return(NULL)
    }
    plot(pred) + labs(title = paste("time *", dimv))
  })
  pl <- Filter(Negate(is.null), pl)
  if (length(pl)) {
    fig <- Reduce(`+`, pl) +
      patchwork::plot_annotation(title = paste("MC ~ time * (Dim.1+Dim.2+Dim.3) -", nm))
    ggsave(file.path("07-img", paste0("marginal_", nm, ".png")), fig, width = 14, height = 5)
    cat("Efectos marginales guardados en 07-img/marginal_", nm, ".png\n", sep = "")
  }
}
