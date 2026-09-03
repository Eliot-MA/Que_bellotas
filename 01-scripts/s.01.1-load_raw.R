###############################################################################
## s.01.1 - Load raw data from the desiccation sensitivity experiment        ##
##                                                                           ##
## Reads the raw CSV exports of both experimental phases, reshapes them from ##
## wide format (one column block per sampling time) to long format, joins    ##
## species and provenance information by acorn ID ranges, and recodes the    ##
## germination records into a clean 0/1 encoding.                            ##
##                                                                           ##
## Inputs  (00-data/, read-only):                                            ##
##   - germination_phase_I.csv     Phase 1 records (acorn IDs 1-1600)        ##
##   - germination_phase_II.csv    Phase 2 records (local IDs, shifted       ##
##                                 +1600 to match the provenance lookup)     ##
##   - id_species_provenances.csv  Species/provenance lookup by ID range     ##
##                                                                           ##
## Output (in memory; exported to disk only by the master script):           ##
##   - df.germ.records             One row per acorn x sampling time         ##
##                                                                           ##
## Units: weights in grams. Datetimes in local time (Europe/Madrid).         ##
###############################################################################

# --- 0. Libraries -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(fuzzyjoin)

# --- 1. Load raw data -------------------------------------------------------

# All three files are semicolon-separated with comma decimals (Spanish locale)
# and carry a UTF-8 BOM.

rD.germ.phase1 <- read.csv2(
  "00-data/germination_phase_I.csv",
  fileEncoding = "UTF-8-BOM",
  check.names  = FALSE,
  stringsAsFactors = FALSE
)

rD.germ.phase2 <- read.csv2(
  "00-data/germination_phase_II.csv",
  fileEncoding = "UTF-8-BOM",
  check.names  = FALSE,
  stringsAsFactors = FALSE
)

rD.id.species <- read.csv2(
  "00-data/id_species_provenances.csv",
  fileEncoding = "UTF-8-BOM",
  check.names  = FALSE,
  stringsAsFactors = FALSE
)

# --- 2. Helper functions ----------------------------------------------------

# Clean column names:
#   - trim surrounding whitespace (e.g. "peso seco ", "Germinada ")
#   - the two trailing "Observaciones" columns are positionally renamed:
#     first = dry-weight observations, second = germination observations
#   - fix the "Germianda" typo present in the phase 2 export
clean_names <- function(df) {
  names(df) <- trimws(names(df))

  obs_idx <- which(names(df) == "Observaciones")
  if (length(obs_idx) == 2) {
    names(df)[obs_idx[1]] <- "observations_dry_weight"
    names(df)[obs_idx[2]] <- "observations_germination"
  }

  names(df)[names(df) == "Germianda"] <- "Germinada"

  df
}

# Reshape one phase from wide to long:
#   - drop metadata columns superseded by the provenance lookup
#   - stack the per-time blocks (fecha/peso/hora/minuto/observaciones tN)
#   - drop rows without a date (empty filler rows created by the time blocks)
#   - tag the phase and shift acorn IDs of phase 2 (+offset) so that IDs are
#     unique across phases and fall inside the provenance lookup ranges
standardise_phase <- function(df, phase_id, id_offset = 0L) {
  df |>
    clean_names() |>
    dplyr::select(-any_of(c("especie", "procedencia", "codigo", "numero_bellota"))) |>
    rename(acorn_id = id_bellota) |>
    pivot_longer(
      cols          = matches("^(fecha|peso|hora|minuto|observaciones) t\\d+$"),
      names_to      = c(".value", "time"),
      names_pattern = "(fecha|peso|hora|minuto|observaciones) t(\\d+)"
    ) |>
    mutate(across(where(is.character), ~ na_if(.x, ""))) |>
    filter(!is.na(fecha)) |>
    mutate(
      phase     = phase_id,
      acorn_id  = as.integer(acorn_id) + id_offset,
      time      = as.character(time),
      datetime  = make_datetime(
        year  = year(dmy(fecha)),
        month = month(dmy(fecha)),
        day   = day(dmy(fecha)),
        hour  = as.integer(hora),
        min   = as.integer(minuto)
      )
    ) |>
    rename(fresh_weight = peso, dry_weight = `peso seco`) |>
    dplyr::select(-fecha, -hora, -minuto)
}

# --- 3. Standardise and combine phases --------------------------------------

# Phase 2 offset: derived from the provenance lookup instead of hardcoded
# (first acorn ID of phase 2 minus one).
phase2_offset <- min(rD.id.species$`Id-bellota-p`[rD.id.species$Fase == 2]) - 1L

df.germ.records <- bind_rows(
  standardise_phase(rD.germ.phase1, phase_id = 1L),
  standardise_phase(rD.germ.phase2, phase_id = 2L, id_offset = phase2_offset)
)

# Safety check: acorn IDs must be unique across phases
stopifnot(
  !any(duplicated(
    df.germ.records |> distinct(acorn_id, phase) |> pull(acorn_id)
  ))
)

# --- 4. Join species and provenance by acorn ID ranges -----------------------

df.germ.records <- df.germ.records |>
  fuzzy_left_join(
    rD.id.species,
    by = c(
      "acorn_id" = "Id-bellota-p",
      "acorn_id" = "Id-bellota-f",
      "phase"    = "Fase"
    ),
    match_fun = list(`>=`, `<=`, `==`)
  ) |>
  rename(species = Especie, provenance = Procedencia) |>
  dplyr::select(-`Id-bellota-p`, -`Id-bellota-f`, -Fase) |>
  relocate(species, provenance, .after = phase)

# --- 5. Recode germination records -------------------------------------------

# 5.1. "Moho" (mould) is not a germination state: move it to the germination
# observations and set the record to missing.
n_moho <- df.germ.records |>
  filter(Germinada %in% "Moho" | Emergida %in% "Moho") |>
  nrow()

df.germ.records <- df.germ.records |>
  mutate(
    moho = Germinada %in% "Moho" | Emergida %in% "Moho",
    observations_germination = case_when(
      moho & is.na(observations_germination) ~ "Moho",
      moho ~ paste(observations_germination, "Moho", sep = "; "),
      TRUE  ~ observations_germination
    ),
    Germinada = if_else(Germinada %in% "Moho", NA_character_, Germinada),
    Emergida  = if_else(Emergida  %in% "Moho", NA_character_, Emergida)
  ) |>
  dplyr::select(-moho)

# 5.2. Binary encoding: "SI" -> 1, "NO" or missing -> 0.
#      Blank cells mean the acorn was simply not recorded as germinated
#      (emerged) on the field sheets, so they are genuine zeros.
df.germ.records <- df.germ.records |>
  mutate(
    Germinada = case_when(Germinada %in% "SI" ~ "1",
                          Germinada %in% "NO" ~ "0",
                          is.na(Germinada)    ~ "0",
                          TRUE ~ Germinada),
    Emergida  = case_when(Emergida %in% "SI" ~ "1",
                          Emergida %in% "NO" ~ "0",
                          is.na(Emergida)    ~ "0",
                          TRUE ~ Emergida)
  )

# 5.3. Ambiguity fix.
#      Acorns with exactly two records (t0 + one later time) have an
#      unreliable germination/emergence entry at t0: on the field sheets the
#      cross could belong to either visit. Set those t0 entries to missing;
#      they are recovered later from the last record of the acorn.
ambiguity_mask <- df.germ.records |>
  group_by(acorn_id) |>
  summarise(
    n_records = n(),
    n_t0      = sum(time == "0"),
    .groups   = "drop"
  )

n_ambiguous <- df.germ.records |>
  left_join(ambiguity_mask, by = "acorn_id") |>
  filter(n_records == 2, n_t0 == 1, time == "0") |>
  nrow()

df.germ.records <- df.germ.records |>
  left_join(ambiguity_mask, by = "acorn_id") |>
  mutate(
    Germinada = if_else(n_records == 2 & n_t0 == 1 & time == "0",
                        NA_character_, Germinada),
    Emergida  = if_else(n_records == 2 & n_t0 == 1 & time == "0",
                        NA_character_, Emergida)
  ) |>
  dplyr::select(-n_records, -n_t0)

# 5.4. Final numeric types and column order
df.germ.records <- df.germ.records |>
  mutate(
    germinated = as.numeric(Germinada),
    emerged    = as.numeric(Emergida)
  ) |>
  dplyr::select(
    acorn_id, time, phase, species, provenance, datetime,
    fresh_weight, dry_weight, germinated, emerged,
    observations_dry_weight, observations_germination
  ) |>
  arrange(acorn_id, as.integer(time))

# --- 6. Quality-control summary ---------------------------------------------

cat("=== s.01.1 load_raw ===\n")
cat("Phase 1:", nrow(df.germ.records[df.germ.records$phase == 1, ]), "records,",
    df.germ.records |> filter(phase == 1) |> distinct(acorn_id) |> nrow(), "acorns\n")
cat("Phase 2:", nrow(df.germ.records[df.germ.records$phase == 2, ]), "records,",
    df.germ.records |> filter(phase == 2) |> distinct(acorn_id) |> nrow(), "acorns",
    "(ID offset applied: +", phase2_offset, ")\n")

sp_summary <- df.germ.records |>
  group_by(species) |>
  summarise(n_acorns = n_distinct(acorn_id), .groups = "drop")
cat("\nAcorns per species:\n")
print(sp_summary, row.names = FALSE)

cat("\n'Moho' entries moved to observations:", n_moho, "\n")
cat("Ambiguous t0 germination entries set to NA:", n_ambiguous, "\n")

n_unmatched <- sum(is.na(df.germ.records$species))
if (n_unmatched > 0) {
  warning(n_unmatched, " records could not be matched to a species/provenance.",
          " Check id_species_provenances.csv ranges.")
} else {
  cat("All records matched to a species/provenance.\n")
}

missing_fw <- sum(is.na(df.germ.records$fresh_weight))
cat("Records with missing fresh weight:", missing_fw, "\n")


