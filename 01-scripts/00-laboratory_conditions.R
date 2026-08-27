###############################################################################
## Laboratory conditions — HOBO datalogger data                              ##
## Temperature, relative humidity, and VPD throughout experimental phases    ##
###############################################################################

# =============================================================================
# STEP 1: Visual exploration — detect anomalies and label experimental phases
# =============================================================================

# --- 0. Libraries -----------------------------------------------------------

library(tidyverse)

# --- 1. Load data -----------------------------------------------------------

# HOBO CSV exports: hourly records, GMT+01:00
# Columns: row number, datetime, temp (°C), RH (%), status flags

hobo19 <- read.csv("00-data/hobo19.csv", stringsAsFactors = FALSE)
hobo18 <- read.csv("00-data/hobo18.csv", stringsAsFactors = FALSE)

# --- 2. Clean and standardise -----------------------------------------------

clean_hobo <- function(df, logger_id) {
  df |>
    dplyr::select(1:4) |>
    dplyr::rename(
      n      = 1,
      date   = 2,
      temp_c = 3,
      rh_pct = 4
    ) |>
    dplyr::mutate(
      date       = mdy_hms(date, tz = "Europe/Madrid"),
      temp_c     = as.numeric(temp_c),
      rh_pct     = as.numeric(rh_pct),
      datalogger = logger_id
    )
}

df_hobo19 <- clean_hobo(hobo19, "HOBO-21961472")
df_hobo18 <- clean_hobo(hobo18, "HOBO-21984651")

# Merge both dataloggers
df_lab <- bind_rows(df_hobo19, df_hobo18)

# Quick summary
cat("=== Data summary ===\n")
cat("HOBO-21961472:", nrow(df_hobo19), "records\n")
cat("  Date range:", as.character(min(df_hobo19$date, na.rm = TRUE)),
    "to", as.character(max(df_hobo19$date, na.rm = TRUE)), "\n")
cat("  Temp:", range(df_hobo19$temp_c, na.rm = TRUE), "°C\n")
cat("  RH:  ", range(df_hobo19$rh_pct, na.rm = TRUE), "%\n\n")
cat("HOBO-21984651:", nrow(df_hobo18), "records\n")
cat("  Date range:", as.character(min(df_hobo18$date, na.rm = TRUE)),
    "to", as.character(max(df_hobo18$date, na.rm = TRUE)), "\n")
cat("  Temp:", range(df_hobo18$temp_c, na.rm = TRUE), "°C\n")
cat("  RH:  ", range(df_hobo18$rh_pct, na.rm = TRUE), "%\n")

# --- 3. Missing values and basic anomalies ----------------------------------

cat("\n=== Missing values ===\n")
cat("HOBO-21961472 - NA temp:", sum(is.na(df_hobo19$temp_c)),
    "| NA rh:", sum(is.na(df_hobo19$rh_pct)), "\n")
cat("HOBO-21984651 - NA temp:", sum(is.na(df_hobo18$temp_c)),
    "| NA rh:", sum(is.na(df_hobo18$rh_pct)), "\n")

cat("\n=== Possible anomalies ===\n")
cat("HOBO-21961472 - Temp < 0°C:", sum(df_hobo19$temp_c < 0, na.rm = TRUE),
    "| Temp > 40°C:", sum(df_hobo19$temp_c > 40, na.rm = TRUE), "\n")
cat("HOBO-21961472 - RH < 10%:", sum(df_hobo19$rh_pct < 10, na.rm = TRUE),
    "| RH > 99%:", sum(df_hobo19$rh_pct > 99, na.rm = TRUE), "\n")
cat("HOBO-21984651 - Temp < 0°C:", sum(df_hobo18$temp_c < 0, na.rm = TRUE),
    "| Temp > 40°C:", sum(df_hobo18$temp_c > 40, na.rm = TRUE), "\n")
cat("HOBO-21984651 - RH < 10%:", sum(df_hobo18$rh_pct < 10, na.rm = TRUE),
    "| RH > 99%:", sum(df_hobo18$rh_pct > 99, na.rm = TRUE), "\n")

# --- 4. Load experimental phases --------------------------------------------

phases <- read.csv("00-data/experimental_phases.csv", stringsAsFactors = FALSE) |>
  mutate(
    start = as.POSIXct(start, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid"),
    end   = as.POSIXct(end,   format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid")
  )

# Label positioning: desiccation at bottom to avoid overlap with germination_I
phases_labels <- phases |>
  mutate(
    x_pos = start + (end - start) / 2,
    y_pos = ifelse(phase == "desiccation", -Inf, Inf),
    vjust = ifelse(phase == "desiccation", 1.5, 1.5)
  )

# Phases with Date for daily plots
phases_date <- phases |>
  mutate(start_date = as.Date(start), end_date = as.Date(end))

phases_labels_date <- phases |>
  mutate(
    x_pos_date = as.Date(start + (end - start) / 2),
    y_pos      = ifelse(phase == "desiccation", -Inf, Inf),
    vjust      = ifelse(phase == "desiccation", 1.5, 1.5)
  )

# Assign colours for consistent use across plots
phase_colours <- c(
  "desiccation"    = "#E69F00",
  "germination_I"  = "#56B4E9",
  "germination_II" = "#009E73",
  "chamber"        = "#CC79A7"
)

# Assign phase to each record based on time windows
df_lab <- df_lab |>
  mutate(phase = NA_character_)

for (i in seq_len(nrow(phases))) {
  idx <- df_lab$date >= phases$start[i] & df_lab$date <= phases$end[i]
  df_lab$phase[idx] <- phases$phase[i]
}

# =============================================================================
# STEP 2: Calculate Vapor Pressure Deficit (VPD)
# =============================================================================

# VPD = esat * (1 - RH/100)
# esat = 0.6108 * exp(17.27 * T / (T + 237.3))  [kPa, Tetens 1930]

df_lab <- df_lab |>
  mutate(
    esat = 0.6108 * exp((17.27 * temp_c) / (temp_c + 237.3)),
    vpd  = esat * (1 - rh_pct / 100)
  )

cat("\n=== VPD summary ===\n")
cat("HOBO-21961472 - VPD range:", range(df_lab$vpd[df_lab$datalogger == "HOBO-21961472"], na.rm = TRUE), "kPa\n")
cat("HOBO-21984651 - VPD range:", range(df_lab$vpd[df_lab$datalogger == "HOBO-21984651"], na.rm = TRUE), "kPa\n")

# =============================================================================
# PLOTS — three core figures with T, RH, and VPD
# =============================================================================

# Common label for faceted variables
df_long <- df_lab |>
  pivot_longer(
    cols      = c(temp_c, rh_pct, vpd),
    names_to  = "variable",
    values_to = "value"
  ) |>
  mutate(
    variable = factor(variable,
      levels = c("temp_c", "rh_pct", "vpd"),
      labels = c("Temperature (°C)", "Relative humidity (%)", "VPD (kPa)")
    )
  )

# --- Plot 1: Distributions -------------------------------------------------

p_dist <- df_long |>
  ggplot(aes(x = value, fill = datalogger)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  labs(
    title = "Distribution of laboratory environmental conditions",
    x = "Value", y = "Count", fill = "Datalogger"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

ggsave("07-img/laboratory_distribution.png", p_dist,
       width = 10, height = 8, dpi = 150)

# --- Plot 2: Hourly evolution with phase labels ----------------------------

p_hourly <- ggplot() +
  geom_rect(
    data = phases,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = phase),
    alpha = 0.15
  ) +
  geom_line(
    data = df_long,
    aes(x = date, y = value, colour = datalogger),
    alpha = 0.5, linewidth = 0.3
  ) +
  geom_label(
    data = phases_labels |> filter(phase != "desiccation"),
    aes(x = x_pos, y = y_pos, label = description, fill = phase),
    vjust = 1.5, size = 3, fontface = "bold", alpha = 0.8
  ) +
  geom_label(
    data = phases_labels |> filter(phase == "desiccation"),
    aes(x = x_pos, y = y_pos, label = description, fill = phase),
    vjust = -0.5, size = 3, fontface = "bold", alpha = 0.8
  ) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = phase_colours, guide = "none") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title = "Hourly evolution of laboratory environmental conditions",
    x = "Date", y = "Value", colour = "Datalogger"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

ggsave("07-img/laboratory_hourly_evolution.png", p_hourly,
       width = 14, height = 9, dpi = 150)

# --- Plot 3: Daily means with SD ribbon and phase labels -------------------

df_daily <- df_lab |>
  mutate(date_day = as.Date(date)) |>
  group_by(date_day, datalogger) |>
  summarise(
    temp_mean = mean(temp_c, na.rm = TRUE),
    temp_sd   = sd(temp_c, na.rm = TRUE),
    rh_mean   = mean(rh_pct, na.rm = TRUE),
    rh_sd     = sd(rh_pct, na.rm = TRUE),
    vpd_mean  = mean(vpd, na.rm = TRUE),
    vpd_sd    = sd(vpd, na.rm = TRUE),
    n_records = n(),
    .groups   = "drop"
  )

df_daily_long <- df_daily |>
  pivot_longer(
    cols      = c(temp_mean, rh_mean, vpd_mean),
    names_to  = "variable",
    values_to = "mean_val"
  ) |>
  mutate(
    sd_val = case_when(
      variable == "temp_mean" ~ temp_sd,
      variable == "rh_mean"   ~ rh_sd,
      variable == "vpd_mean"  ~ vpd_sd
    ),
    variable = factor(variable,
      levels = c("temp_mean", "rh_mean", "vpd_mean"),
      labels = c("Temperature (°C)", "Relative humidity (%)", "VPD (kPa)")
    )
  )

p_daily <- ggplot() +
  geom_rect(
    data = phases_date,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf, fill = phase),
    alpha = 0.15
  ) +
  geom_ribbon(
    data = df_daily_long,
    aes(x = date_day, ymin = mean_val - sd_val, ymax = mean_val + sd_val,
        fill = datalogger),
    alpha = 0.2
  ) +
  geom_line(
    data = df_daily_long,
    aes(x = date_day, y = mean_val, colour = datalogger),
    linewidth = 0.6
  ) +
  geom_label(
    data = phases_labels_date |> filter(phase != "desiccation"),
    aes(x = x_pos_date, y = y_pos, label = description, fill = phase),
    vjust = 1.5, size = 3, fontface = "bold", alpha = 0.8
  ) +
  geom_label(
    data = phases_labels_date |> filter(phase == "desiccation"),
    aes(x = x_pos_date, y = y_pos, label = description, fill = phase),
    vjust = -0.5, size = 3, fontface = "bold", alpha = 0.8
  ) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = phase_colours, guide = "none") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title = "Daily mean of laboratory environmental conditions",
    subtitle = "Ribbon = ±1 SD",
    x = "Date", y = "Daily mean", colour = "Datalogger"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

ggsave("07-img/laboratory_daily_summary.png", p_daily,
       width = 14, height = 9, dpi = 150)

cat("\n=== Plots saved to 07-img/ ===\n")
cat("  - laboratory_distribution.png\n")
cat("  - laboratory_hourly_evolution.png\n")
cat("  - laboratory_daily_summary.png\n")

# =============================================================================
# STEP 3: Statistical tests — datalogger consistency & phase comparison
# =============================================================================

# --- Test 1: Datalogger consistency (replicates) ----------------------------
# Both loggers measured simultaneously → paired comparison
# H0: no difference between dataloggers

cat("\n========================================\n")
cat("TEST 1: Datalogger consistency (paired)\n")
cat("========================================\n")

# Pivot to wide format: one column per datalogger, matched by datetime
df_wide <- df_lab |>
  select(date, datalogger, temp_c, rh_pct, vpd) |>
  mutate(datalogger = str_replace_all(datalogger, "-", "_")) |>
  pivot_wider(
    names_from  = datalogger,
    values_from = c(temp_c, rh_pct, vpd)
  )

# Paired t-tests
test_dlogger_temp <- t.test(
  df_wide$temp_c_HOBO_21961472,
  df_wide$temp_c_HOBO_21984651,
  paired = TRUE
)

test_dlogger_rh <- t.test(
  df_wide$rh_pct_HOBO_21961472,
  df_wide$rh_pct_HOBO_21984651,
  paired = TRUE
)

test_dlogger_vpd <- t.test(
  df_wide$vpd_HOBO_21961472,
  df_wide$vpd_HOBO_21984651,
  paired = TRUE
)

cat("\nTemperature:\n")
print(test_dlogger_temp)
cat("\nRelative humidity:\n")
print(test_dlogger_rh)
cat("\nVPD:\n")
print(test_dlogger_vpd)

# Effect sizes (mean difference)
cat("\nMean differences (HOBO-21961472 − HOBO-21984651):\n")
cat("  Temp:", mean(df_wide$temp_c_HOBO_21961472 - df_wide$temp_c_HOBO_21984651, na.rm = TRUE), "°C\n")
cat("  RH:  ", mean(df_wide$rh_pct_HOBO_21961472 - df_wide$rh_pct_HOBO_21984651, na.rm = TRUE), "%\n")
cat("  VPD: ", mean(df_wide$vpd_HOBO_21961472 - df_wide$vpd_HOBO_21984651, na.rm = TRUE), "kPa\n")

# --- Test 2: Phase comparison (desiccation vs germination I vs II) ----------
# Three groups → one-way ANOVA + post-hoc pairwise Tukey HSD
# H0: no difference between phases

cat("\n========================================\n")
cat("TEST 2: Phase comparison (ANOVA + Tukey)\n")
cat("========================================\n")

# Filter: desiccation + germination I & II (exclude chamber — markedly different)
df_phases <- df_lab |>
  filter(phase %in% c("desiccation", "germination_I", "germination_II")) |>
  mutate(phase = factor(phase, levels = c("desiccation", "germination_I", "germination_II")))

cat("\nRecords per phase:\n")
table(df_phases$phase)

# --- Temperature ---
aov_temp <- aov(temp_c ~ phase, data = df_phases)
cat("\n--- Temperature ---\n")
print(summary(aov_temp))
tukey_temp <- TukeyHSD(aov_temp)
cat("\nTukey HSD pairwise comparisons:\n")
print(tukey_temp)

# --- Relative humidity ---
aov_rh <- aov(rh_pct ~ phase, data = df_phases)
cat("\n--- Relative humidity ---\n")
print(summary(aov_rh))
tukey_rh <- TukeyHSD(aov_rh)
cat("\nTukey HSD pairwise comparisons:\n")
print(tukey_rh)

# --- VPD ---
aov_vpd <- aov(vpd ~ phase, data = df_phases)
cat("\n--- VPD ---\n")
print(summary(aov_vpd))
tukey_vpd <- TukeyHSD(aov_vpd)
cat("\nTukey HSD pairwise comparisons:\n")
print(tukey_vpd)

# =============================================================================
# STEP 4: Summary statistics for the manuscript
# =============================================================================

cat("\n========================================\n")
cat("SUMMARY — for the manuscript\n")
cat("========================================\n")

# Two groups: lab conditions (desiccation + germination I & II) vs chamber
df_stats <- df_lab |>
  filter(!is.na(phase)) |>
  mutate(
    group = ifelse(phase == "chamber", "chamber", "lab")
  ) |>
  group_by(group) |>
  summarise(
    temp_mean = mean(temp_c, na.rm = TRUE),
    temp_se   = sd(temp_c, na.rm = TRUE) / sqrt(n()),
    rh_mean   = mean(rh_pct, na.rm = TRUE),
    rh_se     = sd(rh_pct, na.rm = TRUE) / sqrt(n()),
    vpd_mean  = mean(vpd, na.rm = TRUE),
    vpd_se    = sd(vpd, na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  ) |>
  mutate(
    temp_str = sprintf("%.2f ± %.2f", temp_mean, temp_se),
    rh_str   = sprintf("%.1f ± %.1f", rh_mean, rh_se),
    vpd_str  = sprintf("%.3f ± %.3f", vpd_mean, vpd_se)
  )

cat("\nLab conditions (desiccation + germination I & II):\n")
df_stats |> filter(group == "lab") |> print()
cat("\nChamber (separate):\n")
df_stats |> filter(group == "chamber") |> print()

# Save to CSV
write.csv(df_stats, "00-data/laboratory_conditions_stats.csv", row.names = FALSE)
cat("\nStats saved to 00-data/laboratory_conditions_stats.csv\n")
