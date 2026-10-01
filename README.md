# Que_bellotas

Analyses conducted within the framework of the QueVadis project on the desiccation
rate and desiccation sensitivity of acorns (*Quercus* spp.), and on the relationship
between desiccation dynamics and the functional traits of the acorns.

The project covers two complementary experiments plus laboratory monitoring:

- **Desiccation-rate experiment** (`d.` scripts): repeated weighing of acorns over
  time under controlled laboratory conditions; trait measurements and modelling of
  drying dynamics.
- **Desiccation-sensitivity experiment** (`s.` scripts): germination response of
  acorns exposed to progressively drier states; estimation of the moisture content
  at which germination probability drops to 50% (MC50).
- **Laboratory conditions** (`00-laboratory_conditions.R`): HOBO datalogger
  temperature/RH records, VPD calculation and phase-wise summaries.

## Repository structure

```
Que_bellotas/
├── 00-data/              Data: raw inputs (read-only) and derived outputs
├── 01-scripts/           Analysis scripts (see conventions below)
├── 06-html/              HTML dashboards and model outputs
├── 07-img/               Figures
├── 08-reports/           Quarto reports and rendered PDFs/HTML
└── 09-bib/               Bibliography
```

## Script conventions

### Naming

| Pattern | Meaning |
|--------------------------|---------|
| `00-*.R` | Project-wide utilities shared across experiments |
| `d.X-name.R` | **Master** script of step `X`, desiccation-rate experiment |
| `d.X.Y-name.R` | **Child** script called by master `d.X`, substep `Y` |
| `s.X-name.R` | Master script of step `X`, desiccation-sensitivity experiment |
| `s.X.Y-name.R` | Child script called by master `s.X`, substep `Y` |

### Rules

1. **Masters orchestrate, children compute.** A master `source()`s its children
   sequentially and performs **all** file exports to `00-data/`. Children build
   objects in memory only.
2. **Raw inputs are read-only.** Raw data live in `00-data/` and are never
   overwritten; derived tables are written back to `00-data/` exclusively by
   master scripts.
3. **Object names:** `rD.*` for raw data as loaded, `df.*` for processed data
   frames.
4. **Units:** weights in grams; datetimes in local time (`Europe/Madrid`).
5. **Language:** file names and comments in English.
6. **Paths** are relative to the project root (open the project via
   `Que_bellotas.Rproj`).

## Pipelines

### Desiccation-sensitivity experiment

Goal: build the germination table (dry weight, moisture content, germination and
emergence per acorn and sampling time) needed for the MC50 analyses.

```
s.01-load_sensitivity_exp.R       master: sources children + exports final CSV
├── s.01.1-load_raw.R             raw phase CSVs -> long tidy records        [done]
├── s.01.2-error_correction.R     manual-review patches + drop list, with
│                                 audit trail exported by the master         [done]
├── s.01.3-dry_weight_table.R     calibration table, FW0->DW model per
│                                 species (form + structure selected by
│                                 diagnostics/CV), DW for every acorn        [done]
└── s.01.4-germination_table.R    moisture content + outcomes -> df.analysis [done]
```

Final exports: `00-data/sensitivity_germination_long.csv` (one row per acorn)
and `00-data/error_correction_log.csv` (audit trail of the manual corrections).

s.01.3 additionally writes its own transparency artifacts (agreed exception to
the master-only-export rule): diagnostic plots to `07-img/dw_model_diagnostics/`
and decision/validation tables plus a selection log to
`00-data/tablas_resumen/` (`dw_*` files).

```
s.02-germination_glm.R           germination GLMs: MC x species response,
│                                drying-batch screening, dredge selection,
│                                DHARMa checks; Q. robur excluded (lot-level
│                                viability failure)                     [done]
└── figures                      RO-exclusion justification, phase-stratified
                                 and batch-marginalized predictive curves

s.03-mc50_critical_moisture.R    MC50 per species, marginalized over batches
│                                on the link scale: delta method (primary)
│                                vs stratified-bootstrap percentile CI
│                                (validation)                           [done]
└── figures                      forest plots: dual-CI transparency report
                                 and article figure coloured by bioclimate
```

Both scripts follow the same transparency pattern agreed for s.01.3: each
exports its own tables to `00-data/tablas_resumen/` (`glm_*`, `mc50_*`) and
its own figures to `07-img/germ_glm_diagnostics/` and
`07-img/mc50_estimation/`.

### Desiccation-rate experiment

```
d.01-load_desiccation_exp.R      master: loads d.01.1 + d.01.2, exports
│                                  00-data/desiccation_traits_long.csv
├── d.01.1-load_data.R           raw wide table -> long format
└── d.01.2-derived_variables.R   derived variables (rates, times)

d.02-trait_tables.R              trait summaries per species/provenance
d.03-trait_famd.R                FAMD on acorn traits
d.04.model_species.R             species-level models

d.05.0.model_traits.R            master: trait effects on desiccation rate
├── d.05.1.heterogeneus_effects.R   [exploratory] heterogeneity is real
├── d.05.2.phylo_data.R             [check] phylogeny + covariance matrix
├── d.05.3.reference_model.R        [REFERENCE MODEL] glmmTMB, random slopes
│                                    by species and prov_code, no phylogeny
├── d.05.4_phylo_check.R            [check] brms + phylogenetic covariance;
│                                    not run in the main flow (hours of
│                                    sampling), launch standalone
└── d.05.5.slope_effects_figures.R  trait effect on the slope from the
                                     time_s:Dim.x coefficients, PRE and POST

The reference model is `Moisture_content ~ time_s * (Dim.1 + Dim.2 + Dim.3)
+ (0 + time_s | species) + (1 + time_s | prov_code) + (1 | id_bellota)`, fitted
separately for PRE (t < 94 h) and POST (t > 94 h). Species enter as *random
slopes* rather than fixed per-species coefficients, so the article's question is
answered by one contrast per trait and heterogeneity is described as variance.
The species term carries no intercept on purpose: `prov_code` is nested inside
`species`, so a species intercept is already representable as a provenance
intercept and estimating both leaves a redundant dimension. That is not a
hypothesis, it is what the fit shows — the variant with a species intercept does
not converge in PRE (non-positive-definite Hessian); in POST it fits with a
singular-convergence warning and a finite AIC, and is kept in the comparison.
`d.05.3.reference_model.R` fits that variant as `m.ref_int` for the record,
compares everything by AIC/Akaike weights with non-converged models reported as
excluded, and stops if the reference itself fails to converge.

Phylogeny is a robustness check, not the reference: it carries no detectable
signal at this taxonomic scale (see `08-reports/brms_filogenias_2.qmd`).

The within/between species decomposition was removed from the scripts; it is
kept only in the reports `Heterogeneus_effects_acorn_traits.qmd` and
`d.07_brms_filogenias.qmd`.
```

### Laboratory conditions

```
00-laboratory_conditions.R       HOBO dataloggers, VPD, phase statistics
```
