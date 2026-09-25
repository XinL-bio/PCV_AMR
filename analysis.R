#!/usr/bin/env Rscript
# PCV introduction and pneumococcal AMR analysis.
# Usage and data requirements: see README.md.

required <- c("readxl", "dplyr", "ggplot2", "fixest", "did", "tibble", "patchwork")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(fixest)
  library(patchwork)
})

atlas_file <- Sys.getenv("ATLAS_SP_XLSX", unset = "atlas_sp.xlsx")
pcv_file <- Sys.getenv("PCV_INTRO_XLSX", unset = "PCV_intro.xlsx")
out_dir <- Sys.getenv("PCV_OUTPUT_DIR", unset = "results")
if (!file.exists(atlas_file) || !file.exists(pcv_file)) {
  stop("Provide atlas_sp.xlsx and PCV_intro.xlsx, or set the input path variables.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260723)

check_columns <- function(data, required_names, name) {
  absent <- setdiff(required_names, names(data))
  if (length(absent)) stop(name, " is missing: ", paste(absent, collapse = ", "))
}
check_sir <- function(x, label) {
  invalid <- setdiff(unique(stats::na.omit(as.character(x))),
                     c("Susceptible", "Intermediate", "Resistant"))
  if (length(invalid)) stop(label, " has unexpected categories: ",
                            paste(invalid, collapse = ", "))
}
write_table <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "")
}

atlas <- readxl::read_excel(atlas_file)
pcv <- readxl::read_excel(pcv_file)
drug_columns <- c("Penicillin_I", "Erythromycin_I", "Clindamycin_I",
                  "Levofloxacin_I", "Ceftriaxone_I")
check_columns(atlas, c("Species", "Country", "Year", "Age Group", "Source",
                       drug_columns), "atlas_sp")
check_columns(pcv, c("Country", "Year_all"), "PCV_intro")
pcv <- pcv %>% mutate(Year_all = suppressWarnings(as.numeric(Year_all)))
if (anyNA(atlas$Country) || anyNA(atlas$Year) || anyNA(pcv$Country)) {
  stop("Country and surveillance Year must be present before joining.")
}
if (anyDuplicated(pcv$Country)) stop("PCV_intro has multiple rows per Country.")
if (!is.numeric(atlas$Year)) stop("Atlas Year must be numeric.")
if (any(!is.na(pcv$Year_all) & (pcv$Year_all < 1900 | pcv$Year_all > 2100))) {
  stop("PCV introduction years outside the expected calendar-year range.")
}

sp <- atlas %>% filter(Species == "Streptococcus pneumoniae")
if (nrow(sp) == 0L) stop("No S. pneumoniae isolates in atlas_sp.")
for (nm in drug_columns) check_sir(sp[[nm]], nm)
unmatched <- setdiff(unique(sp$Country), pcv$Country)
if (length(unmatched)) stop("Unmatched country names: ", paste(unmatched, collapse = ", "))
joined <- sp %>% left_join(pcv %>% select(Country, Year_all),
                           by = "Country", relationship = "many-to-one")
if (nrow(joined) != nrow(sp)) stop("PCV country join changed the isolate count.")
if (all(c("Study", "Isolate Id") %in% names(sp))) {
  duplicate_ids <- sp %>% count(Study, `Isolate Id`) %>% filter(n > 1)
  if (nrow(duplicate_ids)) {
    warning(nrow(duplicate_ids), " repeated Study/Isolate Id combinations; ",
            "these records remain in the analysis; review privately.")
  }
}

# Analyze the 68-country study cohort.
country_audit <- joined %>%
  group_by(Country) %>%
  summarise(first_observed = min(Year), last_observed = max(Year),
    introduction = first(Year_all),
    n_isolates = n(), .groups = "drop") %>%
  mutate(status = case_when(
    is.na(introduction) ~ "introduction year missing: cohort set to Inf",
    introduction <= first_observed ~ "already treated at first observation",
    introduction > max(joined$Year) ~ "not treated in study window",
    TRUE ~ "observed before introduction"
  ))
write_table(country_audit, "qc_countries.csv")
analysis_sp <- joined
if (nrow(analysis_sp) != 39029L || n_distinct(analysis_sp$Country) != 68L) {
  stop("Input differs from the study cohort: expected 39,029 isolates / 68 countries.")
}

# MDR: ≥3 non-susceptible classes among ≥3 tested classes.
ns <- function(x) x %in% c("Intermediate", "Resistant")
analysis_sp <- analysis_sp %>%
  mutate(
    beta_ns = case_when(
      is.na(Penicillin_I) & is.na(Ceftriaxone_I) ~ NA,
      ns(Penicillin_I) | ns(Ceftriaxone_I) ~ TRUE,
      TRUE ~ FALSE
    ),
    macro_ns = if_else(is.na(Erythromycin_I), NA, ns(Erythromycin_I)),
    clinda_ns = if_else(is.na(Clindamycin_I), NA, ns(Clindamycin_I)),
    fluoro_ns = if_else(is.na(Levofloxacin_I), NA, ns(Levofloxacin_I))
  ) %>% rowwise() %>%
  mutate(n_classes_tested = sum(!is.na(c(beta_ns, macro_ns, clinda_ns, fluoro_ns))),
         n_classes_ns = sum(c(beta_ns, macro_ns, clinda_ns, fluoro_ns), na.rm = TRUE),
         MDR = if_else(n_classes_tested >= 3L, n_classes_ns >= 3L, NA)) %>%
  ungroup()
write_table(analysis_sp %>% count(n_classes_tested, n_classes_ns),
            "qc_mdr_class_counts.csv")

variables <- c(Penicillin = "Penicillin_I", Macrolide = "Erythromycin_I",
               Clindamycin = "Clindamycin_I", MDR = "MDR")
make_panel <- function(data, outcome) {
  check_columns(data, outcome, "analysis data")
  selected <- data %>% filter(!is.na(.data[[outcome]]))
  if (outcome == "MDR") {
    selected <- selected %>% mutate(non_susceptible = as.integer(.data[[outcome]]))
  } else {
    selected <- selected %>% mutate(non_susceptible = as.integer(ns(.data[[outcome]])))
  }
  panel <- selected %>%
    group_by(Country, Year) %>%
    summarise(n_tested = n(), n_non_susceptible = sum(non_susceptible),
              resistance = n_non_susceptible / n_tested,
              Year_all = first(Year_all), .groups = "drop") %>%
    mutate(cohort = if_else(is.na(Year_all) | Year_all > max(analysis_sp$Year),
                            Inf, as.numeric(Year_all)),
           rel_year = Year - Year_all)
  if (anyDuplicated(panel[c("Country", "Year")])) stop("Duplicate country-year cells.")
  panel
}

fit_sa <- function(panel) {
  feols(resistance ~ sunab(cohort, Year, ref.p = -1) | Country + Year,
        data = panel, weights = ~n_tested, cluster = ~Country)
}
sa_row <- function(fit, panel, label) {
  a <- aggregate(fit, agg = "att")
  data.frame(outcome = label, ATT = unname(a[1, "Estimate"]),
             SE = unname(a[1, "Std. Error"]),
             p_value = unname(a[1, "Pr(>|t|)"]),
             n_countries = n_distinct(panel$Country),
             n_country_years = nrow(panel),
             n_tested = sum(panel$n_tested))
}

panels <- lapply(variables, function(x) make_panel(analysis_sp, x))
write_table(bind_rows(lapply(names(panels), function(x) {
  data.frame(outcome = x, n_tested = sum(panels[[x]]$n_tested),
             n_countries = n_distinct(panels[[x]]$Country),
             n_country_years = nrow(panels[[x]]))
})), "qc_outcome_samples.csv")
models_sa <- lapply(panels, fit_sa)
table_overall <- bind_rows(lapply(names(panels), function(x) {
  sa_row(models_sa[[x]], panels[[x]], x)
}))
table_overall <- table_overall %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE,
         CI_low_pp = 100 * (ATT - 1.96 * SE),
         CI_high_pp = 100 * (ATT + 1.96 * SE))
write_table(table_overall, "table_overall.csv")

age_levels <- c("0 - 17", "18 - 30", "31 - 60", "61+")
age_rows <- list()
age_models <- list()
for (age in age_levels) for (outcome in names(variables)) {
  data_age <- analysis_sp %>% filter(`Age Group` == age)
  panel_age <- make_panel(data_age, variables[[outcome]])
  key <- paste(age, outcome, sep = "__")
  if (nrow(panel_age) < 20L || n_distinct(panel_age$Country) < 5L) {
    warning("Too little data for ", key, "; skipping estimate")
    next
  }
  fit <- tryCatch(fit_sa(panel_age), error = function(e) {
    warning("Failed subgroup ", key, ": ", conditionMessage(e)); NULL
  })
  if (is.null(fit)) next
  age_models[[key]] <- fit
  age_rows[[key]] <- sa_row(fit, panel_age, outcome) %>% mutate(age = age)
}
table_age <- bind_rows(age_rows) %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE,
         CI_low_pp = 100 * (ATT - 1.96 * SE),
         CI_high_pp = 100 * (ATT + 1.96 * SE))
write_table(table_age, "table_age.csv")

# Specimen categories used for the subgroup analysis.
source_groups <- list(Invasive = c("Blood", "CSF"),
                      Respiratory = c("Sputum", "Bronchus", "Bronchoalveolar lavage",
                                      "Endotracheal aspirate", "Trachea"))
write_table(analysis_sp %>% count(Source, sort = TRUE), "qc_source_values.csv")
specimen_rows <- list()
for (group in names(source_groups)) for (outcome in c("Penicillin", "Macrolide")) {
  subgroup <- analysis_sp %>% filter(Source %in% source_groups[[group]])
  panel_sub <- make_panel(subgroup, variables[[outcome]])
  key <- paste(group, outcome, sep = "__")
  if (nrow(panel_sub) < 20L || n_distinct(panel_sub$Country) < 5L) next
  fit <- tryCatch(fit_sa(panel_sub), error = function(e) {
    warning("Failed subgroup ", key, ": ", conditionMessage(e)); NULL
  })
  if (!is.null(fit)) specimen_rows[[key]] <-
    sa_row(fit, panel_sub, outcome) %>% mutate(specimen_group = group)
}
table_specimen <- bind_rows(specimen_rows)
if (nrow(table_specimen)) table_specimen <- table_specimen %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE)
write_table(table_specimen, "table_specimen.csv")

# Baseline cutoffs include countries with future introduction years.
baseline_rows <- list()
for (outcome in c("Penicillin", "Macrolide")) {
  main_panel <- panels[[outcome]]
  baseline <- main_panel %>%
    filter(rel_year >= -5, rel_year <= -1) %>%
    group_by(Country) %>%
    summarise(baseline_rate = weighted.mean(resistance, n_tested),
              n_pre_years = n(), .groups = "drop")
  if (nrow(baseline) < 9L) {
    warning("Too few baseline countries for ", outcome)
    next
  }
  limits <- as.numeric(stats::quantile(baseline$baseline_rate,
                                        probs = c(0.33, 0.67), na.rm = TRUE))
  baseline <- baseline %>%
    mutate(group = case_when(baseline_rate <= limits[1] ~ "Low",
                             baseline_rate <= limits[2] ~ "Middle",
                             TRUE ~ "High"))
  write_table(baseline, paste0("qc_baseline_", outcome, ".csv"))
  grouped_panel <- main_panel %>% left_join(baseline, by = "Country")
  for (group in c("Low", "Middle", "High")) {
    sub <- grouped_panel %>% filter(group == .env$group | !is.finite(cohort))
    if (n_distinct(sub$Country[is.finite(sub$cohort)]) < 3L ||
        n_distinct(sub$Country[!is.finite(sub$cohort)]) < 1L) {
      warning("Insufficient treatment/control support for ", outcome, "/", group)
      next
    }
    fit <- tryCatch(fit_sa(sub), error = function(e) {
      warning("Failed baseline group ", outcome, "/", group, ": ", conditionMessage(e)); NULL
    })
    if (!is.null(fit)) baseline_rows[[paste(outcome, group)]] <-
      sa_row(fit, sub, outcome) %>% mutate(baseline_group = group,
                                           cutoff_low = limits[1],
                                           cutoff_high = limits[2],
                                           n_pre_years_min = min(baseline$n_pre_years[baseline$group == group]))
  }
  if (outcome == "Penicillin") {
    high_trim <- grouped_panel %>%
      filter(group == "High" | !is.finite(cohort), n_tested >= 5L)
    trim_fit <- fit_sa(high_trim)
    baseline_rows[["Penicillin High >=5"]] <-
      sa_row(trim_fit, high_trim, outcome) %>%
      mutate(baseline_group = "High; >=5 isolates/cell",
             cutoff_low = limits[1], cutoff_high = limits[2],
             n_pre_years_min = min(baseline$n_pre_years[baseline$group == "High"]))
  }
}
table_baseline <- bind_rows(baseline_rows)
if (nrow(table_baseline)) table_baseline <- table_baseline %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE)
write_table(table_baseline, "table_baseline_exploratory.csv")

# Income classification snapshot used in this study.
income_lookup <- tibble::tribble(
  ~Country, ~income_group,
  "Argentina", "Upper middle", "Australia", "High",
  "Austria", "High", "Belgium", "High",
  "Brazil", "Upper middle", "Bulgaria", "High",
  "Canada", "High", "Chile", "High",
  "China", "Upper middle", "Colombia", "Upper middle",
  "Costa Rica", "Upper middle", "Croatia", "High",
  "Czech Republic", "High", "Denmark", "High",
  "Dominican Republic", "Upper middle", "Finland", "High",
  "France", "High", "Germany", "High",
  "Greece", "High", "Guatemala", "Upper middle",
  "Honduras", "Lower middle", "Hong Kong", "High",
  "Hungary", "High", "India", "Lower middle",
  "Ireland", "High", "Israel", "High",
  "Italy", "High", "Jamaica", "Upper middle",
  "Japan", "High", "Jordan", "Upper middle",
  "Kenya", "Lower middle", "Korea, South", "High",
  "Kuwait", "High", "Latvia", "High",
  "Lithuania", "High", "Malaysia", "Upper middle",
  "Mauritius", "Upper middle", "Mexico", "Upper middle",
  "Morocco", "Lower middle", "Netherlands", "High",
  "New Zealand", "High", "Norway", "High",
  "Oman", "High", "Pakistan", "Lower middle",
  "Panama", "High", "Philippines", "Upper middle",
  "Poland", "High", "Portugal", "High",
  "Qatar", "High", "Romania", "High",
  "Russia", "Upper middle", "Saudi Arabia", "High",
  "Serbia", "Upper middle", "Singapore", "High",
  "Slovak Republic", "High", "Slovenia", "High",
  "South Africa", "Upper middle", "Spain", "High",
  "Sweden", "High", "Switzerland", "High",
  "Taiwan", "High", "Thailand", "Upper middle",
  "Tunisia", "Lower middle", "Turkey", "Upper middle",
  "Ukraine", "Lower middle", "United Kingdom", "High",
  "Venezuela", "Upper middle", "Vietnam", "Upper middle"
)
if (anyDuplicated(income_lookup$Country) ||
    length(setdiff(unique(analysis_sp$Country), income_lookup$Country))) {
  stop("Country income lookup missing or duplicated countries.")
}
write_table(income_lookup, "qc_income_lookup.csv")
income_rows <- list()
for (outcome in c("Penicillin", "Macrolide")) {
  panel_income <- panels[[outcome]] %>%
    left_join(income_lookup, by = "Country", relationship = "many-to-one")
  for (group in c("High", "Upper middle", "Lower middle")) {
    sub <- panel_income %>% filter(income_group == .env$group | !is.finite(cohort))
    fit <- fit_sa(sub)
    income_rows[[paste(outcome, group)]] <- sa_row(fit, sub, outcome) %>%
      mutate(income_group = group,
             treated_countries = n_distinct(sub$Country[is.finite(sub$cohort)]),
             controls_other_income = n_distinct(sub$Country[
               !is.finite(sub$cohort) & sub$income_group != group]))
  }
}
table_income <- bind_rows(income_rows) %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE)
write_table(table_income, "table_income.csv")

# Adult-program labels are country-level.
adult_vax_countries <- c(
  "Argentina", "Australia", "Austria", "Belgium", "Canada", "Chile",
  "Costa Rica", "Germany", "Spain", "France", "United Kingdom", "Greece",
  "Hong Kong", "Ireland", "Israel", "Italy", "Japan", "Korea, South",
  "Lithuania", "Mexico", "Panama", "Philippines", "Qatar", "Sweden")
if (length(adult_vax_countries) != 24L ||
    length(setdiff(adult_vax_countries, analysis_sp$Country))) {
  stop("Adult-program country mapping differs from the study cohort.")
}
adult_rows <- list()
for (outcome in c("Penicillin", "Macrolide")) {
  panel_adult <- make_panel(analysis_sp %>% filter(`Age Group` == "61+"),
                            variables[[outcome]]) %>%
    mutate(adult_group = if_else(Country %in% adult_vax_countries,
                                 "Adult program", "Pediatric only"))
  for (group in c("Adult program", "Pediatric only")) {
    sub <- panel_adult %>% filter(adult_group == .env$group | !is.finite(cohort))
    fit <- fit_sa(sub)
    adult_rows[[paste(outcome, group)]] <- sa_row(fit, sub, outcome) %>%
      mutate(adult_group = group,
             treated_countries = n_distinct(sub$Country[is.finite(sub$cohort)]))
  }
}
table_adult <- bind_rows(adult_rows) %>%
  mutate(ATT_pp = 100 * ATT, SE_pp = 100 * SE)
write_table(table_adult, "table_adult_program.csv")

# C&S uses the same panels as the primary analysis.
fit_cs <- function(panel) {
  cs_data <- panel %>%
    mutate(country_id = as.integer(factor(Country)),
           gname = if_else(is.finite(cohort), cohort, 0))
  gt <- did::att_gt(yname = "resistance", tname = "Year", idname = "country_id",
               gname = "gname", data = cs_data, weightsname = "n_tested",
               control_group = "notyettreated", allow_unbalanced_panel = TRUE,
               clustervars = "country_id")
  list(gt = gt, simple = did::aggte(gt, type = "simple", na.rm = TRUE),
       dynamic = did::aggte(gt, type = "dynamic", na.rm = TRUE))
}
models_cs <- lapply(names(panels), function(outcome) {
  tryCatch(fit_cs(panels[[outcome]]), error = function(e) {
    warning("C&S failed for ", outcome, ": ", conditionMessage(e)); NULL
  })
})
names(models_cs) <- names(panels)
cs_rows <- bind_rows(lapply(names(models_cs), function(outcome) {
  z <- models_cs[[outcome]]
  if (is.null(z)) return(NULL)
  bind_rows(lapply(c("simple", "dynamic"), function(which) {
    agg <- z[[which]]
    data.frame(outcome = outcome, aggregation = which,
               ATT = agg$overall.att, SE = agg$overall.se,
               ATT_pp = 100 * agg$overall.att,
               SE_pp = 100 * agg$overall.se,
               CI_low_pp = 100 * (agg$overall.att - 1.96 * agg$overall.se),
               CI_high_pp = 100 * (agg$overall.att + 1.96 * agg$overall.se),
               pretrend_Wald_p = if (length(z$gt$Wpval) == 1L) z$gt$Wpval
                                 else NA_real_)
  }))
}))
write_table(cs_rows, "table_cs_sensitivity.csv")

# Figure 1/2 use the study's drug order and significance colors.
drug_order <- c("Penicillin", "Macrolide", "Clindamycin", "MDR")
forest_colors <- c("Not significant" = "grey55",
                   "Significant (p<0.05)" = "firebrick")
forest_data <- function(x) {
  x %>% mutate(
    display_outcome = if_else(outcome == "MDR", "MDR (>=3 classes)", outcome),
    display_outcome = factor(display_outcome,
      levels = rev(c("Penicillin", "Macrolide", "Clindamycin", "MDR (>=3 classes)"))),
    sig = factor(if_else(p_value < 0.05, "Significant (p<0.05)",
                        "Not significant"), levels = names(forest_colors)))
}
forest_plot <- function(x) {
  ggplot(forest_data(x), aes(x = ATT_pp, y = display_outcome, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = CI_low_pp, xmax = CI_high_pp),
                   height = 0.18, linewidth = 0.65) +
    geom_point(size = 3) +
    scale_color_manual(values = forest_colors, drop = FALSE) +
    labs(x = "ATT (percentage points)", y = NULL, colour = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey96", colour = "black"))
}
f1 <- forest_plot(table_overall)
ggsave(file.path(out_dir, "figure1_overall.png"), f1,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")
f2 <- forest_plot(table_age %>% filter(age %in% c("0 - 17", "61+")) %>%
                    mutate(Population = if_else(age == "0 - 17",
                                                "Children (0–17)",
                                                "Older adults (61+)"),
                           Population = factor(Population,
                             levels = c("Children (0–17)", "Older adults (61+)")))) +
  facet_wrap(~Population, ncol = 1)
ggsave(file.path(out_dir, "figure2_age.png"), f2,
       width = 6.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "figure1_2_combined.png"),
       (f1 + f2 + plot_layout(guides = "collect")) & theme(legend.position = "bottom"),
       width = 13, height = 6, dpi = 300, bg = "white")
dynamic_rows <- bind_rows(lapply(names(models_cs), function(outcome) {
  z <- models_cs[[outcome]]
  if (is.null(z)) return(NULL)
  critical <- z$dynamic$crit.val.egt
  if (length(critical) != 1L || !is.finite(critical)) critical <- 1.96
  data.frame(outcome = outcome, event_year = z$dynamic$egt,
             ATT = z$dynamic$att.egt, SE = z$dynamic$se.egt,
             critical_value = critical,
             observed_treated_countries = vapply(z$dynamic$egt, function(e) {
               p <- panels[[outcome]]
               n_distinct(p$Country[is.finite(p$cohort) & p$Year - p$cohort == e])
             }, integer(1)))
}))
write_table(dynamic_rows, "table_dynamic.csv")
if (nrow(dynamic_rows)) {
  # Four plots retain their own axes; ATT remains on the proportion scale.
  event_data <- dynamic_rows %>%
    mutate(lo = ATT - critical_value * SE,
           hi = ATT + critical_value * SE,
           period = factor(if_else(event_year < 0, "Pre-introduction",
                                   "Post-introduction"),
                           levels = c("Pre-introduction", "Post-introduction")))
  panel_titles <- c(
    Penicillin = "A. Penicillin Non-Susceptibility",
    Macrolide = "B. Macrolide Non-Susceptibility",
    Clindamycin = "C. Clindamycin Non-Susceptibility",
    MDR = "D. Multidrug Non-Susceptibility (MDR)")
  event_panel <- function(drug) {
    ggplot(event_data %>% filter(outcome == drug),
           aes(x = event_year, y = ATT, colour = period)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.55) +
      geom_errorbar(aes(ymin = lo, ymax = hi),
                    width = 0, linewidth = 0.6) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Pre-introduction" = "#E87D72",
                                    "Post-introduction" = "#56BCC2")) +
      scale_y_continuous(breaks = seq(-1.0, 0.5, by = 0.5)) +
      coord_cartesian(ylim = c(-1.1, 0.7)) +
      scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
      labs(title = panel_titles[[drug]],
           x = "Years Relative to PCV Introduction",
           y = "ATT (Percentage Points)", colour = NULL) +
      theme_classic(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 13),
            axis.title = element_text(face = "bold"),
            legend.position = "none",
            panel.grid.major.y = element_line(colour = "grey92",
                                              linetype = "dotted"))
  }
  event_plots <- lapply(drug_order, event_panel)
  f3 <- (event_plots[[1]] + event_plots[[2]]) /
        (event_plots[[3]] + event_plots[[4]])
  ggsave(file.path(out_dir, "figure3_dynamic.png"), f3,
         width = 10, height = 8.5, dpi = 300, bg = "white")
}

saveRDS(list(SA = models_sa, age_SA = age_models, CS = models_cs),
        file.path(out_dir, "models.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
