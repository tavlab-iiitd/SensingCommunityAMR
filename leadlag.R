# Lead-lag analysis (prewhitened cross-correlation) and Supplementary Figures 8-11.
#
# For a pair of monthly series x and y: fit an ARIMA model to x (auto.arima), filter y with the
# same model, and cross-correlate the two residual series at lags -12 to +12 months. Missing months
# stay missing (no interpolation) and at least 24 paired months are required. The peak is the lag
# with the largest absolute correlation; p-values are unadjusted (exploratory analysis).
# Negative lag: x leads y. Positive lag: y leads x.
#
#   1. Prescription vs resistance (urine E. coli and K. pneumoniae): x = monthly share of all
#      antibiotic prescriptions, y = monthly % resistant to the same antibiotic.
#   2. Resistance vs resistance (urine E. coli and K. pneumoniae): x and y = monthly % resistant
#      to two different antibiotics.
#
# Input : data/prescription_resistance_monthly.csv, data/ast_clean.csv
# Output: results/tables/leadlag_*.csv, results/figures/supp_fig8-11_*.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forecast)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

MIN_MONTHS <- 24
LAG_MAX <- 12

# ---------------------------------------------------------------- method

on_calendar <- function(residuals, n) {
  # put residuals back on the full monthly grid (auto.arima can drop leading missing months)
  out <- rep(NA_real_, n)
  position <- as.integer(round(as.numeric(time(residuals))))
  keep <- position >= 1 & position <= n
  out[position[keep]] <- as.numeric(residuals)[keep]
  out
}

prewhiten <- function(x, y) {
  if (sum(is.finite(x) & is.finite(y)) < MIN_MONTHS ||
      sd(x[is.finite(x)]) == 0 || sd(y[is.finite(y)]) == 0) return(NULL)
  fit_x <- tryCatch(auto.arima(ts(x)), error = function(e) NULL)
  if (is.null(fit_x)) return(NULL)
  fit_y <- tryCatch(Arima(ts(y), model = fit_x), error = function(e) NULL)   # same filter applied to y
  if (is.null(fit_y)) return(NULL)
  px <- on_calendar(residuals(fit_x), length(x))
  py <- on_calendar(residuals(fit_y), length(y))
  if (sum(is.finite(px) & is.finite(py)) < MIN_MONTHS) return(NULL)
  list(x = px, y = py, arima_x = paste(arimaorder(fit_x), collapse = ","))
}

ccf_by_lag <- function(px, py) {
  # cross-correlation cor(x[t + lag], y[t]); n = residual pairs available at each lag
  cc <- ccf(px, py, lag.max = LAG_MAX, plot = FALSE, na.action = na.pass)
  out <- data.frame(lag = as.integer(cc$lag), ccf = as.numeric(cc$acf))
  out$n <- sapply(out$lag, function(k) {
    t <- max(1, 1 - k):min(length(px), length(px) - k)
    sum(is.finite(px[t + k]) & is.finite(py[t]))
  })
  out$p_value <- 2 * (1 - pnorm(abs(out$ccf) * sqrt(out$n)))
  out
}

lead_lag <- function(x, y, name_x, name_y) {
  pw <- prewhiten(x, y)
  if (is.null(pw)) return(NULL)
  cc <- ccf_by_lag(pw$x, pw$y)
  peak <- cc[which.max(ifelse(is.finite(cc$ccf) & cc$n >= MIN_MONTHS - LAG_MAX, abs(cc$ccf), -Inf)), ]
  data.frame(series_x = name_x, series_y = name_y, lag_months = peak$lag, ccf = peak$ccf,
             p_value = peak$p_value, n_used = peak$n, arima_x = pw$arima_x,
             interpretation = if (peak$lag < 0) paste(name_x, "leads", name_y, "by", -peak$lag, "months")
                              else if (peak$lag > 0) paste(name_y, "leads", name_x, "by", peak$lag, "months")
                              else "contemporaneous")
}

monthly_grid <- function(d) {
  d %>% complete(month = seq(min(month), max(month), by = "month"), antibiotic) %>%
    arrange(antibiotic, month)
}

# ---------------------------------------------------------------- data

merged <- read.csv("data/prescription_resistance_monthly.csv") %>%
  mutate(organism = tolower(trimws(organism_clean)), antibiotic = tolower(trimws(salt_name)),
         month = as.Date(paste0(month_year, "-01")))

ast <- read.csv("data/ast_clean.csv") %>%
  mutate(organism = tolower(trimws(organism_clean)), antibiotic = tolower(trimws(antibiotic_clean)),
         month = as.Date(paste0(month_year, "-01"))) %>%
  filter(test_name == "URINE", nzchar(antibiotic), antibiotic != "na")

aware <- unique(rbind(merged[, c("antibiotic", "aware_2025")], ast[, c("antibiotic", "aware_2025")]))
label <- function(ab) paste0(ab, " [", aware$aware_2025[match(ab, aware$antibiotic)], "]")

compared <- tolower(c(
  "Amikacin", "Amoxicillin/clavulanic acid", "Ampicillin", "Aztreonam", "Cefepime", "Cefixime",
  "Cefotaxime", "Cefoxitin", "Ceftazidime", "Ceftriaxone", "Cefuroxime", "Ciprofloxacin", "Colistin",
  "Ertapenem", "Fosfomycin", "Gentamicin", "Imipenem", "Levofloxacin", "Meropenem", "Nitrofurantoin",
  "Norfloxacin", "Ofloxacin", "Piperacillin", "Piperacillin/tazobactam", "Tetracycline",
  "Tigecycline", "Tobramycin", "Trimethoprim/sulfamethoxazole"))

prescription_panel <- list()
resistance_wide <- list()

for (org in c("escherichia coli", "klebsiella pneumoniae")) {
  slug <- gsub(" ", "_", org)

  # 1. prescription share (x) vs % resistant (y), same antibiotic
  panel <- merged %>% filter(organism == org) %>%
    group_by(month, antibiotic) %>%
    summarise(pct_resistant = ifelse(sum(total_tests) > 0, 100 * sum(RESISTANT) / sum(total_tests), NA_real_),
              exposure = if (any(!is.na(pres_normalized))) pres_normalized[!is.na(pres_normalized)][1] else NA_real_,
              .groups = "drop") %>%
    monthly_grid()
  prescription_panel[[org]] <- panel

  eligible <- panel %>% group_by(antibiotic) %>%
    summarise(months = sum(is.finite(exposure) & is.finite(pct_resistant))) %>%
    filter(months >= MIN_MONTHS) %>% arrange(desc(months))
  results <- bind_rows(lapply(eligible$antibiotic, function(ab) {
    d <- filter(panel, antibiotic == ab)
    out <- lead_lag(d$exposure, d$pct_resistant, paste(ab, "Rx"), paste(ab, "resistance"))
    if (!is.null(out)) cbind(antibiotic = ab, aware_2025 = aware$aware_2025[match(ab, aware$antibiotic)], out)
  })) %>% arrange(desc(abs(ccf)), p_value)
  write.csv(results, paste0("results/tables/leadlag_prescription_vs_resistance_", slug, ".csv"), row.names = FALSE)

  # 2. % resistant to antibiotic A (x) vs antibiotic B (y)
  panel <- ast %>% filter(organism == org) %>%
    group_by(month, antibiotic) %>%
    summarise(tested = sum(SENSITIVE + RESISTANT + INTERMEDIATE),
              pct_resistant = ifelse(tested > 0, 100 * sum(RESISTANT) / tested, NA_real_), .groups = "drop") %>%
    monthly_grid()
  eligible <- panel %>% group_by(antibiotic) %>%
    summarise(months = sum(is.finite(pct_resistant)), tested = sum(tested, na.rm = TRUE)) %>%
    filter(months >= MIN_MONTHS, tested >= 200, antibiotic %in% compared)
  wide <- panel %>% select(month, antibiotic, pct_resistant) %>%
    pivot_wider(names_from = antibiotic, values_from = pct_resistant) %>% arrange(month)
  resistance_wide[[org]] <- wide

  results <- bind_rows(lapply(combn(sort(eligible$antibiotic), 2, simplify = FALSE), function(p) {
    out <- lead_lag(wide[[p[1]]], wide[[p[2]]], p[1], p[2])
    if (!is.null(out)) cbind(antibiotic_1 = p[1], antibiotic_2 = p[2], out)
  })) %>% arrange(desc(abs(ccf)), p_value)
  write.csv(results, paste0("results/tables/leadlag_resistance_vs_resistance_", slug, ".csv"), row.names = FALSE)
}

# ---------------------------------------------------------------- Supplementary Figures 8-11

plot_pair <- function(x, y, title) {
  pw <- prewhiten(x, y)
  plot(ccf(pw$x, pw$y, lag.max = LAG_MAX, plot = FALSE, na.action = na.pass), main = title)
}

ecoli <- resistance_wide[["escherichia coli"]]
png("results/figures/supp_fig8_ecoli_resistance_vs_resistance.png", width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
for (p in list(c("ceftazidime", "trimethoprim/sulfamethoxazole"), c("ampicillin", "ceftriaxone"),
               c("ceftriaxone", "trimethoprim/sulfamethoxazole"), c("ceftriaxone", "tetracycline"))) {
  plot_pair(ecoli[[p[1]]], ecoli[[p[2]]], paste(label(p[1]), "vs", label(p[2])))
}
invisible(dev.off())

prescription_figure <- function(file, org, short, antibiotics) {
  png(file.path("results/figures", file), width = 7 * length(antibiotics), height = 5.4, units = "in", res = 300)
  par(mfrow = c(1, length(antibiotics)))
  for (ab in antibiotics) {
    d <- filter(prescription_panel[[org]], antibiotic == ab)
    plot_pair(d$exposure, d$pct_resistant, paste0(label(ab), " - Pres normalised vs % Resistant (", short, ")"))
  }
  invisible(dev.off())
}
prescription_figure("supp_fig9_kpneumoniae_prescription_vs_resistance.png", "klebsiella pneumoniae", "Kpn",
                    c("ceftazidime", "levofloxacin"))
prescription_figure("supp_fig10_kpneumoniae_prescription_vs_resistance.png", "klebsiella pneumoniae", "Kpn",
                    c("amoxicillin/clavulanic acid", "aztreonam"))
prescription_figure("supp_fig11_ecoli_prescription_vs_resistance.png", "escherichia coli", "E. coli", "cefotaxime")
