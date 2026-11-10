source("R/setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forecast)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

MIN_MONTHS <- 24
LAG_MAX <- 12

on_calendar <- function(residuals, n) {
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

spec <- data.frame(
  figure = c(rep(7, 4), rep(8, 2), rep(9, 2), 10),
  panel = c("a", "b", "c", "d", "a", "b", "a", "b", "a"),
  kind = c(rep("resistance_resistance", 4), rep("prescription_resistance", 5)),
  organism = c(rep("escherichia coli", 4), rep("klebsiella pneumoniae", 4), "escherichia coli"),
  antibiotic_x = c("ceftazidime", "ampicillin", "ceftriaxone", "ceftriaxone", "ceftazidime",
                   "levofloxacin", "amoxicillin/clavulanic acid", "aztreonam", "cefotaxime"),
  antibiotic_y = c("trimethoprim/sulfamethoxazole", "ceftriaxone", "trimethoprim/sulfamethoxazole",
                   "tetracycline", "ceftazidime", "levofloxacin", "amoxicillin/clavulanic acid",
                   "aztreonam", "cefotaxime"))
curves <- inputs <- residual_series <- peaks <- plot_objects <- list()
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]
  if (s$kind == "resistance_resistance") {
    d <- resistance_wide[[s$organism]]
    x <- d[[s$antibiotic_x]]; y <- d[[s$antibiotic_y]]
  } else {
    d <- filter(prescription_panel[[s$organism]], antibiotic == s$antibiotic_x)
    x <- d$exposure; y <- d$pct_resistant
  }
  pw <- prewhiten(x, y)
  if (is.null(pw)) stop("Selected figure pair has insufficient data: ", s$antibiotic_x)
  cc <- ccf_by_lag(pw$x, pw$y)
  obj <- ccf(pw$x, pw$y, lag.max = LAG_MAX, plot = FALSE, na.action = na.pass)
  plot_objects[[i]] <- obj
  cc$original_plot_bound <- qnorm(.975) / sqrt(obj$n.used)
  cc$lag_specific_nominal_bound <- qnorm(.975) / sqrt(cc$n)
  cc$crosses_original_bound <- abs(cc$ccf) > cc$original_plot_bound
  cc$nominal_p_below_05 <- cc$p_value < .05
  curves[[i]] <- cbind(s[rep(1, nrow(cc)), ], cc, row.names = NULL)
  peaks[[i]] <- cbind(s, lead_lag(x, y, s$antibiotic_x, s$antibiotic_y))
  inputs[[i]] <- cbind(s[rep(1, length(x)), ], month = as.character(d$month), x = x, y = y, row.names = NULL)
  residual_series[[i]] <- cbind(s[rep(1, length(x)), ], month = as.character(d$month),
                               residual_x = pw$x, residual_y = pw$y, row.names = NULL)
}
write.csv(bind_rows(curves), "results/tables/leadlag_plotted_lags.csv", row.names = FALSE)
write.csv(bind_rows(peaks), "results/tables/leadlag_plotted_peaks.csv", row.names = FALSE)
write.csv(bind_rows(inputs), "results/checks/leadlag_plotted_inputs.csv", row.names = FALSE)
write.csv(bind_rows(residual_series), "results/checks/leadlag_plotted_residuals.csv", row.names = FALSE)

for (figure in unique(spec$figure)) {
  ix <- which(spec$figure == figure)
  cols <- min(2, length(ix)); rows <- ceiling(length(ix) / cols)
  png(sprintf("results/figures/supp_figure_%02d_leadlag.png", figure),
      width = 7 * cols, height = 5.6 * rows, units = "in", res = 200)
  par(mfrow = c(rows, cols), mar = c(6, 4.5, 5.3, 1.5))
  for (i in ix) {
    s <- spec[i, ]; cc <- curves[[i]]
    title <- if (s$kind == "resistance_resistance") {
      paste0(s$panel, "  ", paste(strwrap(paste(s$antibiotic_x, "vs", s$antibiotic_y),
                                       width = 52), collapse = "\n"), "\n", s$organism)
    } else {
      paste0(s$panel, "  ", s$antibiotic_x,
             "\nPrescription share (x) vs resistance (y)\n", s$organism)
    }
    limits <- range(cc$ccf, -cc$lag_specific_nominal_bound, cc$lag_specific_nominal_bound)
    plot(cc$lag, cc$ccf, type = "h", ylim = limits * 1.05, main = title,
         cex.main = .85, xlab = "Lag (months): negative = x leads y", ylab = "Prewhitened cross-correlation")
    abline(h = 0, col = "grey40")
    abline(h = c(-1, 1) * cc$original_plot_bound[1], col = "blue", lty = 2)
    lines(cc$lag, cc$lag_specific_nominal_bound, col = "#b45309", lty = 3)
    lines(cc$lag, -cc$lag_specific_nominal_bound, col = "#b45309", lty = 3)
    mtext("Blue: paper bound; brown: lag-specific nominal bound (unadjusted)", side = 1, line = 4.4, cex = .64)
  }
  invisible(dev.off())
}
print(bind_rows(peaks)[, c("figure", "panel", "antibiotic_x", "lag_months", "ccf", "p_value")])
writeLines(capture.output(sessionInfo()), "results/checks/R_session_leadlag.txt")
