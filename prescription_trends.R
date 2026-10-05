# Supplementary Tables 4 and 5: prescription trends of the ten most prescribed products in each
# AWaRe group.
#
# For each product, a beta-binomial GLMM (logit link) of its monthly prescriptions out of all
# antibiotic prescriptions (Table 4) or all medicine prescriptions (Table 5) in each state, with
# month as a fixed effect and random intercepts and slopes by state. States with antibiotic activity
# in at least 75% of months are included. If a model does not converge, simpler models are tried in
# order. Trends are average marginal effects in percentage points; p-values are BH-adjusted within
# each AWaRe group.
#
# Input : data/prescriptions_clean.csv, data/State_Month_Prescription_updated.csv
# Output: results/tables/supp_table4_prescription_trends.csv,
#         results/tables/supp_table5_prescription_trends_all_medicines.csv

library(data.table)
library(zoo)
library(glmmTMB)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

rx <- fread("data/prescriptions_clean.csv")[state_included == TRUE]
rx[, `:=`(state = tools::toTitleCase(trimws(patient_rx_state)), ym = as.yearmon(month_year, "%Y-%m"))]

medicines <- fread("data/State_Month_Prescription_updated.csv")[month_year >= "2023-01" & month_year <= "2025-12"]
medicines[, `:=`(state = tools::toTitleCase(trimws(Rx_State)), ym = as.yearmon(month_year, "%Y-%m"))]

antibiotic_totals <- rx[, .(den = sum(Rx_Count)), by = .(state, ym)]
medicine_totals <- medicines[presc > 0, .(den = sum(presc)), by = .(state, ym)]

months <- sort(unique(antibiotic_totals$ym))
coverage <- antibiotic_totals[, .(share = uniqueN(ym) / length(months)), by = state]
states <- coverage[share >= 0.75, state]

top <- rx[aware_2025 %in% c("Access", "Watch", "Reserve", "Not recommended"),
          .(total = sum(Rx_Count)), by = .(aware_2025, antibiotic = product_clean)][order(aware_2025, -total)]
top <- top[, head(.SD, 10), by = aware_2025]
products <- unique(top$antibiotic)

build_panel <- function(totals) {
  panel <- CJ(antibiotic = products, state = sort(states), ym = months, unique = TRUE)
  panel <- merge(panel, totals[state %in% states & ym %in% months], by = c("state", "ym"), all.x = TRUE)
  panel <- panel[!is.na(den) & den > 0]
  counts <- rx[product_clean %in% products, .(rx_ab = sum(Rx_Count)), by = .(antibiotic = product_clean, state, ym)]
  panel <- merge(panel, counts, by = c("antibiotic", "state", "ym"), all.x = TRUE)
  panel[is.na(rx_ab), rx_ab := 0]          # state-months without the product count as zero
  panel[, `:=`(other_rx = den - rx_ab, month_index = as.numeric((ym - min(months)) * 12), share = rx_ab / den)]
  panel
}

control <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))

fit_model <- function(data) {
  slope <- cbind(rx_ab, other_rx) ~ month_index + (1 + month_index | state)
  intercept <- cbind(rx_ab, other_rx) ~ month_index + (1 | state)
  candidates <- list(
    list(slope, betabinomial("logit"), "beta-binomial", "random intercept + random slope by state"),
    list(intercept, betabinomial("logit"), "beta-binomial", "random intercept by state"),
    list(slope, binomial("logit"), "binomial", "random intercept + random slope by state"),
    list(intercept, binomial("logit"), "binomial", "random intercept by state"))
  for (m in candidates) {
    fit <- tryCatch(suppressWarnings(glmmTMB(m[[1]], family = m[[2]], data = data, control = control)),
                    error = function(e) NULL)
    if (!is.null(fit) && isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess) &&
        is.finite(as.numeric(logLik(fit)))) {
      return(list(fit = fit, model_type = m[[3]], random_structure = m[[4]]))
    }
  }
  NULL
}

trend_table <- function(panel) {
  rows <- lapply(products, function(product) {
    data <- panel[antibiotic == product]
    if (nrow(data) < 30 || sum(data$rx_ab) == 0) return(NULL)
    model <- fit_model(data)
    if (is.null(model)) return(NULL)

    coefs <- summary(model$fit)$coefficients$cond
    beta <- coefs["month_index", "Estimate"]
    se <- coefs["month_index", "Std. Error"]
    fixed <- fixef(model$fit)$cond
    p <- plogis(fixed[["(Intercept)"]] + fixed[["month_index"]] * data$month_index)
    pp <- mean(p * (1 - p)) * 100            # logit slope -> percentage points per month

    data.table(antibiotic = product, aware_2025 = top[antibiotic == product, aware_2025][1],
               model_type = model$model_type, random_structure = model$random_structure,
               n_states = uniqueN(data$state), n_months = uniqueN(data$ym),
               mean_share_pct = round(mean(data$share) * 100, 4),
               monthly_ppt_change = round(pp * beta, 5),
               annual_ppt_change = round(pp * beta * 12, 5),
               CI95_low_ppt_per_year = round(pp * (beta - qnorm(0.975) * se) * 12, 5),
               CI95_high_ppt_per_year = round(pp * (beta + qnorm(0.975) * se) * 12, 5),
               p_value = coefs["month_index", "Pr(>|z|)"])
  })
  out <- rbindlist(rows)
  out[, p_value_adj_BH := p.adjust(p_value, method = "BH"), by = aware_2025]
  out[order(aware_2025, -mean_share_pct)]
}

fwrite(trend_table(build_panel(antibiotic_totals)), "results/tables/supp_table4_prescription_trends.csv")
fwrite(trend_table(build_panel(medicine_totals)), "results/tables/supp_table5_prescription_trends_all_medicines.csv")
