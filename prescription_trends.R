source("R/setup.R")
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

fit_prescription_model <- function(data, context) {
  slope <- cbind(rx_ab, other_rx) ~ month_index + (1 + month_index | state)
  intercept <- cbind(rx_ab, other_rx) ~ month_index + (1 | state)
  candidates <- list(
    list(slope, betabinomial("logit"), "beta-binomial", "random intercept + random slope by state"),
    list(intercept, betabinomial("logit"), "beta-binomial", "random intercept by state"),
    list(slope, binomial("logit"), "binomial", "random intercept + random slope by state"),
    list(intercept, binomial("logit"), "binomial", "random intercept by state"))
  for (m in candidates) {
    messages <- character()
    fit <- tryCatch(withCallingHandlers(
      glmmTMB(m[[1]], family = m[[2]], data = data, control = control),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) { messages <<- c(messages, conditionMessage(e)); NULL })
    record_fit("results/checks/prescription_fit_attempts.csv", context, m[[3]], m[[4]], fit,
               paste(messages, collapse = " | "))
    if (!is.null(fit) && isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess) &&
        is.finite(as.numeric(logLik(fit)))) {
      return(list(fit = fit, model_type = m[[3]], random_structure = m[[4]]))
    }
  }
  NULL
}

source("R/monthly_marginal_effect.R")
if (file.exists("results/checks/prescription_fit_attempts.csv")) file.remove("results/checks/prescription_fit_attempts.csv")
prescription_tables <- function(panel, analysis) {
  revised_rows <- state_rows <- list()
  for (product in products) {
    d <- panel[antibiotic == product]
    if (nrow(d) < 30 || sum(d$rx_ab) == 0) next
    stopifnot(all(d$rx_ab >= 0), all(d$other_rx >= 0), all(d$den > 0))
    context <- paste(analysis, product, sep = " | ")
    cat("Fitting", context, "\n"); flush.console()
    model <- fit_prescription_model(d, context)
    if (is.null(model)) stop("All candidate models failed: ", context)
    fit <- model$fit
    effect <- monthly_marginal_effect(fit, d, "month_index")
    group <- top[antibiotic == product, aware_2025][1]
    common <- data.table(antibiotic = product, aware_2025 = group,
      model_type = model$model_type, random_structure = model$random_structure,
      n_states = uniqueN(d$state), n_months = uniqueN(d$ym),
      mean_share_pct = round(100 * mean(d$share), 4))
    revised_rows[[product]] <- cbind(common, data.table(
      monthly_ppt_change = effect$estimate,
      CI95_low_ppt_per_month = effect$low, CI95_high_ppt_per_month = effect$high,
      p_value = effect$p_value))

    state_rows[[product]] <- cbind(data.table(analysis = analysis, antibiotic = product), effect$state_summary)
    save_model(fit, d, paste("prescription", analysis, product))
    last_prescription_model <<- fit
    last_prescription_data <<- d
    cat(sprintf("  %.5f pp/month [%.5f, %.5f]\n", effect$estimate, effect$low, effect$high))
  }
  revised <- rbindlist(revised_rows)
  revised[, p_value_adj_BH := p.adjust(p_value, "BH"), by = aware_2025]
  setorder(revised, aware_2025, -mean_share_pct)
  stem <- if (analysis == "all_antibiotics") "supp_table4_prescription_trends" else
            "supp_table5_prescription_trends_all_medicines"
  fwrite(revised, paste0("results/tables/", stem, "_monthly.csv"))
  fwrite(rbindlist(state_rows), paste0("results/checks/prescription_state_effects_", analysis, ".csv"))
  revised
}

prescription_abx <- prescription_tables(build_panel(antibiotic_totals), "all_antibiotics")
prescription_medicines <- prescription_tables(build_panel(medicine_totals), "all_medicines")
print(prescription_abx)
print(prescription_medicines)
writeLines(capture.output(sessionInfo()), "results/checks/R_session_prescription.txt")
