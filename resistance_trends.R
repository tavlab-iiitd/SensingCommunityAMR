# Supplementary Tables 6-8: resistance trends in urine isolates of E. coli, K. pneumoniae and
# P. aeruginosa.
#
# For each organism-antibiotic pair, a beta-binomial GLMM (logit link) of resistant isolates out of
# isolates tested per state-month, with month as a fixed effect and random intercepts and slopes by
# state. States are included if they have data in at least 75% of the months between the first and
# last observation of the pair and at least two pincodes. If a model fails (error, warning,
# non-positive-definite Hessian), binomial and random-intercept-only models are tried in order.
# Trends are average marginal effects in percentage points per year; p-values are BH-adjusted
# within each organism. Intrinsic resistance pairs (CLSI M100 Appendix B) are not modelled.
#
# Input : data/ast_clean.csv
# Output: results/tables/supp_table6-8_resistance_trends_<organism>.csv

library(data.table)
library(zoo)
library(glmmTMB)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

organisms <- c("Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa")
antibiotics <- c(
  "Ampicillin", "Cefuroxime", "Amoxicillin/Clavulanic acid", "Cefotaxime",
  "Ceftriaxone", "Nitrofurantoin", "Cefoperazone/Sulbactam", "Levofloxacin",
  "Tobramycin", "Ciprofloxacin", "Gentamicin", "Ceftazidime", "Norfloxacin",
  "Imipenem", "Ofloxacin", "Cefepime", "Aztreonam", "Meropenem", "Amikacin",
  "Piperacillin/Tazobactam", "Fosfomycin", "Tigecycline", "Colistin",
  "Cefixime", "Cefoxitin", "Ertapenem", "Piperacillin",
  "Trimethoprim/Sulfamethoxazole", "Tetracycline", "Ticarcillin")
not_modelled <- list(
  "Klebsiella pneumoniae" = c("Ampicillin", "Ticarcillin"),
  "Pseudomonas aeruginosa" = c("Cefuroxime", "Cefotaxime", "Ceftriaxone", "Cefixime", "Cefoxitin",
                               "Ertapenem", "Tetracycline", "Tigecycline", "Trimethoprim/Sulfamethoxazole",
                               "Nitrofurantoin", "Ampicillin", "Amoxicillin/Clavulanic acid",
                               "Fosfomycin"))   # last one: no CLSI interpretive criteria

ast <- fread("data/ast_clean.csv")
ast[, `:=`(state = tools::toTitleCase(trimws(order_state)), pincode = as.character(trimws(Order_Pincode)))]
ast <- ast[test_name == "URINE" & antibiotic_clean %in% antibiotics &
             !is.na(state) & state != "" & !grepl("^unknown$|^na$", state, ignore.case = TRUE) &
             !is.na(pincode) & pincode != "" & !grepl("^unknown$|^na$", pincode, ignore.case = TRUE)]
ast[, ym := as.yearmon(as.Date(end_state_date))]

control <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))

fit_model <- function(d) {
  slope <- cbind(R, Total - R) ~ seq_num + (1 + seq_num | state)
  intercept <- cbind(R, Total - R) ~ seq_num + (1 | state)
  candidates <- list(
    list(slope, betabinomial("logit"), "beta-binomial", "intercept+slope|state"),
    list(slope, binomial("logit"), "binomial", "intercept+slope|state"),
    list(intercept, betabinomial("logit"), "beta-binomial", "intercept|state"),
    list(intercept, binomial("logit"), "binomial", "intercept|state"))
  for (m in candidates) {
    fit <- tryCatch(glmmTMB(m[[1]], family = m[[2]], data = d, control = control),
                    error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fit) && isTRUE(fit$sdr$pdHess) && is.finite(logLik(fit)) &&
        !any(is.na(summary(fit)$coefficients$cond[, "Std. Error"]))) {
      return(list(fit = fit, model_type = m[[3]], random_structure = m[[4]]))
    }
  }
  NULL
}

for (organism in organisms) {
  counts <- ast[organism_clean == organism,
                .(R = sum(RESISTANT), Total = sum(RESISTANT + INTERMEDIATE + SENSITIVE)),
                by = .(antibiotic = antibiotic_clean, state, pincode, ym)][Total > 0]
  rows <- list()

  for (ab in setdiff(antibiotics, not_modelled[[organism]])) {
    m <- counts[antibiotic == ab]
    if (nrow(m) == 0) next

    # states observed in >=75% of months between first and last observation, with >=2 pincodes
    span <- as.integer(round((max(m$ym) - min(m$ym)) * 12)) + 1L
    covered <- m[, .(share = uniqueN(ym) / max(span, 1L)), by = state]
    m <- m[state %in% covered[share >= 0.75, state]]
    m <- m[state %in% m[, .(pincodes = uniqueN(pincode)), by = state][pincodes >= 2, state]]
    if (nrow(m) == 0) next

    d <- m[, .(R = sum(R), Total = sum(Total)), by = .(state, ym)][Total > 0]
    if (uniqueN(d$state) < 3 || nrow(d) < 12) next
    d[, `:=`(seq_num = as.numeric((ym - min(ym)) * 12), state = factor(state))]

    model <- fit_model(d)
    if (is.null(model)) next

    coefs <- summary(model$fit)$coefficients$cond
    beta <- coefs["seq_num", "Estimate"]
    se <- coefs["seq_num", "Std. Error"]
    p <- fitted(model$fit)
    slope_pp <- mean(p * (1 - p))            # x 1200: logit slope -> percentage points per year

    rows[[ab]] <- data.table(
      organism = organism, antibiotic = ab,
      model_type = model$model_type, random_structure = model$random_structure,
      ppt_per_year = beta * slope_pp * 1200,
      ppt_per_year_CI_low = (beta - 1.96 * se) * slope_pp * 1200,
      ppt_per_year_CI_high = (beta + 1.96 * se) * slope_pp * 1200,
      p_value = coefs["seq_num", "Pr(>|z|)"], p_value_adj_BH = NA_real_,
      mean_resistance_pct = 100 * sum(d$R) / sum(d$Total),
      n_states = uniqueN(d$state), n_pincodes = uniqueN(m$pincode), n_state_months = nrow(d))
  }

  out <- rbindlist(rows)
  out[, p_value_adj_BH := p.adjust(p_value, method = "BH")]
  number <- match(organism, organisms) + 5
  fwrite(out, sprintf("results/tables/supp_table%d_resistance_trends_%s.csv", number,
                      gsub(" ", "_", tolower(organism))))
}
