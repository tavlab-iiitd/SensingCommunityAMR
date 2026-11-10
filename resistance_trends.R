source("R/setup.R")
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

fit_resistance_model <- function(d, context) {
  slope <- cbind(R, Total - R) ~ seq_num + (1 + seq_num | state)
  intercept <- cbind(R, Total - R) ~ seq_num + (1 | state)
  candidates <- list(
    list(slope, betabinomial("logit"), "beta-binomial", "intercept+slope|state"),
    list(slope, binomial("logit"), "binomial", "intercept+slope|state"),
    list(intercept, betabinomial("logit"), "beta-binomial", "intercept|state"),
    list(intercept, binomial("logit"), "binomial", "intercept|state"))
  for (m in candidates) {
    message <- ""
    fit <- tryCatch(glmmTMB(m[[1]], family = m[[2]], data = d, control = control),
      error = function(e) { message <<- conditionMessage(e); NULL },
      warning = function(w) { message <<- conditionMessage(w); NULL })
    record_fit("results/checks/resistance_fit_attempts.csv", context, m[[3]], m[[4]], fit, message)
    if (!is.null(fit) && isTRUE(fit$sdr$pdHess) && is.finite(logLik(fit)) &&
        !any(is.na(summary(fit)$coefficients$cond[, "Std. Error"]))) {
      return(list(fit = fit, model_type = m[[3]], random_structure = m[[4]]))
    }
  }
  NULL
}

source("R/monthly_marginal_effect.R")
if (file.exists("results/checks/resistance_fit_attempts.csv")) file.remove("results/checks/resistance_fit_attempts.csv")
resistance_tables <- list()
resistance_eligibility <- list()
for (organism in organisms) {
  counts <- ast[organism_clean == organism,
    .(R = sum(RESISTANT), Total = sum(RESISTANT + INTERMEDIATE + SENSITIVE)),
    by = .(antibiotic = antibiotic_clean, state, pincode, ym)][Total > 0]
  revised_rows <- state_rows <- list()
  for (ab in setdiff(antibiotics, not_modelled[[organism]])) {
    m <- counts[antibiotic == ab]
    if (!nrow(m)) next
    span <- as.integer(round((max(m$ym) - min(m$ym)) * 12)) + 1L
    covered <- m[, .(share = uniqueN(ym) / max(span, 1L)), by = state]
    m <- m[state %in% covered[share >= .75, state]]
    m <- m[state %in% m[, .(pincodes = uniqueN(pincode)), by = state][pincodes >= 2, state]]
    d <- m[, .(R = sum(R), Total = sum(Total)), by = .(state, ym)][Total > 0]
    eligible <- uniqueN(d$state) >= 3 && nrow(d) >= 12
    resistance_eligibility[[paste(organism, ab)]] <- data.table(
      organism = organism, antibiotic = ab, n_states = uniqueN(d$state),
      n_state_months = nrow(d), eligible = eligible)
    if (!eligible) next
    d[, `:=`(seq_num = as.numeric((ym - min(ym)) * 12), state = factor(state))]
    stopifnot(all(d$R >= 0), all(d$Total >= d$R))
    context <- paste(organism, ab, sep = " | ")
    cat("Fitting", context, "\n"); flush.console()
    model <- fit_resistance_model(d, context)
    if (is.null(model)) stop("All candidate models failed: ", context)
    fit <- model$fit
    stopifnot(fit$fit$convergence == 0)
    effect <- monthly_marginal_effect(fit, d, "seq_num")
    common <- data.table(organism = organism, antibiotic = ab,
      model_type = model$model_type, random_structure = model$random_structure)
    details <- data.table(mean_resistance_pct = 100 * sum(d$R) / sum(d$Total),
      n_states = uniqueN(d$state), n_pincodes = uniqueN(m$pincode), n_state_months = nrow(d))
    revised_rows[[ab]] <- cbind(common, data.table(
      ppt_per_month = effect$estimate, ppt_per_month_CI_low = effect$low,
      ppt_per_month_CI_high = effect$high, p_value = effect$p_value), details)

    state_rows[[ab]] <- cbind(data.table(organism = organism, antibiotic = ab), effect$state_summary)
    save_model(fit, d, paste("resistance", organism, ab))
    last_resistance_model <- fit
    last_resistance_data <- d
    cat(sprintf("  %.5f pp/month [%.5f, %.5f]\n", effect$estimate, effect$low, effect$high))
  }
  revised <- rbindlist(revised_rows)
  revised[, p_value_adj_BH := p.adjust(p_value, "BH")]
  slug <- gsub(" ", "_", tolower(organism))
  number <- match(organism, organisms) + 5
  stem <- sprintf("results/tables/supp_table%d_resistance_trends_%s", number, slug)
  fwrite(revised, paste0(stem, "_monthly.csv"))
  fwrite(rbindlist(state_rows), paste0("results/checks/resistance_state_effects_", slug, ".csv"))
  resistance_tables[[organism]] <- revised
  print(revised)
}
fwrite(rbindlist(resistance_eligibility), "results/checks/resistance_eligibility.csv")
writeLines(capture.output(sessionInfo()), "results/checks/R_session_resistance.txt")
