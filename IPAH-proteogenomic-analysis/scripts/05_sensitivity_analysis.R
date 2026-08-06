library(dplyr)

T_Cohort_sen <- read.csv("./data/data_clean/T_overall_cohort.csv", header = TRUE) 

T_Cohort_sen <- T_Cohort_sen %>%
  dplyr::filter(sec_pah_dia == 0) %>%

  dplyr::filter(chd_dia == 0) %>%

  dplyr::filter(hiv_dia == 0) %>%

  dplyr::filter(ctd_dia == 0) %>%

  dplyr::filter(pht_dia == 0) %>%

  dplyr::filter(sch_dia == 0) %>%

  dplyr::filter(hht_dia == 0) %>%

  
  dplyr::filter(lh_dia == 0) %>%

  
  dplyr::filter(copd_dia == 0) %>%

  dplyr::filter(ild_dia == 0) %>%

  dplyr::filter(osa_dia == 0) %>%

  
  dplyr::filter(pe_dia == 0) %>%

  dplyr::mutate(across(any_of(factor_vars), as.factor))    

table(T_Cohort_sen$pah_dia)

T_cox_MR_cis_cross_C1 <-read.csv("./data/Validation/c1_cox_MR_cis_cross.csv", header = TRUE)
T_cox_MR_cis_cross_C2 <-read.csv("./data/Validation/c2_cox_MR_cis_cross.csv", header = TRUE)

sig_prot <- dplyr::bind_rows(T_cox_MR_cis_cross_C1, T_cox_MR_cis_cross_C2)

sig_pro_list <- unique(tolower(as.character(sig_prot$protein_name)))

match_vars <- c("age", "sex", "BMI", "SBP_a", "DBP_a", "smoking01", "alcohol01",
                "hypt_dia", "cad_dia", "hf_dia", "af_dia", 
                "dm_dia", "ckd_dia",
                "ALT", "Creatinine", "CRP", "Glucose", "HbA1c", "LDL", "Triglycerides", "nppb", "ntprobnp",
                "walk_pace", "MET_total", "FVC")

match_by_matchit <- function(df1, df2,
                             covars,

                             ratio = 1,

                             method = "nearest",

                             distance = "logit",

                             caliper = NULL,

                             exact = NULL,

                             replace = FALSE,

                             seed = 1,

                             estimand = "ATT"

) {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Please install MatchIt first: install.packages('MatchIt')")
  }
  stopifnot(all(covars %in% names(df1)), all(covars %in% names(df2)))
  if (!is.null(exact)) stopifnot(all(exact %in% names(df1)), all(exact %in% names(df2)))
  
  set.seed(seed)
  
  df1$.treated <- 1
  df2$.treated <- 0
  dat <- rbind(df1, df2)
  

  fml <- as.formula(paste(".treated ~", paste(covars, collapse = " + ")))
  
  m.out <- MatchIt::matchit(
    formula  = fml,
    data     = dat,
    method   = method,
    distance = distance,
    ratio    = ratio,
    replace  = replace,
    caliper  = caliper,
    exact    = exact,
    estimand = estimand
  )
  
  matched <- MatchIt::match.data(m.out)
  list(
    matchit = m.out,
    matched_data = matched
  )
}

#----DO----
T_Cohort_sen_patients <- T_Cohort_sen %>%
  dplyr::filter(pah_dia == 1)
T_Cohort_sen_control <- T_Cohort_sen %>%
  dplyr::filter(pah_dia == 0)

T_sen_match110 <- match_by_matchit(T_Cohort_sen_patients, T_Cohort_sen_control,
                                       covars = match_vars,
                                       ratio = 10,
                                       method = "nearest",
                                       distance = "logit",
                                       caliper = NULL,
                                       seed = 20260305)

T_sen_match110_matched <- T_sen_match110$matched_data

Table_sen_unmatched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_Cohort_sen)
Table_sen_unmatched_summary <- (print(Table_sen_unmatched, noSpaces = TRUE, smd = TRUE))

Table_sen_unmatched_summary <- data.frame(Variable = rownames(Table_sen_unmatched_summary), Table_sen_unmatched_summary)
write_xlsx(Table_sen_unmatched_summary,  "./data/sensitivity/Table_sen_unmatched.xlsx") 

Table_sen_matched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_sen_match110_matched)
Table_sen_matched_summary <- (print(Table_sen_matched, noSpaces = TRUE, smd = TRUE))

Table_sen_matched_summary <- data.frame(Variable = rownames(Table_sen_matched_summary), Table_sen_matched_summary)
write_xlsx(Table_sen_matched_summary,  "./data/sensitivity/Table_sen_matched.xlsx") 

T_sen_match110_matched_male <- T_sen_match110_matched %>%
  filter(sex == 1)
T_sen_match110_matched_female <- T_sen_match110_matched %>%
  filter(sex == 0)

covs_sen = c("hf_dia")
logit_sensitive <- mclapply(sig_pro_list, function(var) {
  formula_str <- paste0("pah_dia ~ ", var, " + ", paste(covs_sen, collapse = " + "))
  logit_model <- glm(as.formula(formula_str), data = T_sen_match110_matched_female, family = binomial)
  logit_tidy <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = FALSE)
  logit_tidy_var <- logit_tidy %>% filter(term == var)
  logit_tidy_var$variable_tested <- var
  return(logit_tidy_var)
}, mc.cores = 1)

logit_sensitive_all <- bind_rows(logit_sensitive) %>%
  mutate(
    CI_low  = exp(log(estimate) - 1.96 * std.error),
    CI_high = exp(log(estimate) + 1.96 * std.error),
    `OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", estimate, CI_low, CI_high),
    p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

write_xlsx(logit_sensitive_all, "./data/sensitivity/logit_sensitive_all_female_matched.xlsx")

Table_sen_proteins <- CreateTableOne(vars = sig_pro_list, strata = "pah_dia", data = T_sen_match110_matched_female)
Table_sen_proteins_summary <- (print(Table_sen_proteins, noSpaces = TRUE, smd = TRUE))

Table_sen_proteins_summary <- data.frame(Variable = rownames(Table_sen_proteins_summary), Table_sen_proteins_summary)
write_xlsx(Table_sen_proteins_summary,  "./data/sensitivity/Table_sen_female_proteins.xlsx") 
  
