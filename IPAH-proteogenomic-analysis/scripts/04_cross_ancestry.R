library(dplyr)
library(writexl)

T_Cohort1_Black <- read.csv("./data/Cohorts/T_Cohort1_Black.csv", header = TRUE)
T_Cohort_control_Black <- read.csv("./data/Cohorts/T_Cohort_control_Black.csv", header = TRUE)
T_Cohort1_Black <- bind_rows(T_Cohort1_Black, T_Cohort_control_Black) %>%
  dplyr::mutate(Ethnicity = as.character("Black"))

T_Cohort2_Black <- read.csv("./data/Cohorts/T_Cohort2_Black.csv", header = TRUE)  %>%
  dplyr::mutate(Ethnicity = as.character("Black"))

T_Cohort1_Yellow <- read.csv("./data/Cohorts/T_Cohort1_Yellow.csv", header = TRUE)
T_Cohort_control_Yellow <- read.csv("./data/Cohorts/T_Cohort_control_Yellow.csv", header = TRUE)
T_Cohort1_Yellow <- bind_rows(T_Cohort1_Yellow, T_Cohort_control_Yellow) %>%
  dplyr::mutate(Ethnicity = as.character("Yellow"))

T_Cohort2_Yellow <- read.csv("./data/Cohorts/T_Cohort2_Yellow.csv", header = TRUE)

T_Cohort1_Others <- read.csv("./data/Cohorts/T_Cohort1_Others.csv", header = TRUE)
T_Cohort_control_Others <- read.csv("./data/Cohorts/T_Cohort_control_Others.csv", header = TRUE)
T_Cohort1_Others <- bind_rows(T_Cohort1_Others, T_Cohort_control_Others) %>%
  dplyr::mutate(Ethnicity = as.character("Others"))

T_Cohort2_Others <- read.csv("./data/Cohorts/T_Cohort2_Others.csv", header = TRUE) %>%
  dplyr::mutate(Ethnicity = as.character("Black"))

factor_vars <- c("sex", "ethnicity_C","smoking01", "alcohol01", "pah_dia", "sec_pah_dia",
                 "hypt_dia", "cad_dia", "hf_dia", "af_dia", "lh_dia",
                 "copd_dia", "ild_dia", "osa_dia", "pe_dia",
                 "dm_dia", "ckd_dia",
                 "chd_dia", "ctd_dia", "hiv_dia", "pht_dia", "sch_dia", "hht_dia", 
                 "ERA", "PDE5i", "CCB", "transplantation", "all_death", "pah_death", "cvd_death", 
                 "walk_pace", "walk_pcae_2012", "walk_pace_2014", "walk_pace_2019", 
                 "walk_SOB", "walk_SOB_2012", "walk_SOB_2014", "walk_SOB_2019", 
                 "Comp_workload", "Comp_workload_2012", "VO2_max_moedl", "VO2_max_moedl_2012", 
                 "ACCE_data_avail", "ACCE_data_avail2")

#C1
T_Cohort1_nonWhite <- bind_rows(T_Cohort1_Black, T_Cohort1_Yellow, T_Cohort1_Others) %>%
  dplyr::mutate(across(any_of(factor_vars), as.factor)) 

#C2
T_Cohort2_nonWhite <- bind_rows(T_Cohort2_Black, T_Cohort2_Others) %>%
  dplyr::mutate(across(any_of(factor_vars), as.factor)) 

#C1+C2
T_Cohort_nonWhite <- bind_rows(T_Cohort1_Black, T_Cohort1_Yellow, T_Cohort1_Others, T_Cohort2_Black, T_Cohort2_Others) %>%
  dplyr::mutate(across(any_of(factor_vars), as.factor)) 

table(T_Cohort_nonWhite$pah_dia)

T_cox_MR_cis_cross_C1 <-read.csv("./data/Validation/c1_cox_MR_cis_cross.csv", header = TRUE)
T_cox_MR_cis_cross_C2 <-read.csv("./data/Validation/c2_cox_MR_cis_cross.csv", header = TRUE)

sig_prot <- dplyr::bind_rows(T_cox_MR_cis_cross_C1, T_cox_MR_cis_cross_C2)

sig_pro_list <- unique(tolower(as.character(sig_prot$protein_name)))

match_vars <- c("age", "sex", "Ethnicity", "BMI", "SBP_a", "DBP_a", "smoking01", "alcohol01",
                "hypt_dia", "cad_dia", "hf_dia", "af_dia", "lh_dia",
                "copd_dia", "osa_dia", "pe_dia",

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
T_Cohort_nonWhite_patients <- T_Cohort_nonWhite %>%
  dplyr::filter(pah_dia == 1)
T_Cohort_nonWhite_control <- T_Cohort_nonWhite %>%
  dplyr::filter(pah_dia == 0)

T_nonWhite_match11 <- match_by_matchit(T_Cohort_nonWhite_patients, T_Cohort_nonWhite_control,
                                   covars = match_vars,
                                   ratio = 10,
                                   method = "nearest",
                                   distance = "logit",
                                   caliper = NULL,
                                   seed = 20260305)

T_nonWhite_match11_matched <- T_nonWhite_match11$matched_data

Table_nonWhite_unmatched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_Cohort_nonWhite)
Table_nonWhite_unmatched_summary <- (print(Table_nonWhite_unmatched, noSpaces = TRUE, smd = TRUE))

Table_nonWhite_unmatched_summary <- data.frame(Variable = rownames(Table_nonWhite_unmatched_summary), Table_nonWhite_unmatched_summary)
write_xlsx(Table_nonWhite_unmatched_summary,  "./data/Cohorts/Table_nonWhite_unmatched.xlsx") 

Table_nonWhite_matched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_nonWhite_match11_matched)
Table_nonWhite_matched_summary <- (print(Table_nonWhite_matched, noSpaces = TRUE, smd = TRUE))

Table_nonWhite_matched_summary <- data.frame(Variable = rownames(Table_nonWhite_matched_summary), Table_nonWhite_matched_summary)
write_xlsx(Table_nonWhite_matched_summary,  "./data/Cohorts/Table_nonWhite_matched.xlsx") 

#----C1----

T_nonWhite_match11_matched_to_test <- T_Cohort_nonWhite %>%

  dplyr::select(pah_dia, sex, Ethnicity, any_of(sig_pro_list))

covs_nonWhite <- c("sex","Ethnicity")

logit_nonWhite <- mclapply(sig_pro_list, function(var) {
  formula_str <- paste0("pah_dia ~ ", var, " + ", paste(covs_nonWhite, collapse = " + "))
  logit_model <- glm(as.formula(formula_str), data = T_nonWhite_match11_matched_to_test, family = binomial)
  logit_tidy <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = FALSE)
  logit_tidy_var <- logit_tidy %>% filter(term == var)
  logit_tidy_var$variable_tested <- var
  return(logit_tidy_var)
}, mc.cores = 5)

logit_nonWhite_all <- bind_rows(logit_nonWhite) %>%
  mutate(
    CI_low  = exp(log(estimate) - 1.96 * std.error),
    CI_high = exp(log(estimate) + 1.96 * std.error),
    `OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", estimate, CI_low, CI_high),
    p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

write_xlsx(logit_nonWhite_all, "./data/cross_ancestry/logit_nonWhite_all_matched.xlsx")

negative_proteins <- c("apoc1", "lrrn1")

logit_nonWhite <- mclapply(sig_pro_list, function(var) {
  
  formula_str <- paste0(
    "pah_dia ~ ", var, " + ",
    paste(covs_nonWhite, collapse = " + ")
  )
  
  logit_model <- glm(
    as.formula(formula_str),
    data   = T_nonWhite_match11_matched_to_test,
    family = binomial()
  )
  

  logit_tidy_var <- broom::tidy(
    logit_model,
    exponentiate = FALSE,
    conf.int = FALSE
  ) %>%
    filter(term == var) %>%
    mutate(variable_tested = var)
  
  logit_tidy_var
}, mc.cores = 5)

logit_nonWhite_all <- bind_rows(logit_nonWhite) %>%
  mutate(expected_direction = if_else(
    variable_tested %in% negative_proteins, "negative", "positive"),
    z_value = estimate / std.error,
    p_one_sided = if_else(expected_direction == "positive",
      pnorm(z_value, lower.tail = FALSE), 
      pnorm(z_value, lower.tail = TRUE),
    observed_direction = case_when(estimate > 0 ~ "positive", 
                                   estimate < 0 ~ "negative",
                                   TRUE         ~ "null"),
    direction_concordance = if_else(observed_direction == expected_direction, 
                                    "Concordant",
                                    "Discordant"),
    OR      = exp(estimate),
    CI_low  = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    `OR (95%CI)` = sprintf(
      "%.2f (%.2f-%.2f)",
      OR, CI_low, CI_high
    ),
    p_value_one_sided = if_else(p_one_sided < 0.001,
                                "<0.001",
                                sprintf("%.3f", p_one_sided)
                                )
    )
    )

logit_nonWhite_table <- logit_nonWhite_all %>%
  dplyr::select(
    variable_tested,
    OR,
    CI_low,
    CI_high,
    `OR (95%CI)`,
    p_one_sided,
    p_value_one_sided,
    expected_direction,
    observed_direction,
    direction_concordance
  )

Table_nonW_proteins <- CreateTableOne(vars = sig_pro_list, strata = "pah_dia", data = T_Cohort_nonWhite)
Table_nonW_proteins_summary <- (print(Table_nonW_proteins, noSpaces = TRUE, smd = TRUE))

Table_nonW_proteins_summary <- data.frame(Variable = rownames(Table_nonW_proteins_summary), Table_nonW_proteins_summary)
write_xlsx(Table_nonW_proteins_summary,  "./data/cross_ancestry/Table_cross_ancestry_proteins.xlsx") 

T_overall_B <- read.csv("./data/data_clean/T_PAH_HELTH_Cross_Protein_Black.csv", header = TRUE)

T_overall_Y <- read.csv("./data/data_clean/T_PAH_HELTH_Cross_Protein_Yellow.csv", header = TRUE)

T_overall_O <- read.csv("./data/data_clean/T_PAH_HELTH_Cross_Protein_Others.csv", header = TRUE)

T_overall_BYO <- bind_rows(T_overall_B, T_overall_Y, T_overall_O)  %>%
  dplyr::filter(!(pah_dia == 0 & eid %in% eid_list_to_delete)) %>%

  dplyr::filter(sec_pah_dia == 0)
