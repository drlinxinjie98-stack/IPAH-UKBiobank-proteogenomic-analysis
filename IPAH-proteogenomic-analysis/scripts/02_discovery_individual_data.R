library(dplyr)
library(MASS)
library(nnet)
library(survival)
library(broom) 
library(lubridate)
library(MatchIt)
library(optmatch)
library(tableone)
library(parallel)
library(writexl)

T_overall <- read.csv("./data/data_clean/T_PAH_HELTH_Cross_Protein.csv", header = TRUE)

table(T_overall$pah_dia)
T_protein_delBNP <- read.csv("./data/data_clean/T_protein.csv", header = TRUE) %>%
                    dplyr::select(-nppb, -ntprobnp)

eid_list_to_delete <- read.csv("./data/ukb_delete.csv", header = TRUE) %>% 
  filter(p30903_i0 == 1) %>%
  pull(eid)

date_vars <- c("date_join", "date_death", "Date_transplantation", "Date_pah", "Date_sec_pah", 
               "Date_hpyt", "Date_cad", "Date_hf", "Date_af", "Date_lh",
               "Date_copd", "Date_ild", "Date_osa", "Date_pe",
               "Date_dm", "Date_ckd", 
               "Date_chd", "Date_ctd", "Date_hiv", "Date_pht", "Date_sch", "Date_hht")

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

char_vars <- c("eid")

which(colnames(T_overall) == "a1bg"); which(colnames(T_overall) == "zpr1")

T_overall_cohort <- T_overall %>%
  dplyr::filter(!(pah_dia == 0 & eid %in% eid_list_to_delete)) %>%

  dplyr::filter(!(transplantation == 1 & Date_transplantation < date_join)) %>%

  dplyr::mutate(across(any_of(char_vars), as.character)) %>%

  dplyr::mutate(across(where(is.character), ~na_if(., ""))) %>%

  dplyr::mutate(across(any_of(date_vars), as.Date)) %>%

  
  dplyr::mutate(hypt_dia = as.integer(!is.na(Date_hypt) & Date_hypt <= date_join),

                cad_dia  = as.integer(!is.na(Date_cad)  & Date_cad  <= date_join),
                hf_dia   = as.integer(!is.na(Date_hf)   & Date_hf   <= date_join),
                af_dia   = as.integer(!is.na(Date_af)   & Date_af   <= date_join),
                lh_dia   = as.integer(!is.na(Date_lh)   & Date_lh   <= date_join),
                
                
                copd_dia = as.integer(!is.na(Date_copd) & Date_copd <= date_join),
                ild_dia  = as.integer(!is.na(Date_ild)  & Date_ild  <= date_join),
                osa_dia  = as.integer(!is.na(Date_osa)  & Date_osa  <= date_join),
                pe_dia   = as.integer(!is.na(Date_pe)   & Date_pe   <= date_join),
                
                ckd_dia  = as.integer(!is.na(Date_ckd)  & Date_ckd  <= date_join),
                dm_dia   = as.integer(!is.na(Date_dm)   & Date_dm   <= date_join),
                
                ctd_dia  = as.integer(!is.na(Date_ctd)  & Date_ctd  <= date_join),
                hiv_dia  = as.integer(!is.na(Date_hiv)  & Date_hiv  <= date_join),
                pht_dia  = as.integer(!is.na(Date_pht)  & Date_pht  <= date_join),
                sch_dia  = as.integer(!is.na(Date_sch)  & Date_sch  <= date_join),
                hht_dia  = as.integer(!is.na(Date_hht)  & Date_hht  <= date_join)) %>%
  
  dplyr::mutate(all_death  = case_when(transplantation == 1 ~ 1L,
                               TRUE ~ all_death),
                date_death = case_when(transplantation == 1 ~ Date_transplantation,
                               TRUE ~ date_death)) %>%

  dplyr::mutate(Date_pah = case_when(
                    pah_dia == 0 & is.na(Date_pah) ~ as.Date("2024-12-01"),
                    TRUE ~ Date_pah),

         FU_pah_dia = time_length(interval(date_join, Date_pah), "month")) %>%

  dplyr::mutate(date_death = case_when(
                    all_death == 0 & is.na(date_death) ~ as.Date("2024-12-01"),
                    TRUE ~ date_death),

         FU_death   = time_length(interval(date_join, date_death), "month")) %>%

  dplyr::mutate(across(any_of(factor_vars), as.factor)) %>%

  {
    num_vars <- setdiff(names(.), c(char_vars, date_vars, factor_vars))
    mutate(., across(all_of(num_vars), ~as.numeric(as.character(.))))
  } %>%

  dplyr::select(eid, date_join, age, sex, BMI, SBP_a, DBP_a, smoking01, alcohol01, 
                pah_dia, Date_pah, FU_pah_dia, sec_pah_dia,

                hypt_dia, cad_dia, hf_dia, af_dia, lh_dia,

                copd_dia, ild_dia, osa_dia, pe_dia, 
                dm_dia, ckd_dia, 
                chd_dia, ctd_dia, hiv_dia, pht_dia, sch_dia, hht_dia,

                ALT, Creatinine, CRP, Glucose, HbA1c, LDL, Triglycerides, nppb, ntprobnp,

                walk_pace, MET_total, FVC,
                ERA, PDE5i, CCB,
                all_death, date_death, FU_death)

T_overall_cohort <- merge(T_overall_cohort, T_protein_delBNP, by = 'eid', all.x = TRUE, all.y = FALSE)
#write.csv(T_overall_cohort, "./data/data_clean/T_overall_cohort.csv", row.names = FALSE)

table(T_overall_cohort$pah_dia)

#---------------------------------Cohort1---------------------------------------
T_Cohort1 <- subset(T_overall_cohort, pah_dia == 1 & (is.na(Date_pah) | Date_pah > date_join) & 
                   (sec_pah_dia == 0))
table(T_Cohort1$sex)

#write.csv(T_Cohort1, "./data/Cohorts/T_Cohort1.csv", row.names = FALSE)
#write.csv(T_Cohort1, "./data/Cohorts/T_Cohort1_Black.csv", row.names = FALSE)
#write.csv(T_Cohort1, "./data/Cohorts/T_Cohort1_Yellow.csv", row.names = FALSE)
#write.csv(T_Cohort1, "./data/Cohorts/T_Cohort1_Others.csv", row.names = FALSE)

#---------------------------------Cohort2---------------------------------------
T_Cohort2 <- subset(T_overall_cohort, pah_dia == 1 & (is.na(Date_pah) | Date_pah < date_join) & 
                      (sec_pah_dia == 0))
table(T_Cohort2$pah_dia) # (n=141)
sum(T_Cohort2$sex == 0 & T_Cohort2$all_death == 1, na.rm = TRUE)

#write.csv(T_Cohort2, "./data/Cohorts/T_Cohort2.csv", row.names = FALSE)
#write.csv(T_Cohort2, "./data/Cohorts/T_Cohort2_Black.csv", row.names = FALSE)
#write.csv(T_Cohort2, "./data/Cohorts/T_Cohort2_Yellow.csv", row.names = FALSE)
#write.csv(T_Cohort2, "./data/Cohorts/T_Cohort2_Others.csv", row.names = FALSE)

T_Cohort2_dead <- T_Cohort2 %>%
  dplyr::filter(all_death == 1) # (n=60)
  

T_Cohort2_surv <- T_Cohort2 %>%
  dplyr::filter(all_death == 0) # (n=81)

#---------------------------------Cohort_health_control-------------------------
T_Cohort_control <- subset(T_overall_cohort, pah_dia == 0 & sec_pah_dia == 0)

table(T_Cohort_control$sex)

#write.csv(T_Cohort_control, "./data/Cohorts/T_Cohort_control.csv", row.names = FALSE)
#write.csv(T_Cohort_control, "./data/Cohorts/T_Cohort_control_Black.csv", row.names = FALSE)
#write.csv(T_Cohort_control, "./data/Cohorts/T_Cohort_control_Yellow.csv", row.names = FALSE)
#write.csv(T_Cohort_control, "./data/Cohorts/T_Cohort_control_Others.csv", row.names = FALSE)

#================================Discovery======================================

match_vars <- c("age", "sex", "BMI", "SBP_a", "DBP_a", "smoking01", "alcohol01",
                "hypt_dia", "cad_dia", "hf_dia", "af_dia", "lh_dia",
                "copd_dia", "ild_dia", "osa_dia", "pe_dia",
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
#===================================Cohort1=====================================

T_C1_C_match14 <- match_by_matchit(T_Cohort1, T_Cohort_control,
                                   covars = match_vars,
                                   ratio = 4,
                                   method = "nearest",
                                   distance = "logit",
                                   caliper = NULL,
                                   seed = 20260305)

T_C1_C_match14_matched <- T_C1_C_match14$matched_data

write.csv(T_C1_C_match14_matched, "./data/Cohorts/T_C1_C_match14_matched.csv", row.names = FALSE)

T_C1_C_unmatched <- rbind(T_Cohort1, T_Cohort_control)

Table_C1_C_unmatched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_C1_C_unmatched)
Table_C1_C_unmatched_summary <- (print(Table_C1_C_unmatched, noSpaces = TRUE, smd = TRUE))

Table_C1_C_unmatched_summary <- data.frame(Variable = rownames(Table_C1_C_unmatched_summary), Table_C1_C_unmatched_summary)
write_xlsx(Table_C1_C_unmatched_summary,  "./data/Cohorts/Table_C1_C_unmatched.xlsx") 

Table_C1_C_matched <- CreateTableOne(vars = match_vars, strata = "pah_dia", data = T_C1_C_match14_matched)
Table_C1_C_matched_summary <- (print(Table_C1_C_matched, noSpaces = TRUE, smd = TRUE))

Table_C1_C_matched_summary <- data.frame(Variable = rownames(Table_C1_C_matched_summary), Table_C1_C_matched_summary)
write_xlsx(Table_C1_C_matched_summary,  "./data/Cohorts/Table_C1_C_matched.xlsx")

#----cox(fdr p < 0.05)----

T_C1_C_match14_matched <- read.csv('./data/Cohorts/T_C1_C_match14_matched.csv', header = TRUE)

covs_C1 <- c("hf_dia")

which(colnames(T_C1_C_match14_matched) == "a1bg"); which(colnames(T_C1_C_match14_matched) == "zpr1")
key_vars <- names(T_C1_C_match14_matched)[49:2966]

T_C1_C_match14_matched <- T_C1_C_match14_matched %>%
  dplyr::mutate(pah_dia = as.integer(as.character(pah_dia)),
                FU_pah_dia = as.numeric(FU_pah_dia))

cox_c1 <- mclapply(key_vars, function(var) {
  formula_str <- paste0("Surv(FU_pah_dia, pah_dia) ~ ", var, " + ", paste(covs_C1, collapse = " + "))

  cox_model <- coxph(as.formula(formula_str), data = T_C1_C_match14_matched)
  cox_tidy <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  cox_tidy_var <- cox_tidy %>% filter(term == var)
  cox_tidy_var$variable_tested <- var
  return(cox_tidy_var)
}, mc.cores = 5)

cox_c1_all <- bind_rows(cox_c1) %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::rename_with(~ paste0(.x, "_cox"), -term)

write.csv(cox_c1_all, "./data/Discovery/Cohort1/all_cox_c1.csv", row.names = TRUE)

cox_c1_sig <- cox_c1_all %>% 
  dplyr::filter(p_adj_fdr_cox < 0.05) %>%

  dplyr::mutate(direction_cox = case_when(estimate_cox > 1 ~ "positive", 
                                          TRUE ~ "negative")) 

cox_c1_sig_list <- unique(as.character(cox_c1_sig$term))

write.csv(cox_c1_sig, "./data/Discovery/Cohort1/all_cox_c1_sig.csv", row.names = TRUE)

#----logistic(BH p < 0.05)----

logit_c1 <- mclapply(key_vars, function(var) {
  formula_str <- paste0("pah_dia ~ ", var, " + ", paste(covs_C1, collapse = " + "))
  logit_model <- glm(as.formula(formula_str), data = T_C1_C_match14_matched, family = binomial)
  logit_tidy <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = FALSE)
  logit_tidy_var <- logit_tidy %>% filter(term == var)
  logit_tidy_var$variable_tested <- var
  return(logit_tidy_var)
}, mc.cores = 5)

logit_c1_all <- bind_rows(logit_c1) %>% 
  dplyr::rename_with(~ paste0(.x, "_logit"), -term) %>% 
  dplyr::mutate(direction_logit = case_when(estimate_logit > 1 ~ "positive", 
                                            TRUE ~ "negative")) 
write.csv(logit_c1_all, "./data/Discovery/Cohort1/all_logit_c1.csv", row.names = TRUE)

logit_c1_sig <- logit_c1_all %>%
  dplyr::filter(term %in% cox_c1_sig_list) %>%
  dplyr::mutate(p_adj_BH_logit = p.adjust(p.value_logit, method = "BH")) %>%
  dplyr::filter(p_adj_BH_logit < 0.05)

logit_c1_verified_sig_list <- unique(as.character(logit_c1_sig$term))

write.csv(logit_c1_sig, "./data/Discovery/Cohort1/all_logit_c1_sig.csv", row.names = TRUE)

T_discover_C1_overall <- merge(cox_c1_sig, logit_c1_sig, by = 'term', all.x = FALSE, all.y = FALSE) %>%
  dplyr::filter(!is.na(direction_cox), !is.na(direction_logit), direction_cox == direction_logit)

write.csv(T_discover_C1_overall, "./data/Discovery/Cohort1/all_cox_logit_c1_sig.csv", row.names = FALSE)

#===================================Cohort2=====================================

T_C2_match11 <- match_by_matchit(T_Cohort2_dead, T_Cohort2_surv,
                                   covars = match_vars,
                                   ratio = 1,
                                   method = "nearest",
                                   distance = "logit",
                                   caliper = NULL,
                                   seed = 20260305)

T_C2_match11_matched <- T_C2_match11$matched_data

write.csv(T_C2_match11_matched, "./data/Cohorts/T_C2_match11_matched.csv", row.names = FALSE)

Table_C2_unmatched <- CreateTableOne(vars = match_vars, strata = "all_death", data = T_Cohort2)
Table_C2_unmatched_summary <- (print(Table_C2_unmatched, noSpaces = TRUE, smd = TRUE))

Table_C2_unmatched_summary <- data.frame(Variable = rownames(Table_C2_unmatched_summary), Table_C2_unmatched_summary)
write_xlsx(Table_C2_unmatched_summary,  "./data/Cohorts/Table_C2_unmatched.xlsx") 

Table_C2_matched <- CreateTableOne(vars = match_vars, strata = "all_death", data = T_C2_match11_matched)
Table_C2_matched_summary <- (print(Table_C2_matched, noSpaces = TRUE, smd = TRUE))

Table_C2_matched_summary <- data.frame(Variable = rownames(Table_C2_matched_summary), Table_C2_matched_summary)
write_xlsx(Table_C2_matched_summary,  "./data/Cohorts/Table_C2_matched.xlsx")

#----cox(fdr < 0.05)----

T_C2_match11_matched <- read.csv('./data/Cohorts/T_C2_match11_matched.csv', header = TRUE)

covs_C2 <- c("walk_pace", "ntprobnp")

which(colnames(T_C2_match11_matched) == "a1bg"); which(colnames(T_C2_match11_matched) == "zpr1")
key_vars <- names(T_C2_match11_matched)[49:2966]

T_C2_match11_matched <- T_C2_match11_matched %>%
  dplyr::mutate(all_death = as.integer(as.character(all_death)),
                FU_death = as.numeric(FU_death))

cox_c2 <- mclapply(key_vars, function(var) {
  formula_str <- paste0("Surv(FU_death, all_death) ~ ", var, " + ", paste(covs_C2, collapse = " + "))

  cox_model <- coxph(as.formula(formula_str), data = T_C2_match11_matched)
  cox_tidy <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  cox_tidy_var <- cox_tidy %>% filter(term == var)
  cox_tidy_var$variable_tested <- var
  return(cox_tidy_var)
}, mc.cores = 5)

cox_c2_all <- bind_rows(cox_c2) %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::rename_with(~ paste0(.x, "_cox"), -term)

write.csv(cox_c2_all, "./data/Discovery/Cohort2/all_cox_c2.csv", row.names = FALSE)

cox_c2_sig <- cox_c2_all %>% 
  dplyr::filter(p_adj_fdr_cox < 0.05) %>%

  dplyr::mutate(direction_cox = case_when(estimate_cox > 1 ~ "positive", 
                                          TRUE ~ "negative")) 

cox_c2_sig_list <- unique(as.character(cox_c2_sig$term))

write.csv(cox_c2_sig, "./data/Discovery/Cohort2/all_cox_c2_sig.csv", row.names = FALSE)

#----logistic(BH < 0.05)----

logit_c2 <- mclapply(key_vars, function(var) {
  formula_str <- paste0("all_death ~ ", var, " + ", paste(covs_C2, collapse = " + "))
  logit_model <- glm(as.formula(formula_str), data = T_C2_match11_matched, family = binomial)
  logit_tidy <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = FALSE)
  logit_tidy_var <- logit_tidy %>% filter(term == var)
  logit_tidy_var$variable_tested <- var
  return(logit_tidy_var)
}, mc.cores = 5)

logit_c2_all <- bind_rows(logit_c2) %>% 
  dplyr::rename_with(~ paste0(.x, "_logit"), -term) %>% 
  dplyr::mutate(direction_logit = case_when(estimate_logit > 1 ~ "positive", 
                                            TRUE ~ "negative")) 
write.csv(logit_c2_all, "./data/Discovery/Cohort2/all_logit_c2.csv", row.names = FALSE)

logit_c2_sig <- logit_c2_all %>%
  dplyr::filter(term %in% cox_c2_sig_list) %>%
  dplyr::mutate(p_adj_BH_logit = p.adjust(p.value_logit, method = "BH")) %>%
  dplyr::filter(p_adj_BH_logit < 0.05)

logit_c2_verified_sig_list <- unique(as.character(logit_c2_sig$term))

write.csv(logit_c2_sig, "./data/Discovery/Cohort2/all_logit_c2_sig.csv", row.names = FALSE)

T_discover_C2_overall <- merge(cox_c2_sig, logit_c2_sig, by = 'term', all.x = FALSE, all.y = FALSE) %>%
  dplyr::filter(!is.na(direction_cox), !is.na(direction_logit), direction_cox == direction_logit)

write.csv(T_discover_C2_overall, "./data/Discovery/Cohort2/all_cox_logit_c2_sig.csv", row.names = FALSE)

#===================================Cohort_sex==================================

#cohort1_male
T_C1_C_match14_matched_male <- T_C1_C_match14_matched %>%
  dplyr::filter(sex == 1)
table(T_C1_C_match14_matched_male$pah_dia)
#cohort1_female
T_C1_C_match14_matched_female <- T_C1_C_match14_matched %>%
  dplyr::filter(sex == 0)
table(T_C1_C_match14_matched_female$pah_dia)

#cohort2_male51853
T_C2_match11_matched_male <- T_C2_match11_matched %>%
  dplyr::filter(sex == 1)
#cohort2_female
T_C2_match11_matched_female <- T_C2_match11_matched %>%
  dplyr::filter(sex == 0)

run_cox_logit_discovery <- function(
    dat,
    match_vars,

    time_var,          # e.g. "FU_death"
    event_var,         # e.g. "all_death" (0/1)

    protein_range = c("a1bg", "zpr1"),

    protein_vars = NULL,

    mc.cores = 5,

    cox_fdr_cutoff = 0.05,
    logit_bh_cutoff = 0.05,

    prefix1 = "male",
    prefix2 = "C2",
    out_dir = NULL,

    write_all = TRUE,

    write_sig = TRUE

) {
  stopifnot(is.data.frame(dat))
  stopifnot(all(match_vars %in% names(dat)))
  stopifnot(time_var %in% names(dat), event_var %in% names(dat))
  

  if (is.null(protein_vars)) {
    stopifnot(all(protein_range %in% names(dat)))
    i1 <- which(names(dat) == protein_range[1])[1]
    i2 <- which(names(dat) == protein_range[2])[1]
    if (i1 > i2) stop("Invalid protein_range: the start column follows the end column.")
    key_vars <- names(dat)[i1:i2]
  } else {
    stopifnot(all(protein_vars %in% names(dat)))
    key_vars <- protein_vars
  }
  

  dat <- dat |>
    dplyr::mutate(
      .event_num = as.integer(as.character(.data[[event_var]])),
      .time_num  = as.numeric(.data[[time_var]])
    )
  

  cox_list <- parallel::mclapply(key_vars, function(var) {
    tryCatch({
      fml <- as.formula(paste0(
        "survival::Surv(.time_num, .event_num) ~ `", var, "` + ",
        paste0("`", match_vars, "`", collapse = " + ")
      ))
      fit <- survival::coxph(fml, data = dat)
      
      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
        dplyr::filter(.data$term == var | .data$term == paste0("`", var, "`")) |>
        dplyr::mutate(variable_tested = var)
    }, error = function(e) {
      NULL
    })
  }, mc.cores = mc.cores)
  
  cox_all <- dplyr::bind_rows(cox_list) |>
    dplyr::mutate(p_adj_fdr = p.adjust(.data$p.value, method = "fdr")) |>
    dplyr::rename_with(~ paste0(.x, "_cox"), -term)
  
  cox_sig <- cox_all |>
    dplyr::filter(.data$p.value_cox < cox_fdr_cutoff) |> # nominal p < 0.01 as the cut-off 
    dplyr::mutate(direction_cox = dplyr::case_when(
      .data$estimate_cox > 1 ~ "positive",
      TRUE ~ "negative"
    ))
  
  cox_sig_list <- unique(as.character(cox_sig$term))
  

  logit_list <- parallel::mclapply(key_vars, function(var) {
    tryCatch({
      fml <- as.formula(paste0(
        ".event_num ~ `", var, "` + ",
        paste0("`", match_vars, "`", collapse = " + ")
      ))
      fit <- stats::glm(fml, data = dat, family = stats::binomial())
      
      broom::tidy(fit, exponentiate = TRUE, conf.int = FALSE) |>
        dplyr::filter(.data$term == var | .data$term == paste0("`", var, "`")) |>
        dplyr::mutate(variable_tested = var)
    }, error = function(e) {
      NULL
    })
  }, mc.cores = mc.cores)
  
  logit_all <- dplyr::bind_rows(logit_list) |>
    dplyr::rename_with(~ paste0(.x, "_logit"), -term) |>
    dplyr::mutate(direction_logit = dplyr::case_when(
      .data$estimate_logit > 1 ~ "positive",
      TRUE ~ "negative"
    ))
  
  logit_sig <- logit_all |>
    dplyr::filter(.data$term %in% cox_sig_list) |>
    dplyr::mutate(p_adj_BH_logit = p.adjust(.data$p.value_logit, method = "BH")) |>
    dplyr::filter(.data$p.value_logit < logit_bh_cutoff) # nominal p < 0.01 as the cut-off 
  
  logit_verified_sig_list <- unique(as.character(logit_sig$term))
  

  discover_overall <- merge(cox_sig, logit_sig, by = "term", all = FALSE) |>
    dplyr::filter(!is.na(.data$direction_cox), !is.na(.data$direction_logit),
                  .data$direction_cox == .data$direction_logit)
  

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    if (write_all) {
      utils::write.csv(cox_all, file.path(out_dir, paste0(prefix1, "cox_", prefix2, ".csv")), row.names = FALSE)
      utils::write.csv(logit_all, file.path(out_dir, paste0(prefix1, "logit_", prefix2, ".csv")), row.names = FALSE)
    }
    if (write_sig) {
      utils::write.csv(cox_sig, file.path(out_dir, paste0(prefix1, "cox_", prefix2, "_sig.csv")), row.names = FALSE)
      utils::write.csv(logit_sig, file.path(out_dir, paste0(prefix1, "logit_", prefix2, "_sig.csv")), row.names = FALSE)
      utils::write.csv(discover_overall, file.path(out_dir, paste0(prefix1, "cox_logit_", prefix2, "_sig.csv")), row.names = FALSE)
    }
  }
  
  list(
    key_vars = key_vars,
    cox_all = cox_all,
    cox_sig = cox_sig,
    cox_sig_list = cox_sig_list,
    logit_all = logit_all,
    logit_sig = logit_sig,
    logit_verified_sig_list = logit_verified_sig_list,
    discover_overall = discover_overall
  )
}

#C1_male/C1_female/C2_male/C2_female
cox_logit_discovery <- run_cox_logit_discovery(
  dat = T_C2_match11_matched_female,
  match_vars = covs_C2,
  time_var = "FU_death",
  event_var = "all_death",
  protein_range = c("a1bg", "zpr1"),
  mc.cores = 5,
  cox_fdr_cutoff = 0.01,    # (incidence) nominal p < 0.01; (mortality) nominal p < 0.01
  logit_bh_cutoff = 0.01,   
  prefix1 = "female_",
  prefix2 = "c2",
  out_dir = "./data/Discovery/Cohort2"
)

#C1_overall;C1_male;C1_female
T_discover_C1_overall <- read.csv('./data/Discovery/Cohort1/all_cox_logit_c1_sig.csv', header = TRUE) %>%
  mutate(source = "overall")
T_discover_C1_male    <- read.csv('./data/Discovery/Cohort1/male_cox_logit_c1_sig.csv', header = TRUE) %>%
  mutate(source = "male")
T_discover_C1_female  <- read.csv('./data/Discovery/Cohort1/female_cox_logit_c1_sig.csv', header = TRUE) %>%
  mutate(source = "female")

#C2_overall;C2_male;C2_female
T_discover_C2_overall <- read.csv('./data/Discovery/Cohort2/all_cox_logit_c2_sig.csv', header = TRUE) %>%
  mutate(source = "overall")
T_discover_C2_male    <- read.csv('./data/Discovery/Cohort2/male_cox_logit_c2_sig.csv', header = TRUE) %>%
  mutate(source = "male")
T_discover_C2_female  <- read.csv('./data/Discovery/Cohort2/female_cox_logit_c2_sig.csv', header = TRUE) %>%
  mutate(source = "female")

#C1
T_discovery_C1_sigP <- bind_rows(T_discover_C1_overall, T_discover_C1_male, T_discover_C1_female) %>%
  mutate(source = factor(source, levels = c("overall", "male", "female"))) %>%
  arrange(term, source) %>%          
  distinct(term, .keep_all = TRUE)   

write.csv(T_discovery_C1_sigP, "./data/Discovery/Union_P/T_discovery_C1_union_P.csv", row.names = FALSE)

#C2
T_discovery_C2_sigP <- bind_rows(T_discover_C2_overall, T_discover_C2_female) %>%

  mutate(source = factor(source, levels = c("overall", "female"))) %>%
  arrange(term, source) %>%          
  distinct(term, .keep_all = TRUE)  

write.csv(T_discovery_C2_sigP, "./data/Discovery/Union_P/T_discovery_C2_union_P.csv", row.names = FALSE)

T_protein_unimput <- read.csv("./data/Original/protein.csv", header = TRUE)
T_protein_unimput_select <- T_protein_unimput %>%
  dplyr::select(eid, edn1, anxa2, spint1, ltbp2, thbs2, col4a1, qsox1, c7, slamf7, igsf3, lrrn1, nptx1)

T_Cohort1_unimp <- T_Cohort1 %>%
  dplyr::select(1:48) %>%
  dplyr::mutate(eid = as.character(eid)) %>%
  left_join(T_protein_unimput_select %>%
  dplyr::mutate(eid = as.character(eid)),by = "eid")

T_Cohort2_unimp <- T_Cohort2 %>%
  dplyr::select(1:48) %>%
  dplyr::mutate(eid = as.character(eid)) %>%
  left_join(T_protein_unimput_select %>%
              dplyr::mutate(eid = as.character(eid)),by = "eid")

T_Cohort2_unimp_dead <- T_Cohort2_unimp %>%
  dplyr::filter(all_death == 1) # (n=60)

T_Cohort2_unimp_surv <- T_Cohort2_unimp %>%
  dplyr::filter(all_death == 0) # (n=81)

T_Cohortcontrol_unimp <- T_Cohort_control %>%
  dplyr::select(1:48) %>%
  dplyr::mutate(eid = as.character(eid)) %>%
  left_join(T_protein_unimput_select %>%
              dplyr::mutate(eid = as.character(eid)),by = "eid")

#----C1----
T_C1_C_match14_unimp <- match_by_matchit(T_Cohort1_unimp, T_Cohortcontrol_unimp,
                                   covars = match_vars,
                                   ratio = 4,
                                   method = "nearest",
                                   distance = "logit",
                                   caliper = NULL,
                                   seed = 20260305)

T_C1_C_match14_unimp_matched <- T_C1_C_match14_unimp$matched_data

covs_C1 <- c("hf_dia")

key_vars <- c("edn1", "anxa2", "spint1", "ltbp2", "thbs2", "col4a1", "qsox1", "c7", "slamf7", "igsf3", "lrrn1", "nptx1")

T_C1_C_match14_unimp_matched <- T_C1_C_match14_unimp_matched %>%
  dplyr::mutate(
    pah_dia = as.integer(as.character(pah_dia)),
    FU_pah_dia = as.numeric(FU_pah_dia)
  )

cox_c1_unimp <- mclapply(key_vars, function(var) {
  
  formula_str <- paste0(
    "Surv(FU_pah_dia, pah_dia) ~ ",
    var, " + ",
    paste(covs_C1, collapse = " + ")
  )
  
  model_data <- T_C1_C_match14_unimp_matched %>%
    dplyr::select(FU_pah_dia, pah_dia, all_of(var), all_of(covs_C1)) %>%
    tidyr::drop_na()
  
  cox_model <- coxph(
    as.formula(formula_str),
    data = model_data
  )
  
  cox_tidy <- broom::tidy(
    cox_model,
    exponentiate = TRUE,
    conf.int = TRUE
  )
  
  cox_tidy_var <- cox_tidy %>%
    dplyr::filter(term == var) %>%
    dplyr::mutate(
      variable_tested = var,
      n_used = nrow(model_data),
      n_event = sum(model_data$pah_dia == 1, na.rm = TRUE)
    )
  
  return(cox_tidy_var)
  
}, mc.cores = 5)

cox_c1_unimp_all <- bind_rows(cox_c1_unimp) %>%
  dplyr::transmute(
    protein = toupper(term),
    `OR (95%CI)` = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
    pval = ifelse(p.value < 0.001, "< 0.001", sprintf("%.3f", p.value)),
    n_used = n_used,
    n_event = n_event
  )

cox_c1_unimp_all

write.csv(cox_c1_unimp_all,  "./data/Discovery/sensitivity_unimp_results/cox_c1_unimp.csv", row.names = FALSE,  fileEncoding = "UTF-8")
#----C2----
T_C2_C_match11_unimp <- match_by_matchit(T_Cohort2_unimp_dead, T_Cohort2_unimp_surv,
                                         covars = match_vars,
                                         ratio = 1,
                                         method = "nearest",
                                         distance = "logit",
                                         caliper = NULL,
                                         seed = 20260305)

T_C2_C_match11_unimp_matched <- T_C2_C_match11_unimp$matched_data

covs_C2 <- c("walk_pace", "ntprobnp")

key_vars <- c("edn1", "anxa2", "spint1", "ltbp2", "thbs2", "col4a1", "qsox1", "c7", "slamf7", "igsf3", "lrrn1", "nptx1")

T_C2_C_match11_unimp_matched <- T_C2_C_match11_unimp_matched %>%
  dplyr::mutate(
    all_death = as.integer(as.character(all_death)),
    FU_death = as.numeric(FU_death)
  )

cox_c2_unimp <- mclapply(key_vars, function(var) {
  
  formula_str <- paste0(
    "Surv(FU_death, all_death) ~ ",
    var, " + ",
    paste(covs_C2, collapse = " + ")
  )
  
  model_data <- T_C2_C_match11_unimp_matched %>%
    dplyr::select(FU_death, all_death, all_of(var), all_of(covs_C2)) %>%
    tidyr::drop_na()
  
  cox_model <- coxph(
    as.formula(formula_str),
    data = model_data
  )
  
  cox_tidy <- broom::tidy(
    cox_model,
    exponentiate = TRUE,
    conf.int = TRUE
  )
  
  cox_tidy_var <- cox_tidy %>%
    dplyr::filter(term == var) %>%
    dplyr::mutate(
      variable_tested = var,
      n_used = nrow(model_data),
      n_event = sum(model_data$all_death == 1, na.rm = TRUE)
    )
  
  return(cox_tidy_var)
  
}, mc.cores = 5)

cox_c2_unimp_all <- bind_rows(cox_c2_unimp) %>%
  dplyr::transmute(
    protein = toupper(term),
    `OR (95%CI)` = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
    pval = ifelse(p.value < 0.001, "< 0.001", sprintf("%.3f", p.value)),
    n_used = n_used,
    n_event = n_event
  )

cox_c2_unimp_all

write.csv(cox_c2_unimp_all,  "./data/Discovery/sensitivity_unimp_results/cox_c2_unimp.csv", row.names = FALSE,  fileEncoding = "UTF-8")
