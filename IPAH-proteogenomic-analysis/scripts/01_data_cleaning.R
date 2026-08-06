library(ranger)
library(dplyr)
library(MASS)
library(nnet)
library(mice)
library(survival)
library(broom) 
library(lubridate)
library(data.table)
library(stringr)

#===============================================================================
#---------------------------------ICD-10----------------------------------------

T_original_icd10 <- read.csv('./data/Original/icd10.csv', header = TRUE)

date_cols <- grep("^X41280\\.", names(T_original_icd10), value = TRUE)
date_cols <- date_cols[order(as.integer(sub(".*\\.(\\d+)$", "\\1", date_cols)))]

dates_mat <- as.matrix(T_original_icd10[, date_cols])

get_first_date_prefix <- function(code_str, dates_row, prefixes) {
  if (is.na(code_str) || code_str == "") return(as.Date(NA))
  codes <- strsplit(code_str, "\\|")[[1]]
  if (length(codes) == 0) return(as.Date(NA))
  
  n <- min(length(codes), length(dates_row))
  codes <- codes[seq_len(n)]
  dts   <- as.Date(dates_row[seq_len(n)])
  
  hit <- rep(FALSE, n)
  for (p in prefixes) hit <- hit | startsWith(codes, p)
  
  if (!any(hit, na.rm = TRUE)) return(as.Date(NA))
  suppressWarnings(min(dts[hit], na.rm = TRUE))
}

#list
icd_map <- list(
  pah          = c("I270"),  #Primary PAH
  sec_pah      = c("I272"),  #Secondary PAH
  hypt         = c("I1"),    #Hypertension
  cad          = c("I25"),   #Coronary heart disease
  hf           = c("I50"),   #Heart failure
  af           = c("I48"),   #Atrial fibrillation
  lh           = c("I34", "I35", "I05", "I06", "I08", "Q23"), #left heart related disease (nonrheumatic, rheumatic, and congenital)  
  
  copd         = c("J44"),   #COPD
  ild          = c("J849"),  #Interstitial lung disease
  osa          = c("G473"),  #OSAS
  pe           = c("I26"),   #pulmonary embolism
  
  ckd          = c("N18"),   #CKD
  dm           = c("E1"),    #DM
  
  chd          = c("Q200", "Q201", "Q201", "Q202", "Q203", "Q204", "Q205", "Q21", "Q250", "Q262", "Q263"),   #Congenital heart disease
  ctd          = c("L900", "L940", "L93", "M05", "M06", "M08", "M120", "M32", "M33", "M34", "M350", "M351", "M358", "M359"),   #Connective tissue disease
  hiv          = c("B20", "B21", "B22", "B23", "B24", "B973", "Z21"),  #HIV infection
  pht          = c("K766"),   #Portal hypertension
  sch          = c("B65"),   #Schistosoma
  hht          = c("I780")   #Hereditary hemorrhagic telangiectasia

)

#----DO----
code_vec <- T_original_icd10$X41270.0.0
n <- nrow(T_original_icd10)

for (nm in names(icd_map)) {
  pats <- icd_map[[nm]]
  

  pattern_regex <- paste(pats, collapse = "|")
  T_original_icd10[[paste0(nm, "_dia")]] <- as.factor(
    ifelse(grepl(pattern_regex, code_vec, ignore.case = TRUE), 1, 0)
  )
  

  T_original_icd10[[paste0("Date_", nm)]] <- as.Date(unlist(Map(
    f = function(cs, i) get_first_date_prefix(cs, dates_mat[i, ], pats),
    cs = code_vec,
    i  = seq_len(n)
  )))
}

T_overall_icd10 <- T_original_icd10[, c(1, 262:299)]

write.csv(T_overall_icd10, "./data/data_clean/T_overall_icd10.csv", row.names = FALSE)

T_original_cov <- read.csv('./data/Original/data_required_20250814.csv', header = TRUE)
T_medicine_GP <- fread('./data/GP_data/gp_drug_script_cl.csv', quote="\"", fill=TRUE, showProgress=TRUE)
T_supp_PA_0260226 <- read.csv('./data/Original/data_supplied_PA.csv', header =  TRUE)
T_supp_OT_0260303 <- read.csv('./data/Original/data_supplied_operation_date.csv', header = TRUE)

#which(colnames(T_original_cov) == "p20003_i0_a0")
#which(colnames(T_original_cov) == "p20003_i3_a47")
T_medicine_hos <- T_original_cov
dt <- as.data.table(T_medicine_hos)
med_cols <- names(dt)[46:237]

length(med_cols)
head(med_cols)

drug_map <- list(
  ERA_hos        = c(1141186674),

  PDE5i_hos      = c(1141168936, 1141187810),

  CCB_hos        = c(1140879802, 1140879806, 1140861088)

)

for (drug in names(drug_map)) {
  codes <- drug_map[[drug]]
  dt[, (drug) := as.integer(Reduce(`|`, lapply(.SD, `%in%`, table = codes))), .SDcols = med_cols]
}

T_medicine_hos <- as.data.frame(dt)
T_medicine_hos <- T_medicine_hos[, c(1, 344, 345, 346)]

gp <- as.data.table(T_medicine_GP)
gp[, drug_name_low := tolower(drug_name)]

drug_dict <- list(
  ERA_gp   = c("bosentan", "ambrisentan", "macitentan"),
  PDE5i_gp = c("sildenafil", "tadalafil"),
  CCB_gp   = c("nifedipine", "amlodipine", "diltiazem", "verapamil", "felodipine",
             "lercanidipine", "nicardipine", "isradipine", "lacidipine"))

for (cls in names(drug_dict)) {
  keys <- drug_dict[[cls]]
  pat <- paste(keys, collapse = "|")

  gp[, (cls) := as.integer(grepl(pat, drug_name_low))]
}

T_medicine_gp <- gp[, .(
  ERA_gp   = as.integer(any(ERA_gp == 1, na.rm = TRUE)),
  PDE5i_gp = as.integer(any(PDE5i_gp == 1, na.rm = TRUE)),
  CCB_gp   = as.integer(any(CCB_gp == 1, na.rm = TRUE))
), by = eid]

T_medicine <- merge(T_medicine_hos, T_medicine_gp, by = 'eid', all.x = TRUE, all.y = FALSE) %>%
  dplyr::mutate(ERA = case_when(ERA_hos == 1 | ERA_gp == 1 ~ 1,
                                TRUE ~ 0),
                PDE5i = case_when(PDE5i_hos == 1 | PDE5i_gp == 1 ~ 1, 
                                  TRUE ~ 0),
                CCB = case_when(CCB_hos == 1 | CCB_gp == 1 ~ 1,
                                  TRUE ~ 0)) %>%
  dplyr::select(eid, ERA, PDE5i, CCB)
  

write.csv(T_medicine, "./data/data_clean/T_medicine.csv", row.names = FALSE)

T_surgery <- subset(T_original_cov, select = c(eid, p41272))
T_surgery <- merge(T_surgery, T_supp_OT_0260303, by = 'eid', all.x = TRUE, all.y = FALSE)

date_cols <- grep("^p41282_a\\d+$", names(T_surgery), value = TRUE)
date_cols <- date_cols[order(as.integer(sub("^p41282_a", "", date_cols)))]

dates_mat <- as.matrix(T_surgery[, date_cols])

code_vec <- T_surgery$p41272
n <- nrow(T_surgery)
tx_prefix <- c("E53", "K01")

T_surgery$transplantation <- as.factor(ifelse(grepl("E53|K01", code_vec), 1, 0))

T_surgery$Date_transplantation <- as.Date(unlist(Map(
  f = function(cs, i) get_first_date_prefix(cs, dates_mat[i, ], tx_prefix),
  cs = code_vec,
  i  = seq_len(n)
)))

T_surgery <- T_surgery %>%
  dplyr::select(eid, transplantation, Date_transplantation)

write.csv(T_surgery, "./data/data_clean/T_surgery.csv", row.names = FALSE)

#which(colnames(T_death) == "p40000_i0")
#which(colnames(T_death) == "p40002_i0_a14")
death_var <- c("p40001_i0","p40001_i1","p40002_i0_a0","p40002_i0_a1","p40002_i0_a2","p40002_i0_a3",
               "p40002_i0_a4","p40002_i0_a5","p40002_i0_a6","p40002_i0_a7",
               "p40002_i0_a8","p40002_i0_a9","p40002_i0_a10","p40002_i0_a11",
               "p40002_i0_a12","p40002_i0_a13","p40002_i0_a14",
               "p40002_i1_a0","p40002_i1_a1","p40002_i1_a2","p40002_i1_a3",
               "p40002_i1_a4","p40002_i1_a5","p40002_i1_a6","p40002_i1_a7",
               "p40002_i1_a8","p40002_i1_a9","p40002_i1_a10","p40002_i1_a11",
               "p40002_i1_a12","p40002_i1_a13","p40002_i0_a14")
daeth_cvd <- paste0("^(", "I0[0-2]|", "I0[5-9]|", "I1[0-5]|", "I2[0-5]|", "I2[6-8]|", 
                    "I3[0-9]|","I4[0-5]|","I46[0-8]","I4[7-9]|","I5[0-2]|", "I6[0-9]|",
                    "I7[0-7]" ,")")

T_death <- T_original_cov[, c(1, 304:338)]
T_death$all_death <- ifelse(!is.na(T_death$p40001_i0) & T_death$p40001_i0 != "", 1, 0)
T_death$pah_death <- apply(T_death[ , death_var], 1, function(row) {
  if (any(grepl("I270", row, ignore.case = TRUE))) 1 else 0
})
T_death$cvd_death <- apply(T_death[ , death_var], 1, function(row) {
  codes <- row[!is.na(row)]
  if (any(grepl(daeth_cvd, codes, ignore.case = TRUE))) {
    1
  } else {
    0
  }
})

T_death <- T_death %>%
  dplyr::rename(date_death = p40000_i0) %>%
  dplyr::select(eid, date_death, all_death, pah_death, cvd_death)

write.csv(T_death, "./data/data_clean/T_death.csv", row.names = FALSE)

keep_vars <- c(

  "walk_dow","walk_dow_2012","walk_dow_2014","walk_dow_2019",
  "walk_mod","walk_mod_2012","walk_mod_2014","walk_mod_2019",
  "MPA_dow","MPA_dow_2012","MPA_dow_2014","MPA_dow_2019",
  "MPA_mod","MPA_mod_2012","MPA_mod_2014","MPA_mod_2019",
  "VPA_dow","VPA_dow_2012","VPA_dow_2014","VPA_dow_2019",
  "VPA_mod","VPA_mod_2012","VPA_mod_2014","VPA_mod_2019",
  "walk_pace","walk_pace_2012","walk_pace_2014","walk_pace_2019",
  "walk_SOB","walk_SOB_2012","walk_SOB_2014","walk_SOB_2019",
  "MET_total","MET_total_2012","MET_total_2014","MET_total_2019",
  

  "FVC","FVC_2012","FVC_2014","FVC_2019",
  "FEV1","FEV1_2012","FEV1_2014","FEV1_2019",
  "PEF","PEF_2012","PEF_2014","PEF_2019",
  
  "MAX_workload","MAX_workload_2012",
  "Comp_workload","Comp_workload_2012",
  "VO2max_moedl","VO2max_moedl_2012",
  "VO2max","VO2max_2012",
  "ACCE_2013","ACCE_data_avail","ACCE_data_avail2"
)

T_PA <- T_supp_PA_0260226 %>% 
  dplyr::rename(

                walk_dow = p864_i0, walk_dow_2012 = p864_i1, walk_dow_2014 = p864_i2, walk_dow_2019 = p864_i3,
                walk_mod = p874_i0, walk_mod_2012 = p874_i1, walk_mod_2014 = p874_i2, walk_mod_2019 = p874_i3,
                MPA_dow = p884_i0,  MPA_dow_2012 = p884_i1,  MPA_dow_2014 = p884_i2,  MPA_dow_2019 = p884_i3,
                MPA_mod = p894_i0,  MPA_mod_2012 = p894_i1,  MPA_mod_2014 = p894_i2,  MPA_mod_2019 = p894_i3,
                VPA_dow = p904_i0,  VPA_dow_2012 = p904_i1,  VPA_dow_2014 = p904_i2,  VPA_dow_2019 = p904_i3,
                VPA_mod = p914_i0,  VPA_mod_2012 = p914_i1,  VPA_mod_2014 = p914_i2,  VPA_mod_2019 = p914_i3,
                walk_pace = p924_i0, walk_pace_2012 = p924_i1, walk_pace_2014 = p924_i2, walk_pace_2019 = p924_i3,

                walk_SOB = p4717_i0, walk_SOB_2012 = p4717_i1, walk_SOB_2014 = p4717_i2, walk_SOB_2019 = p4717_i3,

                MET_total = p22040_i0, MET_total_2012 = p22040_i1, MET_total_2014 = p22040_i2, MET_total_2019 = p22040_i3,

                FVC_a0      = p3062_i0_a0, FVC_a1      = p3062_i0_a1, FVC_a2      = p3062_i0_a2,
                FVC_2012_a0 = p3062_i1_a0, FVC_2012_a1 = p3062_i1_a1, FVC_2012_a2 = p3062_i1_a2,
                FVC_2014_a0 = p3062_i2_a0, FVC_2014_a1 = p3062_i2_a1, FVC_2014_a2 = p3062_i2_a2,
                FVC_2019_a0 = p3062_i3_a0, FVC_2019_a1 = p3062_i3_a1, FVC_2019_a2 = p3062_i3_a2,
                
                FEV1_a0      = p3063_i0_a0, FEV1_a1      = p3063_i0_a1, FEV1_a2      = p3063_i0_a2,
                FEV1_2012_a0 = p3063_i1_a0, FEV1_2012_a1 = p3063_i1_a1, FEV1_2012_a2 = p3063_i1_a2,
                FEV1_2014_a0 = p3063_i2_a0, FEV1_2014_a1 = p3063_i2_a1, FEV1_2014_a2 = p3063_i2_a2,
                FEV1_2019_a0 = p3063_i3_a0, FEV1_2019_a1 = p3063_i3_a1, FEV1_2019_a2 = p3063_i3_a2,
                
                PEF_a0      = p3064_i0_a0, PEF_a1      = p3064_i0_a1, PEF_a2      = p3064_i0_a2,
                PEF_2012_a0 = p3064_i1_a0, PEF_2012_a1 = p3064_i1_a1, PEF_2012_a2 = p3064_i1_a2,
                PEF_2014_a0 = p3064_i2_a0, PEF_2014_a1 = p3064_i2_a1, PEF_2014_a2 = p3064_i2_a2,
                PEF_2019_a0 = p3064_i3_a0, PEF_2019_a1 = p3064_i3_a1, PEF_2019_a2 = p3064_i3_a2,
                
                MAX_workload  = p6032_i0,  MAX_workload_2012  = p6032_i1,
                Comp_workload = p6020_i0,  Comp_workload_2012 = p6020_i1,
                VO2max_moedl  = p30037_i0, VO2max_moedl_2012  = p30037_i1, 
                VO2max        = p30038_i0, VO2max_2012        = p30038_i1, 
                ACCE_2013     = p90012,    ACCE_data_avail    = p90015, ACCE_data_avail2 = p90016) %>%
  
  dplyr::mutate(FVC = rowMeans(across(c(FVC_a0, FVC_a1, FVC_a2)), na.rm = TRUE),
                 FVC = ifelse(is.nan(FVC), NA_real_, FVC),
                 FVC_2012 = rowMeans(across(c(FVC_2012_a0, FVC_2012_a1, FVC_2012_a2)), na.rm = TRUE),
                 FVC_2012 = ifelse(is.nan(FVC_2012), NA_real_, FVC_2012),
                 FVC_2014 = rowMeans(across(c(FVC_2014_a0, FVC_2014_a1, FVC_2014_a2)), na.rm = TRUE),
                 FVC_2014 = ifelse(is.nan(FVC_2014), NA_real_, FVC_2014),
                 FVC_2019 = rowMeans(across(c(FVC_2019_a0, FVC_2019_a1, FVC_2019_a2)), na.rm = TRUE),
                 FVC_2019 = ifelse(is.nan(FVC_2019), NA_real_, FVC_2019),
                 
                 FEV1 = rowMeans(across(c(FEV1_a0, FEV1_a1, FEV1_a2)), na.rm = TRUE),
                 FEV1 = ifelse(is.nan(FEV1), NA_real_, FEV1),
                 FEV1_2012 = rowMeans(across(c(FEV1_2012_a0, FEV1_2012_a1, FEV1_2012_a2)), na.rm = TRUE),
                 FEV1_2012 = ifelse(is.nan(FEV1_2012), NA_real_, FEV1_2012),
                 FEV1_2014 = rowMeans(across(c(FEV1_2014_a0, FEV1_2014_a1, FEV1_2014_a2)), na.rm = TRUE),
                 FEV1_2014 = ifelse(is.nan(FEV1_2014), NA_real_, FEV1_2014),
                 FEV1_2019 = rowMeans(across(c(FEV1_2019_a0, FEV1_2019_a1, FEV1_2019_a2)), na.rm = TRUE),
                 FEV1_2019 = ifelse(is.nan(FEV1_2019), NA_real_, FEV1_2019),
                 
                 PEF = rowMeans(across(c(PEF_a0, PEF_a1, PEF_a2)), na.rm = TRUE),
                 PEF = ifelse(is.nan(PEF), NA_real_, PEF),
                 PEF_2012 = rowMeans(across(c(PEF_2012_a0, PEF_2012_a1, PEF_2012_a2)), na.rm = TRUE),
                 PEF_2012 = ifelse(is.nan(PEF_2012), NA_real_, PEF_2012),
                 PEF_2014 = rowMeans(across(c(PEF_2014_a0, PEF_2014_a1, PEF_2014_a2)), na.rm = TRUE),
                 PEF_2014 = ifelse(is.nan(PEF_2014), NA_real_, PEF_2014),
                 PEF_2019 = rowMeans(across(c(PEF_2019_a0, PEF_2019_a1, PEF_2019_a2)), na.rm = TRUE),
                 PEF_2019 = ifelse(is.nan(PEF_2019), NA_real_, PEF_2019)) %>%

  dplyr::select(eid, walk_pace, MET_total, FVC) %>%

  dplyr::mutate(walk_pace = case_when(walk_pace %in% c(-7, -3) ~ NA_real_,
                                      TRUE ~ walk_pace),
                walk_pace = factor(walk_pace))

write.csv(T_PA, "./data/data_clean/T_PA.csv", row.names = FALSE)

T_EKG_MRI <- subset(T_original_cov, select =c(eid, p12340_i2, p12340_i3, p12338_i2, p12338_i3, p22332_i2, p22332_i3,
                                                  p22424_i2, p22424_i3, p24103_i2, p24103_i3, p24100_i2, p24100_i3, p24101_i2, p24101_i3, p22420_i2, p22420_i3, p22421_i2, p22421_i3, p22422_i2, p22422_i3,
                                                  p24114_i2, p24114_i3, p24115_i2, p24115_i3, p24116_i2, p24116_i3, p24117_i2, p24117_i3, 
                                                  p31133_i2, p31068_i2, p31067_i2,p31069_i2, p31070_i2)) %>%
  rename(ECG_QRS_2014 = p12340_i2, ECG_QRS_2019 = p12340_i3, ECG_P_2014 = p12338_i2, ECG_P_2019 = p12338_i3, ECG_QTC_2014 = p22332_i2, ECG_QTC_2019 = p22332_i3,
         MRI_co_2014 = p22424_i2, MRI_co_2019 = p22424_i3, MRI_lvef2_2014 = p24103_i2, MRI_lvef2_2019 = p24103_i3, MRI_lvedv2_2014 = p24100_i2, MRI_lvedv2_2019 = p24100_i3, MRI_lvesv2_2014 = p24101_i2, MRI_lvesv2_2019 = p24101_i3,
         MRI_lvef_2014 = p22420_i2, MRI_lvef_2019 = p22420_i3, MRI_lvedv_2014 = p22421_i2, MRI_lvedv_2019 = p22421_i3, MRI_lvesv_2014 = p22422_i2, MRI_lvesv_2019 = p22422_i3,
         MRI_ramv_2014 = p24114_i2, MRI_ramv_2019 = p24114_i3, MRI_ranv_2014 = p24115_i2, MRI_ranv_2019 = p24115_i3, MRI_rasv_2014 = p24116_i2, MRI_rasv_2019 = p24116_i3, MRI_raef_2014 = p24117_i2, MRI_raef_2019 = p24117_i3,
         MRI_rvqs_2014 = p31133_i2, MRI_rvef_2014 = p31068_i2, MRI_rvedv_2014 = p31067_i2, MRI_rvesv_2014 = p31069_i2, MRI_rvsv_2014 = p31070_i2)

write.csv(T_EKG_MRI, "./data/data_clean/T_EKG_MRI.csv", row.names = FALSE)

T_original_prot <- readRDS('./data/Original/protein_imputed_missforest.rds')
T_protein <-  T_original_prot[, c(1, 4:2923)]

write.csv(T_protein, "./data/data_clean/T_protein.csv", row.names = FALSE)

T_original_cov <- read.csv('./data/original/data_required_20250814.csv', header = TRUE)
T_original_smoking <- read.csv('./data/original/data_supplied_smoking.csv', header = TRUE)

T_baseline_cov <- subset(T_original_cov, select = c(eid, p31, p200, p1269_i0, p1279_i0, p4079_i0_a0, p4079_i0_a1, p4080_i0_a0 ,p4080_i0_a1, 
                                                    p20117_i0, p20456, p20457, p21000_i0, p21001_i0, p21022, p30620_i0, p30700_i0, p30710_i0,
                                                    p30740_i0, p30750_i0, p30780_i0, p30870_i0, p30038_i0)) %>% 
  dplyr::rename( sex = p31, date_join = p200, DBP1 = p4079_i0_a0, DBP2 = p4079_i0_a1,
          SBP1 = p4080_i0_a0, SBP2 = p4080_i0_a1, alcohol = p20117_i0, drug_his = p20456, drug_now = p20457, ethnicity = p21000_i0,
          BMI = p21001_i0, age = p21022, ALT = p30620_i0, Creatinine = p30700_i0, CRP = p30710_i0, Glucose = p30740_i0, HbA1c = p30750_i0,
          LDL = p30780_i0, Triglycerides = p30870_i0)
T_baseline_cov$sex <- as.factor(T_baseline_cov$sex)
T_baseline_cov$DBP_a <- (T_baseline_cov$DBP1 + T_baseline_cov$DBP2) / 2
T_baseline_cov$SBP_a <- (T_baseline_cov$SBP1 + T_baseline_cov$SBP2) / 2
T_baseline_cov$alcohol01 <- ifelse(T_baseline_cov$alcohol >= 1, 1, 0)
T_baseline_cov$alcohol01 <- as.factor(T_baseline_cov$alcohol01)
T_baseline_cov$ethnicity <- as.factor(T_baseline_cov$ethnicity)
T_baseline_cov$ethnicity_C <- ifelse(T_baseline_cov$ethnicity %in% c(1, 1001, 1002, 1003), 'White',
                              ifelse(T_baseline_cov$ethnicity %in% c(3, 5, 3001, 3002, 3003, 3004), 'Yellow',
                              ifelse(T_baseline_cov$ethnicity %in% c(4, 4001, 4002, 4003), 'Black',
                              ifelse(T_baseline_cov$ethnicity %in% c(2, 6, 2001, 2002, 2003, 2004), 'Others', NA))))
T_baseline_cov$ethnicity_C <- as.factor(T_baseline_cov$ethnicity_C)
T_baseline_cov <- subset(merge(T_protein, T_baseline_cov, by = "eid", all.x = TRUE, all.y = FALSE))
T_baseline_cov <- subset(T_baseline_cov, select = c(eid, sex, date_join, DBP_a, SBP_a, alcohol01, ethnicity_C,
                                                    BMI, age, ALT, Creatinine, CRP, Glucose, HbA1c,
                                                    LDL, Triglycerides))

T_baseline_smoking <- T_original_smoking %>%
  dplyr::rename(smoking = p20116_i0) %>%
  dplyr::mutate(smoking01 = as.factor(case_when(
                                     smoking %in% c("Previous", "Current") ~ 1,
                                     smoking == "Never" ~ 0,
                                     smoking == "Prefer not to answer" ~ NA_real_,
                                     TRUE ~ NA_real_))) %>%
  dplyr::select(eid, smoking01)

T_baseline_cov_v2 <- merge(T_baseline_cov, T_baseline_smoking, by = 'eid', all.x = TRUE, all.y = FALSE) 
T_baseline_cov_v3 <- merge(T_baseline_cov_v2, T_PA, by = 'eid', all.x = TRUE, all.y = FALSE) %>%
  dplyr::filter(ethnicity_C == 'White')

  
summary(T_baseline_cov_v2$ethnicity_C)

#write.csv(T_baseline_cov_v3, "./data/T_baseline_cov.csv", row.names = FALSE)

meth <- make.method(T_baseline_cov_v3)
meth[] <- "rf"
meth[c("eid", "date_join")] <- ""

pred <- make.predictorMatrix(T_baseline_cov_v3)
pred[, c("eid", "date_join")] <- 0

pred[c("eid", "date_join"), ] <- 0  

imp <- mice(T_baseline_cov_v3, m = 5, method = meth, predictorMatrix = pred, printFlag=FALSE, seed=163)

T_baseline_cov_Ied3 <- complete(imp, action =3)

write.csv(T_baseline_cov_Ied3, "./data/data_clean/T_baseline_cov_Ied3.csv", row.names = FALSE)
#write.csv(T_baseline_cov_Ied3, "./data/data_clean/T_baseline_cov_Ied3_Black.csv", row.names = FALSE)
#write.csv(T_baseline_cov_Ied3, "./data/data_clean/T_baseline_cov_Ied3_Yellow.csv", row.names = FALSE)
#write.csv(T_baseline_cov_Ied3, "./data/data_clean/T_baseline_cov_Ied3_Others.csv", row.names = FALSE)

T_baseline_cov_Ied3 <- read.csv("./data/data_clean/T_baseline_cov_Ied3.csv", header = TRUE)
#T_baseline_cov_Ied3 <- read.csv("./data/data_clean/T_baseline_cov_Ied3_Black.csv", header = TRUE)
#T_baseline_cov_Ied3 <- read.csv("./data/data_clean/T_baseline_cov_Ied3_Yellow.csv", header = TRUE)
#T_baseline_cov_Ied3 <- read.csv("./data/data_clean/T_baseline_cov_Ied3_Others.csv", header = TRUE)

T_overall_icd10 <- read.csv("./data/data_clean/T_overall_icd10.csv", header = TRUE)
T_medicine <- read.csv("./data/data_clean/T_medicine.csv", header = TRUE)
T_surgery <- read.csv("./data/data_clean/T_surgery.csv", header = TRUE)
T_death <- read.csv("./data/data_clean/T_death.csv", header = TRUE)
T_EKG_MRI <- read.csv("./data/data_clean/T_EKG_MRI.csv", header = TRUE)
T_protein <- read.csv("./data/data_clean/T_protein.csv", header = TRUE)

T_overall <- merge(T_baseline_cov_Ied3, T_overall_icd10, by = 'eid', all.x = TRUE, all.y = FALSE)
T_overall <- merge(T_overall, T_medicine, by = 'eid', all.x = TRUE, all.y = FALSE)
T_overall <- merge(T_overall, T_surgery, by = 'eid', all.x = TRUE, all.y = FALSE)
T_overall <- merge(T_overall, T_death, by = 'eid', all.x = TRUE, all.y = FALSE)
T_overall <- merge(T_overall, T_EKG_MRI, by = 'eid', all.x = TRUE, all.y = FALSE)
T_overall <- merge(T_overall, T_protein, by = 'eid', all.x = TRUE, all.y = FALSE)

write.csv(T_overall, "./data/data_clean/T_PAH_HELTH_Cross_Protein.csv", row.names = FALSE)
#write.csv(T_overall, "./data/data_clean/T_PAH_HELTH_Cross_Protein_Black.csv", row.names = FALSE)
#write.csv(T_overall, "./data/data_clean/T_PAH_HELTH_Cross_Protein_Yellow.csv", row.names = FALSE)
#write.csv(T_overall, "./data/data_clean/T_PAH_HELTH_Cross_Protein_Others.csv", row.names = FALSE)

table(T_overall$pah_dia)