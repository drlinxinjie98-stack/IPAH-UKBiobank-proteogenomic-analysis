library(nestedcv)
library(caret)
library(xgboost)
library(dplyr)
library(ggplot2)
library(glmnet)
library(SHAPforxgboost)
library(shapviz)
library(pROC)
library(tidyr)
library(networkD3)
library(ggalluvial)
library(ggsankeyfier)
library(RColorBrewer)
library(viridis)
library(scales)
library(readxl)

C1_premodel <- read.csv("./data/Cohorts/T_Cohort1.csv", header = TRUE)
C2_premodel <- read.csv("./data/Cohorts/T_Cohort2.csv", header = TRUE)
C_control_premodel <- read.csv("./data/Cohorts/T_Cohort_control.csv", header = TRUE)

T_cox_MR_cis_cross_C1 <-read.csv("./data/Validation/c1_cox_MR_cis_cross.csv", header = TRUE)
T_cox_MR_cis_cross_C2 <-read.csv("./data/Validation/c2_cox_MR_cis_cross.csv", header = TRUE)

sig_prot <- dplyr::bind_rows(T_cox_MR_cis_cross_C1, T_cox_MR_cis_cross_C2)

prot_del_acc <- read_excel("./data/cross_ancestry/logit_nonWhite_all_unmatched.xlsx") %>%
  dplyr::filter(p.value < 0.5)
prot_del_acc_list <- unique(tolower(as.character(prot_del_acc$term)))

prot_del_sen <- read_excel("./data/sensitivity/logit_sensitive_all_matched.xlsx") %>%
  filter(term != "rspo1", term != "nos1", term != "shisa5", term != "ccn3")
prot_del_sen_list <- unique(tolower(as.character(prot_del_sen$term)))

#ancestry cross sensitive 
prot_pass_val <- intersect(prot_del_acc_list, prot_del_sen_list)

#cox+MR pass proteins across validated proteins
sig_pro_list <- intersect(unique(tolower(as.character(sig_prot$protein_name))), prot_pass_val)

C1_C_control <- bind_rows(C1_premodel, C_control_premodel)

C_incidence_preM <- C1_C_control %>%
  dplyr::rename_with(~ toupper(.x), .cols = dplyr::any_of(sig_pro_list)) %>%
  dplyr::select(age, sex, BMI, SBP_a, smoking01, alcohol01, ALT, Creatinine, CRP, HbA1c, LDL, Triglycerides, nppb, ntprobnp, 
                walk_pace, MET_total, FVC, 
                pah_dia, Date_pah, FU_pah_dia,
                dplyr::any_of(toupper(sig_pro_list)))

C_incidence_preM_male <- subset(C_incidence_preM, sex == 1)
C_incidence_preM_female <- subset(C_incidence_preM, sex == 0)

C1_C2_control <- bind_rows(C1_premodel, C2_premodel)

C_mortality_preM <- C1_C2_control %>%
  dplyr::rename_with(~ toupper(.x), .cols = dplyr::any_of(sig_pro_list)) %>%
  dplyr::select(age, sex, BMI, SBP_a, smoking01, alcohol01, ALT, Creatinine, CRP, HbA1c, LDL, Triglycerides, nppb, ntprobnp, 
                walk_pace, MET_total, FVC, 
                all_death, date_death, FU_death,
                dplyr::any_of(toupper(sig_pro_list)))

C_mortality_preM_male <- subset(C_mortality_preM, sex == 1)
C_mortality_preM_female<- subset(C_mortality_preM, sex == 0)

table(C1_C2_control$all_death)
sum(C1_C2_control$date_death < C1_C2_control$Date_pah, na.rm = TRUE)

C_preclin_cluster <- read.csv("./data/clustering/Cohort_preclinical_clustered.csv", header = TRUE)
C_clinic_cluster <-  read.csv("./data/clustering/Cohort_clinical_clustered.csv", header = TRUE) 
C_cluster <- bind_rows(C_preclin_cluster, C_clinic_cluster) %>%
  dplyr::rename_with(~ toupper(.x), .cols = dplyr::any_of(sig_pro_list))

sig_pro_list <- toupper(sig_pro_list)

C_mortality_preM_C <- merge(C_mortality_preM, C_cluster, by = c(sig_pro_list, "all_death", "sex"), all.x = TRUE, all.y = FALSE)

C_mortality_preM_DVR <- subset(C_mortality_preM_C, cluster == 1)
C_mortality_preM_RVR <- subset(C_mortality_preM_C, cluster == 2)

run_nested_xgb_batch <- function(data,
                                 feature_sets,
                                 targets,
                                 tune_grid,
                                 outdir = NULL,
                                 seed = 123,
                                 n_outer_folds = 5,
                                 n_inner_folds = 5,
                                 cv_cores = 1,
                                 balance_method = "smote",
                                 metric = "logLoss",
                                 save_files = TRUE) {
  set.seed(seed)
  

  if (!is.null(outdir) && save_files) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  }
  

  lasso_filter <- function(y, x, alpha = 1, nlambda = 100) {
    x_mat <- as.matrix(x)
    y_vec <- ifelse(y == levels(y)[2], 1, 0)
    
    cvfit <- glmnet::cv.glmnet(
      x = x_mat,
      y = y_vec,
      alpha = alpha,
      nlambda = nlambda,
      family = "binomial"
    )
    
    coefs <- coef(cvfit, s = "lambda.min")
    keep <- which(coefs[-1, 1] != 0)
    
    if (length(keep) == 0) keep <- seq_len(ncol(x_mat))
    
    return(keep)
  }
  

  results <- data.frame()
  delong_results <- data.frame()
  data_out <- data
  fits_list <- list()
  shap_rank_list <- list()
  

  pred_store <- list()
  

  roc_data_list <- list()
  

  for (set_name in names(feature_sets)) {
    current_features <- feature_sets[[set_name]]
    
    for (target_name in names(targets)) {
      message("Running: ", set_name, " -> predict ", target_name)
      

      X <- data[, current_features, drop = FALSE]
      y_raw <- targets[[target_name]]
      

      complete_idx <- complete.cases(X) & !is.na(y_raw)
      X_sub <- X[complete_idx, , drop = FALSE]
      y_sub <- y_raw[complete_idx]
      

      y <- factor(y_sub)
      if (length(levels(y)) != 2) {
        stop("Target ", target_name, " is not binary.")
      }
      levels(y) <- c("X0", "X1")
      

      fit <- nestedcv::nestcv.train(
        y = y,
        x = X_sub,
        method = "xgbTree",
        filterFUN = lasso_filter,
        filter_options = list(alpha = 1, nlambda = 200),
        n_outer_folds = n_outer_folds,
        n_inner_folds = n_inner_folds,
        tuneGrid = tune_grid,
        metric = metric,
        balance = balance_method,
        cv.cores = cv_cores,
        savePredictions = "final"
      )
      
      key_name <- paste(target_name, set_name, sep = "__")
      

      fits_list[[key_name]] <- fit
      

      pred <- fit$output
      pred$id_sub <- unlist(fit$outer_folds)
      pred <- pred[order(pred$id_sub), ]
      

      pred_full_class <- rep(NA, nrow(data))
      pred_full_prob  <- rep(NA, nrow(data))
      pred_full_class[which(complete_idx)[pred$id_sub]] <- as.character(pred$predy)
      pred_full_prob[which(complete_idx)[pred$id_sub]]  <- pred$predyp
      
      data_out[[paste0("pred_", target_name, "_", set_name)]] <- pred_full_class
      data_out[[paste0("prob_", target_name, "_", set_name)]] <- pred_full_prob
      
      # ---- AUC + 95% CI ----
      y_true_num <- ifelse(pred$testy == "X1", 1, 0)
      
      roc_obj <- pROC::roc(
        response = y_true_num,
        predictor = pred$predyp,
        ci = TRUE,
        quiet = TRUE
      )
      
      auc_value <- as.numeric(pROC::auc(roc_obj))
      auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
      
      perf <- data.frame(
        set = set_name,
        target = target_name,
        n_features = length(current_features),
        n_complete = sum(complete_idx),
        AUC = auc_value,
        AUC_CI_low = auc_ci[1],
        AUC_CI_mid = auc_ci[2],
        AUC_CI_high = auc_ci[3]
      )
      results <- rbind(results, perf)
      

      pred_store[[key_name]] <- data.frame(
        row_id = which(complete_idx)[pred$id_sub],
        y_true = y_true_num,
        prob = pred$predyp,
        stringsAsFactors = FALSE
      )
      

      roc_data_list[[key_name]] <- data.frame(
        row_id = which(complete_idx)[pred$id_sub],
        target = target_name,
        set = set_name,
        y_true = y_true_num,
        y_true_label = as.character(pred$testy),
        pred_class = as.character(pred$predy),
        pred_prob = pred$predyp,
        stringsAsFactors = FALSE
      )
      
      # ---- SHAP ----
      xgb_model <- fit$final_fit$finalModel
      used_features <- xgb_model$feature_names
      X_mat <- as.matrix(X_sub[, used_features, drop = FALSE])
      
      shap_full <- SHAPforxgboost::shap.values(
        xgb_model = xgb_model,
        X_train = X_mat
      )
      
      shap_rank <- data.frame(
        variable = names(shap_full$mean_shap_score),
        mean_abs_shap = abs(shap_full$mean_shap_score),
        row.names = NULL
      )
      shap_rank <- shap_rank[order(-shap_rank$mean_abs_shap), ]
      
      shap_rank_list[[key_name]] <- shap_rank
      

      if (!is.null(outdir) && save_files) {
        utils::write.csv(
          shap_rank,
          file = file.path(outdir, paste0(target_name, "_SHAP_", set_name, ".csv")),
          row.names = FALSE
        )
      }
      

      top_n <- min(20, nrow(shap_rank))
      top_features <- shap_rank$variable[1:top_n]
      
      shap_long <- SHAPforxgboost::shap.prep(
        xgb_model = xgb_model,
        X_train = X_mat,
        top_n = top_n
      )
      shap_long$variable <- factor(shap_long$variable, levels = top_features)
      
      p <- SHAPforxgboost::shap.plot.summary(shap_long)
      
      if (!is.null(outdir) && save_files) {
        ggplot2::ggsave(
          filename = file.path(outdir, paste0(target_name, "_", set_name, "_shap_sum.png")),
          plot = p,
          width = 6,
          height = 2,
          bg = "white"
        )
      }
    }
  }
  

  for (target_name in names(targets)) {
    target_keys <- names(pred_store)[grepl(paste0("^", target_name, "__"), names(pred_store))]
    
    if (length(target_keys) >= 2) {
      combs <- utils::combn(target_keys, 2, simplify = FALSE)
      
      for (cc in combs) {
        k1 <- cc[1]
        k2 <- cc[2]
        
        d1 <- pred_store[[k1]]
        d2 <- pred_store[[k2]]
        
        dd <- merge(d1, d2, by = c("row_id", "y_true"), suffixes = c("_1", "_2"))
        
        if (nrow(dd) > 1 && length(unique(dd$y_true)) == 2) {
          roc1 <- pROC::roc(dd$y_true, dd$prob_1, quiet = TRUE)
          roc2 <- pROC::roc(dd$y_true, dd$prob_2, quiet = TRUE)
          
          delong_test <- pROC::roc.test(roc1, roc2, method = "delong", paired = TRUE)
          
          set1 <- sub(paste0("^", target_name, "__"), "", k1)
          set2 <- sub(paste0("^", target_name, "__"), "", k2)
          
          delong_row <- data.frame(
            target = target_name,
            model_1 = set1,
            model_2 = set2,
            auc_1 = as.numeric(pROC::auc(roc1)),
            auc_2 = as.numeric(pROC::auc(roc2)),
            auc_diff = as.numeric(pROC::auc(roc1)) - as.numeric(pROC::auc(roc2)),
            p_value = delong_test$p.value,
            method = "paired DeLong test",
            n_shared = nrow(dd)
          )
          
          delong_results <- rbind(delong_results, delong_row)
        }
      }
    }
  }
  

  roc_data_all <- do.call(rbind, roc_data_list)
  rownames(roc_data_all) <- NULL
  

  if (!is.null(outdir) && save_files) {
    utils::write.csv(
      results,
      file = file.path(outdir, "AUC_results.csv"),
      row.names = FALSE
    )
    
    utils::write.csv(
      delong_results,
      file = file.path(outdir, "DeLong_results.csv"),
      row.names = FALSE
    )
    
    utils::write.csv(
      roc_data_all,
      file = file.path(outdir, "ROC_curve_data.csv"),
      row.names = FALSE
    )
  }
  
  return(list(
    results = results,
    delong_results = delong_results,
    roc_data = roc_data_all,
    data_with_preds = data_out,
    fits = fits_list,
    shap_ranks = shap_rank_list
  ))
}

tg <- expand.grid(
  nrounds = 200,
  max_depth = c(3, 5),
  eta = c(0.05, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

#---------------------------------Cmortality_bothsexes------------------------------------

feature_sets <- list(
  proteins = names(C_mortality_preM)[21:32],
  clinics = names(C_mortality_preM)[1:17],
  protein_clinics = c(names(C_mortality_preM)[21:32], names(C_mortality_preM)[1:17])
)

targets <- list(
  outcome = C_mortality_preM$all_death
)

res_C2_bothsexes <- run_nested_xgb_batch(
  data = C_mortality_preM,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_mortality/Bothsexes",
  seed = 321,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 3,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cmortality_male---------------------------------------

feature_sets <- list(
  proteins = names(C_mortality_preM_male)[21:32],
  clinics = names(C_mortality_preM_male)[1:17],
  protein_clinics = c(names(C_mortality_preM_male)[21:32], names(C_mortality_preM_male)[1:17])
)

targets <- list(
  outcome = C_mortality_preM_male$all_death
)

res_C2_male <- run_nested_xgb_batch(
  data = C_mortality_preM_male,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_mortality/Male",
  seed = 321,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 1,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cmortality_female-------------------------------------

feature_sets <- list(
  proteins = names(C_mortality_preM_female)[21:32],
  clinics = names(C_mortality_preM_female)[1:17],
  protein_clinics = c(names(C_mortality_preM_female)[21:32], names(C_mortality_preM_female)[1:17])
)

targets <- list(
  outcome = C_mortality_preM_female$all_death
)

res_C2_female <- run_nested_xgb_batch(
  data = C_mortality_preM_female,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_mortality/Female",
  seed = 321,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 1,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cincidence_bothsexes--------------------------

feature_sets <- list(
  proteins = names(C_incidence_preM)[21:32],
  clinics = names(C_incidence_preM)[1:17],
  protein_clinics = c(names(C_incidence_preM)[21:32], names(C_incidence_preM)[1:17])
)

targets <- list(
  outcome = C_incidence_preM$pah_dia
)

res_C1_bothsexes <- run_nested_xgb_batch(
  data = C_incidence_preM,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_incidence/Bothsexes",
  seed = 123,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 3,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cincidence_male---------------------------------------

feature_sets <- list(
  proteins = names(C_incidence_preM_male)[21:32],
  clinics = names(C_incidence_preM_male)[1:17],
  protein_clinics = c(names(C_incidence_preM_male)[21:32], names(C_incidence_preM_male)[1:17])
)

targets <- list(
  outcome = C_incidence_preM_male$pah_dia
)

res_C1_bothsexes <- run_nested_xgb_batch(
  data = C_incidence_preM_male,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_incidence/Male",
  seed = 123,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 3,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cincidence_female------------------------------

feature_sets <- list(
  proteins = names(C_incidence_preM_female)[21:32],
  clinics = names(C_incidence_preM_female)[1:17],
  protein_clinics = c(names(C_incidence_preM_female)[21:32], names(C_incidence_preM_female)[1:17])
)

targets <- list(
  outcome = C_incidence_preM_female$pah_dia
)

res_C1_bothsexes <- run_nested_xgb_batch(
  data = C_incidence_preM_female,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_incidence/Female",
  seed = 123,
  n_outer_folds = 5,
  n_inner_folds = 5,
  cv_cores = 3,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

#---------------------------------Cincidence_cluster----------------------------
#----DVR----

which(colnames(C_mortality_preM_DVR) == "sex")
which(colnames(C_mortality_preM_DVR) == "FVC")

feature_sets <- list(
  proteins = names(C_mortality_preM_RVR)[1:12],
  clinics = names(C_mortality_preM_RVR)[14:30],
  protein_clinics = c(names(C_mortality_preM_RVR)[1:12], names(C_mortality_preM_RVR)[14:30])
)

targets <- list(
  outcome = C_mortality_preM_RVR$all_death
)

res_C2_bothsexes <- run_nested_xgb_batch(
  data = C_mortality_preM_RVR,
  feature_sets = feature_sets,
  targets = targets,
  tune_grid = tg,
  outdir = "./Results/05PredictiveM/C_RVR/Bothsexes",
  seed = 321,
  n_outer_folds = 3,
  n_inner_folds = 3,
  cv_cores = 3,
  balance_method = "smote",
  metric = "logLoss",
  save_files = TRUE
)

plot_ROC <- function(data_path,
                                    target_name,
                                    clinic_set = "clinic_sets",
                                    protein_set = "protein_sets",
                                    combined_set = "procli_sets",
                                    title_text,
                                    save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(pROC)
  

  df <- read.csv(data_path, header = TRUE)
  

  df_sub <- df %>%
    dplyr::filter(target == target_name) %>%
    dplyr::filter(set %in% c(clinic_set, protein_set, combined_set))
  

  df_clinic  <- df_sub %>% dplyr::filter(set == clinic_set)
  df_protein <- df_sub %>% dplyr::filter(set == protein_set)
  df_procli  <- df_sub %>% dplyr::filter(set == combined_set)
  

  roc_clinic  <- pROC::roc(df_clinic$y_true,  df_clinic$pred_prob, quiet = TRUE)
  roc_protein <- pROC::roc(df_protein$y_true, df_protein$pred_prob, quiet = TRUE)
  roc_procli  <- pROC::roc(df_procli$y_true,  df_procli$pred_prob, quiet = TRUE)
  

  test_protein_vs_clinic   <- pROC::roc.test(roc_protein, roc_clinic, method = "delong", paired = TRUE)
  test_combined_vs_clinic  <- pROC::roc.test(roc_procli, roc_clinic, method = "delong", paired = TRUE)
  test_combined_vs_protein <- pROC::roc.test(roc_procli, roc_protein, method = "delong", paired = TRUE)
  
  auc_results <- data.frame(
    Comparison = c("Proteomic vs Clinical", "Combined vs Clinical", "Combined vs Proteomic"),
    p_value = c(test_protein_vs_clinic$p.value,
                test_combined_vs_clinic$p.value,
                test_combined_vs_protein$p.value)
  ) %>%
    dplyr::mutate(
      p_label = ifelse(
        p_value < 0.001,
        "p < 0.001",
        paste0("p = ", formatC(p_value, digits = 3, format = "f"))
      )
    )
  
  print(auc_results)
  

  roc_df <- dplyr::bind_rows(
    data.frame(
      fpr = 1 - roc_clinic$specificities,
      tpr = roc_clinic$sensitivities,
      model = paste0("Clinical (AUC=", round(pROC::auc(roc_clinic), 3), ")")
    ),
    data.frame(
      fpr = 1 - roc_protein$specificities,
      tpr = roc_protein$sensitivities,
      model = paste0("Proteomic (AUC=", round(pROC::auc(roc_protein), 3), ")")
    ),
    data.frame(
      fpr = 1 - roc_procli$specificities,
      tpr = roc_procli$sensitivities,
      model = paste0("Combined (AUC=", round(pROC::auc(roc_procli), 3), ")")
    )
  )
  

  roc_df_clean <- roc_df %>%
    dplyr::group_by(model, fpr) %>%
    dplyr::summarise(tpr = max(tpr), .groups = "drop")
  

  model_levels <- c(
    paste0("Clinical (AUC=", round(pROC::auc(roc_clinic), 3), ")"),
    paste0("Proteomic (AUC=", round(pROC::auc(roc_protein), 3), ")"),
    paste0("Combined (AUC=", round(pROC::auc(roc_procli), 3), ")")
  )
  
  roc_df_clean$model <- factor(roc_df_clean$model, levels = model_levels)
  model_colors <- c("#D0372C", "#4C7BB3", "#8D539F")
  

  p <- ggplot(roc_df_clean, aes(x = fpr, y = tpr, color = model)) +
    geom_step(linewidth = 1.2, direction = "vh") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray20") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(values = setNames(model_colors, levels(roc_df_clean$model))) +
    labs(title = title_text, x = "1 - Specificity", y = "Sensitivity", color = "Model") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      legend.position = c(0.7, 0.2),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18, color = "black")
    )
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, width = 7, height = 7, plot = p, bg = "white")
  }
  
  return(list(
    plot = p,
    auc_results = auc_results,
    roc_objects = list(
      clinic = roc_clinic,
      protein = roc_protein,
      combined = roc_procli
    )
  ))
}

#----Cmortality_bothsexes----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_mortality/Bothsexes/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Both sexes",
  save_path = "./Results/05PredictiveM/C_mortality/Bothsexes/ROC_curve.pdf"
)

#----Cmortality_male----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_mortality/Male/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Male",
  save_path = "./Results/05PredictiveM/C_mortality/Male/ROC_curve.pdf"
)
#----Cmortality_female----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_mortality/Female/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Female",
  save_path = "./Results/05PredictiveM/C_mortality/Female/ROC_curve.pdf"
)

#----Cincidence_bothsexes----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_incidence/Bothsexes/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Both sexes",
  save_path = "./Results/05PredictiveM/C_incidence/Bothsexes/ROC_curve.pdf"
)

#----Cincidence_male----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_incidence/Male/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Male",
  save_path = "./Results/05PredictiveM/C_incidence/Male/ROC_curve.pdf"
)

#----Cincidence_female----
res_plot <- plot_ROC(
  data_path = "./Results/05PredictiveM/C_incidence/Female/ROC_curve_data.csv",
  target_name = "outcome",
  clinic_set = "clinics",
  protein_set = "proteins",
  combined_set = "protein_clinics",
  title_text = "Female",
  save_path = "./Results/05PredictiveM/C_incidence/Female/ROC_curve.pdf"
)

#---------------------------------Dumbbell plot---------------------------------

plot_shap_dumbbell <- function(dat,
                                sort_by = c("mean", "diff_abs", "male", "female", "input"),
                                male_col = "relative_shap_male",
                                female_col = "relative_shap_female",
                                var_col = "variable",
                                diff_col = "diff_abs",
                                male_label = "Male",
                                female_label = "Female",
                                male_color = "#4C78A8",
                                female_color = "#E76F51",
                                line_color = "black",
                                point_size = 3,
                                line_width = 1,
                                x_lab = "Relative SHAP contribution") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  
  sort_by <- match.arg(sort_by)
  
  plot_dat <- dat %>%
    dplyr::select(all_of(c(var_col, male_col, female_col, diff_col))) %>%
    dplyr::rename(
      variable = all_of(var_col),
      relative_shap_male = all_of(male_col),
      relative_shap_female = all_of(female_col),
      diff_abs = all_of(diff_col)
    ) %>%
    dplyr::mutate(
      avg_shap = (relative_shap_male + relative_shap_female) / 2
    )
  
  if (sort_by == "mean") {
    plot_dat <- plot_dat %>% dplyr::arrange(desc(avg_shap))
  } else if (sort_by == "diff_abs") {
    plot_dat <- plot_dat %>% dplyr::arrange(desc(diff_abs))
  } else if (sort_by == "male") {
    plot_dat <- plot_dat %>% dplyr::arrange(desc(relative_shap_male))
  } else if (sort_by == "female") {
    plot_dat <- plot_dat %>% dplyr::arrange(desc(relative_shap_female))
  }
  
  plot_dat <- plot_dat %>%
    dplyr::mutate(variable = factor(variable, levels = rev(variable)))
  
  point_dat <- plot_dat %>%
    tidyr::pivot_longer(
      cols = c(relative_shap_male, relative_shap_female),
      names_to = "sex",
      values_to = "relative_shap"
    ) %>%
    dplyr::mutate(
      sex = factor(
        sex,
        levels = c("relative_shap_male", "relative_shap_female"),
        labels = c(male_label, female_label)
      )
    )
  
  p <- ggplot() +
    geom_segment(
      data = plot_dat,
      aes(
        x = relative_shap_male,
        xend = relative_shap_female,
        y = variable,
        yend = variable
      ),
      color = line_color,
      linewidth = line_width
    ) +
    geom_point(
      data = point_dat,
      aes(x = relative_shap, y = variable, color = sex),
      size = point_size
    ) +
    scale_color_manual(
      values = setNames(
        c(male_color, female_color),
        c(male_label, female_label)
      )
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.08))
    ) +
    labs(
      x = x_lab,
      y = NULL,
      color = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = point_size))) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 12, face = "plain", color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 13),
      plot.margin = margin(10, 20, 10, 10)
    )
  
  return(p)
}
#----Cmortality----
SHAP_mortality_male   <- read.csv("./Results/05PredictiveM/C_mortality/Male/outcome_SHAP_proteins.csv", header = TRUE) %>%
  dplyr::mutate(relative_shap_male = mean_abs_shap / sum(mean_abs_shap))
SHAP_mortality_female <- read.csv("./Results/05PredictiveM/C_mortality/Female/outcome_SHAP_proteins.csv", header = TRUE) %>%
  dplyr::mutate(relative_shap_female = mean_abs_shap / sum(mean_abs_shap))

SHAP_mortality <- merge(SHAP_mortality_male, SHAP_mortality_female, by = "variable", all.x = TRUE, all.y = TRUE) %>%
  dplyr::mutate(relative_shap_male   = ifelse(is.na(relative_shap_male), 0, relative_shap_male),
                relative_shap_female = ifelse(is.na(relative_shap_female), 0, relative_shap_female)) %>%
  dplyr::mutate(diff_abs = relative_shap_male - relative_shap_female) %>%
  dplyr::select(variable, relative_shap_male, relative_shap_female, diff_abs)

plot_Db_mortality <- plot_shap_dumbbell(dat = SHAP_mortality, sort_by = "mean")
ggsave(filename = "./Results/05PredictiveM/C_mortality/SHAP_dumbbell.pdf", 
       plot = plot_Db_mortality, width = 3,height = 5)

#----Cincidence----
SHAP_incidence_male   <- read.csv("./Results/05PredictiveM/C_incidence/Male/outcome_SHAP_proteins.csv", header = TRUE) %>%
  dplyr::mutate(relative_shap_male = mean_abs_shap / sum(mean_abs_shap))
SHAP_incidence_female <- read.csv("./Results/05PredictiveM/C_incidence/Female/outcome_SHAP_proteins.csv", header = TRUE) %>%
  dplyr::mutate(relative_shap_female = mean_abs_shap / sum(mean_abs_shap))

SHAP_incidence <- merge(SHAP_incidence_male, SHAP_incidence_female, by = "variable", all.x = TRUE, all.y = TRUE) %>%
  dplyr::mutate(relative_shap_male   = ifelse(is.na(relative_shap_male), 0, relative_shap_male),
                relative_shap_female = ifelse(is.na(relative_shap_female), 0, relative_shap_female)) %>%
  dplyr::mutate(diff_abs = relative_shap_male - relative_shap_female) %>%
  dplyr::select(variable, relative_shap_male, relative_shap_female, diff_abs)

plot_Db_incidence <- plot_shap_dumbbell(dat = SHAP_incidence, sort_by = "mean")
ggsave(filename = "./Results/05PredictiveM/C_incidence/SHAP_dumbbell.pdf", 
       plot = plot_Db_incidence, width = 3,height = 5)

#=================================protein*sex===================================
#---------------------------------Cmortality------------------------------------
C_mortality_sex <- C_mortality_preM %>%
  dplyr::mutate(across(any_of(factor_vars), as.factor)) %>%
  dplyr::mutate(all_death = as.numeric(as.character(all_death)))

sex_vars <- c("age", "BMI", "SBP_a", "smoking01", "alcohol01",
                "ALT", "Creatinine", "CRP", "HbA1c", "LDL", "Triglycerides", "nppb", "ntprobnp",
                "walk_pace", "MET_total", "FVC")
sig_pro_list <- toupper(sig_pro_list)

cox_mortality_sex <- mclapply(sig_pro_list, function(var) {
  formula_str <- paste0("Surv(FU_death, all_death) ~ ", var, " * sex + ", paste(sex_vars, collapse = " + "))

  cox_model <- coxph(as.formula(formula_str), data = C_mortality_sex)
  cox_tidy <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  cox_tidy_interaction <- cox_tidy %>%
    dplyr::filter(grepl(":", term)) %>%
    dplyr::filter(grepl(var, term)) %>%
    dplyr::mutate(variable_tested = var)
  return(cox_tidy_interaction)
}, mc.cores = 5)

cox_mortality_sex_all <- bind_rows(cox_mortality_sex)

#---------------------------------Cincidence------------------------------------
C_incidence_sex <- C_incidence_preM %>%
  dplyr::mutate(across(any_of(factor_vars), as.factor)) %>%
  dplyr::mutate(pah_dia = as.numeric(as.character(pah_dia)))

sex_vars <- c("age", "BMI", "SBP_a", "smoking01", "alcohol01",
              "ALT", "Creatinine", "CRP", "HbA1c", "LDL", "Triglycerides", "nppb", "ntprobnp",
              "walk_pace", "MET_total", "FVC")
sig_pro_list <- toupper(sig_pro_list)

cox_incidence_sex <- mclapply(sig_pro_list, function(var) {
  formula_str <- paste0("Surv(FU_pah_dia, pah_dia) ~ ", var, " * sex + ", paste(sex_vars, collapse = " + "))

  cox_model <- coxph(as.formula(formula_str), data = C_incidence_sex)
  cox_tidy <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  cox_tidy_interaction <- cox_tidy %>%
    dplyr::filter(grepl(":", term)) %>%
    dplyr::filter(grepl(var, term)) %>%
    dplyr::mutate(variable_tested = var)
  return(cox_tidy_interaction)
}, mc.cores = 5)

cox_incidence_sex_all <- bind_rows(cox_incidence_sex)

