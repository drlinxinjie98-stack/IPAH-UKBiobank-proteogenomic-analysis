library(dplyr)
library(ggplot2)
library(ggrepel)

MR_UKB_cis <- read.csv('./data/MR_PAH_v2/results_MR_UKB_cis.csv', header = TRUE)
MR_Fenland_cis <- read.csv('./data/MR_PAH_v2/results_MR_Fenland_cis.csv', header = TRUE) %>%
  dplyr::rename(protein_name = geneName)
MR_deCODE_cis <- read.csv('./data/MR_PAH_v2/results_MR_deCODE_cis.csv', header = TRUE) %>%
  dplyr::rename(protein_name = gene_name_Ensembl)

T_discovery_C1_sigP <- read.csv('./data/Discovery/Union_P/T_discovery_C1_union_P.csv', header = TRUE) %>%
  dplyr::rename(protein_name = term) %>%
  dplyr::mutate(protein_name = toupper(protein_name)) %>%
  dplyr::select(protein_name, estimate_cox, conf.low_cox, conf.high_cox, p.value_cox, p_adj_fdr_cox, direction_cox, source)
discovery_C1_sigP_list <- unique(as.character(T_discovery_C1_sigP$protein_name))

T_discovery_C2_sigP <- read.csv('./data/Discovery/Union_P/T_discovery_C2_union_P.csv', header = TRUE) %>%
  dplyr::rename(protein_name = term) %>%
  dplyr::mutate(protein_name = toupper(protein_name)) %>%
  dplyr::select(protein_name, estimate_cox, conf.low_cox, conf.high_cox, p.value_cox, p_adj_fdr_cox, direction_cox, source)
discovery_C2_sigP_list <- unique(as.character(T_discovery_C2_sigP$protein_name))

MR_validation <- function(discovery_sigP_list, MR_df) {
  protein_source_name <- deparse(substitute(MR_df))

  MR_dircon_list <- MR_df %>%
    dplyr::mutate(
      b_sign = dplyr::case_when(
        is.na(b) ~ NA_real_,
        b > 0 ~ 1,
        b < 0 ~ -1,
        TRUE ~ 0
      )
    ) %>%
    dplyr::group_by(exposure, outcome) %>%
    dplyr::filter(dplyr::n_distinct(b_sign[!is.na(b_sign) & b_sign != 0]) <= 1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(protein_name) %>%
    unique() %>%
    as.character()
  

  MR_sig <- MR_df %>%
    dplyr::filter(protein_name %in% discovery_sigP_list) %>%

    dplyr::filter(nsnp >= 3) %>%

    dplyr::filter(method == "Inverse variance weighted") %>%

    dplyr::mutate(pval_mr_adj = p.adjust(pval_mr, method = "BH")) %>%
    dplyr::filter(pval_mr < 0.05) %>%                          # nominal p<0.05
    dplyr::mutate(direction_mr = if_else(or > 1, "positive", "negative")) %>%
    dplyr::mutate(I2 = pmax(0, (Q_ivw - Q_df_ivw) / Q_ivw)) %>%
    dplyr::filter(I2 < 0.5) %>%

    dplyr::filter(pval_pleio > 0.05) %>%

    dplyr::filter(protein_name %in% MR_dircon_list) %>%

    dplyr::mutate(protein_source = protein_source_name) %>%  
    dplyr::select(protein_name, protein_source, exposure, outcome, method, nsnp, pval_mr, pval_mr_adj,
                  or, or_lci95, or_uci95, direction_mr, I2, pval_pleio)
  
  return(MR_sig)
}

#DO
#C1
MR_UKB_cis_C1_sig <- MR_validation(discovery_sigP_list = discovery_C1_sigP_list,
                                   MR_df = MR_UKB_cis)
MR_deCODE_cis_C1_sig <- MR_validation(discovery_sigP_list = discovery_C1_sigP_list,
                                   MR_df = MR_deCODE_cis)
MR_Fenland_cis_C1_sig <- MR_validation(discovery_sigP_list = discovery_C1_sigP_list,
                                   MR_df = MR_Fenland_cis)
MR_allthree_cis_C1_sig <- dplyr::bind_rows(MR_UKB_cis_C1_sig, MR_deCODE_cis_C1_sig, MR_Fenland_cis_C1_sig)
#write.csv(MR_allthree_cis_C1_sig, "./data/MR_PAH_v2/c1_corss_sig_MR.csv", row.names = FALSE)

#C2
MR_UKB_cis_C2_sig <- MR_validation(discovery_sigP_list = discovery_C2_sigP_list,
                                   MR_df = MR_UKB_cis)
MR_deCODE_cis_C2_sig <- MR_validation(discovery_sigP_list = discovery_C2_sigP_list,
                                   MR_df = MR_deCODE_cis)
MR_Fenland_cis_C2_sig <- MR_validation(discovery_sigP_list = discovery_C2_sigP_list,
                                   MR_df = MR_Fenland_cis)
MR_allthree_cis_C2_sig <- dplyr::bind_rows(MR_UKB_cis_C2_sig, MR_deCODE_cis_C2_sig, MR_Fenland_cis_C2_sig)
#write.csv(MR_allthree_cis_C2_sig, "./data/MR_PAH_v2/c2_corss_sig_MR.csv", row.names = FALSE)

T_cox_MR_cis_cross_C1 <- merge(T_discovery_C1_sigP, MR_allthree_cis_C1_sig, by = "protein_name", all.x = FALSE, all.y = TRUE)
write.csv(T_cox_MR_cis_cross_C1, "./data/Validation/c1_cox_MR_cis_cross.csv", row.names = FALSE)

T_cox_MR_cis_cross_C1_unique <- T_cox_MR_cis_cross_C1 %>%
  group_by(protein_name) %>%
  slice_min(order_by = pval_mr, n = 1, with_ties = FALSE) %>%
  ungroup()

T_cox_MR_cis_cross_C2 <- merge(T_discovery_C2_sigP, MR_allthree_cis_C2_sig, by = "protein_name", all.x = FALSE, all.y = TRUE)
write.csv(T_cox_MR_cis_cross_C2, "./data/Validation/c2_cox_MR_cis_cross.csv", row.names = FALSE)

T_cox_MR_cis_cross_C2_unique <- T_cox_MR_cis_cross_C2 %>%
  group_by(protein_name) %>%
  slice_min(order_by = pval_mr, n = 1, with_ties = FALSE) %>%
  ungroup()

#==================================Figure_validation_cross======================

#----C1----
MR_unique_C1 <- dplyr::bind_rows(MR_UKB_cis, MR_deCODE_cis, MR_Fenland_cis) %>%

  dplyr::filter(protein_name %in% discovery_C1_sigP_list) %>%

  dplyr::filter(nsnp >= 3) %>%

  dplyr::filter(method == "Inverse variance weighted") %>%

  group_by(protein_name) %>%
  slice_min(order_by = pval_mr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%

  dplyr::select(protein_name, pval_mr,
                or, or_lci95, or_uci95)
MR_cox_unique_C1 <- merge(MR_unique_C1, T_discovery_C1_sigP, by = 'protein_name',all.x = TRUE, all.y = FALSE) %>%
  dplyr::select(protein_name, pval_mr, or, or_lci95, or_uci95, p.value_cox, p_adj_fdr_cox, estimate_cox, conf.low_cox, conf.high_cox, source)

#----C2----
MR_unique_C2 <- dplyr::bind_rows(MR_UKB_cis, MR_deCODE_cis, MR_Fenland_cis) %>%

  dplyr::filter(protein_name %in% discovery_C2_sigP_list) %>%

  dplyr::filter(nsnp >= 3) %>%

  dplyr::filter(method == "Inverse variance weighted") %>%

  group_by(protein_name) %>%
  slice_min(order_by = pval_mr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%

  dplyr::select(protein_name, pval_mr,
                or, or_lci95, or_uci95)
MR_cox_unique_C2 <- merge(MR_unique_C2, T_discovery_C2_sigP, by = 'protein_name',all.x = TRUE, all.y = FALSE) %>%
  dplyr::select(protein_name, pval_mr, or, or_lci95, or_uci95, p.value_cox, p_adj_fdr_cox, estimate_cox, conf.low_cox, conf.high_cox, source)

#----C1----
scatter_C1 <- MR_cox_unique_C1 %>%
  mutate(s_group = case_when(pval_mr < 0.05 & p_adj_fdr_cox < 0.05 & source == "overall"    ~ "Pass MR (Both sexes)",
                             pval_mr < 0.05 & p.value_cox   < 0.05 & source == "male"       ~ "Pass MR (Male)",
                             pval_mr < 0.05 & p.value_cox   < 0.05 & source == "female"     ~ "Pass MR (Female)",
                             TRUE ~ "Fail MR"))%>%
  mutate(s_group = factor(s_group, levels = c("Fail MR", "Pass MR (Both sexes)", "Pass MR (Male)", "Pass MR (Female)"))) %>%
  arrange(s_group)
  
top_labels <- scatter_C1 %>%
  dplyr::filter(pval_mr < 0.05)

Figure_scat_C1 <- ggplot(scatter_C1, aes(x = log10(estimate_cox), y = log10(or))) +
  
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           fill = "pink3", alpha = 0.1) +
  
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           fill = "pink3", alpha = 0.1) +
  
  geom_errorbarh(aes(xmin = log10(conf.low_cox), xmax = log10(conf.high_cox), color = s_group), 
                 width = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = log10(or_lci95), ymax = log10(or_uci95), color = s_group), 
                width = 0, linewidth = 0.2) +
  geom_point(aes(size = -log10(pval_mr), color = s_group), alpha = 0.9) +
  geom_text_repel(data = top_labels, aes(label = protein_name), size = 6, max.overlaps = 20,
                  box.padding = 0.4, point.padding = 0.1, segment.size = 0.5, 
                  min.segment.length = 0, seed = 123) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("Pass MR (Both sexes)" = "#5B4B8A",
                                "Pass MR (Male)" = "#6B9CC5",
                                "Pass MR (Female)" = "#D26661",
                                "Fail MR" = "gray80")) +
  scale_size(range = c(1,5)) +
  labs(title = " ",
       x = expression(log[10]*"HR in Cox"),
       y = expression(log[10]*"OR in cis-MR"),
       color = "") +
  guides(size = "none") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 14))
ggsave(filename = "./Results/Validation/C1_scar.pdf", width = 10, height = 5, plot = Figure_scat_C1, bg="white")

#----C2----
scatter_C2 <- MR_cox_unique_C2 %>%
  mutate(s_group = case_when(pval_mr < 0.05 & p_adj_fdr_cox < 0.05 & source == "overall"    ~ "Pass MR (Both sexes)",
                             pval_mr < 0.05 & p.value_cox   < 0.05 & source == "male"       ~ "Pass MR (Male)",
                             pval_mr < 0.05 & p.value_cox   < 0.05 & source == "female"     ~ "Pass MR (Female)",
                             TRUE ~ "Fail MR")) %>%
  mutate(s_group = factor(s_group, levels = c("Fail MR", "Pass MR (Both sexes)", "Pass MR (Male)", "Pass MR (Female)"))) %>%
  arrange(s_group)

top_labels <- scatter_C2 %>%
  dplyr::filter(pval_mr < 0.05)

Figure_scat_C2 <- ggplot(scatter_C2, aes(x = log10(estimate_cox), y = log10(or))) +
  
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           fill = "pink3", alpha = 0.1) +
  
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           fill = "pink3", alpha = 0.1) +
  
  geom_errorbarh(aes(xmin = log10(conf.low_cox), xmax = log10(conf.high_cox), color = s_group), 
                 width = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = log10(or_lci95), ymax = log10(or_uci95), color = s_group), 
                width = 0, linewidth = 0.2) +
  geom_point(aes(size = -log10(pval_mr), color = s_group), alpha = 0.9) +
  geom_text_repel(data = top_labels, aes(label = protein_name), size = 6, max.overlaps = 20,
                  box.padding = 0.4, point.padding = 0.1, segment.size = 0.5, 
                  min.segment.length = 0, seed = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("Pass MR (Both sexes)" = "#5B4B8A",
                                "Pass MR (Male)" = "#6B9CC5",
                                "Pass MR (Female)" = "#D26661",
                                "Fail MR" = "gray80")) +
  scale_size(range = c(1,5)) +
  labs(title = " ",
       x = expression(log[10]*"HR in Cox"),
       y = expression(log[10]*"OR in cis-MR"),
       color = "") +
  guides(size = "none") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 14))
ggsave(filename = "./Results/Validation/C2_scar.pdf",width = 10,height = 5, plot = Figure_scat_C2, bg="white")

