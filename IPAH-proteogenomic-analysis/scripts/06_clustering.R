library(cluster)
library(factoextra)
library(survminer)
library(survival)
library(ggplot2)
library(tableone)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(lubridate)
library(readxl)
library(tidyr)
library(purrr)
library(broom)
library(emmeans)

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

#----C1----
T_Cohort1_to_cluster <- read.csv("./data/Cohorts/T_Cohort1.csv", header = TRUE) %>%
  dplyr::mutate(FU_death = FU_death/12) %>%

  dplyr::select(sex, all_death, FU_death, dplyr::any_of(sig_pro_list)) %>%
  dplyr::rename_with(toupper, .cols = 4:15)
  
#----C2----
T_Cohort2_to_cluster <- read.csv("./data/Cohorts/T_Cohort2.csv", header = TRUE) %>%
  dplyr::mutate(FU_death = FU_death/12) %>%

  dplyr::select(sex, all_death, FU_death, dplyr::any_of(sig_pro_list)) %>%
  dplyr::rename_with(toupper, .cols = 4:15)

wss <- sapply(1:5, function(k){
  kmeans(T_Cohort2_to_cluster, centers = k)$tot.withinss
})

avg_sil <- sapply(2:5, function(k){
  ss <- silhouette(kmeans(T_Cohort2_to_cluster, centers=k)$cluster, dist(T_Cohort2_to_cluster))
  mean(ss[, 3])
})

sig_pro_list <- unique(toupper(as.character(sig_prot$protein_name)))

cluster_with_pca <- function(data, vars = sig_pro_list, centers = 2, seed = 123) {
  
  set.seed(seed)
  
  cluster_mat <- data %>%
    dplyr::select(dplyr::any_of(vars)) %>%
    dplyr::rename_with(toupper)
  
  kmeans_result <- kmeans(cluster_mat, centers = centers)
  pca_res <- prcomp(cluster_mat, scale. = TRUE)
  
  pca_df <- data.frame(pca_res$x[, 1:2])
  pca_df$cluster <- factor(kmeans_result$cluster)
  
  clustered_data <- data %>%
    dplyr::bind_cols(pca_df)
  
  cluster_summary <- clustered_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      death_1_n = sum(all_death == 1, na.rm = TRUE),
      death_1_ratio = mean(all_death == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    clustered_data = clustered_data,
    cluster_summary = cluster_summary,
    pca_data = pca_df
  )
}

#DO
T_Cohort1_clustered <- cluster_with_pca(T_Cohort1_to_cluster)
T_Cohort2_clustered <- cluster_with_pca(T_Cohort2_to_cluster)

T_Cohort1_clustered$clustered_data$cluster <- factor(
  dplyr::recode(as.character(T_Cohort1_clustered$clustered_data$cluster), `1` = "2", `2` = "1"),
  levels = c("1", "2"))

T_Cohort1_clustered$pca_data$cluster <- factor(
  dplyr::recode(as.character(T_Cohort1_clustered$pca_data$cluster), `1` = "2", `2` = "1"),
  levels = c("1", "2"))

#write.csv(T_Cohort1_clustered$clustered_data, "./data/clustering/Cohort_preclinical_clustered.csv", row.names = FALSE)
#write.csv(T_Cohort2_clustered$clustered_data, "./data/clustering/Cohort_clinical_clustered.csv", row.names = FALSE)

scatter_plot <- ggplot(T_Cohort2_clustered$pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(title = "Clustering by robust proteins (k-means)",
       x = "Principal Component one",
       y = "Principal Component two") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",
                     labels = c("DVR", "RVR")) + 
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, by = 5)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 3)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
    axis.title.x = element_text(size = 16, face = "plain", color = "black"),
    axis.title.y = element_text(size = 16, face = "plain", color = "black"), 
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, face = "plain", color = "black"),
    legend.text = element_text(size = 14, color = "black"))

pca_df$cluster_bin <- ifelse(pca_df$cluster == 1, 0, 1)
fit <- glm(cluster_bin ~ PC1 + PC2, data = pca_df, family = binomial)
x_range <- seq(min(pca_df$PC1), max(pca_df$PC1), length.out = 200)
y_range <- seq(min(pca_df$PC2), max(pca_df$PC2), length.out = 200)
grid <- expand.grid(PC1 = x_range, PC2 = y_range)
grid$prob <- predict(fit, newdata = grid, type = "response")

#scatter_plot <- scatter_plot +
#  geom_contour(data = grid, aes(x = PC1, y = PC2, z = prob),
#               breaks = 0.5, color = "black", linetype = "dashed", linewidth = 0.5)

ggsave(filename = "./Results/Cluster/V3/C2_cluster_scatter.pdf", 
       width =8,height = 7, plot = scatter_plot, bg="white")

b <- coef(fit)
cat("Decision boundary: PC2 = -(", b[1], "+", b[2], "* PC1 ) /", b[3], "\n")

T_cohort1_clustered_male <- T_Cohort1_clustered$clustered_data %>%
  filter(sex == 1)
T_cohort1_clustered_female <- T_Cohort1_clustered$clustered_data %>%
  filter(sex == 0)

T_cohort2_clustered_male <- T_Cohort2_clustered$clustered_data %>%
  filter(sex == 1)
T_cohort2_clustered_female <- T_Cohort2_clustered$clustered_data %>%
  filter(sex == 0)

fit_KM_plot <- survfit(Surv(FU_death, all_death) ~ cluster, data = T_cohort1_clustered_female)
surv_object <- Surv(time = T_cohort1_clustered_female$FU_death, event = T_cohort1_clustered_female$all_death)
log_rank_test <- survdiff(surv_object ~ cluster, data = T_cohort1_clustered_female)
log_rank_test
pdf("./Results/Cluster/V3/C1_cluster_KM_female.pdf", width = 10, height = 10)
KM_plot <- ggsurvplot(fit_KM_plot, data = T_cohort1_clustered_female, censor = FALSE,
                    xlab = "Years", font.x = c(30, "plain", "black"),
                    ylab = "Freedom from mortality", font.y = c(30, "plain", "black"),
                    break.x.by = 3, break.y.by = 0.25,
                    xlim = c(0, 19), ylim = c(0, 1),
                    font.tickslab = c(30, "plain", "black"),
                    axes.offset = TRUE,
                    conf.int = TRUE, conf.int.alpha = 0.08, conf.int.style = "ribbon",
                    pval = FALSE, pval.method = TRUE, pval.size = 5, pval.coord = c(0, 0.1),
                    palette = c("#D0352B", "#4A7CB3"),
                    risk.table = TRUE, tables.y.text = FALSE, tables.height = 0.2, risk.table.pos = "out",
                    risk.table.fontsize = 8, risk.table.title = "Number at risk",
                    legend = "bottom", legend.title = "",
                    legend.labs = c("DVR", "RVR"),
                    font.legend = c(30, "bold", "black"),
                    tables.theme = theme(
                      plot.title = element_text(size = 24, face = "bold", color = "black", hjust = 0),
                      axis.text.y = element_blank(),
                      legend.position = "none",
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank()
                    ))

KM_plot$plot <- KM_plot$plot + 
  ggtitle("Female") +
  theme(plot.title = element_text(size = 34, face = "bold", hjust = 0.5, color = "black")) + 
  annotate("text", x = 2, y = 0.2, label = "Log-rank P = 0.002", 
           fontface = "bold.italic", size = 8, hjust = 0, vjust = 0, color = "black")
print(KM_plot, newpage = FALSE)
dev.off()

desired_order <- c("cluster1_male", "cluster1_female", "cluster2_male", "cluster2_female")

C_heatmap_data <- T_Cohort1_clustered$clustered_data %>%
  mutate(cluster = if_else(sex == 1,
                           paste0("cluster", cluster, "_male"),
                           paste0("cluster", cluster, "_female")
  ),
  cluster = factor(cluster, levels = desired_order)

  ) %>%
  dplyr::select(4:15, cluster)

cluster_vec <- C_heatmap_data$cluster
row_order <- order(cluster_vec)

scale_minus1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(0, length(x)))

  }
  2 * (x - rng[1]) / (rng[2] - rng[1]) - 1
}

protein_data <- C_heatmap_data[row_order, setdiff(colnames(C_heatmap_data), "cluster")]

desired_protein_order <- c(
  "EDN1", "ANXA2", "SPINT1",
  "LTBP2", "THBS2", "COL4A1", "QSOX1",
  "C7", "SLAMF7",
  "IGSF3", "LRRN1", "NPTX1"
)
protein_data <- protein_data[, intersect(desired_protein_order, colnames(protein_data)), drop = FALSE]

protein_data <- apply(protein_data, 2, scale_minus1_1)
protein_data <- as.matrix(protein_data)

annotation_row <- data.frame(Cluster = cluster_vec[row_order])
rownames(annotation_row) <- rownames(protein_data)

cluster_colors <- list(Cluster = c(cluster1_male = "#D62728", 
                                   cluster1_female = "#E3918F",
                                   cluster2_male = "#1F77B4",
                                   cluster2_female = "#A2BDD8"))

pdf("./Results/04Clustering/C1_Heatmap.pdf", width = 6, height = 12)
Heatmap_plot <- pheatmap(
  protein_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("#6B9CC5", "white", "#D26661"))(100),
  annotation_row = annotation_row,
  annotation_colors = cluster_colors,
  main = "Cluster proteomic signatures stratified by sex (preclinical, z-score)",
  display_numbers = FALSE,
  show_rownames = FALSE,

  fontsize = 9,
  fontsize_col = 10,
  angle_col = 90,
  border_color = NA,
  annotation_legend = FALSE)

print(Heatmap_plot, newpage = FALSE)
dev.off()

group_mean_male <- T_cohort1_clustered_male %>%
  group_by(cluster) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%

  dplyr::mutate(cluster = dplyr::recode(as.character(cluster), `1` = "DVR", `2` = "RVR"),
                cluster = paste0(cluster, " (male)")) %>%
  dplyr::select(cluster, 5:16)
group_mean_female <- T_cohort1_clustered_female %>%
  group_by(cluster) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%

  dplyr::mutate(cluster = dplyr::recode(as.character(cluster), `1` = "DVR", `2` = "RVR"),
                cluster = paste0(cluster, " (female)")) %>%
  dplyr::select(cluster, 5:16)

group_mean <- rbind(group_mean_male, group_mean_female)

group_mean_t <- t(group_mean)
colnames(group_mean_t) <- group_mean_t[1, ]
group_mean_t <- group_mean_t[-1, , drop = FALSE]

#gredient_protein <- group_mean_t %>%
#  as.data.frame() %>%
#  dplyr::filter((cluster1_male > cluster1_female) & ((cluster1_female > cluster2_male)|(cluster1_female > cluster2_female)) |
#                (cluster1_male < cluster1_female) & ((cluster1_female < cluster2_male)|(cluster1_female < cluster2_female))
#                )

protein_values <- group_mean[, -1]

rownames(protein_values) <- group_mean$cluster

normalize_minus1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))

  return(2 * (x - rng[1]) / diff(rng) - 1)
}
protein_values_norm <- as.data.frame(apply(protein_values, 2, normalize_minus1_1))

desired_order <- c("DVR (male)", "DVR (female)", 
                   "RVR (male)", "RVR (female)")
protein_values_norm <- protein_values_norm[desired_order, ]

desired_protein_order <- c("EDN1", "ANXA2", "SPINT1",
                           "LTBP2", "THBS2", "COL4A1", "QSOX1",
                           "C7", "SLAMF7",
                           "IGSF3", "LRRN1", "NPTX1")
protein_values_norm <- protein_values_norm[, intersect(desired_protein_order, colnames(protein_values_norm)), drop = FALSE]

pdf("./Results/Cluster/V3/C1_Heatmap_mean.pdf", width = 15, height = 3)
C2_Heatmap <- pheatmap(
  protein_values_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("#6B9CC5", "white", "#D26661"))(50),
  main = "",
  display_numbers = FALSE,
  fontsize = 16,
  fontsize_row = 14,
  fontsize_col = 14,
  border_color = NA,
  angle_col = 90)
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(AnnotationDbi)

T_original_prot <- readRDS('./data/Original/protein_imputed_missforest.rds')

res <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = '38928',
                             keytype = "ENTREZID",
                             columns = c("SYMBOL"))

background_protein <- data.frame(term = colnames(T_original_prot)[4:2923])
background_protein$term <- toupper(background_protein$term)
background_protein$term[background_protein$term == 'AMY1A_AMY1B_AMY1C'] <- 'AMY1A'
background_protein$term[background_protein$term == 'ARNTL'] <- 'BMAL1'
background_protein$term[background_protein$term == 'BAP18'] <- 'C17orf49'
background_protein$term[background_protein$term == 'BOLA2_BOLA2B'] <- 'BOLA2'
background_protein$term[background_protein$term == 'BTNL10'] <- 'BTNL10'
background_protein$term[background_protein$term == 'C19ORF12'] <- 'C19orf12'
background_protein$term[background_protein$term == 'C2ORF69'] <- 'C2orf69'
background_protein$term[background_protein$term == 'C7ORF50'] <- 'C7orf50'
background_protein$term[background_protein$term == 'C9ORF40'] <- 'C9orf40'
background_protein$term[background_protein$term == 'CERT'] <- 'CERT1'
background_protein$term[background_protein$term == 'CGB3_CGB5_CGB8'] <- 'CGB5'
background_protein$term[background_protein$term == 'CKMT1A_CKMT1B'] <- 'CKMT1A'
background_protein$term[background_protein$term == 'CTAG1A_CTAG1B'] <- 'CTAG1A'
background_protein$term[background_protein$term == 'DEFB103A_DEFB103B'] <- 'DEFB103A'
background_protein$term[background_protein$term == 'DDX58'] <- 'RIGI'
background_protein$term[background_protein$term == 'DEFA1_DEFA1B'] <- 'DEFA1'
background_protein$term[background_protein$term == 'DEFB104A_DEFB104B'] <- 'DEFB104A'
background_protein$term[background_protein$term == 'DEFB4A_DEFB4B'] <- 'DEFB4A'
background_protein$term[background_protein$term == 'EBI3_IL27'] <- 'EBI3'
background_protein$term[background_protein$term == 'ERVV_1'] <- 'ERVV-1'
background_protein$term[background_protein$term == 'FUT3_FUT5'] <- 'FUT3'
background_protein$term[background_protein$term == 'GBA'] <- 'GBA1'
background_protein$term[background_protein$term == 'GPR15L'] <- 'GPR15LG'
background_protein$term[background_protein$term == 'HLA_A'] <- 'HLA-A'
background_protein$term[background_protein$term == 'HLA_DRA'] <- 'HLA-DRA'
background_protein$term[background_protein$term == 'HLA_E'] <- 'HLA-E'
background_protein$term[background_protein$term == 'IL12A_IL12B'] <- 'IL12A'
background_protein$term[background_protein$term == 'LGALS7_LGALS7B'] <- 'LGALS7'
background_protein$term[background_protein$term == 'MICB_MICA'] <- 'MICB'
background_protein$term[background_protein$term == 'LEG1'] <- 'C6orf58'
background_protein$term[background_protein$term == 'MENT'] <- 'C1orf56'
background_protein$term[background_protein$term == 'MYLPF'] <- 'MYL2'
background_protein$term[background_protein$term == 'PALM2'] <- 'PALM2AKAP2'
background_protein$term[background_protein$term == 'SARG'] <- 'C1orf116'
background_protein$term[background_protein$term == 'NTPROBNP'] <- 'NPPB'
background_protein$term[background_protein$term == 'SKIV2L'] <- 'HLP'
background_protein$term[background_protein$term == 'SPACA5_SPACA5B'] <- 'SPACA5'
background_protein$term[background_protein$term == 'WARS'] <- 'WARS1'
background_protein$term[background_protein$term == 'BTNL10'] <- 'BTNL10P'

Overall_gene_ids <- bitr(background_protein$term, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

input_genes <- background_protein$term
mapped_genes <- Overall_gene_ids$SYMBOL
unmapped_genes <- setdiff(input_genes, mapped_genes)
unmapped_genes

Overall_gene_ids <- unique(Overall_gene_ids$ENTREZID)

C1C2_epi_MR_pro_lis_df <- as.data.frame(sig_pro_list)%>%
  dplyr::mutate(term = toupper(sig_pro_list))

C1C2_epi_MR_pro_ids <- bitr(C1C2_epi_MR_pro_lis_df$term, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)

input_genes <- C1C2_epi_MR_pro_lis_df$term
mapped <- bitr(input_genes,
               fromType = "SYMBOL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)
unmapped_genes <- setdiff(input_genes, mapped$SYMBOL)
unmapped_genes

#-GO
C1C2_epi_MR_pro_ego <- enrichGO(gene = C1C2_epi_MR_pro_ids$ENTREZID,
                   universe = Overall_gene_ids,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1) %>% as.data.frame()
#-KEGG
C1C2_epi_MR_pro_ekegg <- enrichKEGG(gene = C1C2_epi_MR_pro_ids$ENTREZID,
                       universe = Overall_gene_ids,
                       organism = 'hsa',
                       pvalueCutoff = 1,
                       qvalueCutoff = 1) %>% as.data.frame() %>%
  dplyr::filter(pvalue < 0.1) 
#-Reactome
C1C2_epi_MR_pro_reactome <- enrichPathway(gene = C1C2_epi_MR_pro_ids$ENTREZID, 
                             universe = Overall_gene_ids,
                             organism = "human", 
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             readable = TRUE) %>% as.data.frame() %>%
  dplyr::filter(pvalue < 0.1) 

#--KEGG--
#ECM-receptor interaction
KEGG_ECM_receptor_interaction <- C1C2_epi_MR_pro_ekegg %>%
  dplyr::filter(Description == 'ECM-receptor interaction') %>%
  dplyr::select(geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, "/")) %>%
  as.list()
#TGF-beta signaling pathway
KEGG_TGF_beta_sig_p <- C1C2_epi_MR_pro_ekegg %>%
  dplyr::filter(Description == 'TGF-beta signaling pathway') %>%
  dplyr::select(geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, "/")) %>%
  as.list()
#Hypertrophic cardiomyopathy
KEGG_HCM <- C1C2_epi_MR_pro_ekegg %>%
  dplyr::filter(Description == 'Hypertrophic cardiomyopathy') %>%
  dplyr::select(geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, "/")) %>%
  as.list()
#--Reactome--

Reactome_EFM_organization <- C1C2_epi_MR_pro_reactome %>%
  dplyr::filter(Description == 'Elastic fibre formation') %>%
  dplyr::select(geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, "/")) %>%

  as.list()
#IGFBP
Reactome_IGFBP_organization <- C1C2_epi_MR_pro_reactome %>%
  dplyr::filter(ID == 'R-HSA-381426') %>%
  dplyr::select(geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, "/")) %>%

  as.list()

KEGG_key_p <- unlist(c(KEGG_ECM_receptor_interaction, KEGG_TGF_beta_sig_p, KEGG_HCM))

KEGG_vector <- as.character(unlist(KEGG_key_p))
KEGG_symbols <- mapIds(org.Hs.eg.db,
                          keys = KEGG_vector,
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first") %>%
  unname()

KEGG_RECT <- unlist(c(KEGG_symbols, Reactome_EFM_organization, Reactome_IGFBP_organization))
KEGG_RECT <- tolower(KEGG_RECT)

C2_orig_heat_v2_mean <- C2_orig_heat_v2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
C2_orig_heat_v2_mean <- C2_orig_heat_v2_mean[, c("cluster", KEGG_RECT), drop = FALSE]

protein_values <- C2_orig_heat_v2_mean[, -1]

rownames(protein_values) <- C2_orig_heat_v2_mean$cluster

normalize_minus1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))

  return(2 * (x - rng[1]) / diff(rng) - 1)
}
protein_values_norm <- as.data.frame(apply(protein_values, 2, normalize_minus1_1))

desired_order <- c("cluster1_male", "cluster1_female", 
                   "cluster2_male", "cluster2_female")
protein_values_norm <- protein_values_norm[desired_order, ]

pdf("C2_Heatmap_KEGG_RECT_v2.pdf", width = 25, height = 5)
C2_Heatmap <- pheatmap(
  protein_values_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("#6B9CC5", "white", "#D26661"))(50),
  main = "",
  display_numbers = FALSE,
  fontsize = 30,
  fontsize_row = 30,
  fontsize_col = 30,
  border_color = NA,
  angle_col = 90)
dev.off()

T_cohort_clustered_sexinter <- T_Cohort1_clustered$clustered_data 
T_cohort_clustered_sexinter <- T_Cohort2_clustered$clustered_data 

options(width = 200)
dat_raw <- T_cohort_clustered_sexinter
protein_cols <- c(
  "ANXA2", "C7", "EDN1", "IGSF3", "LTBP2", "NPTX1",
  "QSOX1", "SLAMF7", "SPINT1", "THBS2", "COL4A1", "LRRN1"
)

stop_if_missing <- function(cols, data) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop("The following variables are missing from the data: ", paste(missing_cols, collapse = ", "))
  }
}

stop_if_missing(c("sex", "cluster", protein_cols), dat_raw)

non_numeric <- protein_cols[!vapply(dat_raw[protein_cols], is.numeric, logical(1))]
if (length(non_numeric)) {
  stop("The following protein variables are not numeric: ", paste(non_numeric, collapse = ", "))
}

dat <- dat_raw %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    cluster = factor(cluster, levels = c(2, 1), labels = c("RVR", "DVR"))
  )

if (any(is.na(dat$sex) & !is.na(dat_raw$sex))) {
  stop("The sex variable contains values other than 0 and 1.")
}
if (any(is.na(dat$cluster) & !is.na(dat_raw$cluster))) {
  stop("The cluster variable contains values other than 1 and 2.")
}

sample_distribution <- dat %>%
  count(sex, cluster, name = "N") %>%
  arrange(sex, cluster)

group_summary <- dat %>%
  select(sex, cluster, all_of(protein_cols)) %>%
  pivot_longer(all_of(protein_cols), names_to = "protein", values_to = "abundance") %>%
  group_by(protein, sex, cluster) %>%
  summarise(
    N = sum(!is.na(abundance)),
    mean = mean(abundance, na.rm = TRUE),
    SD = sd(abundance, na.rm = TRUE),
    median = median(abundance, na.rm = TRUE),
    Q1 = quantile(abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(abundance, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(mean, SD, median, Q1, Q3), ~ round(.x, 4))) %>%
  arrange(protein, sex, cluster)

fit_one_protein <- function(protein_name) {
  analysis_dat <- dat %>%
    select(sex, cluster, abundance = all_of(protein_name)) %>%
    drop_na()
  
  observed_groups <- nrow(distinct(analysis_dat, sex, cluster))
  if (observed_groups < 4) {
    return(list(
      protein = protein_name,
      error = paste0("Only ", observed_groups, " sex-cluster combinations were observed; the interaction cannot be estimated reliably.")
    ))
  }
  
  fit <- lm(abundance ~ cluster * sex, data = analysis_dat)
  
  coefficients <- tidy(fit, conf.int = TRUE, conf.level = 0.95)
  
  interaction <- coefficients %>%
    filter(grepl(":", term)) %>%
    transmute(
      protein = protein_name,
      N = nobs(fit),
      interaction_term = term,
      beta_interaction = estimate,
      standard_error = std.error,
      CI_lower = conf.low,
      CI_upper = conf.high,
      p_interaction = p.value
    )
  
  sex_specific <- emmeans(fit, ~ cluster | sex) %>%
    contrast(method = list("DVR - RVR" = c(-1, 1)), adjust = "none") %>%
    summary(infer = c(TRUE, TRUE), level = 0.95) %>%
    as.data.frame() %>%
    transmute(
      protein = protein_name,
      sex = as.character(sex),
      contrast = as.character(contrast),
      estimate_DVR_minus_RVR = estimate,
      standard_error = SE,
      df = df,
      CI_lower = lower.CL,
      CI_upper = upper.CL,
      p_value = p.value,
      direction = case_when(
        estimate_DVR_minus_RVR > 0 ~ "DVR > RVR",
        estimate_DVR_minus_RVR < 0 ~ "DVR < RVR",
        TRUE ~ "DVR = RVR"
      )
    )
  
  estimated_means <- emmeans(fit, ~ cluster * sex) %>%
    summary(infer = c(TRUE, TRUE), level = 0.95) %>%
    as.data.frame() %>%
    transmute(
      protein = protein_name,
      cluster = as.character(cluster),
      sex = as.character(sex),
      estimated_mean = emmean,
      standard_error = SE,
      df = df,
      CI_lower = lower.CL,
      CI_upper = upper.CL
    )
  
  list(
    protein = protein_name,
    fit = fit,
    interaction = interaction,
    sex_specific = sex_specific,
    estimated_means = estimated_means,
    coefficients = coefficients,
    error = NULL
  )
}

model_results <- map(protein_cols, fit_one_protein)
successful_models <- keep(model_results, ~ is.null(.x$error))

interaction_results <- successful_models %>%
  map_dfr("interaction") %>%
  mutate(
    FDR_interaction = p.adjust(p_interaction, method = "BH"),
    interaction_direction = case_when(
      beta_interaction > 0 ~ "Male DVR-RVR difference is greater than Female",
      beta_interaction < 0 ~ "Male DVR-RVR difference is smaller than Female",
      TRUE ~ "No difference"
    ),
    interaction_significant_P = if_else(p_interaction < 0.05, "Yes", "No"),
    interaction_significant_FDR = if_else(FDR_interaction < 0.05, "Yes", "No"),
    across(
      c(beta_interaction, standard_error, CI_lower, CI_upper, p_interaction, FDR_interaction),
      ~ signif(.x, 4)
    )
  ) %>%
  arrange(p_interaction)

sex_specific_results <- successful_models %>%
  map_dfr("sex_specific") %>%
  group_by(sex) %>%
  mutate(FDR_within_sex = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significant_P = if_else(p_value < 0.05, "Yes", "No"),
    significant_FDR = if_else(FDR_within_sex < 0.05, "Yes", "No"),
    across(
      c(estimate_DVR_minus_RVR, standard_error, CI_lower, CI_upper, p_value, FDR_within_sex),
      ~ signif(.x, 4)
    )
  ) %>%
  arrange(protein, sex)

estimated_mean_results <- successful_models %>%
  map_dfr("estimated_means") %>%
  mutate(across(c(estimated_mean, standard_error, CI_lower, CI_upper), ~ round(.x, 4))) %>%
  arrange(protein, sex, cluster)

model_errors <- keep(model_results, ~ !is.null(.x$error)) %>%
  map_dfr(~ data.frame(protein = .x$protein, error = .x$error))

result_bundle <- list(
  sample_distribution = sample_distribution,
  interaction_results = interaction_results,
  sex_specific_DVR_vs_RVR = sex_specific_results,
  estimated_group_means = estimated_mean_results,
  descriptive_group_summary = group_summary,
  LRRN1_interaction = filter(interaction_results, protein == "LRRN1"),
  LRRN1_sex_specific = filter(sex_specific_results, protein == "LRRN1"),
  LRRN1_estimated_means = filter(estimated_mean_results, protein == "LRRN1"),
  LRRN1_raw_summary = filter(group_summary, protein == "LRRN1"),
  model_errors = model_errors
)

print_section <- function(title, x, subtitle = NULL) {
  cat("\n============================================================\n")
  cat(title, "\n", sep = "")
  if (!is.null(subtitle)) cat(subtitle, "\n", sep = "")
  cat("============================================================\n")
  print(as.data.frame(x), row.names = FALSE)
}

print_section("1. SAMPLE DISTRIBUTION", result_bundle$sample_distribution)
print_section(
  "2. CLUSTER x SEX INTERACTION RESULTS FOR ALL PROTEINS",
  result_bundle$interaction_results,
  "Model: protein ~ cluster * sex; interaction: Male(DVR-RVR) - Female(DVR-RVR)"
)
print_section(
  "3. DVR - RVR DIFFERENCE WITHIN EACH SEX",
  result_bundle$sex_specific_DVR_vs_RVR,
  "estimate > 0: protein levels are higher in DVR than RVR; estimate < 0: protein levels are lower in DVR than RVR"
)
print_section("4. LRRN1: CLUSTER x SEX INTERACTION", result_bundle$LRRN1_interaction)
print_section("5. LRRN1: DVR - RVR DIFFERENCE WITHIN EACH SEX", result_bundle$LRRN1_sex_specific)
print_section("6. LRRN1: MODEL-ESTIMATED MEANS", result_bundle$LRRN1_estimated_means)
print_section("7. LRRN1: RAW DESCRIPTIVE STATISTICS", result_bundle$LRRN1_raw_summary)

if (nrow(result_bundle$model_errors) > 0) {
  print_section("8. MODEL ERRORS", result_bundle$model_errors)
}

cat("\n============================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("All results have been saved in the result_bundle object.\n")
cat("============================================================\n")

