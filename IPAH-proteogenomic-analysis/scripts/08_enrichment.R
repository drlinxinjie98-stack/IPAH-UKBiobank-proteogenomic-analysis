library(ReactomePA)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggprism)
library(ggforce)      
library(tibble)

T_original_prot <- readRDS('./data/Original/protein_imputed_missforest.rds') %>%
  dplyr::select(-c(nppb, ntprobnp))

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

C1C2_sigP_list <- union(discovery_C1_sigP_list, discovery_C2_sigP_list)

#============================cluster_analysis===================================

res <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = '38928',
                             keytype = "ENTREZID",
                             columns = c("SYMBOL"))
#--rename_map--
rename_map <- c(
  "AMY1A_AMY1B_AMY1C" = "AMY1A",
  "ARNTL" = "BMAL1",
  "BAP18" = "C17orf49",
  "BOLA2_BOLA2B" = "BOLA2",
  "BTNL10" = "BTNL10P",
  "C19ORF12" = "C19orf12",
  "C2ORF69" = "C2orf69",
  "C7ORF50" = "C7orf50",
  "C9ORF40" = "C9orf40",
  "CERT" = "CERT1",
  "CGB3_CGB5_CGB8" = "CGB5",
  "CKMT1A_CKMT1B" = "CKMT1A",
  "CTAG1A_CTAG1B" = "CTAG1A",
  "DEFB103A_DEFB103B" = "DEFB103A",
  "DDX58" = "RIGI",
  "DEFA1_DEFA1B" = "DEFA1",
  "DEFB104A_DEFB104B" = "DEFB104A",
  "DEFB4A_DEFB4B" = "DEFB4A",
  "EBI3_IL27" = "EBI3",
  "ERVV_1" = "ERVV-1",
  "FUT3_FUT5" = "FUT3",
  "GBA" = "GBA1",
  "GPR15L" = "GPR15LG",
  "HLA_A" = "HLA-A",
  "HLA_DRA" = "HLA-DRA",
  "HLA_E" = "HLA-E",
  "IL12A_IL12B" = "IL12A",
  "LGALS7_LGALS7B" = "LGALS7",
  "MICB_MICA" = "MICB",
  "LEG1" = "C6orf58",
  "MENT" = "C1orf56",
  "MYLPF" = "MYL2",
  "PALM2" = "PALM2AKAP2",
  "SARG" = "C1orf116",
  "SKIV2L" = "HLP",
  "SPACA5_SPACA5B" = "SPACA5",
  "WARS" = "WARS1"
)
#-rename
background_protein$term <- toupper(background_protein$term)
idx <- background_protein$term %in% names(rename_map)
background_protein$term[idx] <- rename_map[background_protein$term[idx]]

#-match

Overall_gene_ids <- bitr(background_protein$term, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

input_genes <- background_protein$term
mapped_genes <- Overall_gene_ids$SYMBOL
unmapped_genes <- setdiff(input_genes, mapped_genes)
unmapped_genes

Overall_gene_ids <- unique(Overall_gene_ids$ENTREZID)

#----------------------------C1+C2-------------------------------------------------
#--rename--
C1C2_sigP_list <- toupper(C1C2_sigP_list)

idx <- C1C2_sigP_list %in% names(rename_map)
C1C2_sigP_list[idx] <- rename_map[C1C2_sigP_list[idx]]

T_C1C2_union_gene_ids <- bitr(C1C2_sigP_list, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)
T_C1C2_union_gene_ids <- unique(T_C1C2_union_gene_ids$ENTREZID)

#-GO
C1C2_ego <- enrichGO(gene = T_C1C2_union_gene_ids,
                   universe = Overall_gene_ids,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1) %>% as.data.frame()
#-KEGG
options(timeout = 300)
C1C2_ekegg <- enrichKEGG(gene = T_C1C2_union_gene_ids,
                       universe = Overall_gene_ids,
                       organism = 'hsa',
                       pvalueCutoff = 1,
                       qvalueCutoff = 1) %>% as.data.frame()
#-Reactome
C1C2_reactome <- enrichPathway(gene = T_C1C2_union_gene_ids, 
                             universe = Overall_gene_ids,
                             organism = "human", 
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             readable = TRUE) %>% as.data.frame()
#---------------------------------GO--------------------------------------------
GO <- C1C2_ego %>% as_tibble()
GO_GO <- GO %>%
  arrange(p.adjust) %>%  
  top_n(10, wt = -pvalue) %>%  
  ungroup() %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,
                           levels = rev(c('BP', 'CC', 'MF')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

color_GO <- c("BP" = "#8ECFC9", "MF" = "#FFBE7A", "CC" = "#FA7F6F")

pdf("./Results/Discovery/Figure_enrich/C1C2_enrichGO.pdf", width = 7, height = 7)
GO_GO %>% 
  ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_col(aes(y = Description), width = 0.7, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 7) +
  geom_point(aes(x = -1, size = Count), shape = 21) +
  geom_text(aes(x = -1, label = Count, size = 7)) +
  scale_size_continuous(name = "Count", range = c(5, 16)) +
  geom_segment(data = data.frame(x = 0, y = 0, xend = 16, yend = 0),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = color_GO) +  
  labs(title = "GO", x = expression(-log[10]~"adjusted P"), y = NULL) +
  scale_x_continuous(breaks = seq(0, 17, 2), expand = expansion(c(0, 0)), limits = c(-2, 17)) +
  theme_prism() +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "plain"),
    axis.title = element_text(size = 20, face = 'bold'),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

dev.off()
#--------------------------------------KEGG-------------------------------------
KEGG <- C1C2_ekegg %>% as_tibble() %>%
  mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything()) 
KEGG_KEGG <- KEGG %>%
  arrange(p.adjust) %>%  
  top_n(5, wt = -pvalue) %>%  
  ungroup() %>%
  ungroup() %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

pdf("./Results/Discovery/Figure_enrich/C1C2_enrichKEGG.pdf", width = 7, height = 4)

KEGG_KEGG %>% 
  ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_col(aes(y = Description), width = 0.7, alpha = 0.8, fill = "#BEB8DC") +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 7) +
  geom_point(aes(x = -1, size = Count), shape = 21, fill = "#BEB8DC", color = "black") +
  geom_text(aes(x = -1, label = Count, size = 7)) +
  scale_size_continuous(name = "Count", range = c(5, 16)) +
  geom_segment(data = data.frame(x = 0, y = 0, xend = 16, yend = 0),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = "#BEB8DC") +  
  labs(title = "KEGG", x = expression(-log[10]~"adjusted P"), y = NULL) +
  scale_x_continuous(breaks = seq(0, 17, 2), expand = expansion(c(0, 0)), limits = c(-2, 17)) +
  theme_prism() +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "plain"),
    axis.title = element_text(size = 20, face = 'bold'),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
dev.off()
#-----------------------------------Reactome------------------------------------
REACTOME <- C1C2_reactome %>% as_tibble() %>%
  mutate(ONTOLOGY = "Reactome") %>%
  dplyr::select(ONTOLOGY, everything()) 
REACTOME_REACTOME <- REACTOME %>%
  arrange(p.adjust) %>%  
  top_n(10, wt = -pvalue) %>%  
  ungroup() %>%
  ungroup() %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

pdf("./Results/Discovery/Figure_enrich/C1C2_enrichReactome.pdf", width = 7, height = 9)

REACTOME_REACTOME %>% 
  ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_col(aes(y = Description), width = 0.7, alpha = 0.8, fill = "#E7DAD2") +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 7) +
  geom_point(aes(x = -1, size = Count), shape = 21, fill = "#E7DAD2", color = "black") +
  geom_text(aes(x = -1, label = Count, size = 7)) +
  scale_size_continuous(name = "Count", range = c(5, 16)) +
  geom_segment(data = data.frame(x = 0, y = 0, xend = 16, yend = 0),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = '#E7DAD2') +  
  labs(title = "Reactome", x = expression(-log[10]~"adjusted P"), y = NULL) +
  scale_x_continuous(breaks = seq(0, 17, 2), expand = expansion(c(0, 0)), limits = c(-2, 17)) +
  theme_prism() +
  coord_cartesian(clip = "off") +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "plain"),
    axis.title = element_text(size = 20, face = 'bold'),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
dev.off()

