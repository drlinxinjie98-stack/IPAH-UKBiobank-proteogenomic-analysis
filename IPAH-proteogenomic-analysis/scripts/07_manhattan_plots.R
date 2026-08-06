library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(cowplot)

#--C1--
T_C1_all <- read.csv("./data/Discovery/Cohort1/all_cox_c1.csv", header = TRUE)

T_C1_sig <- read.csv("./data/Discovery/Cohort1/all_cox_logit_c1_sig.csv", header = TRUE)

T_C1_male_all <- read.csv("./data/Discovery/Cohort1/male_cox_c1.csv", header = TRUE)
T_C1_male_sig <- read.csv("./data/Discovery/Cohort1/male_cox_logit_c1_sig.csv", header = TRUE)
T_C1_female_all <- read.csv("./data/Discovery/Cohort1/female_cox_c1.csv", header = TRUE)
T_C1_female_sig <- read.csv("./data/Discovery/Cohort1/female_cox_logit_c1_sig.csv", header = TRUE)
#--C2--
T_C2_all <- read.csv("./data/Discovery/Cohort2/all_cox_c2.csv", header = TRUE)
T_C2_sig <- read.csv("./data/Discovery/Cohort2/all_cox_logit_c2_sig.csv", header = TRUE)
T_C2_male_all <- read.csv("./data/Discovery/Cohort2/male_cox_c2.csv", header = TRUE)
T_C2_male_sig <- read.csv("./data/Discovery/Cohort2/male_logit_c2_sig.csv", header = TRUE)
T_C2_female_all <- read.csv("./data/Discovery/Cohort2/female_cox_c2.csv", header = TRUE)
T_C2_female_sig <- read.csv("./data/Discovery/Cohort2/female_cox_logit_c2_sig.csv", header = TRUE)

#===================================C1==========================================
#--overall--
C1_overall_Manh <- T_C1_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C1_sig$term, "validated", "non-validated"),
    direction = case_when(
      p_adj_fdr_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
      p_adj_fdr_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
      p_adj_fdr_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive ",
      p_adj_fdr_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
      TRUE ~ "not significant"),
    set = 'Both sexes')

#--male--
C1_male_Manh <- T_C1_male_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C1_male_sig$term, "validated", "non-validated"),
                direction = case_when(
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
                  TRUE ~ "not significant"),
                set = 'Male')
#--female--
C1_female_Manh <- T_C1_female_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C1_female_sig$term, "validated", "non-validated"),
                direction = case_when(
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
                  TRUE ~ "not significant"),
                set = 'Female')

#--union--
C1_union_Manh <- bind_rows(C1_overall_Manh, C1_male_Manh, C1_female_Manh) %>% 
  dplyr::select(term, estimate_cox, std.error_cox, statistic_cox, p.value_cox, conf.low_cox, conf.high_cox, variable_tested_cox, p_adj_fdr_cox, direction, set) %>%
  dplyr::mutate(term = toupper(term))

yr  <- max(abs(log10((C1_union_Manh$estimate_cox))), na.rm = TRUE)
mid_h <- 0.12 * yr

width_jit <- 0.28

T_C1_Manh <- C1_union_Manh %>%
  mutate(x_num = case_when(set == 'Both sexes' ~ 1, 
                           set == 'Male' ~ 2, 
                           set == 'Female' ~ 3,
                           TRUE ~ NA_real_),
    x_jit = x_num + runif(n(), -width_jit, width_jit))

T_C1_Manh_lab <- T_C1_Manh %>% 
  filter(direction %in% c("validated positive", "validated negative")) %>% 
  group_by(set, direction) %>% 
  arrange(p.value_cox, desc(abs(estimate_cox))) %>% 
  slice_head(n = 5) %>% 
  ungroup()

label_text <- unique(T_C1_Manh$set)
T_C1_Manh_band <- data.frame(
  set = factor(unique(T_C1_Manh$set)),
  x_num = 1:length(unique(T_C1_Manh$set)),
  label = label_text,
  y = 0,
  h = mid_h)

Figure_man_C1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = dplyr::filter(T_C1_Manh, direction == "not significant"),
             aes(x = x_jit, y = log10(estimate_cox)),
             colour = "grey80", size = 1.3, alpha = 0.35) +
  geom_point(data = dplyr::filter(T_C1_Manh, direction != "not significant"),
             aes(x = x_jit, y = log10(estimate_cox), colour = direction),
             size = 1.9, alpha = 0.9) +

  ggrepel::geom_label_repel(
    data = T_C1_Manh_lab,
    aes(x = x_jit, y = log10(estimate_cox), label = term),
    size = 4.5, colour = "black",
    label.size = 0.2, label.padding = unit(0.12, "lines"),
    seed = 2, box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, segment.size = 0.5, segment.color = "black",
    max.overlaps = Inf
  ) +

  geom_tile(data = T_C1_Manh_band,
            aes(x = x_num, y = y, height = 2*h, width = 0.8, fill = set),
            alpha = 0.3, color = NA) +
  geom_text(data = T_C1_Manh_band,
            aes(x = x_num, y = y, label = set),
            size = 5.5, fontface = "bold", color = "grey10") +
  scale_fill_manual(values = c("Both sexes" = "#BCBAB9", "Male" = "#74B3D7", "Female" = "#F28892"), guide = "none") + 
  scale_color_manual(values = c(`validated positive` = "#a73336", `non-validated positive` = "#fe8264", `validated negative` = "#333aab", `non-validated negative` = "#b5dbe6"),
                     name = NULL) +
  scale_x_continuous(breaks = T_C1_Manh_band$x_num, labels = T_C1_Manh_band$set, expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  labs(x = NULL, y = expression(log[10]~"HR"), title = "") +
  theme_half_open() +
  theme(
    panel.grid   = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22), 
    axis.text.x  = element_blank(),

    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y  = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",
    legend.text  = element_text(size = 15),
    legend.key.size = unit(12, "pt"),
    axis.line    = element_line(colour = "black", linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8))
ggsave(filename = "./Results/Discovery/Figure_Manh/C1_Manh.pdf", 
       width = 12,height = 5, plot = Figure_man_C1, bg="white")
#===================================C2==========================================
#--overall--
C2_overall_Manh <- T_C2_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C2_sig$term, "validated", "non-validated"),
                direction = case_when(
                  p_adj_fdr_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
                  p_adj_fdr_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
                  p_adj_fdr_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive ",
                  p_adj_fdr_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
                  TRUE ~ "not significant"),
                set = 'Both sexes')

#--male--
C2_male_Manh <- T_C2_male_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C2_male_sig$term, "validated", "non-validated"),
                direction = case_when(
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
                  TRUE ~ "not significant"),
                set = 'Male')
#--female--
C2_female_Manh <- T_C2_female_all %>%
  dplyr::mutate(dir_con = if_else(term %in% T_C2_female_sig$term, "validated", "non-validated"),
                direction = case_when(
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'validated' ~ "validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'validated' ~ "validated negative",
                  p.value_cox < 0.05 & estimate_cox > 1 & dir_con == 'non-validated' ~ "non-validated positive",
                  p.value_cox < 0.05 & estimate_cox < 1 & dir_con == 'non-validated' ~ "non-validated negative",
                  TRUE ~ "not significant"),
                set = 'Female')

#--union--
C2_union_Manh <- bind_rows(C2_overall_Manh, C2_male_Manh, C2_female_Manh) %>% 
  dplyr::select(term, estimate_cox, std.error_cox, statistic_cox, p.value_cox, conf.low_cox, conf.high_cox, variable_tested_cox, p_adj_fdr_cox, direction, set) %>%
  dplyr::mutate(term = toupper(term))

yr  <- max(abs(log10((C2_union_Manh$estimate))), na.rm = TRUE)
mid_h <- 0.07 * yr

width_jit <- 0.28
T_C2_Manh <- C2_union_Manh %>%
  mutate(x_num = case_when(set == 'Both sexes' ~ 1, 
                           set == 'Male' ~ 2, 
                           set == 'Female' ~ 3, 
                           TRUE ~ NA_real_),
         x_jit = x_num + runif(n(), -width_jit, width_jit))

T_C2_Manh_lab <- T_C2_Manh %>% 
  filter(direction %in% c("validated positive", "validated negative")) %>% 
  group_by(set, direction) %>% 
  arrange(p.value_cox, desc(abs(estimate_cox))) %>% 
  slice_head(n = 5) %>% 
  ungroup()

label_text <- unique(T_C2_Manh$set)

T_C2_Manh_band <- data.frame(
  set = factor(unique(T_C2_Manh$set)),
  x_num = 1:length(unique(T_C2_Manh$set)),
  label = label_text,
  y = 0,
  h = mid_h)

Figure_man_C2 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = dplyr::filter(T_C2_Manh, direction == "not significant"),
             aes(x = x_jit, y = log10(estimate_cox)),
             colour = "grey80", size = 1.3, alpha = 0.35) +
  geom_point(data = dplyr::filter(T_C2_Manh, direction != "not significant"),
             aes(x = x_jit, y = log10(estimate_cox), colour = direction),
             size = 1.9, alpha = 0.9) +

  ggrepel::geom_label_repel(
    data = T_C2_Manh_lab,
    aes(x = x_jit, y = log10(estimate_cox), label = term),
    size = 4.5, colour = "black",
    label.size = 0.2, label.padding = unit(0.1, "lines"),
    seed = 2, box.padding = 0.8, point.padding = 0.1, force = 2,
    min.segment.length = 0, segment.size = 0.5, segment.color = "black",
    max.overlaps = Inf
  ) +

  geom_tile(data = T_C2_Manh_band,
            aes(x = x_num, y = y, height = 2*h, width = 0.8, fill = set),
            alpha = 0.3, color = NA) +
  geom_text(data = T_C2_Manh_band,
            aes(x = x_num, y = y, label = set),
            size = 5.5, fontface = "bold", color = "grey10") +
  scale_fill_manual(values = c("Both sexes" = "#BCBAB9", "Male" = "#74B3D7", "Female" = "#F28892"), guide = "none") + 
  scale_color_manual(values = c(`validated positive` = "#a73336", `non-validated positive` = "#fe8264", `validated negative` = "#333aab", `non-validated negative` = "#b5dbe6"),
                     name = NULL) +
  scale_x_continuous(breaks = T_C2_Manh_band$x_num, labels = T_C2_Manh_band$set, expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  labs(x = NULL, y = expression(log[10]~"HR"), title ="") +
  theme_half_open() +
  theme(
    panel.grid   = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22), 
    axis.text.x  = element_blank(),

    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y  = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",
    legend.text  = element_text(size = 15),
    legend.key.size = unit(12, "pt"),
    axis.line    = element_line(colour = "black", linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8))
ggsave(filename = "./Results/Discovery/Figure_Manh/C2_Manh.pdf", 
       width = 12,height = 5, plot = Figure_man_C2, bg="white")

