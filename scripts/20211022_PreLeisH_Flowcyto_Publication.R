#################################################################
##                    Loading in the packages                   #
#################################################################

suppressMessages(library(here))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))
suppressMessages(library(lmerTest))

##################################################################
##                    Setting global variables                   #
##################################################################

options(ggrepel.max.overlaps = Inf)

text_font <- 'Roboto Condensed'

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#7AC5FF", "#C6FDEC", "#0348A6"))

#################################################################
##                      Loading in the data                     #
#################################################################

full_data <- suppressMessages(read_csv(here("data", "raw_data", "20221112_TimepointDB_PrimaryvChronic.csv"))) %>% filter(., Selected_Final == 'Yes') %>%
  mutate(Time_At_Risk_Years = Time_At_Risk_Days/365) %>% mutate(Relapse_Rate = Nr_Relapses_PreLeisH / Time_At_Risk_Years) %>% relocate(., c("Time_At_Risk_Years", "Relapse_Rate"), .after = Time_At_Risk_Days)

##----------------------------------------------------------------
##                          Boxplot data                         -
##----------------------------------------------------------------

boxplot_data <- full_data %>% mutate(RecodedTimepoint = case_when(
  RecodedTimepoint == 'W4' & Patient_SampleID != '0475' ~ 'EOT',
  RecodedTimepoint == 'W8' ~ 'EOT',
  TRUE ~ RecodedTimepoint
)) %>%
  mutate(RecodedTimepoint = replace_na(RecodedTimepoint, 'D0')) %>% filter(., Group != 'HEC') %>%
  filter(., RecodedTimepoint %in% c('D0'))

# temp_boxplot_EOT <- boxplot_data %>% filter(., RecodedTimepoint %in% c('EOT'))
# temp_boxplot_AL_HIV <- boxplot_data %>% filter(., Group %in% c('HIV', 'AL_HIV'))
# 
# boxplot_data <- bind_rows(temp_boxplot_EOT, temp_boxplot_AL_HIV)
# rm(temp_boxplot_EOT, temp_boxplot_AL_HIV)

boxplot_data <- boxplot_data %>% mutate(Aggregated_Group = case_when(Group %in% c('HIV', 'AL_HIV') ~ Group,
                                                                     Group %in% c('Primary_VL_HIV', 'Chronic_VL_HIV') ~ 'VL_HIV'))

boxplot_data <- boxplot_data %>% mutate_if(is.character, factor)

boxplot_data$Group <- factor(boxplot_data$Group, levels = c('HIV', 'AL_HIV', 'Primary_VL_HIV', 'Chronic_VL_HIV'))
boxplot_data$Aggregated_Group <- factor(boxplot_data$Aggregated_Group, levels = c('HIV', 'AL_HIV', 'VL_HIV'))

##---------------------------------------------------------------
##                        Timeseries data                       -
##---------------------------------------------------------------

timeseries_data <- full_data %>% filter(., Group %in% c('Primary_VL_HIV', 'Chronic_VL_HIV')) %>% filter(., !str_starts(RecodedTimepoint, 'Pre_')) %>% mutate_if(is.character, factor)
timeseries_data$Group <- forcats::fct_drop(timeseries_data$Group)

timeseries_data$Group <- factor(timeseries_data$Group, levels = c('Primary_VL_HIV', 'Chronic_VL_HIV'))
timeseries_data$RecodedTimepoint <- factor(timeseries_data$RecodedTimepoint, levels = c('D0', 'W1', 'W2', 'W3', 'W4', 'W6', 'W8',
                                                                        'Post_M3', 'Post_M4', 'Post_M5', 'Post_M6', 'Post_M7', 'Post_M8', 'Post_M9', 'Post_M10', 'Post_M11', 'Post_M12',
                                                                        'Post_M13', 'Post_M14', 'Post_M16', 'Post_M18', 'Post_M19', 'Post_M22'))

levels(timeseries_data$RecodedTimepoint) <- gsub("Post_", "", levels(timeseries_data$RecodedTimepoint))


##################################################################
##                          Statistics                           #
##################################################################

used_markers <- c('ICK_CD4+', 'ICT_CD3+CD8+_%', 'ICT_CD3-CD56+CD16+_%',
                  "ICK_CD4+_IFNg", 'ICK_CD4+_CXCR3_IFNg', 'ICK_CD4+_IL17', 'ICK_CD4+_FOXP3',
                  "EC_CD8-_PD1+", "EC_CD8-_LAG3+", "EC_CD8-_TIM3+", "ICT_CD3+CD8-_TIGIT", "EC_CD8-_KLRG1+", "EC_CD8-_CD57+",
                  "ICT_CD3+CD8+_IFNg", "ICT_CD3+CD8+_CD107", "ICK_CD4-_IL17",
                  "EC_CD8+_PD1+", "EC_CD8+_LAG3+", "EC_CD8+_TIM3+", "ICT_CD3+CD8+_TIGIT", "EC_CD8+_KLRG1+", "EC_CD8+_CD57+",
                  "ICT_CD3-CD56+CD16+_IFNg", "ICT_CD3-CD56+CD16+_CD107", "ICT_CD3-CD56+CD16+_TIGIT")

##---------------------------------------------------------------
##                        Cross-sectional                       -
##---------------------------------------------------------------

markers <- full_data %>% select(starts_with(c('ICK', 'ICT', 'EC'))) %>% select(!contains(c('Selection', 'Availability'))) %>% colnames(.)

# anova_between_all_groups <- sapply(markers[-43], function(x) {
#   aov_result <- boxplot_data %>% aov(unlist(.[x]) ~ Group, data = .)
#   anova(aov_result)$`Pr(>F)`[1]
# }, simplify = FALSE, USE.NAMES = T)

MWU_csect <- sapply(used_markers, function(x) {
  myformula <- as.formula(paste0("`", x,"`", " ", "~ Aggregated_Group"))
  boxplot_data %>% rstatix::wilcox_test(data = ., formula = myformula) %>% add_xy_position(x = "Aggregated_Group")
}, simplify = F, USE.NAMES = T)

MWU_hiv_v_ALHIV_pvals <- vector(length = length(used_markers))
MWU_hiv_v_VLHIV_pvals <- vector(length = length(used_markers))
MWU_ALHIV_v_VLHIV_pvals <- vector(length = length(used_markers))

for (i in seq_along(MWU_csect)) {
  MWU_hiv_v_ALHIV_pvals[[i]] <- MWU_csect[[i]] %>% filter(group1 == 'HIV', group2 == 'AL_HIV') %>% select(p) %>% pull
  MWU_hiv_v_VLHIV_pvals[[i]] <- MWU_csect[[i]] %>% filter(group1 == 'HIV', group2 == 'VL_HIV') %>% select(p) %>% pull
  MWU_ALHIV_v_VLHIV_pvals[[i]] <- MWU_csect[[i]] %>% filter(group1 == 'AL_HIV', group2 == 'VL_HIV') %>% select(p) %>% pull
}

adjusted_MWU_hiv_v_ALHIV_pvals <- p.adjust(MWU_hiv_v_ALHIV_pvals, method = 'BH')
adjusted_MWU_hiv_v_VLHIV_pvals <- p.adjust(MWU_hiv_v_VLHIV_pvals, method = 'BH')
adjusted_MWU_ALHIV_v_VLHIV_pvals <- p.adjust(MWU_ALHIV_v_VLHIV_pvals, method = 'BH')

for (i in seq_along(MWU_csect)) {
  MWU_csect[[i]] <- MWU_csect[[i]] %>% mutate(adjusted_p = case_when(group1 == 'HIV' & group2 == 'AL_HIV' ~ ifelse(adjusted_MWU_hiv_v_ALHIV_pvals[[i]] < 0.001,
                                                                                                                   "<0.001",
                                                                                                                   as.character(format(round(adjusted_MWU_hiv_v_ALHIV_pvals[[i]],3), nsmall=3))),
                                                                     group1 == 'HIV' & group2 == 'VL_HIV' ~ ifelse(adjusted_MWU_hiv_v_VLHIV_pvals[[i]] < 0.001,
                                                                                                                   "<0.001",
                                                                                                                   as.character(format(round(adjusted_MWU_hiv_v_VLHIV_pvals[[i]],3), nsmall = 3))),
                                                                     group1 == 'AL_HIV' & group2 == 'VL_HIV' ~ ifelse(adjusted_MWU_ALHIV_v_VLHIV_pvals[[i]] < 0.001,
                                                                                                                      "<0.001",
                                                                                                                      as.character(format(round(adjusted_MWU_ALHIV_v_VLHIV_pvals[[i]],3), nsmall = 3)))))
}

MWU_between_primary_chronic <- sapply(used_markers, function(x) {
  boxplot_data_only_primary_chronic <- filter(boxplot_data, Group %in% c('Primary_VL_HIV', 'Chronic_VL_HIV'))
  boxplot_data_only_primary_chronic$Group <- forcats::fct_drop(boxplot_data_only_primary_chronic$Group)
  mwu_primary_chronic <- boxplot_data_only_primary_chronic %>% wilcox.test(unlist(.[x]) ~ Group, data = .)
  mwu_primary_chronic$p.value
}, simplify = FALSE, USE.NAMES = T)

##----------------------------------------------------------------
##                          Longitudinal                         -
##----------------------------------------------------------------

##### GROUP * TIME
longitudinal_mixed_models <- sapply(markers[-43], function(x) {
  lmer_data <- timeseries_data %>% .[x] %>% pull
  lmer_test <- lmer(lmer_data ~ Group*RecodedTimepoint_in_Days + (1|Patient_SampleID), data = timeseries_data)
}, simplify = FALSE, USE.NAMES = TRUE)

# residual_plots <- sapply(longitudinal_mixed_models, function(x) {
#   ggqqplot(residuals(x))
# }, simplify = FALSE, USE.NAMES = TRUE)

longitudinal_anova_of_models <- sapply(longitudinal_mixed_models, function(x) {
  anova(x)
}, simplify = FALSE, USE.NAMES = TRUE)

Longitudinal_Group_pvals <- sapply(longitudinal_anova_of_models, function(x) {
  x$`Pr(>F)`[1]
}, simplify = T, USE.NAMES = T)

Timepoint_pvals <- sapply(longitudinal_anova_of_models, function(x) {
  x$`Pr(>F)`[2]
}, simplify = T, USE.NAMES = T)

Longitudinal_Interaction_pvals <- sapply(longitudinal_anova_of_models, function(x) {
  x$`Pr(>F)`[3]
}, simplify = T, USE.NAMES = T)

Longitudinal_significant_groups <- Longitudinal_Group_pvals[Longitudinal_Group_pvals < 0.05]
Longitudinal_significant_interactions <- Longitudinal_Interaction_pvals[Longitudinal_Interaction_pvals < 0.05]

##### GROUP * RELAPSE RATE
relapse_rate_mixed_models <- sapply(markers[-43], function(x) {
  lmer_data <- timeseries_data %>% .[x] %>% pull
  lmer_test <- lmer(lmer_data ~ Group*Relapse_Rate + (1|Patient_SampleID), data = timeseries_data)
}, simplify = FALSE, USE.NAMES = TRUE)

relapse_rate_anova_of_models <- sapply(relapse_rate_mixed_models, function(x) {
  anova(x)
}, simplify = FALSE, USE.NAMES = TRUE)

relapse_rate_Group_pvals <- sapply(relapse_rate_anova_of_models, function(x) {
  x$`Pr(>F)`[1]
}, simplify = T, USE.NAMES = T)

relapse_rate_pvals <- sapply(relapse_rate_anova_of_models, function(x) {
  x$`Pr(>F)`[2]
}, simplify = T, USE.NAMES = T)

relapse_rate_Interaction_pvals <- sapply(relapse_rate_anova_of_models, function(x) {
  x$`Pr(>F)`[3]
}, simplify = T, USE.NAMES = T)

relapse_rate_significant_groups <- relapse_rate_Group_pvals[relapse_rate_Group_pvals < 0.05]
relapse_rate_significant_interactions <- relapse_rate_Interaction_pvals[relapse_rate_Interaction_pvals < 0.05]

#################################################################
##                  Global plotting functions                   #
#################################################################

Cross_sectional_boxplot_generator <- function(marker_list, step_size_list) {
  
  boxplot_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    title <- names(marker_list)[[x]]
    plot_title <- paste0("%", sub(".*[%] *(.*?) *[(].*", "\\1", title), " at D0") #regex to capture everything between % and (
    y_axis_title <- paste0("%",sub(".*[(] *(.*?) *[)].*", "\\1", title)) #regex to capture everything between ( and )
    
    temp_data <- boxplot_data
    temp_data <- temp_data %>% filter(., Group %in% c('Primary_VL_HIV', 'Chronic_VL_HIV'))
    temp_data$Group <- forcats::fct_drop(temp_data$Group)
    temp_data$Aggregated_Group <- temp_data$Group
    
    combined_data <- bind_rows(boxplot_data, temp_data)
    
    combined_data$Aggregated_Group <- factor(combined_data$Aggregated_Group,
                                             levels = c('HEC', 'HIV', 'AL_HIV', 'VL_HIV', 'Primary_VL_HIV', 'Chronic_VL_HIV'))
    
    my_comparisons <- list(c('HIV', 'AL_HIV'), c('HIV', 'VL_HIV'), c('AL_HIV', 'VL_HIV'))
    
    overall_boxplot <- combined_data %>% filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = Aggregated_Group, y = .data[[tmp]], fill = Aggregated_Group)) +
      ggdist::stat_halfeye(adjust = .2, 
                           width = .6, 
                           justification = -.2, 
                           .width = 0, 
                           point_colour = NA,
                           aes(alpha = 0.5)) + 
      geom_boxplot(width = .175, 
                   outlier.color = NA) +
      #stat_compare_means(comparisons = my_comparisons, size = 7, step.increase = 0.275) +
      stat_pvalue_manual(MWU_csect[[tmp]], hide.ns=F, label="adjusted_p", label.size = 8, step.increase = step_size_list[x], tip.length = 0.025, inherit.aes = F) +
      geom_signif(textsize = 8, comparisons = list(c('Primary_VL_HIV', 'Chronic_VL_HIV')), annotations = format(MWU_between_primary_chronic[[tmp]], digits = 1, nsmall = 3)) +
      geom_point(position = position_nudge(x = -0.15)) +
      geom_vline(aes(xintercept=3.5), linetype="dashed",size=0.7) +
      ggtitle(label = plot_title) +
      labs(y = y_axis_title) +
      theme_bw() +
      scale_fill_manual(values = c("#00A087FF", "#3C5488FF", "#F39B7FFF", "#4DBBD5FF", "#E64B35FF")) +
      scale_x_discrete(expand = c(0.1,0.1))
    
    # overall_boxplot <- overall_boxplot + coord_cartesian(ylim = c(round(min(combined_data[[tmp]], na.rm = T),0), 
    #                                                               ifelse(round(max(combined_data[[tmp]], na.rm = T),0)*1.4 > 100 & round(max(combined_data[[tmp]], na.rm = T),0)*1.4 < 150,
    #                                                                      150,
    #                                                                      (round(max(combined_data[[tmp]], na.rm = T),0)*2)+0.5)))
    
    if (round(max(combined_data[[tmp]], na.rm = T)) <= 5) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 5, by = 1))
    } else if (round(max(combined_data[[tmp]], na.rm = T)) >= 5 && round(max(combined_data[[tmp]], na.rm = T)) <= 10) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 10, by = 2))
    } else if (round(max(combined_data[[tmp]], na.rm = T)) >= 10 && round(max(combined_data[[tmp]], na.rm = T)) <= 20) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 20, by = 4))
    } else if (round(max(combined_data[[tmp]], na.rm = T)) >= 20 && round(max(combined_data[[tmp]], na.rm = T)) <= 50) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 50, by = 10))
    } else if (round(max(combined_data[[tmp]], na.rm = T)) >= 50 && round(max(combined_data[[tmp]], na.rm = T)) <= 100) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 100, by = 20))
    }
    
    overall_boxplot <- overall_boxplot +
      theme(axis.text.y = element_text(size = 26, family = text_font, color = 'black'),
            axis.title.y = element_text(size = 26, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 30, family = text_font, face = 'bold', hjust = 0.5),
            legend.position = 'none',
            panel.grid.major.x = element_blank()) + ggpubr::rotate_x_text(angle = 40)
    
  })
  
  names(boxplot_list) <- names(marker_list)
  return(boxplot_list)
  
}

Longitudinal_timeseries_generator <- function(marker_list) {
  
  timeseries_marker_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    title <- names(marker_list)[[x]]
    ylab_simplified <- str_split(title, ' ', simplify = T)
    
    timeseries_data %>% 
      filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = RecodedTimepoint, y = .data[[tmp]])) +
      geom_line(aes(group = Patient_SampleID, linetype = Group, colour = Group), lwd = 1.5) +
      geom_point(aes(colour = Group, shape = Relapse, size = Relapse)) +
      labs(y = ifelse(length(ylab_simplified) > 2, paste(ylab_simplified[[1]], ylab_simplified[[2]]), ylab_simplified[[1]])) +
      ggtitle(paste("Temporal dynamics of", sub(".*[%] *(.*?) *[(].*", "\\1", title)),
              subtitle = paste0('\u0394 ', 'Group p = ', format(Longitudinal_Group_pvals[tmp], digits = 1, nsmall = 3))) +
      scale_color_manual(values = c("#4DBBD5FF", "#E64B35FF")) +
      scale_size_manual(values = c(4,8)) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 26, family = text_font, color = 'black'),
            axis.title.y = element_text(size = 26, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_blank(),
            plot.subtitle = element_text(size = 28, family = text_font, hjust = 0.5),
            plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
            legend.title = element_text(size = 26, family = text_font, face = 'bold'),
            legend.text = element_text(size = 26, family = text_font),
            legend.key.size = unit(2, "lines"),
            legend.position="bottom") + ggpubr::rotate_x_text() +
      guides(linetype = guide_legend(order = 1, title = NULL, override.aes = list(lwd = c(2,2), colour = c("#4DBBD5FF", "#E64B35FF"))),
             colour = 'none')
    
    # stat_smooth(data = subset(timeseries_data, Group %in% "Primary_VL_HIV"),
    #             method = "lm", se = T, size = 3) +
    # stat_smooth(data = subset(timeseries_data, Group %in% "Chronic_VL_HIV"),
    #             method = "lm", se = T, size = 3) +
    #geom_text(aes(label = Patient_SampleID)) +
    #stat_summary(aes(group = Group), fun = mean, geom = 'line', size = 3) +
    # geom_point(aes(colour = Group)) + geom_line(color = 'gray') +
  })
  
  names(timeseries_marker_list) <- names(marker_list)
  return(timeseries_marker_list)
}

Longitudinal_stat_model_generator <-  function(marker_list) {
 
  stat_model_marker_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    title <- names(marker_list)[[x]]
    ylab_simplified <- str_split(title, ' ', simplify = T)
    
    timeseries_data %>% 
      filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = Relapse_Rate, y = .data[[tmp]])) +
      geom_point(aes(colour = Group)) +
      stat_smooth(aes(colour = Group, fill = Group), se = T, method='lm', formula = y ~ x, size = 2.5) +
      labs(y = ifelse(length(ylab_simplified) > 2, paste(ylab_simplified[[1]], ylab_simplified[[2]]), ylab_simplified[[1]])) +
      ggtitle(paste0("%", sub(".*[%] *(.*?) *[(].*", "\\1", title), " in function of the relapse rate"),
              subtitle = paste0('\u0394 ', 'Group p = ', format(relapse_rate_Interaction_pvals[tmp], digits = 1, nsmall = 3))) +
      scale_color_manual(name = '', values = c("#4DBBD5FF", "#E64B35FF")) +
      scale_fill_manual(name = '', values = c("#4DBBD5FF", "#E64B35FF")) +
      scale_size_manual(values = c(4,8)) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 26, family = text_font),
            axis.title.y = element_text(size = 26, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_text(color='black', size=26, family=text_font),
            plot.subtitle = element_text(size = 28, family = text_font, hjust = 0.5),
            plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
            legend.title = element_text(size = 26, family = text_font, face = 'bold'),
            legend.text = element_text(size = 26, family = text_font),
            legend.position="bottom")
    
    # stat_smooth(data = subset(timeseries_data, Group %in% "Primary_VL_HIV"),
    #             method = "lm", se = T, size = 3) +
    # stat_smooth(data = subset(timeseries_data, Group %in% "Chronic_VL_HIV"),
    #             method = "lm", se = T, size = 3) +
    #geom_text(aes(label = Patient_SampleID)) +
    #stat_summary(aes(group = Group), fun = mean, geom = 'line', size = 3) +
    # geom_point(aes(colour = Group)) + geom_line(color = 'gray') +
  })
  
  names(stat_model_marker_list) <- names(marker_list)
  return(stat_model_marker_list)
   
}

##################################################################
##                      Cellular composition                     #
##################################################################

cell_marker_list <- list('%CD4+ T cells (CD3+CD4+)' = 'ICK_CD4+',
                         '%CD8+ T cells (CD3+CD8+)' = 'ICT_CD3+CD8+_%',
                         '%NK cells (CD3-CD56+CD16+)' = 'ICT_CD3-CD56+CD16+_%')

cross_sectional_composition <- Cross_sectional_boxplot_generator(cell_marker_list, step_size_list = c(0.175, 0.2, 0.15))

longitudinal_composition <- Longitudinal_timeseries_generator(cell_marker_list)

stat_composition <- Longitudinal_stat_model_generator(cell_marker_list)


layout <- "
AABBBCC
DDEEEFF
GGHHHII
"

composition_arranged <- wrap_elements(full = cross_sectional_composition$`%CD4+ T cells (CD3+CD4+)` + coord_cartesian(ylim = c(0, 120))) + wrap_elements(full = longitudinal_composition$`%CD4+ T cells (CD3+CD4+)`) + wrap_elements(full = stat_composition$`%CD4+ T cells (CD3+CD4+)`) +
  wrap_elements(full = cross_sectional_composition$`%CD8+ T cells (CD3+CD8+)` + coord_cartesian(ylim = c(0, 145))) + wrap_elements(full = longitudinal_composition$`%CD8+ T cells (CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) + wrap_elements(full = stat_composition$`%CD8+ T cells (CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = cross_sectional_composition$`%NK cells (CD3-CD56+CD16+)` + coord_cartesian(ylim = c(0, 75))) + wrap_elements(full = longitudinal_composition$`%NK cells (CD3-CD56+CD16+)`) + wrap_elements(full = stat_composition$`%NK cells (CD3-CD56+CD16+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold', vjust = 2))

# for (i in seq_along(composition_arranged$patches$layout$design)) {
#   if (i > 2) {
#     composition_arranged[[i]] <- composition_arranged[[i]] + plot_layout(tag_level = 'new')
#   }
# }

# for (i in seq_along(composition_arranged$patches$plots)) {
#   if (i %% 2 == 1) {
#     composition_arranged[[i+1]] <- composition_arranged[[i+1]] + plot_layout(tag_level = 'new')
#   }
# }

ggsave(here("analyses", "final_output", 'Figure2_Flowcyto_composition.pdf'), composition_arranged, dev = cairo_pdf, height = 17.5, width = 32.5)
ggsave(here("analyses", "final_output", 'Figure2_Flowcyto_composition.png'), composition_arranged, type = 'cairo-png', height = 17.5, width = 32.5)

ggsave(here("analyses", "interim_output", 'Figure2_Flowcyto_composition.png'), composition_arranged, type = 'cairo-png', height = 17.5, width = 32.5)

##################################################################
##                        CD4 Exhaustion                         #
##################################################################


CD4_exhaustion_list <- list('%PD1+ CD8- T cells (PD1+/CD3+CD8-)' = "EC_CD8-_PD1+",
                            '%LAG3+ CD8- T cells (LAG3+/CD3+CD8-)' = "EC_CD8-_LAG3+",
                            '%TIM3+ CD8- T cells (TIM3+/CD3+CD8-)' = "EC_CD8-_TIM3+",
                            '%TIGIT+ CD8- T cells (TIGIT+/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT")

cross_sectional_CD4_exhaustion <- Cross_sectional_boxplot_generator(CD4_exhaustion_list, step_size_list = c(0.15, 0.035, 0.05, 0.125))
longitudinal_CD4_exhaustion <- Longitudinal_timeseries_generator(CD4_exhaustion_list)
stat_CD4_exhaustion <- Longitudinal_stat_model_generator(CD4_exhaustion_list)

layout <- "
AABBBCC
DDEEEFF
GGHHHII
JJKKKLL
"

arranged_CD4_exhaustion <- wrap_elements(full = cross_sectional_CD4_exhaustion$`%PD1+ CD8- T cells (PD1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,147.5))) + wrap_elements(full = longitudinal_CD4_exhaustion$`%PD1+ CD8- T cells (PD1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) + wrap_elements(full = stat_CD4_exhaustion$`%PD1+ CD8- T cells (PD1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = cross_sectional_CD4_exhaustion$`%LAG3+ CD8- T cells (LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,110))) + wrap_elements(full = longitudinal_CD4_exhaustion$`%LAG3+ CD8- T cells (LAG3+/CD3+CD8-)`) + wrap_elements(full = stat_CD4_exhaustion$`%LAG3+ CD8- T cells (LAG3+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_CD4_exhaustion$`%TIM3+ CD8- T cells (TIM3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,32.5))) + wrap_elements(full = longitudinal_CD4_exhaustion$`%TIM3+ CD8- T cells (TIM3+/CD3+CD8-)`) + wrap_elements(full = stat_CD4_exhaustion$`%TIM3+ CD8- T cells (TIM3+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_CD4_exhaustion$`%TIGIT+ CD8- T cells (TIGIT+/CD3+CD8-)` + coord_cartesian(ylim = c(0,145))) + wrap_elements(full = longitudinal_CD4_exhaustion$`%TIGIT+ CD8- T cells (TIGIT+/CD3+CD8-)`) + wrap_elements(full = stat_CD4_exhaustion$`%TIGIT+ CD8- T cells (TIGIT+/CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "Figure4_Flowcyto_CD4Exhaustion.pdf"), arranged_CD4_exhaustion, dev = cairo_pdf, height = 22.5, width = 35)
ggsave(here("analyses", "final_output", "Figure4_Flowcyto_CD4Exhaustion.png"), arranged_CD4_exhaustion, type = 'cairo-png', height = 22.5, width = 35)

ggsave(here("analyses", "interim_output", "Figure4_Flowcyto_CD4Exhaustion.png"), arranged_CD4_exhaustion, type = 'cairo-png', height = 22.5, width = 35)

##################################################################
##                        CD8 Exhaustion                         #
##################################################################

CD8_exhaustion_list <- list('%PD1+ CD8+ T cells (PD1+/CD3+CD8+)' = "EC_CD8+_PD1+",
                            '%LAG3+ CD8+ T cells (LAG3+/CD3+CD8+)' = "EC_CD8+_LAG3+",
                            '%TIM3+ CD8+ T cells (TIM3+/CD3+CD8+)' = "EC_CD8+_TIM3+",
                            '%TIGIT+ CD8+ T cells (TIGIT+/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT")

cross_sectional_CD8_exhaustion <- Cross_sectional_boxplot_generator(CD8_exhaustion_list, step_size = c(0.15, 0.045, 0.125, 0.2))
longitudinal_CD8_exhaustion <- Longitudinal_timeseries_generator(CD8_exhaustion_list)
stat_CD8_exhaustion <- Longitudinal_stat_model_generator(CD8_exhaustion_list)

layout <- "
AABBBCC
DDEEEFF
GGHHHII
JJKKKLL
"

arranged_CD8_exhaustion <- wrap_elements(full = cross_sectional_CD8_exhaustion$`%PD1+ CD8+ T cells (PD1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,165))) + wrap_elements(full = longitudinal_CD8_exhaustion$`%PD1+ CD8+ T cells (PD1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) + wrap_elements(full = stat_CD8_exhaustion$`%PD1+ CD8+ T cells (PD1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = cross_sectional_CD8_exhaustion$`%LAG3+ CD8+ T cells (LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,127.5))) + wrap_elements(full = longitudinal_CD8_exhaustion$`%LAG3+ CD8+ T cells (LAG3+/CD3+CD8+)`) + wrap_elements(full = stat_CD8_exhaustion$`%LAG3+ CD8+ T cells (LAG3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_CD8_exhaustion$`%TIM3+ CD8+ T cells (TIM3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,35))) + wrap_elements(full = longitudinal_CD8_exhaustion$`%TIM3+ CD8+ T cells (TIM3+/CD3+CD8+)`) + wrap_elements(full = stat_CD8_exhaustion$`%TIM3+ CD8+ T cells (TIM3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_CD8_exhaustion$`%TIGIT+ CD8+ T cells (TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,160))) + wrap_elements(full = longitudinal_CD8_exhaustion$`%TIGIT+ CD8+ T cells (TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) + wrap_elements(full = stat_CD8_exhaustion$`%TIGIT+ CD8+ T cells (TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "Figure5_Flowcyto_CD8Exhaustion.pdf"), arranged_CD8_exhaustion, dev = cairo_pdf, height = 22.5, width = 35)
ggsave(here("analyses", "final_output", "Figure5_Flowcyto_CD8Exhaustion.png"), arranged_CD8_exhaustion, type = 'cairo-png', height = 22.5, width = 35)

ggsave(here("analyses", "interim_output", "Figure5_Flowcyto_CD8Exhaustion.png"), arranged_CD8_exhaustion, type = 'cairo-png', height = 22.5, width = 35)


#################################################################
##                        Supplementary                         #
#################################################################

##----------------------------------------------------------------
##                        CD4 Functional                         -
##----------------------------------------------------------------

CD4_functional_list <- list('%IFNg+ CD4+ T cells (IFNg+/CD3+CD4+)' = "ICK_CD4+_IFNg",
                            #'%Th1 cells (IFNg+/CXCR3+/CD3+CD4+)' = 'ICK_CD4+_CXCR3_IFNg', #rather for supplementary
                            '%IL17+ CD4+ T cells (IL17+/CD3+CD4+)' = 'ICK_CD4+_IL17',
                            '%Tregs (FOXP3+/CD3+CD4+)' = 'ICK_CD4+_FOXP3')

cross_sectional_CD4_functional <- Cross_sectional_boxplot_generator(CD4_functional_list, step_size_list = c(0.175, 0.175, 0.225))
longitudinal_CD4_functional <- Longitudinal_timeseries_generator(CD4_functional_list)
stat_CD4_functional <- Longitudinal_stat_model_generator(CD4_functional_list)

layout <- "
AABBBCC
DDEEEFF
GGHHHII
"

arranged_CD4_functional <- wrap_elements(full = cross_sectional_CD4_functional$`%IFNg+ CD4+ T cells (IFNg+/CD3+CD4+)` + coord_cartesian(ylim = c(0,142.5))) + wrap_elements(full = longitudinal_CD4_functional$`%IFNg+ CD4+ T cells (IFNg+/CD3+CD4+)`) + wrap_elements(full = stat_CD4_functional$`%IFNg+ CD4+ T cells (IFNg+/CD3+CD4+)` + coord_cartesian(ylim = c(0,80))) +
  wrap_elements(full = cross_sectional_CD4_functional$`%IL17+ CD4+ T cells (IL17+/CD3+CD4+)` + coord_cartesian(ylim = c(0,11))) + wrap_elements(full = longitudinal_CD4_functional$`%IL17+ CD4+ T cells (IL17+/CD3+CD4+)`) + wrap_elements(full = stat_CD4_functional$`%IL17+ CD4+ T cells (IL17+/CD3+CD4+)`) +
  wrap_elements(full = cross_sectional_CD4_functional$`%Tregs (FOXP3+/CD3+CD4+)`+ coord_cartesian(ylim = c(0,19))) + wrap_elements(full = longitudinal_CD4_functional$`%Tregs (FOXP3+/CD3+CD4+)`) + wrap_elements(full = stat_CD4_functional$`%Tregs (FOXP3+/CD3+CD4+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold', vjust = 2))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure1_Flowcyto_CD4Functional.pdf"), arranged_CD4_functional, dev = cairo_pdf, height = 17.5, width = 34)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure1_Flowcyto_CD4Functional.png"), arranged_CD4_functional, type = 'cairo-png', height = 17.5, width = 34)

ggsave(here("analyses", "interim_output", "SuppFigure1_Flowcyto_CD4Functional.png"), arranged_CD4_functional, type = 'cairo-png', height = 17.5, width = 34)

##---------------------------------------------------------------
##                              Th1                             -
##---------------------------------------------------------------

# Supp_list_Th1 <- c('%Th1 cells (IFNg+/CXCR3+/CD3+CD4+)' = 'ICK_CD4+_CXCR3_IFNg')
# 
# cross_sectional_Supp_Th1 <- Cross_sectional_boxplot_generator(Supp_list_Th1, step_size_list = c(0.1))
# longitudinal_Supp_Th1 <- Longitudinal_timeseries_generator(Supp_list_Th1)
# stat_Supp_Th1 <- Longitudinal_stat_model_generator(Supp_list_Th1)
# 
# arranged_Supp_Th1 <- wrap_elements(full = cross_sectional_Supp_Th1$`%Th1 cells (IFNg+/CXCR3+/CD3+CD4+)` + coord_cartesian(ylim = c(0,107.5))) +
#   longitudinal_Supp_Th1$`%Th1 cells (IFNg+/CXCR3+/CD3+CD4+)` +
#   stat_Supp_Th1$`%Th1 cells (IFNg+/CXCR3+/CD3+CD4+)` +
#   plot_layout(design = c("AABBBCC")) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))
# 
# ggsave(here("analyses", "final_output", "supplemental", "SuppFigure2_Flowcyto_Th1.pdf"), arranged_Supp_Th1, dev = cairo_pdf, height = 7.5, width = 30)
# ggsave(here("analyses", "final_output", SuppFigure2_Flowcyto_Th1.png"), arranged_Supp_Th1, type = 'cairo-png', height = 7.5, width = 30)

##----------------------------------------------------------------
##                        CD8 Functional                         -
##----------------------------------------------------------------

CD8_functional_list <- list('%IFNg+ CD8+ T cells (IFNg+/CD3+CD8+)' = "ICT_CD3+CD8+_IFNg",
                            '%IL17+ CD8+ T cells (IL17+/CD3+CD4-)' = "ICK_CD4-_IL17",
                            '%CD107+ CD8+ T cells (CD107+/CD3+CD8+)' = "ICT_CD3+CD8+_CD107")


cross_sectional_CD8_functional <- Cross_sectional_boxplot_generator(CD8_functional_list, step_size_list = c(0.175, 0.075, 0.175))
longitudinal_CD8_functional <- Longitudinal_timeseries_generator(CD8_functional_list)
stat_CD8_functional <- Longitudinal_stat_model_generator(CD8_functional_list)


layout <- "
AABBBCC
DDEEEFF
GGHHHII
"

arranged_CD8_functional <- wrap_elements(full = cross_sectional_CD8_functional$`%IFNg+ CD8+ T cells (IFNg+/CD3+CD8+)` + coord_cartesian(ylim = c(0,122.5))) + wrap_elements(full = longitudinal_CD8_functional$`%IFNg+ CD8+ T cells (IFNg+/CD3+CD8+)`) + wrap_elements(full = stat_CD8_functional$`%IFNg+ CD8+ T cells (IFNg+/CD3+CD8+)` + coord_cartesian(ylim = c(0,80))) +
  wrap_elements(full = cross_sectional_CD8_functional$`%IL17+ CD8+ T cells (IL17+/CD3+CD4-)`  + coord_cartesian(ylim = c(0,57.5))) + wrap_elements(full = longitudinal_CD8_functional$`%IL17+ CD8+ T cells (IL17+/CD3+CD4-)`) + wrap_elements(full = stat_CD8_functional$`%IL17+ CD8+ T cells (IL17+/CD3+CD4-)`) +
  wrap_elements(full = cross_sectional_CD8_functional$`%CD107+ CD8+ T cells (CD107+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) + wrap_elements(full = longitudinal_CD8_functional$`%CD107+ CD8+ T cells (CD107+/CD3+CD8+)`) + wrap_elements(full = stat_CD8_functional$`%CD107+ CD8+ T cells (CD107+/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental",  "SuppFigure4_Flowcyto_CD8Functional.pdf"), arranged_CD8_functional, dev = cairo_pdf, height = 17.5, width = 35)
ggsave(here("analyses", "final_output", "supplemental",  "SuppFigure4_Flowcyto_CD8Functional.png"), arranged_CD8_functional, type = 'cairo-png', height = 17.5, width = 35)

ggsave(here("analyses", "interim_output", "SuppFigure4_Flowcyto_CD8Functional.png"), arranged_CD8_functional, type = 'cairo-png', height = 17.5, width = 35)

##----------------------------------------------------------------
##                          Senescence                           -
##----------------------------------------------------------------

Supp_list_Senescence <- c('%KLRG1+ CD8- T cells (KLRG1+/CD3+CD8-)' = "EC_CD8-_KLRG1+",
                          '%CD57+ CD8- T cells (CD57+/CD3+CD8-)' = "EC_CD8-_CD57+",
                          '%KLRG1+ CD8+ T cells (KLRG1+/CD3+CD8+)' = "EC_CD8+_KLRG1+",
                          '%CD57+ CD8+ T cells (CD57+/CD3+CD8+)' = "EC_CD8+_CD57+")

cross_sectional_Supp_Senescence <- Cross_sectional_boxplot_generator(Supp_list_Senescence, step_size_list = c(0.15, 0.125, 0.2, 0.225))
longitudinal_Supp_Senescence <- Longitudinal_timeseries_generator(Supp_list_Senescence)
stat_Supp_Senescence <- Longitudinal_stat_model_generator(Supp_list_Senescence)

layout <- "
AABBBCC
DDEEEFF
GGHHHII
JJKKKLL
"

arranged_Supp_Senescence <- wrap_elements(full = cross_sectional_Supp_Senescence$`%KLRG1+ CD8- T cells (KLRG1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,127.5))) + wrap_elements(full = longitudinal_Supp_Senescence$`%KLRG1+ CD8- T cells (KLRG1+/CD3+CD8-)`) + wrap_elements(full = stat_Supp_Senescence$`%KLRG1+ CD8- T cells (KLRG1+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Supp_Senescence$`%CD57+ CD8- T cells (CD57+/CD3+CD8-)` + coord_cartesian(ylim = c(0,105))) + wrap_elements(full = longitudinal_Supp_Senescence$`%CD57+ CD8- T cells (CD57+/CD3+CD8-)`) + wrap_elements(full = stat_Supp_Senescence$`%CD57+ CD8- T cells (CD57+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Supp_Senescence$`%KLRG1+ CD8+ T cells (KLRG1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,150))) + wrap_elements(full = longitudinal_Supp_Senescence$`%KLRG1+ CD8+ T cells (KLRG1+/CD3+CD8+)`) + wrap_elements(full = stat_Supp_Senescence$`%KLRG1+ CD8+ T cells (KLRG1+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Senescence$`%CD57+ CD8+ T cells (CD57+/CD3+CD8+)` + coord_cartesian(ylim = c(0,140))) + wrap_elements(full = longitudinal_Supp_Senescence$`%CD57+ CD8+ T cells (CD57+/CD3+CD8+)`) + wrap_elements(full = stat_Supp_Senescence$`%CD57+ CD8+ T cells (CD57+/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure3_Flowcyto_Senescence.pdf"), arranged_Supp_Senescence, dev = cairo_pdf, height = 22.5, width = 35)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure3_Flowcyto_Senescence.png"), arranged_Supp_Senescence, type = 'cairo-png', height = 22.5, width = 35)

ggsave(here("analyses", "interim_output", "SuppFigure3_Flowcyto_Senescence.png"), arranged_Supp_Senescence, type = 'cairo-png', height = 22.5, width = 35)

##----------------------------------------------------------------
##                              NK                               -
##----------------------------------------------------------------

Supp_list_NK <- c("%IFNg+ NK cells (IFNg+/CD3-CD56+CD16+)" = "ICT_CD3-CD56+CD16+_IFNg",
                  "%CD107+ NK cells (CD107+/CD3-CD56+CD16+)" = "ICT_CD3-CD56+CD16+_CD107",
                  "%TIGIT+ NK cells (TIGIT+/CD3-CD56+CD16+)" = "ICT_CD3-CD56+CD16+_TIGIT")

cross_sectional_Supp_NK <- Cross_sectional_boxplot_generator(Supp_list_NK, step_size = c(0.135, 0.135, 0.135))
longitudinal_Supp_NK <- Longitudinal_timeseries_generator(Supp_list_NK)
stat_Supp_NK <- Longitudinal_stat_model_generator(Supp_list_NK)

layout <- "
AABBBCC
DDEEEFF
GGHHHII
"

arranged_Supp_NK <- wrap_elements(full = cross_sectional_Supp_NK$`%IFNg+ NK cells (IFNg+/CD3-CD56+CD16+)` + coord_cartesian(ylim = c(0,102.5))) + wrap_elements(full = longitudinal_Supp_NK$`%IFNg+ NK cells (IFNg+/CD3-CD56+CD16+)`) + wrap_elements(full = stat_Supp_NK$`%IFNg+ NK cells (IFNg+/CD3-CD56+CD16+)`) +
  wrap_elements(full = cross_sectional_Supp_NK$`%CD107+ NK cells (CD107+/CD3-CD56+CD16+)` + coord_cartesian(ylim = c(0,85))) + wrap_elements(full = longitudinal_Supp_NK$`%CD107+ NK cells (CD107+/CD3-CD56+CD16+)`) + wrap_elements(full = stat_Supp_NK$`%CD107+ NK cells (CD107+/CD3-CD56+CD16+)`) +
  wrap_elements(full = cross_sectional_Supp_NK$`%TIGIT+ NK cells (TIGIT+/CD3-CD56+CD16+)` + coord_cartesian(ylim = c(0,137.5))) + wrap_elements(full = longitudinal_Supp_NK$`%TIGIT+ NK cells (TIGIT+/CD3-CD56+CD16+)`) + wrap_elements(full = stat_Supp_NK$`%TIGIT+ NK cells (TIGIT+/CD3-CD56+CD16+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure7_Flowcyto_NK.pdf"), arranged_Supp_NK, dev = cairo_pdf, height = 20, width = 35)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure7_Flowcyto_NK.png"), arranged_Supp_NK, type = 'cairo-png', height = 20, width = 35)

ggsave(here("analyses", "interim_output", "SuppFigure7_Flowcyto_NK.png"), arranged_Supp_NK, type = 'cairo-png', height = 20, width = 35)
