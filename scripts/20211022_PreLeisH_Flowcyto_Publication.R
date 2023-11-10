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

full_data <- suppressMessages(read_csv(here("data", "raw_data", "20231011_Regated_TimepointDB_NonChronicvChronic.csv"))) %>% filter(., Selected_Final == 'Yes') %>%
  mutate(Time_At_Risk_Years = Time_At_Risk_Days/365) %>% mutate(Relapse_Rate = Nr_Relapses_PreLeisH / Time_At_Risk_Years) %>% relocate(., c("Time_At_Risk_Years", "Relapse_Rate"), .after = Time_At_Risk_Days)

full_data <- full_data %>% mutate(ICT_Selection = case_when(
  ICT_Live < 25 ~ "No",
  ICT_Live >= 25 ~ "Yes",
  is.na(ICT_Live) ~ "No"
),
EC_Selection = case_when(
  EC_Live < 25 ~ 'No',
  EC_Live >= 25 ~ 'Yes',
  is.na(EC_Live) ~ 'No'
))

##---------------------------------------------------------------
##                        Boxplot data D0                       -
##---------------------------------------------------------------

boxplot_data_D0 <- full_data %>% 
  mutate(RecodedTimepoint = case_when(
    RecodedTimepoint == 'W4' & Patient_SampleID != '0475' ~ 'EOT',
    RecodedTimepoint == 'W8' ~ 'EOT',
    TRUE ~ RecodedTimepoint)) %>%
  mutate(RecodedTimepoint = replace_na(RecodedTimepoint, 'D0')) %>% filter(., Group != 'HEC') %>%
  filter(., RecodedTimepoint %in% c('D0'))

boxplot_data_D0 <- boxplot_data_D0 %>% mutate(Aggregated_Group = case_when(Group %in% c('HIV', 'AL_HIV') ~ Group,
                                                                           Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV') ~ 'VL_HIV'))

boxplot_data_D0 <- boxplot_data_D0 %>% mutate_if(is.character, factor)

boxplot_data_D0$Group <- factor(boxplot_data_D0$Group, levels = c('HIV', 'AL_HIV', 'NonChronic_VLHIV', 'Chronic_VLHIV'))
boxplot_data_D0$Aggregated_Group <- factor(boxplot_data_D0$Aggregated_Group, levels = c('HIV', 'AL_HIV', 'VL_HIV'))

##----------------------------------------------------------------
##                        Boxplot data EOT                       -
##----------------------------------------------------------------

boxplot_data_EOT <- full_data %>% 
  mutate(RecodedTimepoint = case_when(
    RecodedTimepoint == 'W4' & Patient_SampleID != '0475' ~ 'EOT',
    RecodedTimepoint == 'W8' ~ 'EOT',
    TRUE ~ RecodedTimepoint)) %>% 
  filter(., RecodedTimepoint %in% c("EOT"))

boxplot_data_EOT <- boxplot_data_EOT %>% mutate_if(is.character, factor)
boxplot_data_EOT$Group <- factor(boxplot_data_EOT$Group, levels = c("NonChronic_VLHIV", "Chronic_VLHIV"))

##---------------------------------------------------------------
##                        Timeseries data                       -
##---------------------------------------------------------------

timeseries_data <- full_data %>% filter(., Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(., !str_starts(RecodedTimepoint, 'Pre_')) %>% mutate_if(is.character, factor)
timeseries_data$Group <- forcats::fct_drop(timeseries_data$Group)

timeseries_data$Group <- factor(timeseries_data$Group, levels = c('NonChronic_VLHIV', 'Chronic_VLHIV'))
timeseries_data$RecodedTimepoint <- factor(timeseries_data$RecodedTimepoint, levels = c('D0', 'W1', 'W2', 'W3', 'W4', 'W6', 'W8',
                                                                                        'Post_M3', 'Post_M4', 'Post_M5', 'Post_M6', 'Post_M7', 'Post_M8', 'Post_M9', 'Post_M10', 'Post_M11', 'Post_M12',
                                                                                        'Post_M13', 'Post_M14', 'Post_M16', 'Post_M18', 'Post_M19', 'Post_M22'))

levels(timeseries_data$RecodedTimepoint) <- gsub("Post_", "", levels(timeseries_data$RecodedTimepoint))
levels(timeseries_data$RecodedTimepoint) <- gsub("W", "TxW", levels(timeseries_data$RecodedTimepoint))

##################################################################
##                          Statistics                           #
##################################################################

markers <- full_data %>% select(starts_with(c('ICK', 'ICT', 'EC'))) %>% select(!contains(c('Selection', 'Availability'))) %>% colnames(.)

used_markers <- c("EC_CD8+", "EC_CD8-", #"ICK_CD4+", #composition
                  #"ICT_CD3+CD8+_IFNg", "ICT_CD3+CD8-_IFNg", "ICT_CD3+CD8+_CD107", "ICT_CD3+CD8-_CD107", #ICT functional prop
                  "ICT_CD3+CD8+_IFNg_MFI", "ICT_CD3+CD8-_IFNg_MFI",  "ICT_CD3+CD8+_CD107_MFI", #ICT functional MFI
                  "EC_CD8-_PD1+", "ICT_CD3+CD8-_TIGIT+", "EC_CD8-_LAG3+", "EC_CD8-_TIM3+", #CD8- Exhaustion
                  "EC_CD8+_PD1+", "ICT_CD3+CD8+_TIGIT+", "EC_CD8+_LAG3+", "EC_CD8+_TIM3+",  #CD8+ Exhaustion
                  "EC_CD8-_PD1+_MFI", "EC_CD8-_LAG3+_MFI", "EC_CD8-_TIM3+_MFI", "ICT_CD3+CD8-_TIGIT_MFI", #CD8- Exhaustion MFI
                  "EC_CD8+_PD1+_MFI", "EC_CD8+_LAG3+_MFI", "EC_CD8+_TIM3+_MFI", "ICT_CD3+CD8+_TIGIT_MFI", #CD8+ Exhaustion MFI
                  "EC_CD8-_PD1+LAG3+", "EC_CD8+_PD1+LAG3+", "EC_CD8-_PD1+TIM3+", "EC_CD8+_PD1+TIM3+", "EC_CD8-_TIM3+LAG3+", "EC_CD8+_TIM3+LAG3+", #EC Co-exhaustion
                  "ICT_CD3+CD8-_TIGIT+_IFNg+", "ICT_CD3+CD8-_TIGIT-_IFNg+", #"ICT_CD3+CD8-_TIGIT+_IFNg_MFI", "ICT_CD3+CD8-_TIGIT-_IFNg_MFI", #CD8- IFNG OUT OF TIGIT+ AND TIGIT- POPULATION
                  "ICT_CD3+CD8+_TIGIT+_IFNg+", "ICT_CD3+CD8+_TIGIT-_IFNg+", #"ICT_CD3+CD8+_TIGIT+_IFNg_MFI", "ICT_CD3+CD8+_TIGIT-_IFNg_MFI", #CD8+ IFNG OUT OF TIGIT+ AND TIGIT- POPULATION
                  "ICT_CD3+CD8-_TIGIT+IFNg-", "ICT_CD3+CD8-_TIGIT-IFNg+", #CD8- TIGIT CO-EXPRESSION
                  "ICT_CD3+CD8+_TIGIT+IFNg-", "ICT_CD3+CD8+_TIGIT-IFNg+", #"ICT_CD3+CD8+_TIGIT+CD107+", #CD8+ TIGIT CO-EXPRESSION
                  # "ICT_CD3+CD8-_TIGIT+IFNg+", "ICT_CD3+CD8-_TIGIT+IFNg-", "ICT_CD3+CD8-_TIGIT-IFNg+",  "ICT_CD3+CD8-_TIGIT-IFNg-", "ICT_CD3+CD8-_TIGIT+CD107+", #CD8- TIGIT CO-EXPRESSION
                  # "ICT_CD3+CD8+_TIGIT+IFNg+", "ICT_CD3+CD8+_TIGIT+IFNg-", "ICT_CD3+CD8+_TIGIT-IFNg+",  "ICT_CD3+CD8+_TIGIT-IFNg-", "ICT_CD3+CD8+_TIGIT+CD107+", #CD8+ TIGIT CO-EXPRESSION
                  "EC_CD8-_KLRG1+_MFI", "EC_CD8-_CD57+_MFI", "EC_CD8+_KLRG1+_MFI", "EC_CD8+_CD57+_MFI") #EC Senescence

##---------------------------------------------------------------
##                        Cross-sectional                       -
##---------------------------------------------------------------

MWU_csect <- sapply(used_markers[-c(28:35)], function(x) {
  myformula <- as.formula(paste0("`", x,"`", " ", "~ Aggregated_Group"))
  if (str_split(x, "_", simplify = T)[,1] == "ICK") {
    boxplot_data_D0 %>% filter(ICK_Selection == "Yes") %>% rstatix::wilcox_test(data = ., formula = myformula) %>% add_xy_position(x = "Aggregated_Group")
  } else if (str_split(x, "_", simplify = T)[,1] == "ICT") {
    boxplot_data_D0 %>% filter(ICT_Selection == "Yes") %>% rstatix::wilcox_test(data = ., formula = myformula) %>% add_xy_position(x = "Aggregated_Group")
  } else if (str_split(x, "_", simplify = T)[,1] == "EC") {
    boxplot_data_D0 %>% filter(EC_Selection == "Yes") %>% rstatix::wilcox_test(data = ., formula = myformula) %>% add_xy_position(x = "Aggregated_Group")
  }
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

MWU_between_primary_chronic_D0 <- sapply(used_markers, function(x) {
  if (str_split(x, "_", simplify = T)[,1] == "ICK") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_D0, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(ICK_Selection == "Yes")
  } else if (str_split(x, "_", simplify = T)[,1] == "ICT") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_D0, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(ICT_Selection == "Yes")
  } else if (str_split(x, "_", simplify = T)[,1] == "EC") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_D0, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(EC_Selection == "Yes")
  }
  boxplot_data_only_primary_chronic$Group <- forcats::fct_drop(boxplot_data_only_primary_chronic$Group)
  mwu_primary_chronic <- boxplot_data_only_primary_chronic %>% wilcox.test(unlist(.[x]) ~ Group, data = .)
  mwu_primary_chronic$p.value
}, simplify = FALSE, USE.NAMES = T)

MWU_between_primary_chronic_EOT <- sapply(used_markers, function(x) {
  if (str_split(x, "_", simplify = T)[,1] == "ICK") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_EOT, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(ICK_Selection == "Yes")
  } else if (str_split(x, "_", simplify = T)[,1] == "ICT") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_EOT, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(ICT_Selection == "Yes")
  } else if (str_split(x, "_", simplify = T)[,1] == "EC") {
    boxplot_data_only_primary_chronic <- filter(boxplot_data_EOT, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV')) %>% filter(EC_Selection == "Yes")
  }
  mwu_primary_chronic <- boxplot_data_only_primary_chronic %>% wilcox.test(unlist(.[x]) ~ Group, data = .)
  mwu_primary_chronic$p.value
}, simplify = FALSE, USE.NAMES = T)

##----------------------------------------------------------------
##                          Longitudinal                         -
##----------------------------------------------------------------

##### GROUP * TIME
longitudinal_mixed_models <- sapply(used_markers, function(x) {
  if (str_split(x, "_", simplify = T)[,1] == "ICK") {
    lmer_data <- timeseries_data %>% filter(ICK_Selection == "Yes")
    lmer_variable <- lmer_data %>% .[x] %>% pull
  } else if (str_split(x, "_", simplify = T)[,1] == "ICT") {
    lmer_data <- timeseries_data %>% filter(ICT_Selection == "Yes")
    lmer_variable <- lmer_data %>% .[x] %>% pull
  } else if (str_split(x, "_", simplify = T)[,1] == "EC") {
    lmer_data <- timeseries_data %>% filter(EC_Selection == "Yes")
    lmer_variable <- lmer_data %>% .[x] %>% pull
  }
  
  lmer_test <- lmer(lmer_variable ~ Group*RecodedTimepoint_in_Days + (1|Patient_SampleID), data = lmer_data)
}, simplify = FALSE, USE.NAMES = TRUE)

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

#################################################################
##                  Global plotting functions                   #
#################################################################

Cross_sectional_boxplot_generator <- function(marker_list, step_size_list) {
  
  boxplot_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    
    title <- names(marker_list)[[x]]
    plot_title <- paste0(sub("*(.*?) *[(].*", "\\1", title), " between groups") #regex to capture everything between % and (
    y_axis_title <- paste0(sub(".*[(] *(.*?) *[)].*", "\\1", title)) #regex to capture everything between ( and )
    
    if (str_split(tmp, "_", simplify = T)[,1] == "ICK") {
      temp_data <- boxplot_data_D0 %>% filter(ICK_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "ICT") {
      temp_data <- boxplot_data_D0 %>% filter(ICT_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "EC") {
      temp_data <- boxplot_data_D0 %>% filter(EC_Selection == "Yes")
    }
    temp_data$Aggregated_Group <- factor(temp_data$Aggregated_Group,
                                         levels = c('HIV', 'AL_HIV', 'VL_HIV'))
    
    my_comparisons <- list(c('HIV', 'AL_HIV'), c('HIV', 'VL_HIV'), c('AL_HIV', 'VL_HIV'))
    
    overall_boxplot <- temp_data %>% filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = Aggregated_Group, y = .data[[tmp]], fill = Aggregated_Group)) +
      geom_boxplot(width = .5, 
                   outlier.color = NA) +
      stat_pvalue_manual(MWU_csect[[tmp]], hide.ns=F, label="adjusted_p", label.size = 8, step.increase = step_size_list[x], tip.length = 0.025, inherit.aes = F) +
      geom_point(position = position_jitterdodge()) +
      ggtitle(label = plot_title) +
      labs(y = y_axis_title) +
      theme_bw() +
      scale_fill_manual(values = c("#00A087FF", "#3C5488FF", "#F39B7FFF", "#4DBBD5FF", "#E64B35FF")) +
      scale_x_discrete(expand = c(0,0))
    
    if (round(max(temp_data[[tmp]], na.rm = T)) <= 5) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 5, by = 1), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 5 && round(max(temp_data[[tmp]], na.rm = T)) <= 10) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 10, by = 2), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 10 && round(max(temp_data[[tmp]], na.rm = T)) <= 20) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 20, by = 4), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 20 && round(max(temp_data[[tmp]], na.rm = T)) <= 50) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 50, by = 10), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 50 && round(max(temp_data[[tmp]], na.rm = T)) <= 100) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 100, by = 20), expand = c(0,0))
    }
    
    overall_boxplot <- overall_boxplot +
      theme(axis.text.y = element_text(size = 24, family = text_font, color = 'black'),
            axis.title.y = element_text(size = 24, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 28, family = text_font, face = 'bold', hjust = 0.5),
            legend.position = 'none',
            panel.grid.major.x = element_blank(),
            plot.margin = margin(0,1,2,1, unit = 'cm'))
    
  })
  
  names(boxplot_list) <- names(marker_list)
  return(boxplot_list)
  
}


ncVL_vs_cVL_boxplot_generator <- function(marker_list, step_size_list) {
  
  boxplot_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    title <- names(marker_list)[[x]]
    plot_title <- paste0(sub("*(.*?) *[(].*", "\\1", title), " in VL-HIV at D0 and EOT") #regex to capture everything between % and (
    y_axis_title <- paste0(sub(".*[(] *(.*?) *[)].*", "\\1", title)) #regex to capture everything between ( and )
    
    if (str_split(tmp, "_", simplify = T)[,1] == "ICK") {
      temp_data_D0 <- boxplot_data_D0 %>% filter(ICK_Selection == "Yes")
      temp_data_EOT <- boxplot_data_EOT %>% filter(ICK_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "ICT") {
      temp_data_D0 <- boxplot_data_D0 %>% filter(ICT_Selection == "Yes")
      temp_data_EOT <- boxplot_data_EOT %>% filter(ICT_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "EC") {
      temp_data_D0 <- boxplot_data_D0 %>% filter(EC_Selection == "Yes")
      temp_data_EOT <- boxplot_data_EOT %>% filter(EC_Selection == "Yes")
    }
    
    temp_data_D0 <- filter(temp_data_D0, Group %in% c('NonChronic_VLHIV', 'Chronic_VLHIV'))
    temp_data_D0$Group <- forcats::fct_drop(temp_data_D0$Group)
    
    levels(temp_data_D0$Group) <- c("ncVL-HIV D0", "cVL-HIV D0")
    levels(temp_data_EOT$Group) <- c("ncVL-HIV EOT", "cVL-HIV EOT")
    
    temp_data <- bind_rows(temp_data_D0, temp_data_EOT)
    
    overall_boxplot <- temp_data %>% filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = Group, y = .data[[tmp]], fill = Group)) +
      geom_boxplot(width = .5, 
                   outlier.color = NA) +
      geom_signif(textsize = 8, comparisons = list(c('ncVL-HIV D0', 'cVL-HIV D0')), annotations = format(MWU_between_primary_chronic_D0[[tmp]], digits = 1, nsmall = 3)) +
      geom_signif(textsize = 8, comparisons = list(c('ncVL-HIV EOT', 'cVL-HIV EOT')), annotations = format(MWU_between_primary_chronic_EOT[[tmp]], digits = 1, nsmall = 3)) +
      geom_point(position = position_jitterdodge()) +
      ggtitle(label = plot_title) +
      labs(y = y_axis_title) +
      theme_bw() +
      scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF", "#4DBBD5FF", "#E64B35FF")) +
      scale_x_discrete(expand = c(0,0))
    
    if (round(max(temp_data[[tmp]], na.rm = T)) <= 5) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 5, by = 1), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 5 && round(max(temp_data[[tmp]], na.rm = T)) <= 10) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 10, by = 2), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 10 && round(max(temp_data[[tmp]], na.rm = T)) <= 20) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 20, by = 4), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 20 && round(max(temp_data[[tmp]], na.rm = T)) <= 50) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 50, by = 10), expand = c(0,0))
    } else if (round(max(temp_data[[tmp]], na.rm = T)) >= 50 && round(max(temp_data[[tmp]], na.rm = T)) <= 100) {
      overall_boxplot <- overall_boxplot + scale_y_continuous(breaks=seq(0, 100, by = 20), expand = c(0,0))
    }
    
    overall_boxplot <- overall_boxplot +
      theme(axis.text.y = element_text(size = 24, family = text_font, color = 'black'),
            axis.title.y = element_text(size = 24, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 28, family = text_font, face = 'bold', hjust = 0.5),
            legend.position = 'none',
            panel.grid.major.x = element_blank(),
            plot.margin = margin(0,1,2,1, unit = 'cm')) #+ ggpubr::rotate_x_text(angle = 40)
    
  })
  
  names(boxplot_list) <- names(marker_list)
  return(boxplot_list)
  
}

Longitudinal_timeseries_generator <- function(marker_list) {
  
  timeseries_marker_list <- lapply(seq_along(marker_list), function(x) {
    
    tmp <- marker_list[[x]]
    title <- names(marker_list)[[x]]
    ylab_simplified <- paste0(sub(".*[(] *(.*?) *[)].*", "\\1", title))
    
    if (str_split(tmp, "_", simplify = T)[,1] == "ICK") {
      temp_data <- timeseries_data %>% filter(ICK_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "ICT") {
      temp_data <- timeseries_data %>% filter(ICT_Selection == "Yes")
    } else if (str_split(tmp, "_", simplify = T)[,1] == "EC") {
      temp_data <- timeseries_data %>% filter(EC_Selection == "Yes")
    }
    
    levels(temp_data$Group) <- c("NonChronic VL-HIV", "Chronic VL-HIV")
    
    temp_data %>% 
      filter(., !is.na(.data[[tmp]])) %>%
      ggplot(aes(x = RecodedTimepoint, y = .data[[tmp]])) +
      geom_line(aes(group = Patient_SampleID, linetype = Group, colour = Group), lwd = 1.5) +
      geom_point(aes(colour = Group, shape = Relapse, size = Relapse)) +
      labs(y = ylab_simplified) +
      ggtitle(paste("Temporal dynamics of", sub("*(.*?) *[(].*", "\\1", title)),
              subtitle = paste0('\u0394 ', 'Group p = ', format(Longitudinal_Group_pvals[tmp], digits = 1, nsmall = 3))) +
      scale_color_manual(values = c("#4DBBD5FF", "#E64B35FF")) +
      scale_size_manual(values = c(4,8)) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 24, family = text_font, color = 'black'),
            axis.title.y = element_text(size = 24, family = text_font),
            axis.text.x = element_text(color='black', size=26, family=text_font),
            axis.title.x = element_blank(),
            plot.subtitle = element_text(size = 26, family = text_font, hjust = 0.5),
            plot.title = element_text(size = 28, family = text_font, face = 'bold', hjust = 0.5),
            legend.title = element_text(size = 26, family = text_font, face = 'bold'),
            legend.text = element_text(size = 26, family = text_font),
            legend.key.size = unit(2, "lines"),
            legend.position="bottom",
            plot.margin = margin(0,0,0,1, unit = 'cm')) + ggpubr::rotate_x_text() +
      guides(linetype = guide_legend(order = 1, title = NULL, override.aes = list(lwd = c(2,2), colour = c("#4DBBD5FF", "#E64B35FF"))),
             colour = 'none')
  })
  
  names(timeseries_marker_list) <- names(marker_list)
  return(timeseries_marker_list)
}

##################################################################
##                      Cellular composition                     #
##################################################################

composition_marker_list <- list('%CD8+ T cells (%CD3+CD8+)' = 'EC_CD8+',
                                '%CD8- T cells (%CD3+CD8-)' = 'EC_CD8-')

cross_sectional_composition <- Cross_sectional_boxplot_generator(composition_marker_list, step_size_list = c(0.175, 0.175))
ncVL_vs_cVL_composition <- ncVL_vs_cVL_boxplot_generator(composition_marker_list, step_size_list = c(0.25, 0.25))
longitudinal_composition <- Longitudinal_timeseries_generator(composition_marker_list)

layout <- "
AABBBCCCC
DDEEEFFFF
"

composition_arranged <- wrap_elements(full = cross_sectional_composition$`%CD8+ T cells (%CD3+CD8+)` + coord_cartesian(ylim = c(0, 135))) +
  wrap_elements(full = ncVL_vs_cVL_composition$`%CD8+ T cells (%CD3+CD8+)` + coord_cartesian(ylim = c(0, 100))) +
  wrap_elements(full = longitudinal_composition$`%CD8+ T cells (%CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_composition$`%CD8- T cells (%CD3+CD8-)` + coord_cartesian(ylim = c(0, 120))) +
  wrap_elements(full = ncVL_vs_cVL_composition$`%CD8- T cells (%CD3+CD8-)` + coord_cartesian(ylim = c(0, 90))) +
  wrap_elements(full = longitudinal_composition$`%CD8- T cells (%CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold', vjust = 2))

ggsave(here("analyses", "final_output", 'Figure2_Flowcyto_composition.pdf'), composition_arranged, dev = cairo_pdf, height = 15, width = 30)
ggsave(here("analyses", "final_output", 'Figure2_Flowcyto_composition.png'), composition_arranged, type = 'cairo-png', height = 15, width = 30)

##################################################################
##                      Mono-inhibition CD8+                     #
##################################################################

Mono_inhibition <- list('%PD1+ CD8+ T cells (%PD1+/CD3+CD8+)' = "EC_CD8+_PD1+", #0.175
                        'MFI of PD1 on CD8+ T cells (MFI PD1/CD3+CD8+)' = "EC_CD8+_PD1+_MFI", #0.1
                        '%TIGIT+ CD8+ T cells (%TIGIT+/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT+", #0.25
                        'MFI of TIGIT on CD8+ T cells (MFI TIGIT/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT_MFI") #0.175

cross_sectional_Mono_inhibition <- Cross_sectional_boxplot_generator(Mono_inhibition, step_size_list = c(0.175, 0.1, 0.25, 0.175))
ncVL_vs_cVL_Mono_inhibition <- ncVL_vs_cVL_boxplot_generator(Mono_inhibition, step_size_list = c(0.15, 0.15, 0.15, 0.15))
longitudinal_Mono_inhibition <- Longitudinal_timeseries_generator(Mono_inhibition)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Mono_inhibition <- wrap_elements(full = cross_sectional_Mono_inhibition$`%PD1+ CD8+ T cells (%PD1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,155))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition$`%PD1+ CD8+ T cells (%PD1+/CD3+CD8+)` + coord_cartesian(ylim = c(0,110))) +
  wrap_elements(full = longitudinal_Mono_inhibition$`%PD1+ CD8+ T cells (%PD1+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Mono_inhibition$`MFI of PD1 on CD8+ T cells (MFI PD1/CD3+CD8+)` + coord_cartesian(ylim = c(0,4000))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition$`MFI of PD1 on CD8+ T cells (MFI PD1/CD3+CD8+)` + coord_cartesian(ylim = c(0,2750))) +
  wrap_elements(full = longitudinal_Mono_inhibition$`MFI of PD1 on CD8+ T cells (MFI PD1/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Mono_inhibition$`%TIGIT+ CD8+ T cells (%TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,145))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition$`%TIGIT+ CD8+ T cells (%TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = longitudinal_Mono_inhibition$`%TIGIT+ CD8+ T cells (%TIGIT+/CD3+CD8+)`) + 
  wrap_elements(full = cross_sectional_Mono_inhibition$`MFI of TIGIT on CD8+ T cells (MFI TIGIT/CD3+CD8+)` + coord_cartesian(ylim = c(0,2250))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition$`MFI of TIGIT on CD8+ T cells (MFI TIGIT/CD3+CD8+)` + coord_cartesian(ylim = c(0,1750))) +
  wrap_elements(full = longitudinal_Mono_inhibition$`MFI of TIGIT on CD8+ T cells (MFI TIGIT/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "Figure4_Flowcyto_MonoInhibition.pdf"), arranged_Mono_inhibition, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "Figure4_Flowcyto_MonoInhibition.png"), arranged_Mono_inhibition, type = 'cairo-png', height = 22.5, width = 37.5)

#################################################################
##                    CD8+ TIGIT Functional                     #
#################################################################

CD8pos_TIGIT_functional <- list('%IFNg+TIGIT- CD8+ T cells (%IFNg+TIGIT-/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT-IFNg+",
                                '%IFNg-TIGIT+ CD8+ T cells (%IFNg-TIGIT+/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT+IFNg-",
                                '%IFNg+ out of TIGIT+ CD8+ T cells (%IFNg+/TIGIT+/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT+_IFNg+",
                                '%IFNg+ out of TIGIT- CD8+ T cells (%IFNg+/TIGIT-/CD3+CD8+)' = "ICT_CD3+CD8+_TIGIT-_IFNg+")

ncVL_vs_cVL_CD8pos_TIGIT_functional <- ncVL_vs_cVL_boxplot_generator(CD8pos_TIGIT_functional, step_size = c(0.05, 0.05, 0.05, 0.05))
longitudinal_CD8pos_TIGIT_functional <- Longitudinal_timeseries_generator(CD8pos_TIGIT_functional)

layout <- "
AABBB
CCDDD
EEFFF
GGHHH
"

arranged_CD8pos_TIGIT <- wrap_elements(full = ncVL_vs_cVL_CD8pos_TIGIT_functional$`%IFNg+TIGIT- CD8+ T cells (%IFNg+TIGIT-/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8pos_TIGIT_functional$`%IFNg+TIGIT- CD8+ T cells (%IFNg+TIGIT-/CD3+CD8+)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8pos_TIGIT_functional$`%IFNg-TIGIT+ CD8+ T cells (%IFNg-TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8pos_TIGIT_functional$`%IFNg-TIGIT+ CD8+ T cells (%IFNg-TIGIT+/CD3+CD8+)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8pos_TIGIT_functional$`%IFNg+ out of TIGIT+ CD8+ T cells (%IFNg+/TIGIT+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8pos_TIGIT_functional$`%IFNg+ out of TIGIT+ CD8+ T cells (%IFNg+/TIGIT+/CD3+CD8+)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8pos_TIGIT_functional$`%IFNg+ out of TIGIT- CD8+ T cells (%IFNg+/TIGIT-/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8pos_TIGIT_functional$`%IFNg+ out of TIGIT- CD8+ T cells (%IFNg+/TIGIT-/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "Figure5_Flowcyto_CD8Pos_TIGIT.pdf"), arranged_CD8pos_TIGIT, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "Figure5_Flowcyto_CD8Pos_TIGIT.png"), arranged_CD8pos_TIGIT, type = 'cairo-png', height = 22.5, width = 37.5)

#################################################################
##                        Supplementary                         #
#################################################################

##----------------------------------------------------------------
##                      Mono-inhibition CD8-                     -
##----------------------------------------------------------------

Mono_inhibition_CD8neg <- list('%PD1+ CD8- T cells (%PD1+/CD3+CD8-)' = "EC_CD8-_PD1+", #0.15
                               'MFI of PD1 on CD8- T cells (MFI PD1/CD3+CD8-)' = "EC_CD8-_PD1+_MFI", #0.075
                               '%TIGIT+ CD8- T cells (%TIGIT+/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT+", #0.12
                               'MFI of TIGIT on CD8- T cells (MFI TIGIT/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT_MFI") #0.175

cross_sectional_Mono_inhibition_CD8neg <- Cross_sectional_boxplot_generator(Mono_inhibition_CD8neg, step_size_list = c(0.15, 0.075, 0.12, 0.175))
ncVL_vs_cVL_Mono_inhibition_CD8neg <- ncVL_vs_cVL_boxplot_generator(Mono_inhibition_CD8neg, step_size_list = c(0.15, 0.15, 0.15, 0.15))
longitudinal_Mono_inhibition_CD8neg <- Longitudinal_timeseries_generator(Mono_inhibition_CD8neg)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Mono_inhibition_CD8neg <- wrap_elements(full = cross_sectional_Mono_inhibition_CD8neg$`%PD1+ CD8- T cells (%PD1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,147.5))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition_CD8neg$`%PD1+ CD8- T cells (%PD1+/CD3+CD8-)` + coord_cartesian(ylim = c(0,110))) +
  wrap_elements(full = longitudinal_Mono_inhibition_CD8neg$`%PD1+ CD8- T cells (%PD1+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Mono_inhibition_CD8neg$`MFI of PD1 on CD8- T cells (MFI PD1/CD3+CD8-)` + coord_cartesian(ylim = c(0,4500))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition_CD8neg$`MFI of PD1 on CD8- T cells (MFI PD1/CD3+CD8-)` + coord_cartesian(ylim = c(0,4250))) +
  wrap_elements(full = longitudinal_Mono_inhibition_CD8neg$`MFI of PD1 on CD8- T cells (MFI PD1/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Mono_inhibition_CD8neg$`%TIGIT+ CD8- T cells (%TIGIT+/CD3+CD8-)` + coord_cartesian(ylim = c(0,110))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition_CD8neg$`%TIGIT+ CD8- T cells (%TIGIT+/CD3+CD8-)` + coord_cartesian(ylim = c(0,90))) +
  wrap_elements(full = longitudinal_Mono_inhibition_CD8neg$`%TIGIT+ CD8- T cells (%TIGIT+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Mono_inhibition_CD8neg$`MFI of TIGIT on CD8- T cells (MFI TIGIT/CD3+CD8-)` + coord_cartesian(ylim = c(0,2250))) +
  wrap_elements(full = ncVL_vs_cVL_Mono_inhibition_CD8neg$`MFI of TIGIT on CD8- T cells (MFI TIGIT/CD3+CD8-)` + coord_cartesian(ylim = c(0,2000))) +
  wrap_elements(full = longitudinal_Mono_inhibition_CD8neg$`MFI of TIGIT on CD8- T cells (MFI TIGIT/CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output",  "supplemental", "SuppFigure03_Flowcyto_MonoInhibition_CD8neg.pdf"), arranged_Mono_inhibition_CD8neg, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output",  "supplemental", "SuppFigure03_Flowcyto_MonoInhibition_CD8neg.png"), arranged_Mono_inhibition_CD8neg, type = 'cairo-png', height = 22.5, width = 37.5)

##................................................................
##              Mono-inhibition Other Markers CD8+               .
##................................................................

Supp_Mono_inhibition_OM_CD8pos <- list('%LAG3+ CD8+ T cells (%LAG3+/CD3+CD8+)' = "EC_CD8+_LAG3+",
                                       'MFI of LAG3 on CD8+ T cells (MFI LAG3/CD3+CD8+)' = "EC_CD8+_LAG3+_MFI", #0.15
                                       '%TIM3+ CD8+ T cells (%TIM3+/CD3+CD8+)' = "EC_CD8+_TIM3+",
                                       'MFI of TIM3 on CD8+ T cells (MFI TIM3/CD3+CD8+)' = "EC_CD8+_TIM3+_MFI") #0.25

cross_sectional_Supp_Mono_inhibition_OM_CD8pos <- Cross_sectional_boxplot_generator(Supp_Mono_inhibition_OM_CD8pos, step_size = c(0.025, 0.15, 0.05, 0.25))
ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8pos <- ncVL_vs_cVL_boxplot_generator(Supp_Mono_inhibition_OM_CD8pos, step_size = c(0.15, 0.15, 0.15, 0.15))
longitudinal_Supp_Mono_inhibition_OM_CD8pos <- Longitudinal_timeseries_generator(Supp_Mono_inhibition_OM_CD8pos)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Supp_Mono_inhibition_CD8pos <- wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8pos$`%LAG3+ CD8+ T cells (%LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,120))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8pos$`%LAG3+ CD8+ T cells (%LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8pos$`%LAG3+ CD8+ T cells (%LAG3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8pos$`MFI of LAG3 on CD8+ T cells (MFI LAG3/CD3+CD8+)` + coord_cartesian(ylim = c(0,2250))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8pos$`MFI of LAG3 on CD8+ T cells (MFI LAG3/CD3+CD8+)` + coord_cartesian(ylim = c(0,2000))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8pos$`MFI of LAG3 on CD8+ T cells (MFI LAG3/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8pos$`%TIM3+ CD8+ T cells (%TIM3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,50))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8pos$`%TIM3+ CD8+ T cells (%TIM3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,40))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8pos$`%TIM3+ CD8+ T cells (%TIM3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8pos$`MFI of TIM3 on CD8+ T cells (MFI TIM3/CD3+CD8+)` + coord_cartesian(ylim = c(0,1500))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8pos$`MFI of TIM3 on CD8+ T cells (MFI TIM3/CD3+CD8+)` + coord_cartesian(ylim = c(0,1000))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8pos$`MFI of TIM3 on CD8+ T cells (MFI TIM3/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure04_Flowcyto_MonoInhibitionOM_CD8pos.pdf"), arranged_Supp_Mono_inhibition_CD8pos, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure04_Flowcyto_MonoInhibitionOM_CD8pos.png"), arranged_Supp_Mono_inhibition_CD8pos, type = 'cairo-png', height = 22.5, width = 37.5)

##................................................................
##              Mono-inhibition Other Markers CD8-               .
##................................................................

Supp_Mono_inhibition_OM_CD8neg <- list('%LAG3+ CD8- T cells (%LAG3+/CD3+CD8-)' = "EC_CD8-_LAG3+",
                                       'MFI of LAG3 on CD8- T cells (MFI LAG3/CD3+CD8-)' = "EC_CD8-_LAG3+_MFI",#0.15
                                       '%TIM3+ CD8- T cells (%TIM3+/CD3+CD8-)' = "EC_CD8-_TIM3+",
                                       'MFI of TIM3 on CD8- T cells (MFI TIM3/CD3+CD8-)' = "EC_CD8-_TIM3+_MFI") #0.25

cross_sectional_Supp_Mono_inhibition_OM_CD8neg <- Cross_sectional_boxplot_generator(Supp_Mono_inhibition_OM_CD8neg, step_size = c(0.025, 0.15, 0.075, 0.25))
ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8neg <- ncVL_vs_cVL_boxplot_generator(Supp_Mono_inhibition_OM_CD8neg, step_size = c(0.15, 0.15, 0.15, 0.15))
longitudinal_Supp_Mono_inhibition_OM_CD8neg <- Longitudinal_timeseries_generator(Supp_Mono_inhibition_OM_CD8neg)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Supp_Mono_inhibition_CD8neg <-
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8neg$`%LAG3+ CD8- T cells (%LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,120))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8neg$`%LAG3+ CD8- T cells (%LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,90))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8neg$`%LAG3+ CD8- T cells (%LAG3+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8neg$`MFI of LAG3 on CD8- T cells (MFI LAG3/CD3+CD8-)` + coord_cartesian(ylim = c(0,1500))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8neg$`MFI of LAG3 on CD8- T cells (MFI LAG3/CD3+CD8-)` + coord_cartesian(ylim = c(0,1500))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8neg$`MFI of LAG3 on CD8- T cells (MFI LAG3/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8neg$`%TIM3+ CD8- T cells (%TIM3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,30))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8neg$`%TIM3+ CD8- T cells (%TIM3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,22.5))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8neg$`%TIM3+ CD8- T cells (%TIM3+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Supp_Mono_inhibition_OM_CD8neg$`MFI of TIM3 on CD8- T cells (MFI TIM3/CD3+CD8-)` + coord_cartesian(ylim = c(0,1250))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Mono_inhibition_OM_CD8neg$`MFI of TIM3 on CD8- T cells (MFI TIM3/CD3+CD8-)` + coord_cartesian(ylim = c(0,1500))) +
  wrap_elements(full = longitudinal_Supp_Mono_inhibition_OM_CD8neg$`MFI of TIM3 on CD8- T cells (MFI TIM3/CD3+CD8-)`) +
  
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure05_Flowcyto_MonoInhibitionOM_CD8neg.pdf"), arranged_Supp_Mono_inhibition_CD8neg, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure05_Flowcyto_MonoInhibitionOM_CD8neg.png"), arranged_Supp_Mono_inhibition_CD8neg, type = 'cairo-png', height = 22.5, width = 37.5)

##---------------------------------------------------------------
##                        Poly-inhibition                       -
##---------------------------------------------------------------

Poly_inhibition <- list('%PD1+LAG3+ CD8+ T cells (PD1+LAG3+/CD3+CD8+)' = "EC_CD8+_PD1+LAG3+",
                        '%PD1+LAG3+ CD8- T cells (PD1+LAG3+/CD3+CD8-)' = "EC_CD8-_PD1+LAG3+",
                        '%PD1+TIM3+ CD8+ T cells (PD1+TIM3+/CD3+CD8+)' = "EC_CD8+_PD1+TIM3+",
                        '%PD1+TIM3+ CD8- T cells (PD1+TIM3+/CD3+CD8-)' = "EC_CD8-_PD1+TIM3+")

cross_sectional_Poly_inhibition <- Cross_sectional_boxplot_generator(Poly_inhibition, step_size_list = c(0.025, 0.025, 0.025, 0.01))
ncVL_vs_cVL_Poly_inhibition <- ncVL_vs_cVL_boxplot_generator(Poly_inhibition, step_size_list = c(0.1, 0.1, 0.1, 0.1))
longitudinal_Poly_inhibition <- Longitudinal_timeseries_generator(Poly_inhibition)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Poly_inhibition <- wrap_elements(full = cross_sectional_Poly_inhibition$`%PD1+LAG3+ CD8+ T cells (PD1+LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = ncVL_vs_cVL_Poly_inhibition$`%PD1+LAG3+ CD8+ T cells (PD1+LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,70))) +
  wrap_elements(full = longitudinal_Poly_inhibition$`%PD1+LAG3+ CD8+ T cells (PD1+LAG3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Poly_inhibition$`%PD1+LAG3+ CD8- T cells (PD1+LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,80))) +
  wrap_elements(full = ncVL_vs_cVL_Poly_inhibition$`%PD1+LAG3+ CD8- T cells (PD1+LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,60))) +
  wrap_elements(full = longitudinal_Poly_inhibition$`%PD1+LAG3+ CD8- T cells (PD1+LAG3+/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_Poly_inhibition$`%PD1+TIM3+ CD8+ T cells (PD1+TIM3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,27.5))) +
  wrap_elements(full = ncVL_vs_cVL_Poly_inhibition$`%PD1+TIM3+ CD8+ T cells (PD1+TIM3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,20))) +
  wrap_elements(full = longitudinal_Poly_inhibition$`%PD1+TIM3+ CD8+ T cells (PD1+TIM3+/CD3+CD8+)`) + 
  wrap_elements(full = cross_sectional_Poly_inhibition$`%PD1+TIM3+ CD8- T cells (PD1+TIM3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,30))) +
  wrap_elements(full = ncVL_vs_cVL_Poly_inhibition$`%PD1+TIM3+ CD8- T cells (PD1+TIM3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,22.5))) +
  wrap_elements(full = longitudinal_Poly_inhibition$`%PD1+TIM3+ CD8- T cells (PD1+TIM3+/CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure06_Flowcyto_PolyInhibition.pdf"), arranged_Poly_inhibition, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure06_Flowcyto_PolyInhibition.png"), arranged_Poly_inhibition, type = 'cairo-png', height = 22.5, width = 37.5)

##...............................................................
##                Poly-inhibition Other Markers                 .
##...............................................................

Supp_Poly_inhibition <- list('%TIM3+LAG3+ CD8+ T cells (TIM3+LAG3+/CD3+CD8+)' = "EC_CD8+_TIM3+LAG3+",
                             '%TIM3+LAG3+ CD8- T cells (TIM3+LAG3+/CD3+CD8-)' = "EC_CD8-_TIM3+LAG3+")

cross_sectional_Supp_Poly_inhibition <- Cross_sectional_boxplot_generator(Supp_Poly_inhibition, step_size = c(0.001, 0.005))
ncVL_vs_cVL_Supp_Poly_inhibition <- ncVL_vs_cVL_boxplot_generator(Supp_Poly_inhibition, step_size = c(0.1, 0.1))
longitudinal_Supp_Poly_inhibition <- Longitudinal_timeseries_generator(Supp_Poly_inhibition)

layout <- "
AABBBCCCC
DDEEEFFFF
"

arranged_Supp_Poly_inhibition <- wrap_elements(full = cross_sectional_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8+ T cells (TIM3+LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,45))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8+ T cells (TIM3+LAG3+/CD3+CD8+)` + coord_cartesian(ylim = c(0,35))) +
  wrap_elements(full = longitudinal_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8+ T cells (TIM3+LAG3+/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8- T cells (TIM3+LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,22.5))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8- T cells (TIM3+LAG3+/CD3+CD8-)` + coord_cartesian(ylim = c(0,20))) +
  wrap_elements(full = longitudinal_Supp_Poly_inhibition$`%TIM3+LAG3+ CD8- T cells (TIM3+LAG3+/CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure07_Flowcyto_PolyInhibitionOM.pdf"), arranged_Supp_Poly_inhibition, dev = cairo_pdf, height = 15, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure07_Flowcyto_PolyInhibitionOM.png"), arranged_Supp_Poly_inhibition, type = 'cairo-png', height = 15, width = 37.5)

##----------------------------------------------------------------
##                          Senescence                           -
##----------------------------------------------------------------

Supp_list_Senescence <- c('MFI of KLRG1 on CD8+ T cells (MFI KLRG1/CD3+CD8+)' = "EC_CD8+_KLRG1+_MFI",
                          'MFI of KLRG1 on CD8- T cells (MFI KLRG1/CD3+CD8-)' = "EC_CD8-_KLRG1+_MFI",
                          'MFI of CD57 on CD8+ T cells (MFI CD57/CD3+CD8+)' = "EC_CD8+_CD57+_MFI",
                          'MFI of CD57 on CD8- T cells (MFI CD57/CD3+CD8-)' = "EC_CD8-_CD57+_MFI")

cross_sectional_Supp_Senescence <- Cross_sectional_boxplot_generator(Supp_list_Senescence, step_size_list = c(0.1, 0.1, 0.1, 0.15))
ncVL_vs_cVL_Supp_Senescence <- ncVL_vs_cVL_boxplot_generator(Supp_list_Senescence, step_size_list = c(0.1, 0.1, 0.1, 0.1))
longitudinal_Supp_Senescence <- Longitudinal_timeseries_generator(Supp_list_Senescence)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
JJKKKLLLL
"

arranged_Supp_Senescence <- 
  wrap_elements(full = cross_sectional_Supp_Senescence$`MFI of KLRG1 on CD8+ T cells (MFI KLRG1/CD3+CD8+)` + coord_cartesian(ylim = c(0,25000))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Senescence$`MFI of KLRG1 on CD8+ T cells (MFI KLRG1/CD3+CD8+)` + coord_cartesian(ylim = c(0,17500))) +
  wrap_elements(full = longitudinal_Supp_Senescence$`MFI of KLRG1 on CD8+ T cells (MFI KLRG1/CD3+CD8+)`) + 
  wrap_elements(full = cross_sectional_Supp_Senescence$`MFI of KLRG1 on CD8- T cells (MFI KLRG1/CD3+CD8-)` + coord_cartesian(ylim = c(0,20000))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Senescence$`MFI of KLRG1 on CD8- T cells (MFI KLRG1/CD3+CD8-)` + coord_cartesian(ylim = c(0,15000))) +
  wrap_elements(full = longitudinal_Supp_Senescence$`MFI of KLRG1 on CD8- T cells (MFI KLRG1/CD3+CD8-)`) + 
  wrap_elements(full = cross_sectional_Supp_Senescence$`MFI of CD57 on CD8+ T cells (MFI CD57/CD3+CD8+)` + coord_cartesian(ylim = c(0,65000))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_Senescence$`MFI of CD57 on CD8+ T cells (MFI CD57/CD3+CD8+)` + coord_cartesian(ylim = c(0,20000))) +
  wrap_elements(full = longitudinal_Supp_Senescence$`MFI of CD57 on CD8+ T cells (MFI CD57/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_Supp_Senescence$`MFI of CD57 on CD8- T cells (MFI CD57/CD3+CD8-)` + coord_cartesian(ylim = c(0,20000))) + 
  wrap_elements(full = ncVL_vs_cVL_Supp_Senescence$`MFI of CD57 on CD8- T cells (MFI CD57/CD3+CD8-)` + coord_cartesian(ylim = c(0,15000))) + 
  wrap_elements(full = longitudinal_Supp_Senescence$`MFI of CD57 on CD8- T cells (MFI CD57/CD3+CD8-)`) + 
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure08_Flowcyto_Senescence.pdf"), arranged_Supp_Senescence, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure08_Flowcyto_Senescence.png"), arranged_Supp_Senescence, type = 'cairo-png', height = 22.5, width = 37.5)

##----------------------------------------------------------------
##                        ICT_Functional                         -
##----------------------------------------------------------------

Supp_ICT_functional <- list('MFI of IFNg on CD8+ T cells (MFI IFNg/CD3+CD8+)' = "ICT_CD3+CD8+_IFNg_MFI",
                            'MFI of IFNg on CD8- T cells (MFI IFNg/CD3+CD8-)' = "ICT_CD3+CD8-_IFNg_MFI",
                            'MFI of CD107 on CD8+ T cells (MFI CD107/CD3+CD8+)' = "ICT_CD3+CD8+_CD107_MFI")

cross_sectional_ICT_functional <- Cross_sectional_boxplot_generator(Supp_ICT_functional, step_size_list = c(0.25, 0.15, 0.175))
ncVL_vs_cVL_ICT_functional <- ncVL_vs_cVL_boxplot_generator(Supp_ICT_functional, step_size_list = c(0.1, 0.1, 0.1))
longitudinal_ICT_functional <- Longitudinal_timeseries_generator(Supp_ICT_functional)

layout <- "
AABBBCCCC
DDEEEFFFF
GGHHHIIII
"

arranged_Supp_ICT_functional <- wrap_elements(full = cross_sectional_ICT_functional$`MFI of IFNg on CD8+ T cells (MFI IFNg/CD3+CD8+)` + coord_cartesian(ylim = c(0,12500))) +
  wrap_elements(full = ncVL_vs_cVL_ICT_functional$`MFI of IFNg on CD8+ T cells (MFI IFNg/CD3+CD8+)` + coord_cartesian(ylim = c(0,8000))) +
  wrap_elements(full = longitudinal_ICT_functional$`MFI of IFNg on CD8+ T cells (MFI IFNg/CD3+CD8+)`) +
  wrap_elements(full = cross_sectional_ICT_functional$`MFI of IFNg on CD8- T cells (MFI IFNg/CD3+CD8-)` + coord_cartesian(ylim = c(0,22500))) +
  wrap_elements(full = ncVL_vs_cVL_ICT_functional$`MFI of IFNg on CD8- T cells (MFI IFNg/CD3+CD8-)` + coord_cartesian(ylim = c(0,17500))) +
  wrap_elements(full = longitudinal_ICT_functional$`MFI of IFNg on CD8- T cells (MFI IFNg/CD3+CD8-)`) +
  wrap_elements(full = cross_sectional_ICT_functional$`MFI of CD107 on CD8+ T cells (MFI CD107/CD3+CD8+)` + coord_cartesian(ylim = c(0,8000))) +
  wrap_elements(full = ncVL_vs_cVL_ICT_functional$`MFI of CD107 on CD8+ T cells (MFI CD107/CD3+CD8+)` + coord_cartesian(ylim = c(0,5500))) +
  wrap_elements(full = longitudinal_ICT_functional$`MFI of CD107 on CD8+ T cells (MFI CD107/CD3+CD8+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold', vjust = 2))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure09_Flowcyto_ICTFunc.pdf"), arranged_Supp_ICT_functional, dev = cairo_pdf, height = 15, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure09_Flowcyto_ICTFunc.png"), arranged_Supp_ICT_functional, type = 'cairo-png', height = 15, width = 37.5)

##---------------------------------------------------------------
##                    CD8- TIGIT Functional                     -
##---------------------------------------------------------------

CD8neg_TIGIT_functional <- list('%IFNg+TIGIT- CD8- T cells (%IFNg+TIGIT-/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT-IFNg+",
                                '%IFNg-TIGIT+ CD8- T cells (%IFNg-TIGIT+/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT+IFNg-",
                                '%IFNg+ out of TIGIT+ CD8- T cells (%IFNg+/TIGIT+/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT+_IFNg+",
                                '%IFNg+ out of TIGIT- CD8- T cells (%IFNg+/TIGIT-/CD3+CD8-)' = "ICT_CD3+CD8-_TIGIT-_IFNg+")

ncVL_vs_cVL_CD8neg_TIGIT_functional <- ncVL_vs_cVL_boxplot_generator(CD8neg_TIGIT_functional, step_size = c(0.05, 0.05, 0.05, 0.05))
longitudinal_CD8neg_TIGIT_functional <- Longitudinal_timeseries_generator(CD8neg_TIGIT_functional)

layout <- "
AABBB
CCDDD
EEFFF
GGHHH
"

arranged_CD8neg_TIGIT <- wrap_elements(full = ncVL_vs_cVL_CD8neg_TIGIT_functional$`%IFNg+TIGIT- CD8- T cells (%IFNg+TIGIT-/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8neg_TIGIT_functional$`%IFNg+TIGIT- CD8- T cells (%IFNg+TIGIT-/CD3+CD8-)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8neg_TIGIT_functional$`%IFNg-TIGIT+ CD8- T cells (%IFNg-TIGIT+/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8neg_TIGIT_functional$`%IFNg-TIGIT+ CD8- T cells (%IFNg-TIGIT+/CD3+CD8-)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8neg_TIGIT_functional$`%IFNg+ out of TIGIT+ CD8- T cells (%IFNg+/TIGIT+/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8neg_TIGIT_functional$`%IFNg+ out of TIGIT+ CD8- T cells (%IFNg+/TIGIT+/CD3+CD8-)`) +
  wrap_elements(full = ncVL_vs_cVL_CD8neg_TIGIT_functional$`%IFNg+ out of TIGIT- CD8- T cells (%IFNg+/TIGIT-/CD3+CD8-)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full =  longitudinal_CD8neg_TIGIT_functional$`%IFNg+ out of TIGIT- CD8- T cells (%IFNg+/TIGIT-/CD3+CD8-)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure10_Flowcyto_CD8Neg_TIGIT.pdf"), arranged_CD8neg_TIGIT, dev = cairo_pdf, height = 22.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure10_Flowcyto_CD8Neg_TIGIT.png"), arranged_CD8neg_TIGIT, type = 'cairo-png', height = 22.5, width = 37.5)


##---------------------------------------------------------------
##                        ICK CD4 Proxy                         -
##---------------------------------------------------------------

Supp_CD4_Proxy <- list('%CD4+ T cells (%CD3+CD4+)' = "ICK_CD4+")

cross_sectional_Supp_CD4_Proxy <- Cross_sectional_boxplot_generator(Supp_CD4_Proxy, step_size = c(0.05))
ncVL_vs_cVL_Supp_CD4_Proxy <- ncVL_vs_cVL_boxplot_generator(Supp_CD4_Proxy, step_size = c(0.1))
longitudinal_Supp_CD4_Proxy <- Longitudinal_timeseries_generator(Supp_CD4_Proxy)

layout <- "
AABBBCCCC
"

arranged_Supp_Poly_CD4_Proxy <- wrap_elements(full = cross_sectional_Supp_CD4_Proxy$`%CD4+ T cells (%CD3+CD4+)` + coord_cartesian(ylim = c(0,100))) +
  wrap_elements(full = ncVL_vs_cVL_Supp_CD4_Proxy$`%CD4+ T cells (%CD3+CD4+)` + coord_cartesian(ylim = c(0,80))) +
  wrap_elements(full = longitudinal_Supp_CD4_Proxy$`%CD4+ T cells (%CD3+CD4+)`) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure15_Flowcyto_CD4Proxy.pdf"), arranged_Supp_Poly_CD4_Proxy, dev = cairo_pdf, height = 7.5, width = 37.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure15_Flowcyto_CD4Proxy.png"), arranged_Supp_Poly_CD4_Proxy, type = 'cairo-png', height = 7.5, width = 37.5)

