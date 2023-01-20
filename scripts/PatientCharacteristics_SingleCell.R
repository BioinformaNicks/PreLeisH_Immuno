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
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))

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

characteristics_data <- suppressMessages(read_csv(here("data", "meta_data", "Patient_Characteristics", "20220823_SingleCell_PrimaryVsChronicCharacteristics.csv")))

characteristics_data$Patient_ID <- as.factor(characteristics_data$Patient_ID)
characteristics_data$Group <- as.factor(characteristics_data$Group)
characteristics_data$Age <- as.numeric(characteristics_data$Age)
characteristics_data$Sex <- as.factor(characteristics_data$Sex)
#changing occupational data to fewer categories without modifying raw data
characteristics_data <- characteristics_data %>% mutate(Occupation = case_when(Occupation %in% c("Sex worker", "Unemployed") ~ "Other",
                                                                               TRUE ~ Occupation))
  
characteristics_data$Occupation <- as.factor(characteristics_data$Occupation)
characteristics_data$BMI <- as.numeric(characteristics_data$BMI)
characteristics_data$Treatment_Regimen <- as.factor(characteristics_data$Treatment_Regimen)
characteristics_data$VL_History <- as.factor(characteristics_data$VL_History)
characteristics_data$Primary_VL <- as.factor(characteristics_data$Primary_VL)
characteristics_data$Months_From_Prev_VL <- as.numeric(characteristics_data$Months_From_Prev_VL)
characteristics_data$CD4_Baseline <- as.numeric(characteristics_data$CD4_Baseline)
characteristics_data$CD4_d0 <- as.numeric(characteristics_data$CD4_d0)
characteristics_data$CD4_Post_M6 <- as.numeric(characteristics_data$CD4_Post_M6)
characteristics_data$CD4_Post_M12 <- as.numeric(characteristics_data$CD4_Post_M12)
characteristics_data$Microscope_Confirmed <- as.factor(characteristics_data$Microscope_Confirmed)
characteristics_data$ART_Baseline <- as.factor(characteristics_data$ART_Baseline)
characteristics_data$VL_Episodes_Baseline <- as.numeric(characteristics_data$VL_Episodes_Baseline)
characteristics_data$Concomitant_Disease <- as.factor(characteristics_data$Concomitant_Disease)

#Create VL-HIV group consisting of Primary + Chronic
characteristics_data <- characteristics_data %>% mutate(Aggregated_Group = ifelse(characteristics_data$Group %in% c('Chronic_VL_HIV', 'Primary_VL_HIV'), 'Symptomatic_VL_HIV', as.character(characteristics_data$Group))) %>%
  relocate(Aggregated_Group, .after = Group)
characteristics_data$Aggregated_Group <- as.factor(characteristics_data$Aggregated_Group)

# HIV and AL-HIV baseline should be compared across D0 VL-HIV, so mutate some baseline to D0 variables
characteristics_data <- characteristics_data %>% mutate(CD4_d0 = case_when(Group %in% c("HIV", "Asymptomatic_HIV") ~ CD4_Baseline,
                                                                           TRUE ~ CD4_d0))
characteristics_data <- characteristics_data %>% mutate(Lymphocytes_D0 = case_when(Group %in% c("HIV", "Asymptomatic_HIV") ~ Lymphocytes_Baseline,
                                                                                   TRUE ~ Lymphocytes_D0))
characteristics_data <- characteristics_data %>% mutate(Platelet_Count_d0 = case_when(Group %in% c("HIV", "Asymptomatic_HIV") ~ Platelet_Count_Baseline,
                                                                                      TRUE ~ Platelet_Count_d0))
characteristics_data <- characteristics_data %>% mutate(Hemoglobin_Count_d0 = case_when(Group %in% c("HIV", "Asymptomatic_HIV") ~ Hemoglobin_Count_Baseline,
                                                                                        TRUE ~ Hemoglobin_Count_d0))

#Create primary vs chronic dataframe for primary vs chronic testing

primary_v_chronic_data <- characteristics_data %>% filter(Group %in% c('Primary_VL_HIV', 'Chronic_VL_HIV'))
primary_v_chronic_data$Group <- forcats::fct_drop(primary_v_chronic_data$Group)

##################################################################
##                        Total Socio-Demo                       #
##################################################################

Total_Age <- summary(characteristics_data$Age)
Total_Male <- characteristics_data %>% filter(Sex == 'Male') %>% nrow()
Total_BMI <- summary(characteristics_data$BMI)
Total_Literate <- paste0(sum(table(characteristics_data$Literate, characteristics_data$Aggregated_Group)['Yes',]),
                         " (",
                         (sum(table(characteristics_data$Literate, characteristics_data$Aggregated_Group)['Yes',] / sum(table(characteristics_data$Literate, characteristics_data$Aggregated_Group))))*100, ")")
Total_Occupation <- cbind(rowSums(table(characteristics_data$Occupation, characteristics_data$Aggregated_Group)),
                          rowSums(table(characteristics_data$Occupation, characteristics_data$Aggregated_Group))/sum(table(characteristics_data$Occupation, characteristics_data$Aggregated_Group))*100)

#################################################################
##                    All groups Socio-Demo                     #
#################################################################

all_groups <- c("HIV", "Asymptomatic_HIV", "Symptomatic_VL_HIV")

##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

All_groups_Age <- sapply(all_groups, function(x) {
  summary_age <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Age) %>% unlist(.) %>% summary(.)
  paste0(round(summary_age['Median'],1), ' (', round(summary_age['1st Qu.'],1), '-', round(summary_age['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Male <- sapply(all_groups, function(x) {
  num_male <- characteristics_data %>% filter(Aggregated_Group == x) %>% filter(Sex == 'Male') %>% nrow()
  total_patients <- characteristics_data %>% filter(Aggregated_Group == x) %>% nrow()
  paste0(num_male, ' (', round((num_male / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_BMI <- sapply(all_groups, function(x) {
  summary_BMI <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(BMI) %>% unlist(.) %>% summary(.)
  paste0(round(summary_BMI['Median'],1), ' (', round(summary_BMI['1st Qu.'],1), '-', round(summary_BMI['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Literate <- cbind(table(characteristics_data$Literate, characteristics_data$Aggregated_Group)['Yes',], 
                             (table(characteristics_data$Literate, characteristics_data$Aggregated_Group)['Yes',] / colSums(table(characteristics_data$Literate, characteristics_data$Aggregated_Group)) * 100))

All_groups_Occupation <- table(characteristics_data$Occupation, characteristics_data$Aggregated_Group)

#################################################################
##                Primary vs Chronic Socio-Demo                 #
#################################################################

##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

primary_v_chronic_Age <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_age <- primary_v_chronic_data %>% filter(Group == x) %>% select(Age) %>% unlist(.) %>% summary(.)
  paste0(round(summary_age['Median'],1), ' (', round(summary_age['1st Qu.'],1), '-', round(summary_age['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Male <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  num_male <- primary_v_chronic_data %>% filter(Group == x) %>% filter(Sex == 'Male') %>% nrow()
  total_patients <- primary_v_chronic_data %>% filter(Group == x) %>% nrow()
  paste0(num_male, ' (', round((num_male / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_BMI <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_BMI <- primary_v_chronic_data %>% filter(Group == x) %>% select(BMI) %>% unlist(.) %>% summary(.)
  paste0(round(summary_BMI['Median'],1), ' (', round(summary_BMI['1st Qu.'],1), '-', round(summary_BMI['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

##---------------------------------------------------------------
##                      Tables for % data                       -
##---------------------------------------------------------------

primary_v_chronic_literate <- cbind(table(primary_v_chronic_data$Literate, primary_v_chronic_data$Group)['Yes',], 
      (table(primary_v_chronic_data$Literate, primary_v_chronic_data$Group)['Yes',] / colSums(table(primary_v_chronic_data$Literate, primary_v_chronic_data$Group)) * 100))

primary_v_chronic_occupation <- table(primary_v_chronic_data$Occupation, primary_v_chronic_data$Group)

##################################################################
##                        Total Clinical                         #
##################################################################

Total_VL_History <- table(characteristics_data$VL_History)
Total_Past_VL_Episodes <- summary(characteristics_data$VL_Episodes_Baseline)
Total_Months_Since_Past_VL <- summary(characteristics_data$Months_From_Prev_VL)
Total_Microscope_Confirmed <- table(characteristics_data$Microscope_Confirmed)
Total_Parasite_Grading <- table(characteristics_data$Parasite_Grading)
Total_VL_Treatment <- table(characteristics_data$Treatment_Regimen)
Total_ART <- table(characteristics_data$ART_Baseline)
Total_Concomitant_Disease <- table(characteristics_data$Concomitant_Disease)

Total_CD4_D0 <- summary(characteristics_data$CD4_d0)
Total_Lymphocytes_D0 <- summary(characteristics_data$Lymphocytes_D0)
Total_Platelet_D0 <- summary(characteristics_data$Platelet_Count_d0)
Total_HB_D0 <- summary(characteristics_data$Hemoglobin_Count_d0)

Total_Lymphocytes_EOT <- summary(characteristics_data$Lymphocytes_EOT)
Total_Platelet_EOT <- summary(characteristics_data$Platelet_Count_EOT)
Total_HB_EOT <- summary(characteristics_data$Hemoglobin_Count_EOT)

Total_CD4_Post_M6 <- summary(characteristics_data$CD4_Post_M6)
Total_Lymphocytes_Post_M6 <- summary(characteristics_data$Lymphocytes_Post_M6)
Total_Platelet_Post_M6 <- summary(characteristics_data$Platelet_Count_Post_M6)
Total_HB_Post_M6 <- summary(characteristics_data$Hemoglobin_Count_Post_M6)

#################################################################
##                      All groups Clinical                     #
#################################################################

##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

All_groups_VL_History <- sapply(all_groups, function(x) {
  num_VL_History <- characteristics_data %>% filter(Aggregated_Group == x) %>% filter(VL_History == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(Aggregated_Group == x) %>% nrow()
  paste0(num_VL_History, ' (', round((num_VL_History / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_episodes<- sapply(all_groups, function(x) {
  summary_episodes <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(VL_Episodes_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_episodes['Median'],1), ' (', round(summary_episodes['1st Qu.'],1), '-', round(summary_episodes['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Months_Prev_VL <- sapply(all_groups, function(x) {
  summary_Months_Prev_VL <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Months_From_Prev_VL) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Months_Prev_VL['Median'],1), ' (', round(summary_Months_Prev_VL['1st Qu.'],1), '-', round(summary_Months_Prev_VL['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Parasite_Grading <- table(characteristics_data$Parasite_Grading, characteristics_data$Aggregated_Group)

All_groups_ART <- sapply(all_groups, function(x) {
  num_ART <- characteristics_data %>% filter(Aggregated_Group == x) %>% filter(ART_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(Aggregated_Group == x) %>% nrow()
  paste0(num_ART, ' (', round((num_ART / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Concomitant_Diseases <- sapply(all_groups, function(x) {
  num_Concomitant_Disease <- characteristics_data %>% filter(Aggregated_Group == x) %>% filter(Concomitant_Disease == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(Aggregated_Group == x) %>% nrow()
  paste0(num_Concomitant_Disease, ' (', round((num_Concomitant_Disease / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)


All_groups_CD4_D0 <- sapply(all_groups, function(x) {
  summary_CD4 <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(CD4_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_CD4['Median'],0), ' (', round(summary_CD4['1st Qu.'],0), '-', round(summary_CD4['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Lymphocytes_D0 <- sapply(all_groups, function(x) {
  summary_Lymphocytes <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Lymphocytes_D0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Platelets_D0 <- sapply(all_groups, function(x) {
  summary_Platelets <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Platelet_Count_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Hemoglobin_D0 <- sapply(all_groups, function(x) {
  summary_Hemoglobin <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Hemoglobin_Count_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Lymphocytes_EOT <- sapply(all_groups, function(x) {
  summary_Lymphocytes <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Lymphocytes_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Platelets_EOT <- sapply(all_groups, function(x) {
  summary_Platelets <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Platelet_Count_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Hemoglobin_EOT <- sapply(all_groups, function(x) {
  summary_Hemoglobin <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Hemoglobin_Count_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_CD4_Post_M6 <- sapply(all_groups, function(x) {
  summary_CD4 <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(CD4_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_CD4['Median'],0), ' (', round(summary_CD4['1st Qu.'],0), '-', round(summary_CD4['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Lymphocytes_Post_M6 <- sapply(all_groups, function(x) {
  summary_Lymphocytes <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Lymphocytes_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Platelets_Post_M6 <- sapply(all_groups, function(x) {
  summary_Platelets <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Platelet_Count_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Hemoglobin_Post_M6 <- sapply(all_groups, function(x) {
  summary_Hemoglobin <- characteristics_data %>% filter(Aggregated_Group == x) %>% select(Hemoglobin_Count_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)


#################################################################
##                  Primary vs Chronic Clinical                 #
#################################################################

##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

primary_v_chronic_VL_History <- cbind(table(primary_v_chronic_data$VL_History, primary_v_chronic_data$Group)['Yes',], 
                                     (table(primary_v_chronic_data$VL_History, primary_v_chronic_data$Group)['Yes',] / colSums(table(primary_v_chronic_data$VL_History, primary_v_chronic_data$Group)) * 100))

primary_v_chronic_episodes <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_episodes <- primary_v_chronic_data %>% filter(Group == x) %>% select(VL_Episodes_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_episodes['Median'],1), ' (', round(summary_episodes['1st Qu.'],1), '-', round(summary_episodes['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Months_Prev_VL <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Months_Prev_VL <- primary_v_chronic_data %>% filter(Group == x) %>% select(Months_From_Prev_VL) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Months_Prev_VL['Median'],1), ' (', round(summary_Months_Prev_VL['1st Qu.'],1), '-', round(summary_Months_Prev_VL['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Microscope_Confirmed <- cbind(table(primary_v_chronic_data$Microscope_Confirmed, primary_v_chronic_data$Group)['Yes',], 
                                                (table(primary_v_chronic_data$Microscope_Confirmed, primary_v_chronic_data$Group)['Yes',] / colSums(table(primary_v_chronic_data$Microscope_Confirmed, primary_v_chronic_data$Group)) * 100))

primary_v_chronic_data <- primary_v_chronic_data %>% mutate(Parasite_Grading = case_when(Parasite_Grading %in% c('Spleen +1', 'Spleen +2', 'Spleen +3', 'Bone Marrow +1',
                                                                                                                 'Bone Marrow +2', 'Bone Marrow +3') ~ '+1 to +3',
                                                                                         Parasite_Grading %in% c('Spleen +4', 'Spleen +5', 'Spleen +6',
                                                                                                                 'Bone Marrow +4', 'Bone Marrow +5', 'Bone Marrow +6') ~ '+4 to +6'))
primary_v_chronic_Parasite_Grading <- table(primary_v_chronic_data$Parasite_Grading, primary_v_chronic_data$Group)

primary_v_chronic_treatment_regimen <- table(primary_v_chronic_data$Treatment_Regimen, primary_v_chronic_data$Group)

primary_v_chronic_ART <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  num_ART <- primary_v_chronic_data %>% filter(Group == x) %>% filter(ART_Baseline == 'Yes') %>% nrow()
  total_patients <- primary_v_chronic_data %>% filter(Group == x) %>% nrow()
  paste0(num_ART, ' (', round((num_ART / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Concomitant_Disease <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  num_Concomitant_Disease <- primary_v_chronic_data %>% filter(Group == x) %>% filter(Concomitant_Disease == 'Yes') %>% nrow()
  total_patients <- primary_v_chronic_data %>% filter(Group == x) %>% nrow()
  paste0(num_Concomitant_Disease, ' (', round((num_Concomitant_Disease / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_CD4_D0 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_CD4 <- primary_v_chronic_data %>% filter(Group == x) %>% select(CD4_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_CD4['Median'],0), ' (', round(summary_CD4['1st Qu.'],0), '-', round(summary_CD4['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Lymphocytes_D0 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Lymphocytes <- primary_v_chronic_data %>% filter(Group == x) %>% select(Lymphocytes_D0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Platelets_D0 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Platelets <- primary_v_chronic_data %>% filter(Group == x) %>% select(Platelet_Count_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Hemoglobin_D0 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Hemoglobin <- primary_v_chronic_data %>% filter(Group == x) %>% select(Hemoglobin_Count_d0) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Lymphocytes_EOT <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Lymphocytes <- primary_v_chronic_data %>% filter(Group == x) %>% select(Lymphocytes_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Platelets_EOT <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Platelets <- primary_v_chronic_data %>% filter(Group == x) %>% select(Platelet_Count_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Hemoglobin_EOT <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Hemoglobin <- primary_v_chronic_data %>% filter(Group == x) %>% select(Hemoglobin_Count_EOT) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_CD4_Post_M6 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_CD4 <- primary_v_chronic_data %>% filter(Group == x) %>% select(CD4_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_CD4['Median'],0), ' (', round(summary_CD4['1st Qu.'],0), '-', round(summary_CD4['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Lymphocytes_Post_M6 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Lymphocytes <- primary_v_chronic_data %>% filter(Group == x) %>% select(Lymphocytes_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Platelets_Post_M6 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Platelets <- primary_v_chronic_data %>% filter(Group == x) %>% select(Platelet_Count_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

primary_v_chronic_Hemoglobin_Post_M6 <- sapply(levels(primary_v_chronic_data$Group), function(x) {
  summary_Hemoglobin <- primary_v_chronic_data %>% filter(Group == x) %>% select(Hemoglobin_Count_Post_M6) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

