#################################################################
##                    Loading in the packages                   #
#################################################################

suppressMessages(library(here))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(forcats))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(ggtext))
suppressMessages(library(Cairo))

##################################################################
##                    Setting global variables                   #
##################################################################

text_font <- 'Roboto Condensed'

here::i_am("scripts/PatientTimeline.R")

#################################################################
##                        Loading in data                       #
#################################################################

dat <- read_delim(here("data", "meta_data", "Patient_Characteristics", "20220712_PatientTimeline.csv"), col_types = cols(.default = "n"), delim = '\t')
dat$pt_id <- as.factor(dat$pt_id)
dat$timepoint <- as.factor(dat$timepoint)
dat$sample <- as.factor(dat$sample)
dat$treatment <- as.factor(dat$treatment)
dat$relapse <- as.factor(dat$relapse)
dat$death <- as.factor(dat$death)

#dat <- dat %>% filter(., !pt_id %in% c('395', '470', '524'))
dat$pt_id <- fct_drop(dat$pt_id)

primary_patients <- c("102", "114", "117", "305", "421", "517", "524")
dat <- dat %>% mutate(group = case_when(pt_id %in% primary_patients ~ 'primary',
                                      TRUE ~ 'chronic'))

# relapsed_patients <- dat %>% filter(., relapse == 1) %>% pull(pt_id) %>% as.character() %>% unique()
# dat <- dat %>% mutate(group = case_when(pt_id %in% relapsed_patients ~ 'relapse',
#                                         TRUE ~ 'cured'))

dat$month <- as.factor(round(as.numeric(levels(dat$timepoint))[dat$timepoint]/30.41667, digit = 1))
dat$month <- factor(paste0(dat$month, 'M'), levels = paste0(sort(as.numeric(levels(dat$month))), 'M'))
dat$timepoint <- factor(dat$timepoint, levels = sort(as.numeric(levels(dat$timepoint))))
dat <- dat %>% 
  group_by(pt_id) %>% 
  mutate(diff = max(as.numeric(levels(timepoint))[timepoint]) - min(as.numeric(levels(timepoint))[timepoint]))

levels(dat$timepoint) <- paste0(levels(dat$timepoint), 'd')

levels(dat$timepoint)[levels(dat$timepoint) %in% c('-593d')] <- '-20M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-533d', '-538d')] <- '-18M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-476d')] <- '-15.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-449d')] <- '-15M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-425d')] <- '-14M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-362d')] <- '-12M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-311d')] <- '-10M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-249d')] <- '-8M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-228d')] <- '-7.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-195d', '-198d')] <- '-6.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-175d', '-185d')] <- '-6M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-161d', '-163d')] <- '-5.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-145d', '-147d', '-149d')] <- '-5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-136d')] <- '-4.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-125d', '-128d')] <- '-4M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-103d', '-106d', '-108d')] <- '-3.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-83d', '-84d', '-89d', '-91d')] <- '-3M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-71d', '-72d', '-79d', '-80d')] <- '-2.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-57d', '-59d', '-60d')] <- '-2M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-42d', '-43d')] <- '-1.5M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-29d', '-32d', '-37d')] <- '-1M'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-13d', '-15d')] <- '-2W'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-10d')] <- '-1W'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('-20d')] <- '-3W'

levels(dat$timepoint)[levels(dat$timepoint)=='0d'] <- 'D0'
levels(dat$timepoint)[levels(dat$timepoint)=='7d'] <- 'W1'
levels(dat$timepoint)[levels(dat$timepoint)=='14d'] <- 'W2'
levels(dat$timepoint)[levels(dat$timepoint)=='21d'] <- 'W3'
levels(dat$timepoint)[levels(dat$timepoint)=='30d'] <- 'EOT'
levels(dat$timepoint)[levels(dat$timepoint)=='34d'] <- 'W5'
levels(dat$timepoint)[levels(dat$timepoint)=='41d'] <- 'W6'
levels(dat$timepoint)[levels(dat$timepoint)=='48d'] <- 'W7'
levels(dat$timepoint)[levels(dat$timepoint)=='55d'] <- 'W8'

levels(dat$timepoint)[levels(dat$timepoint) %in% c('84d', '96d')] <- 'M3'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('105d', '106d')] <- 'M3.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('121d', '122d', '123d')] <- 'M4'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('135d')] <- 'M4.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('145d', '149d')] <- 'M5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('161d', '173d')] <- 'M5.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('177d', '181d')] <- 'M6'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('198d', '202d')] <- 'M6.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('213d', '216d', '217d', '220d')] <- 'M7'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('231d')] <- 'M7.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('238d', '242d', '244d', '249d')] <- 'M8'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('252d', '256d', '258d')] <- 'M8.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('275d', '279d')] <- 'M9'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('299d', '301d', '304d')] <- 'M10'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('312d', '318d')] <- 'M10.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('337d', '339d')] <- 'M11'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('355d', '356d')] <- 'M11.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('362d', '364d')] <- 'M12'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('378d', '384d')] <- 'M12.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('394d', '395d')] <- 'M13'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('430d')] <- 'M14'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('467d', '469d')] <- 'M15.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('497d')] <- 'M16.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('533d')] <- 'M17.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('540d')] <- 'M18'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('557d')] <- 'M18.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('628d')] <- 'M20.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('665d', '667d')] <- 'M22'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('701d')] <- 'M23'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('717d', '718d')] <- 'M23.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('742d')] <- 'M24'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('812d')] <- 'M26.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('851d')] <- 'M28'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('870d')] <- 'M28.5'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('909d')] <- 'M30'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('947d')] <- 'M31'
levels(dat$timepoint)[levels(dat$timepoint) %in% c('1041d')] <- 'M34'

#dat$pt_id <- factor(dat$pt_id, levels = c('210', '114', '403', '517', '305', '483', '117', '478', '170', '122', '235', '475', '153', '434', '415', '315', '421', '314', '145', '104', '102'))

#TO ADD D0 AS EPISODE
dat <- dat %>% mutate(relapse = case_when(timepoint %in% 'D0' ~ "1",
                                          TRUE ~ as.character(relapse)))
dat$relapse <- as.factor(dat$relapse)

#TO REORDER ALONG GROUPS, PRIMARY V CHRONIC
primary_pt_only <- filter(dat, pt_id %in% primary_patients) %>% mutate(pt_id = fct_drop(pt_id))
primary_pt_only <- fct_reorder(primary_pt_only$pt_id, primary_pt_only$diff, .desc = F)
dat$pt_id <- fct_relevel(fct_reorder(dat$pt_id, dat$diff, .fun = max, .desc = F), levels(primary_pt_only))

#TO REORDER ALONG GROUPS, CURED V RELAPSE
# cured_pt_only <- filter(dat, !pt_id %in% relapsed_patients) %>% mutate(pt_id = fct_drop(pt_id))
# cured_pt_only <- fct_reorder(cured_pt_only$pt_id, cured_pt_only$diff, .desc = F)
# dat$pt_id <- fct_relevel(fct_reorder(dat$pt_id, dat$diff, .fun = max, .desc = F), levels(cured_pt_only))

#Create lookup-table to convert to a named vector so we can paste in the nr. episodes baseline into the y-axis of the plot
episodes_baseline_lookup <- distinct(dat, pt_id, .keep_all = T) %>% select(c(pt_id, episodes_baseline))
episodes_baseline_named_vec <- setNames(episodes_baseline_lookup$episodes_baseline, as.character(episodes_baseline_lookup$pt_id))

##################################################################
##                            Graphing                           #
##################################################################

timeline <- dat %>%
  ggplot(aes(x = timepoint,
             y = pt_id,
             group = pt_id)) +
  geom_line(aes(col = treatment), lwd = 1) +
  geom_point(aes(shape = 'Visit'), size = 4, stroke = 4) +
  geom_vline(xintercept = 'D0', size = 2, colour = "#F39B7FFF") +
  geom_point(data=subset(dat, sample == 1), aes(x = timepoint, y = pt_id, shape = 'Sample taken'), col = "black", size = 7, stroke = 1.25) +
  geom_point(data=subset(dat, relapse == 1), aes(x = timepoint, y = pt_id, shape = 'Episode'), col = "#E64B35FF", size = 10, stroke = 2) +
  geom_point(data=subset(dat, death == 1), aes(x = timepoint, y = pt_id, shape = 'Death'), col = "black", size = 10, stroke = 2) +
  labs(x = 'Timepoint', y = 'Patient ID') +
  ggtitle("Prior VL episodes") +
  scale_shape_manual(name = '', values = c('Visit' = 108, 'Sample taken' = 19, 'Episode' = 5, 'Death' = 4, 'D0' = 108)) +
  scale_color_manual(name = 'VL treatment', values=c("black", "#4DBBD5FF"), labels = c('No', 'Yes')) +
  scale_x_discrete(labels = c('-20M', '-12M', '-9M', '-6M', '-3M', '-1M', 'D0', 'W1', 'W2', 'W3', 'EOT', 'M3', 'M6', 'M9', 'M12', 'M18', 'M24', 'M30'),
                   breaks = c('-20M', '-12M', '-9M', '-6M', '-3M', '-1M', 'D0', 'W1', 'W2', 'W3', 'EOT', 'M3', 'M6', 'M9', 'M12', 'M18', 'M24', 'M30')) + 
  guides(shape = guide_legend(order = 2,
                              override.aes=list(size=c(4,7,10,10,4), col = c('black', 'black', "#E64B35FF", 'black', "#F39B7FFF"), stroke = c(4,1.25,2,2,20))),
         colour = guide_legend(order = 1,
                               title.position = 'top', title.hjust = 0.5,
                               override.aes=list(shape = c('-', '-'), lwd = c(1,1)))) +
  # scale_y_discrete(labels = c("*<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' /><br>",
  #                              rep("<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' /><br>", 2),
  #                              "*<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' /><br>",
  #                              rep("<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' /><br>", 2),
  #                              rep("<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' /><br>", 3),
  #                              "*<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' /><br>",
  #                              rep("<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' /><br>",10),
  #                              "*<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' /><br>")) +
  scale_y_discrete(labels = paste0(c("<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                                     "*<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              "†<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              "<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              "*<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              "<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              "†<img src='./scripts/resources/cured_figure_timeline.png'width='32.5' />",
                              rep("<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' />", 4),
                              "*<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' />",
                              rep("<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' />",11),
                              "*<img src='./scripts/resources/relapse_figure_timeline.png'width='32.5' />"),
                              unname(episodes_baseline_named_vec[levels(dat$pt_id)]))) +
  
  # # scale_y_discrete(labels = c("<img src='./scripts/resources/cured_figure_timeline.png'width='35' /><br>",
  #                             "*<img src='./scripts/resources/cured_figure_timeline.png'width='35' /><br>",
  #                             rep("<img src='./scripts/resources/cured_figure_timeline.png'width='35' /><br>", 3),
  #                             "*<img src='./scripts/resources/cured_figure_timeline.png'width='35' /><br>",
  #                             "<img src='./scripts/resources/relapse_figure_timeline.png'width='35' /><br>",
  #                             "*<img src='./scripts/resources/relapse_figure_timeline.png'width='35' /><br>",
  #                             rep("<img src='./scripts/resources/relapse_figure_timeline.png'width='35' /><br>", 12),
  #                             "*<img src='./scripts/resources/relapse_figure_timeline.png'width='35' /><br>")) +
  theme_classic() +
  theme(axis.text.x = element_text(family = text_font, size = 24, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(size = 36, face = 'bold', family = text_font),
        legend.background = element_rect(color='black', linetype='solid'),
        legend.title = element_text(size=26, face='bold', family = text_font),
        legend.text = element_text(size=22, family = text_font),
        legend.position="top",
        legend.box = 'vertical',
        axis.title = element_text(family = text_font, size = 26),
        plot.title = element_text(size=26, face='bold', family = text_font, vjust = -22.5),
        plot.title.position = 'plot') + ggpubr::rotate_x_text(angle = 60)

# facet_grid(~ month, space = 'free_x', scales = 'free_x', switch = 'x')
#   # remove facet spacing on x-direction
#   theme(panel.spacing.x = unit(0,"line")) +
#   # switch the facet strip label to outside
#   # remove background color
#   theme(strip.placement = 'outside',
#         strip.background.x = element_blank())

# # 
ggsave(here('Analyses', 'final_output', 'Figure2_PatientTimeline.pdf'), timeline, height = 15, width = 25, device = cairo_pdf)
ggsave(here('Analyses', 'final_output', 'Figure2_PatientTimeline.png'), timeline, height = 15, width = 25, type = 'cairo-png')
# 
# # 
# ggsave('PatientTimeline.eps', timeline, height = 12.5, width = 20, device = cairo_ps)

##----------------------------------------------------------------
##        Adding patient flowchart and timeline together         -
##----------------------------------------------------------------

patient_flowchart <- image_read_pdf(here("analyses", "final_output", "20221113_PreLeisH_Flowchart.pdf"), density = 300)
patient_flowchart <- image_ggplot(patient_flowchart, interpolate = T)

layout <- "
AAAA
AAAA
BBBB
BBBB
BBBB
"

composition_arranged <- wrap_elements(full = patient_flowchart) + wrap_elements(full = timeline) +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold', vjust = 2))

ggsave(here('Analyses', 'final_output', 'Figure1_PatientTimeline_Flowchart.pdf'), composition_arranged, height = 25, width = 25, device = cairo_pdf)
ggsave(here('Analyses', 'final_output', 'Figure1_PatientTimeline_Flowchart.png'), composition_arranged, height = 25, width = 25, type = 'cairo-png')
