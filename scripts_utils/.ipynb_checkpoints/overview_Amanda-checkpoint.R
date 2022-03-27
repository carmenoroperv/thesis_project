# Author:   Amanda Frydendahl
# Date:     20-01-2022
# Purpose:  plot clinical follow-up for full C2i cohort

library(dplyr)
library(ggpubr)
library(ggplot2)
library(data.table)

############ LOAD DATA ###########
setwd(dir = "O:\\HE_colongruppen/C2i project/Klinisk info/trimmed_files/")
c2i_follow_up <- read.csv("2022-01-20_c2i_follow-up.csv")
c2i_FU <- read.csv("2022-01-20_c2i_FU_length.csv")
c2i_intervention_chemo <- read.csv("2022-01-20_c2i_intervention_chemo.csv")
c2i_intervention_other <- read.csv("2022-01-20_c2i_intervention_other.csv")
c2i_act <- read.csv("2022-01-20_c2i_act.csv")
c2i_image <- read.csv("2022-01-20_c2i_image.csv")
c2i_image_relapse <- read.csv("2022-01-20_c2i_image_relapse.csv")

#test data
setwd(dir = "O:\\HE_colongruppen/C2i project/Results/Phase II/")
c2i_data <- read.csv("Sample_calls/2022-02-23_c2i_phaseII_131pts_high_accuracy.csv", sep = ";")

############ trim data #########
#dataframe with sorted pts. ID's
pts_sort <- c2i_follow_up[order(c2i_follow_up$is_relapse, c2i_follow_up$last_ctscan),] %>% 
  .$biobankID %>% as.data.frame()
pts_sort$pts_order <- 1:146

setnames(x = pts_sort,
         old = ".",
         new = "biobankID")

#add "pts_order" to all data frames
c2i_follow_up <- merge(c2i_follow_up, pts_sort, "biobankID")
c2i_intervention_chemo <- merge(c2i_intervention_chemo, pts_sort, "biobankID")
c2i_intervention_other <- merge(c2i_intervention_other, pts_sort, "biobankID")
c2i_act <- merge(c2i_act, pts_sort, "biobankID")
c2i_image <- merge(c2i_image, pts_sort, "biobankID")
c2i_image_relapse <- merge(c2i_image_relapse, pts_sort, "biobankID")
c2i_data <- merge(c2i_data, pts_sort, "biobankID")

#add op-time to sample_data
c2i_data <- merge(c2i_data, c2i_follow_up, "biobankID") %>% select(., c(biobankID, C2.Test, sampleID, Tumor.Fraction, 
                                                             Coverage, sample_date, pts_order.x, op_date))
#add sample timepoints (months since OP)
c2i_data$sample_timepoint_months <- 
  (as.Date(c2i_data$sample_date, "%d-%m-%Y") - as.Date(c2i_data$op_date, "%Y-%m-%d"))/(365.25/12) %>% as.numeric()

#vector with biobank ID's in the same order as sorted data frame. Used for re-labeling of y-axis
biobankID_order <- c2i_follow_up$biobankID


# set all pre OP samples to timepoint 0
#c2i_data$sample_timepoint_months[c2i_data$sample_timepoint_months < 0 ] <- 0

c2i_follow_up <- c2i_follow_up[order(c2i_follow_up$pts_order),]

c2i_data <- c2i_data[order(c2i_data$pts_order),]

labels <- c2i_data %>% select(biobankID, pts_order.x) %>% unique()

pts_list <- c2i_data$biobankID %>% unique()

############ plot data ######
overview <- ggplot(c2i_follow_up[c2i_follow_up$biobankID %in% pts_list,], aes(x = last_ctscan/(365.25/12), y = as.factor(pts_order))) + 
  #add follow-up length
  geom_linerange(aes(xmin = 0, 
                     xmax = last_ctscan/(365.25/12), 
                     y = as.factor(pts_order))) +
  #add baseline
  geom_vline(xintercept = 0) +
  #add act treatment
  geom_segment(data = c2i_act[c2i_act$biobankID %in% pts_list,], aes(x = act_start_months, y = as.factor(pts_order), 
                                   xend = act_end_months, yend = as.factor(pts_order), 
                                   alpha = 0.1), color = "pink", size = 1.2) +
  #add intervention (chemo)
  geom_segment(data = c2i_intervention_chemo[c2i_intervention_chemo$biobankID %in% pts_list,], aes(x = intervention_start_months, y = as.factor(pts_order),
                                            xend = intervention_end_months, yend = as.factor(pts_order), 
                                            alpha = 0.1), color = "steelblue2", size = 1.2) +
  #add intervention (other)
  geom_point(data = c2i_intervention_other[c2i_intervention_other$biobankID %in% pts_list,], aes(x = intervention_start_months, y = as.factor(pts_order),
                                                color = "intervention"), color = "steelblue2", size = 0.8, shape = 8) +
  #add scans
  geom_point(data = c2i_image[c2i_image$biobankID %in% pts_list,], aes(x = time_since_op_days/(365.25/12),
                                  y = as.factor(pts_order),
                                  color = as.factor(imaging_result)), shape = 6, size = 0.8) +
  #add relapse
  geom_point(data = c2i_image_relapse[c2i_image_relapse$biobankID %in% pts_list,], aes(x = time_since_op_days/(365.25/12),
                                   y = as.factor(pts_order)), color = "red", shape = 18, size = 1.5) +
  #add blood samples
  geom_point(data = c2i_data, aes(x = sample_timepoint_months, y = as.factor(pts_order.x), 
                                   fill = as.factor(C2.Test), stroke = 0.3),
                                   shape = 21, size = 1) +
  scale_color_manual(values = c("0" = "grey", "2" = "yellowgreen", "3" = "forestgreen", 
                                "4" = "pink", "5" = "orange", "9" = "purple", 
                                "black", "black")) +
  scale_fill_manual(values = c("darkgoldenrod1", "black", "white")) +
  #define theme
  scale_y_discrete(labels = labels$biobankID) +
  scale_x_continuous(breaks = seq(0, 80, by = 4)) +
  #xlim(0,60) +
  labs(
    title = "Phase II, n = 131 \n High spec cut-off", 
    y = "", 
    x = "Months since OP"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.text=element_text(size=4),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

overview
#export plots
ggexport(overview, filename = "O:\\HE_colongruppen/C2i project/Results/Phase II/Plots/round1_130pts/overview_131pts_high_spec_with_scans.pdf")
