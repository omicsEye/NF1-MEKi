library(tidyr)
library(dplyr)
library(reshape2)
library(deepath)
library(omicsArt)
library(cowplot)
#source('~/Documents/omicsEye/omicsArt/R/utils.R')
#setting the working directory
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")

number_of_sig_to_keep <- 30
sig_threshold <- 0.1

###############################################################################################

#                                       Tweedieverse                                          #

###############################################################################################
######################
# G 11, 9, 10, 8vs 7 #
######################
## read metabolites
metabolites_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G11_G9_G10_G8_G7/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
#G11
# group 7 is the reference group
metabolites_score_data_G11 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                     metabolites_results$value== "5_4 - 6 weeks of MEKi",]
rownames(metabolites_score_data_G11) <- metabolites_score_data_G11$feature
#G9
metabolites_score_data_G9 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "4_3 weeks of MEKi",]
rownames(metabolites_score_data_G9) <- metabolites_score_data_G9$feature
#G10
metabolites_score_data_G10 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                     metabolites_results$value== "3_4 - 6 weeks of Vehicle",]
rownames(metabolites_score_data_G10) <- metabolites_score_data_G10$feature
#G8
metabolites_score_data_G8 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "2_3 weeks of Vehicle",]
rownames(metabolites_score_data_G8) <- metabolites_score_data_G8$feature

#order by Group 11
order_sig <- rownames(metabolites_score_data_G11)[1:number_of_sig_to_keep]
metabolites_score_data_G11 <- metabolites_score_data_G11[order_sig,]
metabolites_score_data_G11<- metabolites_score_data_G11[order(metabolites_score_data_G11$coef, decreasing=F),]
order_sig <- rownames(metabolites_score_data_G11)
metabolites_score_data_G11 <- within(metabolites_score_data_G11,
                                     feature <- factor(feature,
                                                       levels=order_sig))
#grab metabolites from Group 11 for Groups 9, 10, 8
#group 9
metabolites_score_data_G9 <- metabolites_score_data_G9[rownames(metabolites_score_data_G11),]

metabolites_score_data_G9 <- within(metabolites_score_data_G9,
                                    feature <- factor(feature,
                                                      levels=order_sig))
# group 10
metabolites_score_data_G10 <- metabolites_score_data_G10[rownames(metabolites_score_data_G11),]

metabolites_score_data_G10 <- within(metabolites_score_data_G10,
                                     feature <- factor(feature,
                                                       levels=order_sig))
# Group 8
metabolites_score_data_G8 <- metabolites_score_data_G8[rownames(metabolites_score_data_G11),]

metabolites_score_data_G8 <- within(metabolites_score_data_G8,
                                    feature <- factor(feature,
                                                      levels=order_sig))

metabolites_score_data_G11_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G11, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = '', y_label = '')
metabolites_score_data_G9_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G9, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = '', y_label = '')
metabolites_score_data_G10_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G10, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = '', y_label = '')
metabolites_score_data_G8_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G8, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = '', y_label = '')
##################################
# G6, 5, 4, 3 vs 1: Control Mice #
##################################
## read metabolites
metabolites_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G6_G4_G5_G3_G1/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
#Define Groups 
#G6
metabolites_score_data_G6 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "5_4 - 6 weeks of MEKi",]
rownames(metabolites_score_data_G6) <- metabolites_score_data_G6$feature
#G4
metabolites_score_data_G4 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "4_3 weeks of MEKi",]
rownames(metabolites_score_data_G4) <- metabolites_score_data_G4$feature
#G5
metabolites_score_data_G5 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "3_4 - 6 weeks of Vehicle",]
rownames(metabolites_score_data_G5) <- metabolites_score_data_G5$feature
#G3
metabolites_score_data_G3 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "2_3 weeks of Vehicle",]
rownames(metabolites_score_data_G3) <- metabolites_score_data_G3$feature

#Order by Group 6
order_sig <- rownames(metabolites_score_data_G6)[1:number_of_sig_to_keep]
metabolites_score_data_G6 <- metabolites_score_data_G6[order_sig,]
metabolites_score_data_G6<- metabolites_score_data_G6[order(metabolites_score_data_G6$coef, decreasing=F),]
order_sig <- rownames(metabolites_score_data_G6)
metabolites_score_data_G6 <- within(metabolites_score_data_G6,
                                    feature <- factor(feature,
                                                      levels=order_sig))
#Filter G4
metabolites_score_data_G4 <- metabolites_score_data_G4[rownames(metabolites_score_data_G6),]
metabolites_score_data_G4 <- within(metabolites_score_data_G4,
                                    feature <- factor(feature,
                                                      levels=order_sig))
#Filter G5
metabolites_score_data_G5 <- metabolites_score_data_G5[rownames(metabolites_score_data_G6),]
metabolites_score_data_G5 <- within(metabolites_score_data_G5,
                                    feature <- factor(feature,
                                                      levels=order_sig))
#Filter G3
metabolites_score_data_G3 <- metabolites_score_data_G3[rownames(metabolites_score_data_G6),]
metabolites_score_data_G3 <- within(metabolites_score_data_G3,
                                    feature <- factor(feature,
                                                      levels=order_sig))

metabolites_score_data_G6_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G6, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
metabolites_score_data_G4_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G4, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
metabolites_score_data_G5_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G5, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
metabolites_score_data_G3_diff_bar_tweedie <- diff_bar_plot(metabolites_score_data_G3, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                    fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
###############################################################################################

########################################## Plots ##############################################

###############################################################################################

#################
# NF1 Mice only #
#################
fig3_metabolites_G8_G9_G10_G11_G7_tweedie <- plot_grid(
  metabolites_score_data_G8_diff_bar_tweedie,
  metabolites_score_data_G9_diff_bar_tweedie, 
  metabolites_score_data_G10_diff_bar_tweedie, 
  metabolites_score_data_G11_diff_bar_tweedie + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.line.y = element_blank()),
  ncol = 4,
  byrow = TRUE
) ;fig3_metabolites_G8_G9_G10_G11_G7_tweedie

#####################
# Control Mice only #
#####################
fig3_metabolites_G6_G5_G4_G3_G1_tweedie <- plot_grid(
  metabolites_score_data_G3_diff_bar_tweedie, 
  metabolites_score_data_G4_diff_bar_tweedie,
  metabolites_score_data_G5_diff_bar_tweedie, 
  metabolites_score_data_G6_diff_bar_tweedie + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.line.y = element_blank()),
  ncol = 4,
  byrow = TRUE
) ;fig3_metabolites_G6_G5_G4_G3_G1_tweedie
#################
# Altogether #
#################
# ggdraw plots 
box_association_1 <- readRDS("analysis/Tweedieverse_MOUSE_CDT_G11_G9_G10_G8_G7/figures/Treatment_gg_associations.RDS")
box_association_2 <- readRDS("analysis/Tweedieverse_MOUSE_CDT_G6_G4_G5_G3_G1/figures/Treatment_gg_associations.RDS")

## do plots
fig3_metabolites_tweedie <- ggdraw() +
  draw_plot(metabolites_score_data_G8_diff_bar_tweedie,
            x = 0, y = .52, width = .4, height = .45) +
  draw_plot(metabolites_score_data_G9_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_blank(),
                                                       axis.ticks.y = element_blank(),
                                                       axis.line.y = element_blank()),
            x = .4, y = .52, width = .16, height = .45) +
  draw_plot(metabolites_score_data_G10_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.line.y = element_blank()),
            x = .56, y = .52, width = .22, height = .45) +
  draw_plot(metabolites_score_data_G11_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.line.y = element_blank()),
            x = .78, y = .52, width = .22, height = .45) +
  draw_plot(metabolites_score_data_G3_diff_bar_tweedie,
            x = 0, y = 0.01, width = .3, height = .45) +
  draw_plot(metabolites_score_data_G4_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_blank(),
                                                       axis.ticks.y = element_blank(),
                                                       axis.line.y = element_blank()),
            x = .3, y = 0.01, width = .15, height = .45) +
  draw_plot(metabolites_score_data_G5_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_blank(),
                                                       axis.ticks.y = element_blank(),
                                                       axis.line.y = element_blank()),
            x = .45, y = 0.01, width = .275, height = .45) +
  draw_plot(metabolites_score_data_G6_diff_bar_tweedie + theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_blank(),
                                                       axis.ticks.y = element_blank(),
                                                       axis.line.y = element_blank()),
            x = .725, y = 0.01, width = .275, height = .45) +
  draw_plot(box_association_1[[14]]
            + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45)+
  draw_plot_label((label = c("NF1 Vehicle 3 weeks", "NF1 MEKi 3 weeks", "NF1 Vehicle P.I. 4-6 weeks", "NF1 MEKi P.I.4-6 weeks",
                             "Control Vehicle 3 weeks", "Control MEKi 3 weeks", "Control Vehicle P.I. 4-6 weeks", "Control MEKi P.I. 4-6 weeks")),
                  size = 7,x = c(.15, .35, .55, .75, .15, .35, .55, .75), 
                  y = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5))
fig3_metabolites_tweedie

#save files to folder
ggsave(filename = 'manuscript/Figures/figure_3/fig3_metabolites_tweedie.pdf', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_3/fig3_metabolites_tweedie.png', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)







# #To do: #choose metabolites based on 4-6 week MEKi for both rows
# #first row: NF-1 mice: 3 weeks vehicle,3 weeks MEKi, 4-6 weeks vehicle, 4-6 weeks MEKi
# #second row: control mice: 3 weeks vehicle,3 weeks MEKi, 4-6 weeks vehicle, 4-6 weeks MEKi
# ###############################################################################################
# 
# #                                        Maaslin                                              #
# 
# ###############################################################################################
# ######################
# # G 11, 9, 10, 8vs 7 #
# ######################
# ## read metabolites
# metabolites_results <- read.delim(
#   "analysis/Maaslin2_MOUSE_CDT_G11_G9_G10_G8_G7/all_results.tsv",
#   sep = '\t',
#   header = T,
#   fill = F,
#   comment.char = "" ,
#   check.names = F,
#   #row.names = NA
# )
# #G11
# metabolites_score_data_G11 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                    metabolites_results$value== "5_4 - 6 weeks of MEKi",]
# rownames(metabolites_score_data_G11) <- metabolites_score_data_G11$feature
# #G9
# metabolites_score_data_G9 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "4_3 weeks of MEKi",]
# rownames(metabolites_score_data_G9) <- metabolites_score_data_G9$feature
# #G10
# metabolites_score_data_G10 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "3_4 - 6 weeks of Vehicle",]
# rownames(metabolites_score_data_G10) <- metabolites_score_data_G10$feature
# #G8
# metabolites_score_data_G8 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "2_3 weeks of Vehicle",]
# rownames(metabolites_score_data_G8) <- metabolites_score_data_G8$feature
# 
# #order by Group 11
# order_sig <- rownames(metabolites_score_data_G11)[1:number_of_sig_to_keep]
# metabolites_score_data_G11 <- metabolites_score_data_G11[order_sig,]
# metabolites_score_data_G11<- metabolites_score_data_G11[order(metabolites_score_data_G11$coef, decreasing=F),]
# order_sig <- rownames(metabolites_score_data_G11)
# metabolites_score_data_G11 <- within(metabolites_score_data_G11,
#                                         feature <- factor(feature,
#                                                           levels=order_sig))
# #grab metabolites from Group 11 for Groups 9, 10, 8
# #group 9
# metabolites_score_data_G9 <- metabolites_score_data_G9[rownames(metabolites_score_data_G11),]
# 
# metabolites_score_data_G9 <- within(metabolites_score_data_G9,
#                                             feature <- factor(feature,
#                                                               levels=order_sig))
# # group 10
# metabolites_score_data_G10 <- metabolites_score_data_G10[rownames(metabolites_score_data_G11),]
# 
# metabolites_score_data_G10 <- within(metabolites_score_data_G10,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# # Group 8
# metabolites_score_data_G8 <- metabolites_score_data_G8[rownames(metabolites_score_data_G11),]
# 
# metabolites_score_data_G8 <- within(metabolites_score_data_G8,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# 
# metabolites_score_data_G11_diff_bar <- diff_bar_plot(metabolites_score_data_G11, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                   fdr ="qval", orderby = NA, x_label = '', y_label = '')
# metabolites_score_data_G9_diff_bar <- diff_bar_plot(metabolites_score_data_G9, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                       fdr ="qval", orderby = NA, x_label = '', y_label = '')
# metabolites_score_data_G10_diff_bar <- diff_bar_plot(metabolites_score_data_G10, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                      fdr ="qval", orderby = NA, x_label = '', y_label = '')
# metabolites_score_data_G8_diff_bar <- diff_bar_plot(metabolites_score_data_G8, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                     fdr ="qval", orderby = NA, x_label = '', y_label = '')
# ##################################
# # G6, 5, 4, 3 vs 1: Control Mice #
# ##################################
# ## read metabolites
# metabolites_results <- read.delim(
#   "analysis/Maaslin2_MOUSE_CDT_G6_G4_G5_G3_G1/all_results.tsv",
#   sep = '\t',
#   header = T,
#   fill = F,
#   comment.char = "" ,
#   check.names = F,
#   #row.names = NA
# )
# #Define Groups 
# #G6
# metabolites_score_data_G6 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "5_4 - 6 weeks of MEKi",]
# rownames(metabolites_score_data_G6) <- metabolites_score_data_G6$feature
# #G4
# metabolites_score_data_G4 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "4_3 weeks of MEKi",]
# rownames(metabolites_score_data_G4) <- metabolites_score_data_G4$feature
# #G5
# metabolites_score_data_G5 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "3_4 - 6 weeks of Vehicle",]
# rownames(metabolites_score_data_G5) <- metabolites_score_data_G5$feature
# #G3
# metabolites_score_data_G3 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
#                                                     metabolites_results$value== "2_3 weeks of Vehicle",]
# rownames(metabolites_score_data_G3) <- metabolites_score_data_G3$feature
# 
# #Order by Group 6
# order_sig <- rownames(metabolites_score_data_G6)[1:number_of_sig_to_keep]
# metabolites_score_data_G6 <- metabolites_score_data_G6[order_sig,]
# metabolites_score_data_G6<- metabolites_score_data_G6[order(metabolites_score_data_G6$coef, decreasing=F),]
# order_sig <- rownames(metabolites_score_data_G6)
# metabolites_score_data_G6 <- within(metabolites_score_data_G6,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# #Filter G4
# metabolites_score_data_G4 <- metabolites_score_data_G4[rownames(metabolites_score_data_G6),]
# metabolites_score_data_G4 <- within(metabolites_score_data_G4,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# #Filter G5
# metabolites_score_data_G5 <- metabolites_score_data_G5[rownames(metabolites_score_data_G6),]
# metabolites_score_data_G5 <- within(metabolites_score_data_G5,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# #Filter G3
# metabolites_score_data_G3 <- metabolites_score_data_G3[rownames(metabolites_score_data_G6),]
# metabolites_score_data_G3 <- within(metabolites_score_data_G3,
#                                     feature <- factor(feature,
#                                                       levels=order_sig))
# 
# metabolites_score_data_G6_diff_bar <- diff_bar_plot(metabolites_score_data_G6, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
# metabolites_score_data_G4_diff_bar <- diff_bar_plot(metabolites_score_data_G4, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
# metabolites_score_data_G5_diff_bar <- diff_bar_plot(metabolites_score_data_G5, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
# metabolites_score_data_G3_diff_bar <- diff_bar_plot(metabolites_score_data_G3, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
#                                                     fdr ="qval", orderby = NA, x_label = 'Effect Size', y_label = '')
# ###############################################################################################
# ###############################################################################################
# ########################################## Plots ##############################################
# ###############################################################################################
# ###############################################################################################
# 
# #################
# # NF1 Mice only #
# #################
# fig3_metabolites_G8_G9_G10_G11_G7 <- plot_grid(
#   metabolites_score_data_G8_diff_bar,metabolites_score_data_G9_diff_bar, 
#   metabolites_score_data_G10_diff_bar, metabolites_score_data_G11_diff_bar + theme(axis.title.y = element_blank(),
#                                                                                    axis.text.y = element_blank(),
#                                                                                    axis.ticks.y = element_blank(),
#                                                                                    axis.line.y = element_blank()),
#   ncol = 4,
#   byrow = TRUE
# ) ;fig3_metabolites_G8_G9_G10_G11_G7
# 
# #####################
# # Control Mice only #
# #####################
# fig3_metabolites_G6_G5_G4_G3_G1 <- plot_grid(
#   metabolites_score_data_G3_diff_bar, metabolites_score_data_G4_diff_bar,
#   metabolites_score_data_G5_diff_bar, metabolites_score_data_G6_diff_bar+ 
#     theme(axis.title.y = element_blank(), 
#           axis.text.y = element_blank(),
#           axis.line.y = element_blank()),
#   ncol = 4,
#   byrow = TRUE
# ) ;fig3_metabolites_G6_G5_G4_G3_G1
# #################
# # Altogether #
# #################
# # ggdraw plots 
# fig3_metabolites <- ggdraw() +
#   draw_plot(metabolites_score_data_G8_diff_bar,
#             x = 0, y = .52, width = .4, height = .45) +
#   draw_plot(metabolites_score_data_G9_diff_bar + theme(axis.title.y = element_blank(),
#                                                          axis.text.y = element_blank(),
#                                                          axis.ticks.y = element_blank(),
#                                                          axis.line.y = element_blank()),
#             x = .4, y = .52, width = .16, height = .45) +
#   draw_plot(metabolites_score_data_G10_diff_bar + theme(axis.title.y = element_blank(),
#                                                         axis.text.y = element_blank(),
#                                                         axis.ticks.y = element_blank(),
#                                                         axis.line.y = element_blank()),
#             x = .56, y = .52, width = .22, height = .45) +
#   draw_plot(metabolites_score_data_G11_diff_bar + theme(axis.title.y = element_blank(),
#                                                         axis.text.y = element_blank(),
#                                                         axis.ticks.y = element_blank(),
#                                                         axis.line.y = element_blank()),
#             x = .78, y = .52, width = .22, height = .45) +
#   draw_plot(metabolites_score_data_G3_diff_bar,
#             x = 0, y = 0.01, width = .3, height = .45) +
#  draw_plot(metabolites_score_data_G4_diff_bar + theme(axis.title.y = element_blank(),
#                                                        axis.text.y = element_blank(),
#                                                        axis.ticks.y = element_blank(),
#                                                        axis.line.y = element_blank()),
#            x = .3, y = 0.01, width = .15, height = .45) +
#   draw_plot(metabolites_score_data_G5_diff_bar + theme(axis.title.y = element_blank(),
#                                                         axis.text.y = element_blank(),
#                                                         axis.ticks.y = element_blank(),
#                                                         axis.line.y = element_blank()),
#             x = .45, y = 0.01, width = .275, height = .45) +
#    draw_plot(metabolites_score_data_G6_diff_bar + theme(axis.title.y = element_blank(),
#                                                         axis.text.y = element_blank(),
#                                                         axis.ticks.y = element_blank(),
#                                                         axis.line.y = element_blank()),
#              x = .725, y = 0.01, width = .275, height = .45) +
#   draw_plot_label((label = c("NF1 Vehicle 3 weeks", "NF1 MEKi 3 weeks", "NF1 Vehicle P.I. 4-6 weeks", "NF1 MEKi P.I.4-6 weeks",
#                              "Control Vehicle 3 weeks", "Control MEKi 3 weeks", "Control Vehicle P.I. 4-6 weeks", "Control MEKi P.I. 4-6 weeks")),
#                   size = 7,x = c(.15, .35, .55, .75, .15, .35, .55, .75), 
#                   y = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5))
# fig3_metabolites
# 
# #save files to folder
# ggsave(filename = 'manuscript/Figures/figure_3/fig3_metabolites.pdf', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)
# ggsave(filename = 'manuscript/Figures/figure_3/fig3_metabolites.png', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)
# 
# #####################
# ## with plot grid ##
# ####################
# 
# # fig3_metabolites_altogether <- plot_grid(
# #   metabolites_score_data_G8_diff_bar + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                              axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G9_diff_bar + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                               axis.ticks.y = element_blank(),axis.line.y = element_blank()),
# #   metabolites_score_data_G10_diff_bar+ theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                              axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G11_diff_bar+ theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                              axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G3_diff_bar+ theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                             axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G4_diff_bar+ theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                             axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G5_diff_bar+ theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                             axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   metabolites_score_data_G6_diff_bar + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
# #                                              axis.ticks.y = element_blank(), axis.line.y = element_blank()),
# #   ncol = 4, nrow=2, byrow = TRUE)
# # fig3_metabolites_altogether
# # #Save plots in Box
# # ggsave(filename = 'manuscript/Figures/figure_3/fig3_metabolites_altogether.pdf', 
# #        plot=fig3_metabolites_altogether, width = 300, height = 110, units = "mm", dpi = 350)
# #ggsave(filename = '/Users/rah/Box/NF1_MEKi/Figures/Figure_2/fig2_metabolites.png', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)
# 
# # fig3_metabolites_G6_G5_G4_G3_G1 <- plot_grid(
# #   metabolites_score_data_G3_diff_bar, metabolites_score_data_G4_diff_bar,
# #   metabolites_score_data_G5_diff_bar, metabolites_score_data_G6_diff_bar+  theme(axis.title.y = element_blank(),
# #                                                                                  axis.text.y = element_blank(),
# #                                                                                  axis.ticks.y = element_blank(),
# #                                                                                  axis.line.y = element_blank()),
# #   ncol = 4,
# #   byrow = TRUE
# # ) ;fig3_metabolites_G6_G5_G4_G3_G1