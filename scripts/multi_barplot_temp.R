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

## read metabolites
metabolites_results <- read.delim(
  "analysis/Maaslin2_MOUSE_CDT_G11_G9_G10_G8_G7/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_score_data_G11 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                   metabolites_results$value== "5_4 - 6 weeks of MEKi",]
rownames(metabolites_score_data_G11) <- metabolites_score_data_G11$feature


metabolites_score_data_G9 <- metabolites_results [metabolites_results$metadata=="Treatment" & 
                                                    metabolites_results$value== "4_3 weeks of MEKi",]
rownames(metabolites_score_data_G9) <- metabolites_score_data_G9$feature



order_sig <- rownames(metabolites_score_data_G11)[1:number_of_sig_to_keep]
metabolites_score_data_G11 <- metabolites_score_data_G11[order_sig,]
metabolites_score_data_G11<- metabolites_score_data_G11[order(metabolites_score_data_G11$coef, decreasing=F),]
order_sig <- rownames(metabolites_score_data_G11)
metabolites_score_data_G11 <- within(metabolites_score_data_G11,
                                        feature <- factor(feature,
                                                          levels=order_sig))
metabolites_score_data_G9 <- metabolites_score_data_G9[rownames(metabolites_score_data_G11),]

metabolites_score_data_G9 <- within(metabolites_score_data_G9,
                                            feature <- factor(feature,
                                                              levels=order_sig))


metabolites_score_data_G11_diff_bar <- diff_bar_plot(metabolites_score_data_G11, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
metabolites_score_data_G9_diff_bar <- diff_bar_plot(metabolites_score_data_G9, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')

## read association
#box_association <- readRDS("analysis/Maaslin2_MOUSE_CDT_G11_G9_G10_G8_G7/figures/Treatment_gg_associations.RDS")
## do plots

fig3_metabolites <- plot_grid(
  metabolites_score_data_G9_diff_bar, metabolites_score_data_G11_diff_bar+  theme(axis.title.y = element_blank(),
                                                                                   axis.text.y = element_blank(),
                                                                                   axis.ticks.y = element_blank(),
                                                                                   axis.line.y = element_blank()),
  ncol = 2,
  byrow = TRUE
) 
fig2_metabolites <- ggdraw() +
  draw_plot(metabolites_severe_temp_diff_bar,
            x = 0, y = .47, width = .55, height = .53) +
  draw_plot(metabolites_non_severe_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .55, y = .47, width = .225, height = .53) +
  draw_plot(metabolites_non_covid_temp_diff_bar + theme(axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.line.y = element_blank()),
            x = .775, y = .47, width = .225, height = .53) +
  draw_plot(box_association[[14]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[28]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[83]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[37]] + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .75, y = 0, width = .25, height = .45) +

  draw_plot_label((label = c("a",  "Human", "Mice 2 Hours", "Mice 12 hours", "b", "c", "d", "e")),
                  size = 7,x = c(0, .25, .5, .76, 0, .25, .5, .75), y = c(1, 1, 1, 1, 0.47, 0.47, 0.47, 0.47))
fig2_metabolites

ggsave(filename = '/Users/rah/Box/NF1_MEKi/Figures/Figure_2/fig2_metabolites.pdf', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)
ggsave(filename = '/Users/rah/Box/NF1_MEKi/Figures/Figure_2/fig2_metabolites.png', plot=fig2_metabolites, width = 183, height = 110, units = "mm", dpi = 350)

