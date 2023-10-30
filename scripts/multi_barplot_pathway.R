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

####################################################################################################

#                                        Tweedieverse                                              #

####################################################################################################
#####################
# Humans: G17 vs 16 #
#####################
## read pathway
#G17 vs G16
pathway_results_human <- read.delim(
  "analysis/enrichment_analysis/tweedie_deepath_SUB_PATHWAYS_HD4/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
##################
# Mouse G9 vs G7 #
##################
pathway_results_mouse_G9G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G9_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
##################
# Mouse G11 vs G7 #
##################
pathway_results_mouse_G11G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G11_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
#intersection between pathway names between human and mouse data. Then only order based on humans. 
# or plot enrichment score
#pathway_score_data_human <- pathway_results_human[pathway_results_human$fdr <= sig_threshold,]
# pathway_score_data_human <- pathway_results_human[,c('pathway', 'pval','fdr','set_enrichment_score')]
pathways <- intersect(pathway_results_human$pathway,pathway_results_mouse_G9G7$pathway)

#select human data
pathway_score_data_human <- pathway_results_human
rownames(pathway_score_data_human) <- pathway_score_data_human$pathway
#choose pathways common to both mouse and human data
pathway_score_data_human <- pathway_score_data_human[pathways,]
colnames(pathway_score_data_human) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')

#select mouse G9G7 data
pathway_score_data_mouse_G9G7 <- pathway_results_mouse_G9G7
rownames(pathway_score_data_mouse_G9G7) <- pathway_results_mouse_G9G7$pathway
pathway_score_data_mouse_G9G7 <- pathway_score_data_mouse_G9G7[pathways,]
colnames(pathway_score_data_mouse_G9G7) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')

#select mouse G11G7 data
pathway_score_data_mouse_G11G7 <- pathway_results_mouse_G11G7
rownames(pathway_score_data_mouse_G11G7) <- pathway_score_data_mouse_G11G7$pathway
pathway_score_data_mouse_G11G7 <- pathway_score_data_mouse_G11G7[pathways,]
colnames(pathway_score_data_mouse_G11G7) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')


#order by Humans
# order_sig <- rownames(pathway_score_data_human)[1:number_of_sig_to_keep]
order_sig <- rownames(pathway_score_data_human)
pathway_score_data_human <- pathway_score_data_human[order(pathway_score_data_human$coef, decreasing=F),]
order_sig <- rownames(pathway_score_data_human)[1:number_of_sig_to_keep]
pathway_score_data_human <- pathway_score_data_human[order_sig,]
# order_sig <- rownames(pathway_score_data_human)
# pathway_score_data_human <- within(pathway_score_data_human,
#                                         feature <- factor(pathway,
#                                                           levels=order_sig))
pathway_score_data_human$feature <- factor(pathway_score_data_human$pathway, levels=order_sig)
pathway_score_data_human <- pathway_score_data_human[,c(9,2,3,4,5,6,7,8)]


#grab pathway data from Humans for Mouse Data (G9G7)
pathway_score_data_mouse_G9G7 <- pathway_score_data_mouse_G9G7[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G9G7$feature <- factor(pathway_score_data_mouse_G9G7$pathway, levels=order_sig)
pathway_score_data_mouse_G9G7 <- pathway_score_data_mouse_G9G7[,c(9,2,3,4,5,6,7,8)]
#grab pathway data from Humans for Mouse Data (G11G7)
pathway_score_data_mouse_G11G7 <- pathway_score_data_mouse_G11G7[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G11G7$feature <- factor(pathway_score_data_mouse_G11G7$pathway, levels=order_sig)
pathway_score_data_mouse_G11G7 <- pathway_score_data_mouse_G11G7[,c(9,2,3,4,5,6,7,8)]

pathway_score_data_human_diff_bar <- diff_bar_plot(pathway_score_data_human, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="fdr", orderby = NA, x_label = '', y_label = '')
pathway_score_data_G9G7_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G9G7, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="fdr", orderby = NA, x_label = '', y_label = '')
pathway_score_data_G11G7_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G11G7, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="fdr", orderby = NA, x_label = '', y_label = '')
###############################################################################################

########################################## Plots ##############################################

###############################################################################################

#################
# NF1 Mice only #
#################
# fig4_pathwayenrichment <- plot_grid(
#   pathway_score_data_human_diff_bar, pathway_score_data_G9G7_diff_bar + 
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.line.y = element_blank()),
#   ncol = 2,
#   byrow = TRUE
# ) ;fig4_pathwayenrichment

#################
# Altogether #
#################
# ggdraw plots 
fig4_pathway <- ggdraw() +
  draw_plot(pathway_score_data_human_diff_bar,
            x = 0, y = 0.06, width = .41, height = 0.86) +
  draw_plot(pathway_score_data_G9G7_diff_bar + theme(axis.title.y = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.line.y = element_blank()),
            x = .41, y = 0.06, width = .32, height = 0.86) +
  draw_plot(pathway_score_data_G11G7_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .73, y = 0.06, width = .26, height = 0.86) +
  draw_plot_label((label = c("Humans", "G9 vs G7", "G11 vs G7")),
                  size = 7,x = c(.25, .5, .75), 
                  y = c(.95, .95, .95)) + 
  draw_plot_label((label = c('Enrichment Score')),
                  size = 7,x = 0.5, 
                  y = 0.05) 
fig4_pathway


#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/fig4_pathway.pdf', plot=fig4_pathway, width = 300, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/fig4_pathway.png', plot=fig4_pathway, width = 300, height = 110, units = "mm", dpi = 350)

##################################################################################
# Control Mice
########
# G4G1#
########
pathway_results_mouse_G4G1 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G4_G1/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

########
# G6G1#
########
pathway_results_mouse_G6G1 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G6_G1/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
# pathways <- intersect(pathway_results_human$pathway,pathway_results_mouse_G4G1$pathway)
#select mouse G4G1 data
pathway_score_data_mouse_G4G1 <- pathway_results_mouse_G4G1
rownames(pathway_score_data_mouse_G4G1) <- pathway_score_data_mouse_G4G1$pathway
pathway_score_data_mouse_G4G1 <- pathway_score_data_mouse_G4G1[pathways,]
colnames(pathway_score_data_mouse_G4G1) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')
#select mouse G6G1 data
pathway_score_data_mouse_G6G1 <- pathway_results_mouse_G6G1
rownames(pathway_score_data_mouse_G6G1) <- pathway_score_data_mouse_G6G1$pathway
pathway_score_data_mouse_G6G1 <- pathway_score_data_mouse_G6G1[pathways,]
colnames(pathway_score_data_mouse_G6G1) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')


#grab pathway data from Humans for Mouse Data (G9G7)
pathway_score_data_mouse_G4G1 <- pathway_score_data_mouse_G4G1[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G4G1$feature <- factor(pathway_score_data_mouse_G4G1$pathway, levels=order_sig)
pathway_score_data_mouse_G4G1 <- pathway_score_data_mouse_G4G1[,c(9,2,3,4,5,6,7,8)]
#grab pathway data from Humans for Mouse Data (G11G7)
pathway_score_data_mouse_G6G1 <- pathway_score_data_mouse_G6G1[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G6G1$feature <- factor(pathway_score_data_mouse_G6G1$pathway, levels=order_sig)
pathway_score_data_mouse_G6G1 <- pathway_score_data_mouse_G6G1[,c(9,2,3,4,5,6,7,8)]

#define barplot objects
pathway_score_data_G4G1_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G4G1, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="fdr", orderby = NA, x_label = '', y_label = '')
pathway_score_data_G6G1_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G6G1, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                   fdr ="fdr", orderby = NA, x_label = '', y_label = '')

# ggdraw plots 
fig4b_pathway <- ggdraw() +
  draw_plot(pathway_score_data_human_diff_bar,
            x = 0, y = 0.06, width = .44, height = 0.86) +
  draw_plot(pathway_score_data_G4G1_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .44, y = 0.06, width = .27, height = 0.86) +
  draw_plot(pathway_score_data_G6G1_diff_bar + theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank(),
                                                      axis.line.y = element_blank()),
            x = .71, y = 0.06, width = .27, height = 0.86) +
  draw_plot_label((label = c("Humans", "G4 vs G1", "G6 vs G1")),
                  size = 7,x = c(.25, .5, .75), 
                  y = c(.95, .95, .95)) + 
  draw_plot_label((label = c('Enrichment Score')),
                  size = 7,x = 0.5, 
                  y = 0.05) 
fig4b_pathway


#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/fig4b_control_pathway.pdf', plot=fig4b_pathway, width = 300, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/fig4b_control_pathway.png', plot=fig4b_pathway, width = 300, height = 110, units = "mm", dpi = 350)

#do this for vehicle? 8 vs 7 10 vs 7
##################
# Mouse G8 vs G7 #
##################
pathway_results_mouse_G8G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G8_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
##################
# Mouse G10 vs G7 #
##################
pathway_results_mouse_G10G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G10_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
#intersection between pathway names between human and mouse data. Then only order based on humans. 
# or plot enrichment score
#pathway_score_data_human <- pathway_results_human[pathway_results_human$fdr <= sig_threshold,]
# pathway_score_data_human <- pathway_results_human[,c('pathway', 'pval','fdr','set_enrichment_score')]
pathways <- intersect(pathway_results_human$pathway,pathway_results_mouse_G8G7$pathway)

#select human data
pathway_score_data_human <- pathway_results_human
rownames(pathway_score_data_human) <- pathway_score_data_human$pathway
#choose pathways common to both mouse and human data
pathway_score_data_human <- pathway_score_data_human[pathways,]
colnames(pathway_score_data_human) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')

#select mouse G8G7 data
pathway_score_data_mouse_G8G7 <- pathway_results_mouse_G8G7
rownames(pathway_score_data_mouse_G8G7) <- pathway_results_mouse_G8G7$pathway
pathway_score_data_mouse_G8G7 <- pathway_score_data_mouse_G8G7[pathways,]
colnames(pathway_score_data_mouse_G8G7) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')

#select mouse G10G7 data
pathway_score_data_mouse_G10G7 <- pathway_results_mouse_G10G7
rownames(pathway_score_data_mouse_G10G7) <- pathway_score_data_mouse_G10G7$pathway
pathway_score_data_mouse_G10G7 <- pathway_score_data_mouse_G10G7[pathways,]
colnames(pathway_score_data_mouse_G10G7) <- c('',"pathway","pathway_members","pval","fdr","n","N",'coef')


#order by Humans
# order_sig <- rownames(pathway_score_data_human)[1:number_of_sig_to_keep]
order_sig <- rownames(pathway_score_data_human)
pathway_score_data_human <- pathway_score_data_human[order(pathway_score_data_human$coef, decreasing=F),]
order_sig <- rownames(pathway_score_data_human)[1:number_of_sig_to_keep]
pathway_score_data_human <- pathway_score_data_human[order_sig,]
# order_sig <- rownames(pathway_score_data_human)
# pathway_score_data_human <- within(pathway_score_data_human,
#                                         feature <- factor(pathway,
#                                                           levels=order_sig))
pathway_score_data_human$feature <- factor(pathway_score_data_human$pathway, levels=order_sig)
pathway_score_data_human <- pathway_score_data_human[,c(9,2,3,4,5,6,7,8)]


#grab pathway data from Humans for Mouse Data (G8G7)
pathway_score_data_mouse_G8G7 <- pathway_score_data_mouse_G8G7[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G8G7$feature <- factor(pathway_score_data_mouse_G8G7$pathway, levels=order_sig)
pathway_score_data_mouse_G8G7 <- pathway_score_data_mouse_G8G7[,c(9,2,3,4,5,6,7,8)]
#grab pathway data from Humans for Mouse Data (G10G7)
pathway_score_data_mouse_G10G7 <- pathway_score_data_mouse_G10G7[rownames(pathway_score_data_human),]

pathway_score_data_mouse_G10G7$feature <- factor(pathway_score_data_mouse_G10G7$pathway, levels=order_sig)
pathway_score_data_mouse_G10G7 <- pathway_score_data_mouse_G10G7[,c(9,2,3,4,5,6,7,8)]

pathway_score_data_human_diff_bar <- diff_bar_plot(pathway_score_data_human, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                   fdr ="fdr", orderby = NA, x_label = '', y_label = '')
pathway_score_data_G8G7_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G8G7, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                  fdr ="fdr", orderby = NA, x_label = '', y_label = '')
pathway_score_data_G10G7_diff_bar <- diff_bar_plot(pathway_score_data_mouse_G10G7, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                   fdr ="fdr", orderby = NA, x_label = '', y_label = '')
# ggdraw plots 
fig4c_pathway <- ggdraw() +
  draw_plot(pathway_score_data_human_diff_bar,
            x = 0, y = 0.06, width = .41, height = 0.86) +
  draw_plot(pathway_score_data_G8G7_diff_bar + theme(axis.title.y = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks.y = element_blank(),
                                                     axis.line.y = element_blank()),
            x = .41, y = 0.06, width = .32, height = 0.86) +
  draw_plot(pathway_score_data_G10G7_diff_bar + theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank(),
                                                      axis.line.y = element_blank()),
            x = .73, y = 0.06, width = .26, height = 0.86) +
  draw_plot_label((label = c("Humans", "G8 vs G7", "G10 vs G7")),
                  size = 7,x = c(.25, .5, .75), 
                  y = c(.95, .95, .95)) + 
  draw_plot_label((label = c('Enrichment Score')),
                  size = 7,x = 0.5, 
                  y = 0.05) 
fig4c_pathway


#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/fig4c_pathway.pdf', plot=fig4c_pathway, width = 300, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/fig4c_pathway.png', plot=fig4c_pathway, width = 300, height = 110, units = "mm", dpi = 350)

