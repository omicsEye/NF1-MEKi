library(tidyverse)
library(reshape2)
# library(deepath)
library(omicsArt)
library(cowplot)
#source('~/Documents/omicsEye/omicsArt/R/utils.R')
#setting the working directory
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")

number_of_sig_to_keep <- 20
sig_threshold <- 0.1

## use this@!!!!!
#omePath_SUB_PATHWAYS_HUMAN_CLP_CDT_G17_G16
#omePath_SUB_PATHWAYS_HD4_G17_G16

####
# Compare this to mice omePath_SUB_PATHWAYS_HD4_G17_G16
##
human <- read.delim(
  "analysis/enrichment_analysis/tweedie_deepath_SUB_PATHWAYS_HD4/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

## only human data
# human <- read.delim(
#   "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HUMAN_CLP_CDT_G17_G16/enrichment_stats.tsv",
#   sep = '\t',
#   header = T,
#   fill = F,
#   comment.char = "" ,
#   check.names = F,
#   row.names = 1
# )

G9G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G9_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F, row.names = 1
)
G11G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G11_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
G8G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G8_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

G10G7 <- read.delim(
  "analysis/enrichment_analysis/deepath_SUB_PATHWAYS_Mouse_G10_G7/enrichment_stats.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)


human <- human %>% 
  top_n(n = number_of_sig_to_keep,
        wt = set_enrichment_score) %>% 
  arrange(set_enrichment_score)

order_sig = human$pathway



df_main = human %>% 
  mutate(significant = ifelse(fdr<sig_threshold,
                              'Significant','In-significant'),
         gr = 'Human'
  ) %>% 
  bind_rows(
    G9G7%>% 
      filter(pathway %in% order_sig) %>% 
      mutate(significant = ifelse(fdr<sig_threshold,
                                  'Significant','In-significant'),
             gr = 'G9'
      )
  ) %>% 
  bind_rows(
    G11G7%>% 
      filter(pathway %in% order_sig) %>% 
      mutate(significant = ifelse(fdr<sig_threshold,
                                  'Significant','In-significant'),
             gr = 'G11'
      )
  )%>% 
  bind_rows(
    G10G7%>% 
      filter(pathway %in% order_sig) %>% 
      mutate(significant = ifelse(fdr<sig_threshold,
                                  'Significant','In-significant'),
             gr = 'G10'
      )
  )%>% 
  bind_rows(
    G8G7%>% 
      filter(pathway %in% order_sig) %>% 
      mutate(significant = ifelse(fdr<sig_threshold,
                                  'Significant','In-significant'),
             gr = 'G8'
      )
  )

df_main$comparison = ifelse(df_main$gr %in% c('G8', 'G9'), '3 weeks PI',
                            ifelse(df_main$gr=='Human', 'Human', '4-6 weeks PI'))


df_main$comparison = factor(df_main$comparison,
                            levels = c('Human', '3 weeks PI', '4-6 weeks PI'))
df_main$pathway = factor(df_main$pathway, levels = order_sig)
df_main$gr = factor(df_main$gr, 
                    levels = c('Human', 'G8', 'G9', 'G10', 'G11'))

df_main$treatment = ifelse(df_main$gr %in% c('Human','G9','G11'), 
                           'MEKi',
                           'Vehicle')
fig4 = ggplot(df_main,
                     aes(x = set_enrichment_score,
                         y = pathway,
                         color = significant, 
                         shape = treatment))+
  geom_segment(aes(x=0, xend=set_enrichment_score,
                   y=pathway, yend=pathway),
               color="grey", size = 0.2) +
  scale_shape_manual(values=c(1,2))+
  geom_point(aes(size = -log(fdr)))+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Enrichment score')+
  # ggtitle(label = '3 weeks')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'Insignificant'= '#AA9868'))+
  facet_wrap(~comparison, ncol = 3)+
  # facet_grid(.~time_point + comparison)
  omicsArt::theme_omicsEye()+
  guides(size = guide_legend(nrow= 1, title.position = 'top',
                             title = '-Log(pvalue)',order = 1),
         color = guide_legend(nrow= 1, title.position = 'top',
                              title = 'Statistical Significance',
                              override.aes = list(size = 6), order = 2),
         shape = guide_legend(nrow= 1, title.position = 'top',
                              title = 'Treatment',
                              override.aes = list(size = 6), order = 3)
  )+
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.text  = element_text(size = 8),
        panel.spacing = unit(.5, "lines")
  )

fig4
#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/fig4.pdf', 
       plot=fig4, width = 300, height = 115, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/fig4.png', 
       plot=fig4, width = 300, height = 115, units = "mm", dpi = 350)

###
####
new_list = list()
'omePath_SUB_PATHWAYS_HUMAN_CLP_CDT_G17_G16'
list_of_metabolites = unique(df_main$feature)
list_of_metabolites = make.names(paste(list_of_metabolites))
plots_tw = readRDS('analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HUMAN_CLP_CDT_G17_G16/figures/gg_enrichment_rank.RDS')
plot_names = c("PC Ester", "Sphingomyelin", 
               "DAG Ester",
               "TAG Ester")
legend <- get_legend(fig4)

#
## read metabolites
metabolites_score_data_Human <- read.delim(
  "analysis/Tweedieverse_HUMAN_HD4_G17_G16/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)
drop_ind = which(metabolites_score_data_Human$feature %in% c('pregabalin', 'ranitidine'))
metabolites_score_data_Human = metabolites_score_data_Human[-drop_ind,]
rownames(metabolites_score_data_Human) <- metabolites_score_data_Human$feature
metabolites_score_data_Human = metabolites_score_data_Human %>% 
  arrange((abs(coef))) %>% 
  top_n(wt = abs(coef), 15) %>%
  mutate(significant = ifelse(qval<sig_threshold, 'Significant','Insignificant'))
metabolites_score_data_Human$feature = factor(metabolites_score_data_Human$feature, 
                                              levels = metabolites_score_data_Human$feature)

human_meta_fig = ggplot(metabolites_score_data_Human,
       aes(x = coef, y = feature, color = significant))+
  geom_segment( aes(x=0, xend=coef,
                    y=feature, yend=feature),
                color="grey", size = 0.2) +
  geom_point(size =0.9)+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Coefficient')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'Insignificant'= '#AA9868'))+
  guides(color = guide_legend(nrow= 1, title.position = 'top',
                              title = 'Statistical Significance',
                              override.aes = list(size = 6))
  )+
  ggtitle('Human metabolites')+
  omicsArt::theme_omicsEye()+
  theme(legend.title = element_text(size = 8),
        legend.position="None",
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 10, face = 'bold'),
        legend.text  = element_text(size = 8),
        panel.spacing = unit(0.15, "lines")
  )
box_association <- readRDS("analysis/Tweedieverse_HUMAN_HD4_G17_G16/figures/Treatment_gg_associations.RDS")

box_association[[14]]+
  scale_x_discrete(labels=c("16 (n=10)","17 (n=10)"))+
  scale_fill_manual(values = c('#4e79a7', '#f28e2b'))

box_association_clp <- readRDS("analysis/Tweedieverse_HUMAN_CLP_CDT_G17_G16/figures/Treatment_gg_associations.RDS")
box_association_clp[[5]]
box_association_clp[[7]]
box_association_clp[[8]]
box_association_clp[[12]]

box_association[[14]]+
  scale_x_discrete(labels=c("16 (n=10)","17 (n=10)"))+
  scale_fill_manual(values = c('#4e79a7', '#f28e2b'))

###
human_meta_fig_com <- ggdraw() +
  draw_plot(human_meta_fig,
            x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(box_association_clp[[5]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("No treatment (n=10)",
                                        "2 months of treatment (n=10)"))+ 
              scale_x_discrete(labels = NULL)+
              theme(
                axis.title.x = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.y = element_text(size = 6),
                axis.text.y = element_text(size = 5)),
            x = 0.5, y = 0.65, width = .25, height = .35) +
  draw_plot(box_association_clp[[7]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("No treatment (n=10)",
                                        "2 months of treatment (n=10)"))+
              scale_x_discrete(labels = NULL)+
              theme(
                axis.title.x = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.y = element_text(size = 6),
                axis.text.y = element_text(size = 5)), 
            x = .75, y = 0.65, width = .25, height = .35) +
  draw_plot(box_association_clp[[12]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("No treatment (n=10)",
                                        "2 months of treatment (n=10)"))+
              theme(
                axis.title.x = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.y = element_text(size = 6),
                axis.text.y = element_text(size = 5)),
            x = 0.5, y = 0, width = .25, height = .65) +
  draw_plot(box_association_clp[[8]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("No treatment (n=10)",
                                        "2 months of treatment (n=10)"))+
              theme(
                axis.title.x = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.y = element_text(size = 6),
                axis.text.y = element_text(size = 5)), 
            x = .75, y = 0, width = .25, height = .65) 
  # draw_plot_label((label = c("a","b", "c", "d", "e")),
  #                 size = 7,x = c(0, 0, .25, .5, .75),
  #                 y = c(1,0.47, 0.47, 0.47, 0.47))


#




fig4_comb <- ggdraw() +
  draw_plot(fig4+
              theme(legend.position='none'),
            x = 0, y = .1, width = 1, height = .55)+
draw_plot(legend,
          x = 0.1, y = .0, width = .5, height = .1)
for(i in 1:4){
  fig4_comb <- fig4_comb +
    draw_plot(plots_tw[[plot_names[i]]]+
                theme(plot.title = element_text(size = 7),
                  axis.title.x = element_text(size = 6),
                  axis.text.x = element_text(size = 6),
                  axis.title.y = element_text(size = 6),
                  axis.text.y = element_text(size = 6)),
              x = (i-1)/4, y = 0.65, width = .24, height = .33)
  
}
fig4_comb = fig4_comb+
  draw_plot_label((label = c("a","b", "c", "d", "e")),
                  size = 9,x = c(0, 0, .25, .5, .75), 
                  y = c(1,.33, .33, .33, .33))
fig4_comb


comb <- ggdraw() +
  draw_plot(human_meta_fig_com,
            x = 0, y = .6, width = 1, height = .4)+
  draw_plot(fig4_comb,
            x = 0, y = 0, width = 1, height = .6)+
  draw_plot_label((label = c("a","b", "c", "d")),
                  size = 9,x = c(0, 0.5, .0, 0), 
                  y = c(1,1, .6, .4))

#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/comb1.pdf', 
       plot=comb, width = 7.2, height = 9, dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/comb.png', 
       plot=comb, width = 7.2, height = 9.5, dpi = 350)

#save files to folder
ggsave(filename = 'manuscript/Figures/figure_4/fig4_comb.pdf', 
       plot=fig4_comb, width = 7.2, height = 6, dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_4/fig4_comb.png', 
       plot=fig4_comb, width = 7.2, height = 6, dpi = 350)


# order_sig <- rownames(pathway_score_data_human)
# pathway_score_data_human <- within(pathway_score_data_human,
#                                         feature <- factor(pathway,
#                                                           levels=order_sig))
pathway_score_data_human$feature <- factor(pathway_score_data_human$pathway,
                                           levels=order_sig)
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

