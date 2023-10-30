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

# read the plots and then pass to channel!
# filter based on the metabolites in figure 3!
# plots_tw = readRDS('analysis/Tweedieverse_MOUSE_HD4_G11_G9_G10_G8_G7_make_names/figures/Treatment_gg_associations.RDS')
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



# prepare data for visualization

g11 = metabolites_score_data_G11 %>% 
  arrange(desc(abs(coef))) %>% 
  top_n(wt = abs(coef), 20) %>%
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G11') 

g10 = metabolites_score_data_G10 %>% 
  filter(feature %in% g11$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G10') 
g5 = metabolites_score_data_G5 %>% 
  filter(feature %in% g11$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G5') 

g6 = metabolites_score_data_G6 %>% 
  filter(feature %in% g11$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G6') 

df_4_6 = g11 %>% 
  bind_rows(g10) %>% 
  bind_rows(g5) %>% 
  bind_rows(g6)


df_4_6$feature = factor(df_4_6$feature, levels = rev(g11$feature))
plot_4_6 = ggplot(df_4_6,
                  aes(x = coef, y = feature, color = significant, shape = gr))+
  geom_segment(aes(x=0, xend=coef,
                    y=feature, yend=feature),
                color="grey", size = 0.2) +
  geom_point(size =0.9)+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Coefficient')+
  ggtitle(label = '4-6 weeks')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'In-significant'= '#AA9868'))+
  omicsArt::theme_omicsEye()+
  theme(legend.title = element_blank(),
        legend.position="bottom")+
  guides(color='none')
plot_4_6


g9 = metabolites_score_data_G9 %>% 
  arrange(desc(abs(coef))) %>% 
  top_n(wt = abs(coef), 20) %>%
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G9') 

g8 = metabolites_score_data_G8 %>% 
  filter(feature %in% g9$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G8') 
g3 = metabolites_score_data_G3 %>% 
  filter(feature %in% g9$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G3') 

g4 = metabolites_score_data_G4 %>% 
  filter(feature %in% g9$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'G4') 

df_3 = g9 %>% 
  bind_rows(g8) %>% 
  bind_rows(g3) %>% 
  bind_rows(g4)


df_3$feature = factor(df_3$feature, levels = rev(g9$feature))
plot_3 = ggplot(df_3,
                  aes(x = coef, y = feature, color = significant, shape = gr))+
  geom_segment(aes(x=0, xend=coef,
                   y=feature, yend=feature),
               color="grey", size = 0.2) +
  scale_shape_manual(values=c(4:7))+
  geom_point(size =0.9)+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Coefficient')+
  ggtitle(label = '3 weeks')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'In-significant'= '#AA9868'))+
  omicsArt::theme_omicsEye()+
  theme(legend.title = element_blank(),
        legend.position="bottom")
plot_3

fig3_metabolites = plot_grid(plot_3, plot_4_6,
                                     labels = c('a', 'b'))

#save files to folder
ggsave(filename = 'manuscript/Figures/figure_3/fig3_n.pdf', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_3/fig3_n.png', plot=fig3_metabolites, width = 300, height = 110, units = "mm", dpi = 350)


### paired comparisons
for(i in c(3, 4, 5, 6, 8, 9, 10, 11)){
  temp_df = eval(parse(text = paste0('metabolites_score_data_G', i)))
  temp_df = temp_df %>%
    filter(feature !="dimethyl sulfone") 
    
  if(i == 3){
    tmp1 = temp_df[order(abs(temp_df$coef), decreasing = TRUE), 'feature'][1:5]
    tmp2 = temp_df[temp_df$qval<sig_threshold, 'feature']
    if(!is.null(nrow(tmp2))){
      if(nrow(tmp2)<10){
        tmp = tmp2
      }else{
        tmp = tmp2[1:10]
      }
    }else{
      tmp = tmp1
    }
    
  }else{
    tmp1 = temp_df[order(abs(temp_df$coef), decreasing = TRUE), 'feature'][1:5]
    tmp2 = temp_df[temp_df$qval<sig_threshold, 'feature']
    
    if(!is.null(nrow(tmp2))){
      if(nrow(tmp2)<10){
        tmp = union(tmp, tmp2)
      }else{
        tmp = union(tmp, tmp2[1:10])
      }
    }else{
      tmp = union(tmp, tmp1)
    }
  }
}
must_add_metabolites = c('methionine sulfone',
                         'gamma-glutamyl-epsilon-lysine',
                         'galactonate',
                         'picolinate')
for(metabolite in must_add_metabolites){
  if(!metabolite %in% tmp){
    tmp = c(tmp, metabolite)
  }
}

for(i in c(3, 4, 5, 6, 8, 9, 10, 11)){
  temp_df = eval(parse(text = paste0('metabolites_score_data_G', i)))
  if(i == 3){
    df_main = temp_df %>% 
      filter(feature %in% tmp) %>%
      mutate(significant = ifelse(qval<sig_threshold,
                                  'Significant','In-significant'),
             gr = paste0('G',i),
             # comparison = ifelse(i %in% c(3, 4), 'G3 and G4',
             #                     ifelse(i %in% c(8, 9), 'G8 and G9',
             #                            ifelse(i %in% c(5, 6), 'G5 and G6',
             #                                   'G10 and G11'))))
             comparison = ifelse(i %in% c(3, 4), 'Three_w_control',
                                 ifelse(i %in% c(8, 9), 'Three_w_NF1',
                                        ifelse(i %in% c(5, 6), 'four_w_control',
                                               'four_w_NF1'))))
  }else{
    df_main = df_main %>% 
      bind_rows(
        temp_df %>% 
          filter(feature %in% tmp) %>%
          mutate(significant = ifelse(qval<sig_threshold,
                                      'Significant','Insignificant'),
                 gr = paste0('G',i),
                 # comparison = ifelse(i %in% c(3, 4), 'G3 and G4',
                 #                     ifelse(i %in% c(8, 9), 'G8 and G9',
                 #                            ifelse(i %in% c(5, 6), 'G5 and G6',
                 #                                   'G10 and G11'))))
                 comparison = ifelse(i %in% c(3, 4), 'Three_w_control',
                                     ifelse(i %in% c(8, 9), 'Three_w_NF1',
                                            ifelse(i %in% c(5, 6), 'four_w_control',
                                                   'four_w_NF1'))))
      )
  }
}

# df_main$comparison = factor(df_main$comparison,
#                             levels = c('G3 and G4', 'G8 and G9',
#                                        'G5 and G6', 'G10 and G11'))
df_main$comparison = factor(df_main$comparison,
                            levels = c('Three_w_control', 'Three_w_NF1',
                                       'four_w_control', 'four_w_NF1'))

df_main$treatment = ifelse(df_main$gr %in% c('G4','G9', 'G6', 'G11'), 
                       'MEKi',
                       'Vehicle')
# confirm the order of metabolites based on G3
fact_levels = metabolites_score_data_G3$feature[order(abs(metabolites_score_data_G3$coef),
                                                      decreasing = T)]
fact_levels = fact_levels[fact_levels %in% tmp]

change_name = "3-carboxy-4-methyl-5-pentyl-2-furanpropionate (3-CMPFP)"
fact_levels[fact_levels==change_name] = "3-CMPFP"

df_main$feature[df_main$feature == change_name] = "3-CMPFP"
df_main$feature = factor(df_main$feature, levels = fact_levels)
df_main$gr = factor(df_main$gr, 
                    levels = c('G3', 'G4', 'G5',
                               'G6', 'G7', 'G8',
                               'G9', 'G10', 'G11'))
fig3_4plots = ggplot(df_main,
       aes(x = coef, y = feature,
           color = significant, 
           shape = treatment))+
  geom_segment(aes(x=0, xend=coef,
                   y=feature, yend=feature),
               color="grey", size = 0.2) +
  scale_shape_manual(values=c(1,2))+
  geom_point(aes(size = -log(qval)))+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Coefficient')+
  # ggtitle(label = '3 weeks')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'Insignificant'= '#AA9868'))+
  facet_wrap(~comparison, ncol = 4)+
  # facet_grid(.~time_point + comparison)
  omicsArt::theme_omicsEye()+
  guides(size = guide_legend(nrow= 1, title.position = 'top',
                              title = '-Log(pvalue)'),
         color = guide_legend(nrow= 1, title.position = 'top',
                              title = 'Statistical Significance',
                              override.aes = list(size = 6)),
         shape = guide_legend(nrow= 1, title.position = 'top',
                             title = 'Treatment',
                              override.aes = list(size = 6))
         )+
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 10, face = 'bold'),
        legend.text  = element_text(size = 8),
        panel.spacing = unit(0.15, "lines")
        )

fig3_4plots
#save files to folder
ggsave(filename = 'manuscript/Figures/figure_3/fig3_4plots.pdf', 
       plot=fig3_4plots, width = 300, height = 115, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/Figures/figure_3/fig3_4plots.png', 
       plot=fig3_4plots, width = 300, height = 115, units = "mm", dpi = 350)


####
new_list = list()
list_of_metabolites = unique(df_main$feature)
list_of_metabolites = make.names(paste(list_of_metabolites))
plots_tw = readRDS('analysis/Tweedieverse_MOUSE_HD4_G11_G9_G10_G8_G7_make_names/figures/Treatment_gg_associations.RDS')
for(i in 1:length(plots_tw)){
  y_axis = plots_tw[[i]]$labels$y
  if(y_axis %in% list_of_metabolites){
    new_list[[y_axis]] = plots_tw[[i]]
  }
}


pdf("manuscript/Figures/figure_3/top_plots.pdf",onefile = TRUE)
for(i in 1:length(new_list)){
  print(new_list[[i]])
}
dev.off()
saveRDS(new_list, 'analysis/Tweedieverse_MOUSE_HD4_G11_G9_G10_G8_G7_make_names/figures/top_plots.RDS')

metabolite_list = c('methionine.sulfone', 'picolinate',
                    'galactonate', 'hypoxanthine')

# The palette with grey:
cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                        "#F0E442", "#0072B2", 
                        "#D55E00", "#CC79A7", "#999999")

library(cowplot)                       

fig3_comb <- ggdraw() +
  draw_plot(fig3_4plots,
            x = 0, y = .35, width = 1, height = .65)
for(i in 1:4){
  fig3_comb <- fig3_comb +
    draw_plot(new_list[[metabolite_list[i]]]+
    scale_fill_manual(values = cbPalette[1:5])+
    scale_x_discrete(labels=c("G7 (n=6)", "G11 (n=7)",
                              "G9 (n=5)","G8 (n=6)",
                              "G10 (n=6)"))+
    theme(
      axis.title.x = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text.y = element_text(size = 6)),
    x = (i-1)/4, y = 0, width = .25, height = .33)
    
}

fig3_comb = fig3_comb+
  draw_plot_label((label = c("a","b", "c", "d", "e")),
                size = 9,x = c(0, 0, .25, .5, .75), 
                y = c(1,.33, .33, .33, .33))

#save files to folder
ggsave(filename = 'manuscript/Figures/figure_3/fig3_comb.pdf', 
       plot=fig3_comb, width = 7.2, height = 6, dpi = 300)
ggsave(filename = 'manuscript/Figures/figure_3/fig3_comb.png', 
       plot=fig3_comb, width = 7.2, height = 6,dpi = 300)

