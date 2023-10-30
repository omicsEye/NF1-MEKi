library(tidyr)
library(dplyr)
library(reshape2)
library(deepath)
library(omicsArt)
library(cowplot)
source('~/Documents/omicsEye/omicsArt/R/utils.R')
#setting the working directory
# setwd("~/Box/NF1_MEKi/")
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")
number_of_sig_to_keep <- 20
sig_threshold <- 0.1

## read metabolites
metabolites_Tweedieverse <- read.delim(
  "analysis/Maaslin2_individualcomparisons/Maaslin2_HUMAN_HD4_G17_G16/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_score_data_Human <- metabolites_Tweedieverse #[metabolites_Tweedieverse$metadata=="Treatment" & 
                                                        #   metabolites_Tweedieverse$value== 17 & 
                                                         #  metabolites_Tweedieverse$percent.zero,]
rownames(metabolites_score_data_Human) <- metabolites_score_data_Human$feature

## read metabolites
metabolites_Tweedieverse <- read.delim(
  "analysis/Maaslin2_individualcomparisons/Maaslin2_MOUSE_CDT_G15_G14/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_Tweedieverse_Mice_12hrs <- metabolites_Tweedieverse#[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="15" ,]
rownames(metabolites_Tweedieverse_Mice_12hrs) <- metabolites_Tweedieverse_Mice_12hrs$feature

## read metabolites
metabolites_Tweedieverse <- read.delim(
  "analysis/Maaslin2_individualcomparisons/Maaslin2_MOUSE_CDT_G13_G12/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_Tweedieverse_Mice_2hrs <- metabolites_Tweedieverse#[metabolites_Tweedieverse$metadata=="Group" & metabolites_Tweedieverse$value=="13" ,]
rownames(metabolites_Tweedieverse_Mice_2hrs) <- metabolites_Tweedieverse_Mice_2hrs$feature



# use score_data_severe is reference
metabolites_score_data_Human <- metabolites_score_data_Human[intersect(rownames(metabolites_score_data_Human),
                                                                       intersect(rownames(metabolites_Tweedieverse_Mice_12hrs),
                                                                                 rownames(metabolites_Tweedieverse_Mice_2hrs))),]
order_sig <- rownames(metabolites_score_data_Human)[1:number_of_sig_to_keep]
metabolites_score_data_Human <- metabolites_score_data_Human[order_sig,]
metabolites_score_data_Human<- metabolites_score_data_Human[order(metabolites_score_data_Human$coef, decreasing=F),]
order_sig <- rownames(metabolites_score_data_Human)
metabolites_score_data_Human <- within(metabolites_score_data_Human,
                                        feature <- factor(feature,
                                                          levels=order_sig))
metabolites_Tweedieverse_Mice_12hrs <- metabolites_Tweedieverse_Mice_12hrs[rownames(metabolites_score_data_Human),]
metabolites_Tweedieverse_Mice_12hrs <- within(metabolites_Tweedieverse_Mice_12hrs,
                                            feature <- factor(feature,
                                                              levels=order_sig))

#metabolites_Tweedieverse_Mice_12hrs[is.na(metabolites_Tweedieverse_Mice_12hrs)] <- 0 
#metabolites_Tweedieverse_Mice_12hrs$feature <- rownames(metabolites_Tweedieverse_Mice_12hrs)

metabolites_Tweedieverse_Mice_2hrs <- metabolites_Tweedieverse_Mice_2hrs[rownames(metabolites_score_data_Human),]
metabolites_Tweedieverse_Mice_2hrs <- within(metabolites_Tweedieverse_Mice_2hrs,
                                           feature <- factor(feature,
                                                             levels=order_sig))
#metabolites_Tweedieverse_Mice_2hrs[is.na(metabolites_Tweedieverse_Mice_2hrs)] <- 0 
#metabolites_Tweedieverse_Mice_2hrs$feature <- rownames(metabolites_Tweedieverse_Mice_2hrs) 

metabolites_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_Human,
                                                  threshold = sig_threshold,
                                                  pvalue_col = "pval", 
                                                  method = "none",
                                                  fdr ="qval", orderby = NA,
                                                  x_label = 'Coefficient',
                                                  y_label = '')




metabolites_non_severe_temp_diff_bar <- diff_bar_plot(metabolites_Tweedieverse_Mice_12hrs, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                      fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')
metabolites_non_covid_temp_diff_bar <- diff_bar_plot(metabolites_Tweedieverse_Mice_2hrs, threshold = sig_threshold, pvalue_col = "pval",  method = "none",
                                                     fdr ="qval", orderby = NA, x_label = 'Coefficient', y_label = '')

df_main = metabolites_score_data_Human %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'Human') 
m_new = metabolites_Tweedieverse_Mice_12hrs %>% 
  mutate(significant = ifelse(qval<sig_threshold,
                              'Significant','In-significant'),
         gr = 'Mice 12 hours')
df_main = df_main %>% rbind(m_new) %>% 
  bind_rows(metabolites_Tweedieverse_Mice_2hrs %>% 
              mutate(significant = ifelse(qval<sig_threshold,
                                          'Significant','In-significant'),
                     gr = 'Mice 2 hours')
  )

top_plot = ggplot(df_main,
       aes(x = coef, y = feature, color = significant))+
  geom_segment( aes(x=0, xend=coef,
                    y=feature, yend=feature),
                color="grey", size = 0.2) +
  geom_point(size =0.9)+
  geom_vline(xintercept = 0, 
             color = "gray", size=0.3)+
  ylab('')+xlab('Coefficient')+
  scale_color_manual(values = c('Significant'= '#033C5A',
                                'In-significant'= '#AA9868'))+
  facet_grid(~factor(gr, 
                     levels=c('Human', 
                              'Mice 2 hours',
                              'Mice 12 hours')))+
  omicsArt::theme_omicsEye()+
  theme(strip.text.x = element_text(size = 6,
                                    face = 'bold'),
        legend.title = element_blank())
## read association
box_association <- readRDS("analysis/Tweedieverse_HUMAN_HD4_G17_G16/figures/Treatment_gg_associations.RDS")
## do plots

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

box_association[[14]]+
  scale_x_discrete(labels=c("16 (n=10)","17 (n=10)"))+
  scale_fill_manual(values = c('#4e79a7', '#f28e2b'))
###
fig2_metabolites_n <- ggdraw() +
  draw_plot(top_plot,
            x = 0, y = .47, width = 1, height = .53) +
  draw_plot(box_association[[14]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("G16 (n=10)","G17 (n=10)"))+ 
              theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = 0, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[28]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("G16 (n=10)","G17 (n=10)"))+ 
              theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), 
    x = .25, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[83]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("G16 (n=10)","G17 (n=10)"))+
              theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .5, y = 0, width = .25, height = .45) +
  draw_plot(box_association[[37]] +
              scale_fill_manual(values = c('#4e79a7', '#f28e2b'))+
              scale_x_discrete(labels=c("G16 (n=10)","G17 (n=10)"))+
              theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), 
    x = .75, y = 0, width = .25, height = .45) +
  
  draw_plot_label((label = c("a","b", "c", "d", "e")),
                  size = 7,x = c(0, 0, .25, .5, .75), y = c(1,0.47, 0.47, 0.47, 0.47))

###
fig2_metabolites

ggsave(filename = './manuscript/Figures/Figure_2/fig2_metabolites_n.pdf',
       plot=fig2_metabolites_n, width = 183, height = 110, units = "mm", dpi = 350)
ggsave(filename = './manuscript/Figures/Figure_2/fig2_metabolites_n.png',
       plot=fig2_metabolites_n, width = 183, height = 110, units = "mm", dpi = 350)

