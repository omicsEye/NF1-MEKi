library(tidyverse)
library(reshape2)
library(deepath)
library(omicsArt)
library(cowplot)
source('~/Documents/omicsEye/omicsArt/R/utils.R')
#setting the working directory
# setwd("~/Box/NF1_MEKi/")
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")
number_of_sig_to_keep <- 35
sig_threshold <- 0.1

g15_13 <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G15_G13/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

g14_12 <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G14_G12/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

rownames(g15_13) <- g15_13$feature

# select lowest qvalues
g15_13 = g15_13 %>% 
  top_n(wt = (-qval), 15) %>%
  arrange((abs(coef))) %>% 
  mutate(significant = ifelse(qval<sig_threshold, 'Significant','Insignificant'),
         gr = 'MEKi')
g14 = g14_12 %>% 
  filter(feature %in% g15_13$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold, 'Significant','Insignificant'),
         gr = 'Vehicle') 

df = g15_13 %>% bind_rows(g14)
df$feature = factor(df$feature, levels = g15_13$feature)


plot_1_below = ggplot(df,
                      aes(x = coef, y = feature,
                          color = significant, 
                          shape = gr))+
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
  ggtitle('Immediate effects of MEKi treatment')+
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
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 10, face = 'bold'),
        legend.text  = element_text(size = 8),
        panel.spacing = unit(0.15, "lines")
  )
plot_1_below

plots_tw = readRDS('analysis/Tweedieverse_MOUSE_CDT_G15_G13/figures/Group_gg_associations.RDS')
legend <- get_legend(plot_1_below)

fig1_below = ggdraw() +
  draw_plot(plot_1_below+theme(legend.position='none'),
            x = -0.02, y = 0.1, width = 0.6, height = .9)+
  draw_plot( legend,
             x = 0.1, y = 0, width = .5, height = .1)+
  draw_plot(plots_tw[[1]]+
                         scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
                         scale_x_discrete(labels=c("G13 (n=11)", "G15 (n=7)"))+
                         theme(
                           axis.title.x = element_text(size = 7),
                           axis.text.x = element_text(size = 7),
                           axis.title.y = element_text(size = 7),
                           axis.text.y = element_text(size = 6)),
            x = 0.6, y = 0.5, width = 0.2, height = .4
  )+
  draw_plot(plots_tw[[2]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G13 (n=11)", "G15 (n=7)"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.8, y = 0.5, width = 0.2, height = .4
  )+
  draw_plot(plots_tw[[3]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G13 (n=11)", "G15 (n=7)"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.6, y = 0.1, width = 0.2, height = .4
  )+
  draw_plot(plots_tw[[4]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G13 (n=11)", "G15 (n=7)"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.8, y = 0.1, width = 0.2, height = .4
  )
# +
#   draw_plot_label((label = c("a","b", "c", "d", "e")),
#                   size = 7,x = c(0, 0, .25, .5, .75),
#                   y = c(1,0.47, 0.47, 0.47, 0.47))

fig1_below

#save files to folder
ggsave(filename = 'manuscript/Figures/Fig1_overview/fig1_below.pdf',
       plot=fig1_below, width = 7.2, height = 4,
       units = "in", dpi = 350)
ggsave(filename = 'manuscript/Figures/Fig1_overview/fig1_below.png',
       plot=fig1_below, width = 7.2, height = 4, 
       units = "in", dpi = 350)

