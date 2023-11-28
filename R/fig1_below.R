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

g13_12 <- read.delim(
  "analysis/Tweedieverse_MOUSE_HD4_G13_G12/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

g15_14 <- read.delim(
  "analysis/Tweedieverse_MOUSE_HD4_G15_G14/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
)

rownames(g13_12) <- g13_12$feature

# select lowest qvalues
g13_12 = g13_12 %>% 
  filter(!feature %in% c('inosine')) %>%
  top_n(wt = (-qval), 15) %>%
  arrange((abs(coef))) %>% 
  mutate(significant = ifelse(qval<sig_threshold, 'Significant','Insignificant'),
         gr = '2 hours')
g14 = g15_14 %>% 
  filter(feature %in% g13_12$feature) %>% 
  mutate(significant = ifelse(qval<sig_threshold, 'Significant','Insignificant'),
         gr = '12 hours') 

df = g13_12 %>% bind_rows(g14)
df$feature = factor(df$feature, levels = g13_12$feature)


plot_1_below = ggplot(df,
                      aes(x = coef, y = feature,
                          color = factor(significant), 
                          shape = factor(gr))) +
  geom_segment(aes(x=0, xend=coef,
                   y=feature, yend=feature),
               color="grey", size = 0.2) +
  scale_shape_manual(values=c(0, 6))+
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
                              title = 'Time point',
                              override.aes = list(size = 6))
  )+
  theme(legend.title = element_text(size = 6),
        legend.position="bottom",
        title = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text = element_text(size = 8, face = 'bold'),
        legend.text  = element_text(size = 6),
        panel.spacing = unit(0.15, "lines")
  )
plot_1_below

# plots_tw = readRDS('analysis/Tweedieverse_MOUSE_HD4_G13_G12/figures/Group_gg_associations.RDS')
plots_box_g13_g12 = readRDS('analysis/Tweedieverse_MOUSE_HD4_G13_G12/figures/Treatment_gg_associations.RDS')
plots_box_g15_g14 = readRDS('analysis/Tweedieverse_MOUSE_HD4_G15_G14/figures/Treatment_gg_associations.RDS')
legend <- get_legend(plot_1_below)

# Serotonin G15 v. G14
g13_g12_plots_idx = c(14, 10)
g15_g14_plots_idx = c(21, 3)

fig1_below = ggdraw() +
  draw_plot(plot_1_below+theme(legend.position='none'),
            x = -0.02, y = 0.1, width = 0.5, height = .9)+
  draw_plot( legend,scale = 0.1,
             x = 0, y = 0, width = .01, height = .1)+
  draw_plot(plots_box_g13_g12[[g13_g12_plots_idx[1]]]+
                         scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
                         scale_x_discrete(labels=c("G12", "G13"))+
                         theme(
                           axis.title.x = element_text(size = 7),
                           axis.text.x = element_text(size = 7),
                           axis.title.y = element_text(size = 7),
                           axis.text.y = element_text(size = 6)),
            x = 0.48, y = 0.55, width = 0.26, height = .45
  )+
  draw_plot(plots_box_g13_g12[[g13_g12_plots_idx[2]]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G12", "G13"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.74, y = 0.55, width = 0.26, height = .45
  )+
  draw_plot(plots_box_g15_g14[[g15_g14_plots_idx[1]]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G14", "G15"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.48, y = 0.1, width = 0.26, height = .45
  )+
  draw_plot(plots_box_g15_g14[[g15_g14_plots_idx[2]]]+
              scale_fill_manual(values = c('#00e1d9', '#5e001f'))+
              scale_x_discrete(labels=c("G14", "G15"))+
              theme(
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                axis.text.y = element_text(size = 6)),
            x = 0.74, y = 0.1, width = 0.26, height = .45
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

