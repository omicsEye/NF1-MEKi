library(tidyr)
library(dplyr)
library(reshape2)
library(deepath)
library(omicsArt)
library(ggExtra)
library(ggrepel)

source('~/Documents/omicsEye/omicsArt/R/utils.R')
#setting the working directory
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi/")

sig_threshold <- 0.1

## read metabolites
metabolites_Tweedieverse <- read.delim(
  "analysis/Tweedieverse_HUMAN_HD4_G17_G16/all_results.tsv",
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
  "analysis/Tweedieverse_MOUSE_CDT_G15_G14/all_results.tsv",
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
  "analysis/Tweedieverse_MOUSE_CDT_G13_G12/all_results.tsv",
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
number_of_sig_to_keep <- dim(metabolites_score_data_Human)[1]
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

metabolites_Tweedieverse_Mice_2hrs <- metabolites_Tweedieverse_Mice_2hrs[rownames(metabolites_score_data_Human),]
metabolites_Tweedieverse_Mice_2hrs <- within(metabolites_Tweedieverse_Mice_2hrs,
                                           feature <- factor(feature,
                                                             levels=order_sig))
#colnames(metabolites_score_data_Human) <- paste0("Human ", colnames(metabolites_score_data_Human))
colnames(metabolites_Tweedieverse_Mice_2hrs) <- paste0("Mice 2 Hrs ", colnames(metabolites_Tweedieverse_Mice_2hrs))
colnames(metabolites_Tweedieverse_Mice_12hrs) <- paste0("Mice 12 Hrs ", colnames(metabolites_Tweedieverse_Mice_12hrs))
combined_data <- cbind(metabolites_score_data_Human, metabolites_Tweedieverse_Mice_2hrs, metabolites_Tweedieverse_Mice_12hrs)
upper_lim <- max(combined_data$coef, na.rm = T) 
lower_lim <- min(combined_data$coef, na.rm = T) 
scatter_plot_1 <- ggplot(data=combined_data,aes(`Mice 2 Hrs coef`,
                                                `Mice 12 Hrs coef`,
                                                label = feature)) +
  xlab("Effect size of mice 2 hourse of MEKi treatment") + 
  ylab("Effect size of mice 12 hourse of MEKi treatment")+
  stat_density2d( aes(fill = 'grey', alpha = .75), geom='polygon')+
  geom_point( aes(), fill = 'darkolivegreen4', color = 'darkolivegreen4',
              alpha = .5, shape = 21, size = 1, stroke = 0.2) +
  stat_smooth(method = "glm", col="#FDE725FF", fill="#FDE725FF", size =0.5)+
  xlim(lower_lim, upper_lim) +
  ylim(lower_lim, upper_lim)+
  ggplot2::geom_vline(xintercept=0, col='red', size = .1) +
  ggplot2::geom_hline(yintercept = 0, color = "red",size = 0.1) +
  guides(alpha='none', fill='none')+labs("")+
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(size = 1.5, force = 1,  box.padding = unit(.1, "lines"),
                   point.padding = unit(0.05, "lines"))+
  theme_omicsEye()
scatter_plot_1_1 <- ggMarginal(scatter_plot_1, xparams = list(fill = "#481568FF", col = "#440154FF", alpha=.3), #margins = 'x',
                               yparams = list(fill = "#29AF7FFF", col = "#20A386FF", alpha=.3))
scatter_plot_1_1
scatter_plot_2 <- ggplot(data=combined_data,aes(`coef`, `Mice 2 Hrs coef`, label = feature)) +
  xlab("Effect size of human MEKi treatment") + ylab("Effect size of mice 2 hourse of MEKi treatment") +
  stat_density2d( aes(fill = 'grey', alpha = .75), geom='polygon')+
  geom_point( aes(), fill = 'darkolivegreen4', color = 'darkolivegreen4', alpha = .5, shape = 21, size = 1, stroke = 0.2) +
  stat_smooth(method = "glm", col="#FDE725FF", fill="#FDE725FF", size =0.5)+
  xlim(lower_lim, upper_lim) +
  ylim(lower_lim, upper_lim)+
  ggplot2::geom_vline(xintercept=0, col='red', size = .1) +
  ggplot2::geom_hline(yintercept = 0, color = "red",size = 0.1) +
  guides(alpha='none', fill='none')+labs("")+
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(size = 1.5, force = 1,  box.padding = unit(.1, "lines"),
                  point.padding = unit(0.05, "lines"))+
  theme_omicsEye()
scatter_plot_2_2 <- ggMarginal(scatter_plot_2, xparams = list(fill = "#481568FF", col = "#440154FF", alpha=.3), #margins = 'x',
                               yparams = list(fill = "#29AF7FFF", col = "#20A386FF", alpha=.3))

scatter_plot_3 <- ggplot(data=combined_data,aes(`coef`, `Mice 12 Hrs coef`, label = feature)) +
  xlab("Effect size of human MEKi treatment") + ylab("Effect size of mice 12 hourse of MEKi treatment") +
  stat_density2d( aes(fill = 'grey', alpha = .75), geom='polygon')+
  geom_point( aes(), fill = 'darkolivegreen4', color = 'darkolivegreen4', alpha = .5, shape = 21, size = 1, stroke = 0.2) +
  stat_smooth(method = "glm", col="#FDE725FF", fill="#FDE725FF", size =0.5)+
  xlim(lower_lim, upper_lim ) +
  ylim(lower_lim, upper_lim)+
  ggplot2::geom_vline(xintercept=0, col='red', size = .1) +
  ggplot2::geom_hline(yintercept = 0, color = "red",size = 0.1) +
  guides(alpha='none', fill='none')+labs("")+
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(size = 1.5, force = 1.5,  box.padding = unit(.1, "lines"),
                  point.padding = unit(0.05, "lines"))+
  theme_omicsEye()
scatter_plot_3_3 <- ggMarginal(scatter_plot_3, xparams = list(fill = "#481568FF", col = "#440154FF", alpha=.3), #margins = 'x',
                               yparams = list(fill = "#29AF7FFF", col = "#20A386FF", alpha=.3))


fig_scatters <-   plot_grid(scatter_plot_1_1 ,
                            scatter_plot_2_2 ,
                            scatter_plot_3_3,
                            #rel_widths = c(1, 1, 1,1),
                            #labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', "i"),
                            label_size = 7, ncol = 3, nrow = 1)
fig_scatters
ggsave(filename='fig_scatters_manuscript_format.pdf',
       plot=fig_scatters, width = 7.2, height = 2.5, units = "in", dpi = 350)
ggsave(filename='fig_scatters_presentation_format.pdf',
       plot=fig_scatters, width = 15, height = 5, units = "in", dpi = 350)



#####
# clustering

km_dat = combined_data[, c("coef", "Mice 2 Hrs coef", "Mice 12 Hrs coef")]

km_mat = abs(km_dat)
km_mat = apply(km_mat, 2, FUN = function(x){(x- min(x))/(max(x)-min(x))})

wss = c()
for(i in 2:30){
km = kmeans(km_mat, centers = i, iter.max = 500, nstart = 50)
wss = append(wss, km$tot.withinss)
}
plot(2:30, wss, type = 'b')

km = kmeans(km_mat, centers = 6, iter.max = 500, nstart = 50)

km$centers

library(plotly)
km_dat_plot = km_dat
colnames(km_dat_plot) = c('human', 'mice2', 'mice12')
km_dat_plot$cluster = factor(km$cluster)
fig <- plot_ly(km_dat_plot,
               x = ~human, 
               y = ~mice2,
               z = ~mice12,
               color = ~cluster
)
fig
