library(ggplot2)
library(omicsArt)

setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")
mice_meki_Tweedieverse <- read.delim(file='analysis/Tweedieverse_MOUSE_HD4_G6_G4_G5_G3_G1/all_results.tsv',
                   sep = '\t',
                   header = T,
                   fill = F,
                   comment.char = "" ,
                   check.names = F,)
mice_control_Tweedieverse <- read.delim(file='analysis/Tweedieverse_MOUSE_HD4_G11_G9_G10_G8_G7_make_names/all_results.tsv',
                                   sep = '\t',
                                   header = T,
                                   fill = F,
                                   comment.char = "" ,
                                   check.names = F,)
human_Tweedieverse <- read.delim(file='analysis/Tweedieverse_HUMAN_HD4_G17_G16/all_results.tsv',
                                 sep = '\t',
                                 header = T,
                                 fill = F,
                                 comment.char = "" ,
                                 check.names = F,)

p <- ggplot(mice_meki_Tweedieverse,aes(x=coef, color=value))+
      geom_density() +
  theme_set(theme_omicsEye()); p
p1 <- ggplot(mice_control_Tweedieverse,aes(x=coef, color=value))+
  geom_density()+
  theme_set(theme_omicsEye()); p1
p2 <- ggplot(human_Tweedieverse,aes(x=coef, color=value))+
  geom_density()+
  theme_set(theme_omicsEye()); p2

