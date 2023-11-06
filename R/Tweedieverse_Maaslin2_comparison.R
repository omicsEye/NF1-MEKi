library(tidyverse)



qval_threshold <- 0.1
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")

tweedie_outputs = list.files(path = "./analysis/", pattern = 'Tweedieverse')


for(i in 1:length(tweedi_outputs)){
  analysis = unlist(strsplit(tweedi_outputs[i], 'Tweedieverse_'))[2]
  if(dir.exists(paste('analysis/Maaslin2_', analysis, sep = ""))){
    tweedie_dir = paste('analysis/Tweedieverse_', analysis, sep = "")
    maaslin_dir = paste('analysis/Maaslin2_', analysis, sep = "")
    
  
print(tweedie_dir)
print(maaslin_dir)

maaslin_result = read.delim(paste(maaslin_dir, '/all_results.tsv', sep = ''),
                            sep = '\t', header = T, fill = F, 
                            comment.char = "" , check.names = F)


tweedie_result = read.delim(paste(tweedie_dir, '/all_results.tsv', sep = ''),
                            sep = '\t', header = T, fill = F, 
                            comment.char = "" , check.names = F
                            )
tweedie_result$feature = make.names(tweedie_result$feature)


# combining the results of the two packages
results = maaslin_result %>%
  select(feature, metadata,value, coef, qval) %>%
  rename('coef_mas' = 'coef',
         'qval_mas' = 'qval') %>%
  inner_join(
    tweedie_result %>%
      select(feature, metadata,value, coef, qval) %>%
      rename('coef_tw' = 'coef',
             'qval_tw' = 'qval')
  ) %>%
  mutate(group =
           ifelse(qval_mas<=qval_threshold & qval_tw<=qval_threshold,"Both",
                  ifelse(qval_mas<=qval_threshold & qval_tw>qval_threshold,"Maaslin2 only",
                         ifelse(qval_mas>qval_threshold & qval_tw<=qval_threshold,"Tweedieverse only",
                                'None')))
  )



cr = cor(-log(results$qval_mas), 
         -log(results$qval_tw), 
         use = 'pairwise.complete.obs')
cbp <- c('#e69f00', '#0072b2', '#009e73','#808080')
names(cbp) = c("Both", "Maaslin2 only", "Tweedieverse only", 'None')
txt_y = max(-log(results$qval_tw))*.8
txt_x = max(-log(results$qval_mas))*.05

p = ggplot(results, aes(x = -log(qval_mas), y = -log(qval_tw)))+
  geom_point(alpha = 0.7, size = 1, aes(color = group))+
  xlab(bquote(-log[e](`Maaslin2 qvalue`)))+
  ylab(bquote(-log[e](`Tweedieverse qvalue`)))+
  geom_smooth(method='lm', formula= y~x, size = .5)+
  ggtitle(label = analysis,
          subtitle = paste('Correlation =', round(cr,2)))+
  # annotate("text", txt_x, txt_y, size = 2,
  #          label = paste('Correlation =', round(cr,2)))+
  scale_color_manual(values = cbp, drop = TRUE)+
  omicsArt::theme_omicsEye()

ggsave(filename = paste('analysis/Tweedieverse_Maaslin2_comparison/',analysis,
                        '.pdf', sep = ''),device = 'pdf',
                        p, width = 3, height = 3, units = 'in')

  }
}
