library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(Maaslin2)
library(readxl)
library(omePath)
library(omicsArt)
library(Tweedieverse)
fdr_threshold <- 1.00
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")

# Mice_Row2
loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata$TREATMENT <-  as.character(metadata$TREATMENT)
pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_MOUSE CDT', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_MOUSE CDT', method = 'tsne')

metadata_2 <- metadata[metadata$`BOX NUMBER` == 'Box2',]
data2 <- data[rownames(metadata_2),]
pcoa_plots <- ordplots(data2, metadata_2, output = 'analysis/', outputname = 'pcoa_MOUSE CDT Box1', method = 'pcoa')

#tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_MOUSE CDT', method = 'tsne')

#############################################################################################################
##############
## Combined ##
##############

##################
## Control Mice ##
##################
# Group 6 and 4, 5 and 3 vs 1
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 4, 6,3,5), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==3] <- "2_3 weeks of Vehicle"
metadata_2[metadata_2==5] <- "3_4 - 6 weeks of Vehicle"
metadata_2[metadata_2==4] <- "4_3 weeks of MEKi"
metadata_2[metadata_2==6] <- "5_4 - 6 weeks of MEKi"
data2 <- data[rownames(metadata_2),]
pcoa_plots <- ordplots(data2, metadata_2, output = 'analysis/', outputname = 'pcoa_MOUSE CDT Control', method = 'pcoa')
tsne_plots <- ordplots(data2, metadata_2, output = 'analysis/', outputname = 'tsne_MOUSE CDT Control', method = 'tsne')
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_HD4_G6_G4_G5_G3_G1',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,1_Baseline_of_HGFAP'
)
# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_MOUSE_CDT_G6_G4_G5_G3_G1',
#          reference='1_Baseline_of_HGFAP',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')

##############
## NFP Mice ##
##############
# Group 11 and 9, 10 and 8 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(7, 8,9,10,11), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==7] <- "NF1_Baseline_of_HGFAP" #Nf1-/- and G_2 is Nf1+/-
metadata_2[metadata_2==8] <- "3 weeks of Vehicle"
metadata_2[metadata_2==10] <- "4 - 6 weeks of Vehicle"
metadata_2[metadata_2==9] <- "3 weeks of MEKi"
metadata_2[metadata_2==11] <- "4 - 6 weeks of MEKi"
levels(metadata_2$Treatment) <- c("Baseline_of_HGFAP", "3_weeks_of_Vehicle", "4_6 weeks_of_Vehicle", "3_weeks_of_MEKi", "4_6 weeks_of_MEKi")
colnames(data) <- make.names(colnames(data))
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           '~/Downloads/Tweedieverse_MOUSE_HD4_G11_G9_G10_G8_G7_make_names',
                           max_significance = 0.2,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,Baseline_of_HGFAP'
)
# Group 11 and 9, 10 and 8 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(9,11), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==9] <- "3_weeks_of_MEKi"
metadata_2[metadata_2==11] <- "4 - 6 weeks of MEKi"
levels(metadata_2$Treatment) <- c("3_weeks_of_MEKi", "4_6 weeks_of_MEKi")
#colnames(data) <- make.names(colnames(data))
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           '~/Downloads/Tweedieverse_MOUSE_HD4_G11_G9',
                           max_significance = 0.2,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,3_weeks_of_MEKi'
)
# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_MOUSE_CDT_G11_G9_G10_G8_G7',
#          reference='Treatment:1',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')


###############################
## Heterozygous - Control Mice ##
###############################
# All Control mice vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 4, 6,3,5), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==3] <- "2_3 weeks of Vehicle"
metadata_2[metadata_2==5] <- "3_4 - 6 weeks of Vehicle"
metadata_2[metadata_2==4] <- "4_3 weeks of MEKi"
metadata_2[metadata_2==6] <- "5_4 - 6 weeks of MEKi"
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G2v_',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,1_Baseline_of_HGFAP'
)
# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_MOUSE_CDT_G6_G4_G5_G3_G2',
#          reference='1_Baseline_of_HGFAP',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')
#All NF1 Mice vs Group 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2,8,9,10,11), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==8] <- "2_3 weeks of Vehicle"
metadata_2[metadata_2==10] <- "3_4 - 6 weeks of Vehicle"
metadata_2[metadata_2==9] <- "4_3 weeks of MEKi"
metadata_2[metadata_2==11] <- "5_4 - 6 weeks of MEKi"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G4_G1',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G11_G9_G10_G8_G2',
         reference='Treatment:1',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

##################
## Control Mice ##
##################
# Group 1 vs 4
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 4), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==4] <- "3 weeks of treatment"
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G4_G1',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,1_Baseline_of_HGFAP'
)
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G4_G1',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

#################
#Group 1 vs 6
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 6), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==6] <- "4-6 weeks post-treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G6_G1',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G6_G1',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
#################
#Group 3 vs 1
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 3), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==3] <- "3 weeks post-treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G11_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G3_G1',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
#################
#Group 5 vs 1
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 5), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==5] <- "4-6 weeks post-treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G5_G1',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G5_G1',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

#############################################################################################################
##############
## NF1 Mice ##
##############
# Group 9 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(7, 9), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==7] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==9] <- "3 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G9_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G9_G7',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
##############
# Group 11 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(7, 11), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==7] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==11] <- "4-6 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G11_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G11_G7',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

# Group 11 vs 9
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(9, 11), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==9] <- "3 weeks of treatment"
metadata_2[metadata_2==11] <- "4-6 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G11_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G11_G9',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

#################
#Group 8 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(7, 8), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==7] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==8] <- "3 weeks post-treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G11_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G8_G7',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
#################
#Group 10 vs 7
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(7, 10), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==7] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==10] <- "4-6 weeks post-treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G10_G7',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G10_G7',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

#############################################################################################################
#################################################
## Comparisons to Heterozygous Mice (Group 2) ##
#################################################
# Group 4 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 4), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==4] <- "3 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G4_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G4_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 3 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 3), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==3] <- "3 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G3_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G3_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 9 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 9), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==9] <- "3 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G9_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G9_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 8 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 8), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==8] <- "3 weeks of treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G8_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = 'Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G8_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 6 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 6), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==6] <- "4-6_weeks_of_Treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G6_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = '4-6_Weeks_of_Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G6_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 5 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 5), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==5] <- "4-6_weeks_of_Treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G5_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = '4-6_Weeks_of_Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G5_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

###############
# Group 10 vs 2
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(2, 10), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==2] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==10] <- "4-6_weeks_of_Treatment"
# Tweedieverse::Tweedieverse(data,
#                            metadata_2,
#                            'analysis/Tweedieverse_MOUSE_CDT_G10_G2',
#                            max_significance = 0.1,
#                            plot_heatmap = T,
#                            plot_scatter = T,
#                            standardize = F,
#                            reference = '4-6_Weeks_of_Treatment,1_Baseline_of_HGFAP'
# )
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G10_G2',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
