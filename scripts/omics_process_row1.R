library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
# library(Maaslin2)
library(readxl)
library(omePath)
library(omicsArt)
library(dplyr)

fdr_threshold <- 1.00
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")

loaded_data <- omicsArt::load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX',
                                        type = 'known', sheet = 1, ID = 'Metabolite')

#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE CDT 12-5-19.XLSX'

data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata <- metadata[,c("GROUP NUMBER", "TREATMENT"),drop=F]
colnames(metadata) <- c("Treatment", "Group")
metadata[metadata==12] <- "2 hours after Vehicle"
metadata[metadata==13] <- "2 hours after MEKi"

variances = apply(data, 2, var)
scaled_data = scale(data[, variances>0])
pca = prcomp(scaled_data)
summary(pca)
dim(pca$x)
pcoa_plots <- ordplots(data, metadata, output = 'analysis/',
                       outputname = 'pcoa_box2_CLP CDT',
                       method = 'pcoa')

ordplots(scaled_data, metadata, output = 'analysis/',
         outputname = 'pcoa_box2_CLP CDT',
         method = 'pcoa')

tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_box2_CLP CDT', method = 'tsne')
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("No_treatment", "2 months of treatment"), c( "Treatment"), drop = F]
#metadata_2[metadata_2==16] <- "1_No_treatment"
#metadata_2[metadata_2==17] <- "2 months of treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_HUMAN_CLP_CDT_G17_G16',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)
# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_HUMAN_CLP_CDT_G17_G16',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')



###### HD4 ######
loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE CDT 12-5-19.XLSX'

data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata <- metadata[,c("GROUP NUMBER", "TREATMENT"),drop=F]
metadata[metadata==16] <- "No treatment"
metadata[metadata==17] <- "2 months of treatment"
metadata$TREATMENT <-  as.character(metadata$TREATMENT)
pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_box2_HD4', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_box2_HD4', method = 'tsne')
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("No treatment", "2 months of treatment"), c( "TREATMENT"), drop = F]
colnames(metadata_2) <- c( "Treatment")
#metadata_2[metadata_2==16] <- "1_No_treatment"
#metadata_2[metadata_2==17] <- "2 months of treatment"

# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_HUMAN_HD4_G17_G16',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_HUMAN_HD4_G17_G16',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference="1_No_treatment"
)



# Mice CDT that changed to HD4

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
data_G7 <- data
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata <- metadata[,c("GROUP NUMBER", "TREATMENT"),drop=F]
metadata[metadata==7] <- "NF1_Baseline_of_HGFAP"

metadata[metadata==1] <- "Control Baseline_of_HGFAP"
metadata[metadata==2] <- "Mice NF1.Het"
metadata[metadata==3] <- "2 hours after Vehicle"

metadata[metadata==4] <- "3 weeks of treatment"
metadata[metadata==6] <- "4-6 weeks post-treatment"
metadata[metadata==7] <- "NF1 Baseline_of_HGFAP"
metadata[metadata==12] <- "2 hours after Vehicle"
metadata[metadata==13] <- "2 hours after MEKi"
metadata[metadata==14] <- "12 hours after Vehicle"
metadata[metadata==15] <- "12 hours after MEKi"
metadata$TREATMENT <-  as.character(metadata$TREATMENT)

pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_MOUSE HD4', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_MOUSE Hd4', method = 'tsne')
#metadata_2 <- metadata[metadata$`BOX NUMBER` == 'Box2',]
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("2 hours after Vehicle", "2 hours after MEKi"),"GROUP NUMBER", drop = F]
colnames(metadata_2) <- "Treatment"
#metadata_2[metadata_2==12] <- "1_2 hours after Vehicle"
#metadata_2[metadata_2==13] <- "2 hours after MEKi"
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_HD4_G13_G12',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = "Treatment,2 hours after Vehicle"
)
# Maaslin2(data,
#          metadata_2,
#          'analysis/Maaslin2_MOUSE_CDT_G13_G12',
#          min_abundance = 0,
#          min_prevalence = 0,
#          max_significance = 0.1,
#          #analysis_method = 'CPLM',
#          #random_effects = c("GROUP NUMBER"),
#          standardize = FALSE,
#          #transform = 'LOG',
#          normalization = 'NONE')

metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("12 hours after MEKi", "12 hours after Vehicle"), c( "TREATMENT"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G15_G14',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G15_G14',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("12 hours after MEKi", "2 hours after MEKi"), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G15_G13',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("12 hours after Vehicle", "2 hours after Vehicle"), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G14_G12',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)

metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("1_Baseline_of_HGFAP", "3 weeks of treatment"), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G4_G1',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,Baseline_of_HGFAP'
)
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c("Baseline_of_HGFAP", "4-6 weeks post-treatment"), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_MOUSE_CDT_G6_G1',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,Baseline_of_HGFAP'
)

#group 7 vs 16


metadata_2_G7 <- metadata[metadata$`GROUP NUMBER` %in% c("NF1_Baseline_of_HGFAP"), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2_G7) <- "Treatment"
data_G7 <- data_G7[row.names(metadata_2_G7),]
###### Group 16 HD4 ######
loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')


data_G16 <- loaded_data$data
data_G16 <- numeric_dataframe(data_G16)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata <- metadata[,c("GROUP NUMBER", "TREATMENT"),drop=F]
metadata[metadata==16] <- "No treatment"
metadata[metadata==17] <- "2 months of treatment"
metadata$TREATMENT <-  as.character(metadata$TREATMENT)
metadata_2_G16 <- metadata[metadata$`GROUP NUMBER` %in% c("No treatment"), c( "TREATMENT"), drop = F]
colnames(metadata_2_G16) <- c( "Treatment")

metadata_2 <- rbind(metadata_2_G7, metadata_2_G16)
comm_met <- intersect(colnames(data_G16), colnames(data_G7))
data_G16 <- data_G16[, comm_met]
data_G7 <- data_G7[, comm_met]
data <- rbind(data_G16, data_G7)
Tweedieverse::Tweedieverse(data,
                           metadata_2,
                           'analysis/Tweedieverse_HUMAN_CLP_CDT_G7_G16',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,Baseline_of_HGFAP'
)
