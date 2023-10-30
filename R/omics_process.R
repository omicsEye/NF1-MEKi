
#require(gdata)
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
library(deepath)
library(omicsArt)
fdr_threshold <- 1.00
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")


loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX',
                                       type = 'known', sheet = 1, ID = 'Metabolite')
#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX'
#Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE CDT 12-5-19.XLSX'

data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_box2_CLP CDT', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_box2_CLP CDT', method = 'tsne')
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(16, 17), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==16] <- "Human_NF1"
metadata_2[metadata_2==17] <- "Human_NF1_MEKi"
Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_HUMAN_CLP_CDT_G17_G16',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = "Treatment,Human_NF1"
)



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
metadata$TREATMENT <-  as.character(metadata$TREATMENT)
pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_box2_HD4', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_box2_HD4', method = 'tsne')
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(16, 17), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_HUMAN_HD4_G17_G16',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_HUMAN_HD4_G17_G16',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)



# Mice
loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE CDT 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata
metadata$TREATMENT <-  as.character(metadata$TREATMENT)
pcoa_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'pcoa_MOUSE CDT', method = 'pcoa')
tsne_plots <- ordplots(data, metadata, output = 'analysis/', outputname = 'tsne_MOUSE CDT', method = 'tsne')
#metadata_2 <- metadata[metadata$`BOX NUMBER` == 'Box2',]
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(13, 12), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_MOUSE_CDT_G13_G12',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)
Maaslin2(data,
         metadata_2,
         'analysis/Maaslin2_MOUSE_CDT_G13_G12',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         #transform = 'LOG',
         normalization = 'NONE')

metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(15, 14), c( "GROUP NUMBER"), drop = F]
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
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(15, 13), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_MOUSE_CDT_G15_G13',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(14, 12), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"

Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_MOUSE_CDT_G14_G12',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F
)

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
metadata_2 <- metadata[metadata$`GROUP NUMBER` %in% c(1, 6), c( "GROUP NUMBER"), drop = F]
colnames(metadata_2) <- "Treatment"
metadata_2[metadata_2==1] <- "1_Baseline_of_HGFAP"
metadata_2[metadata_2==6] <- "4-6 weeks post-treatment"
Tweedieverse::Tweedieverse(data, 
                           metadata_2, 
                           'analysis/Tweedieverse_MOUSE_CDT_G6_G1',
                           max_significance = 0.1,
                           plot_heatmap = T,
                           plot_scatter = T,
                           standardize = F,
                           reference = 'Treatment,1_Baseline_of_HGFAP'
)
################Metabolomics data#############################
#####set up#########################

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX',
                                       type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

output_path <- 'analysis/'
output = paste(output_path,'/','enricment_analysis', sep = '')
dir.create(file.path(output), showWarnings = TRUE)
output = paste(output,'/','HUMAN CLP CDT_SUPERPathway', sep = '')
dir.create(file.path(output), showWarnings = TRUE)

tryCatch({
  path_results <- deepath::deepath(
                 input_data = stats_table,
                 input_metadata = NA,
                 meta <- 'TREATMENT',
                 case_label <- case_group,
                 control_label <- control_group,
                 output = output,
                 score_col = 'logFC',
                 pval_threshold = .05,
                 fdr_threshold = NA,
                 Pathway.Subject = NA,#'Metabolic',
                 do_plot = TRUE,
                 mapper_file = features_info,
                 method = 'wilcox',
                 min_member = 2,
                 pathway_col = "SUPER PATHWAY",
                 feature_col = "Metabolite")
}, error = function(err) {
  #invisible(dev.off())
  logging::logerror("Unable to do deepath or its plots!!!")
  logging::logerror(err)
})





