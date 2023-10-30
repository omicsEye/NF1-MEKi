library(omePath)
library(omicsArt)
#library(omicsPattern)
fdr_threshold <- 1.00
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")
#source('~/Documents/omicsEye/omicsPattern/R/utils.R')
#source('~/Documents/omicsEye/omicsPattern/R/pcl_utils.R')

# #Maaslin2
# biomarker_results <- read.delim(
#   "analysis/Maaslin2_HUMAN_CLP_CDT_G17_G16/all_results.tsv",
#   sep = '\t',
#   header = T,
#   fill = F,
#   comment.char = "" ,
#   check.names = F,
#   row.names = 1
# )
# score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="2 months of treatment",]
# 
# 
# loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX',
#                          type = 'known', sheet = 1, ID = 'Metabolite')
# data <- loaded_data$data
# data <- numeric_dataframe(data)
# features_info <- loaded_data$feature_metadata
# 
# #SUB and SUPER PATHWAYS
# omePath_result <- omePath::omePath(
#   input_data = score_data_filtered,
#   output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_HUMAN_CLP_CDT",
#   score_col = 'coef',
#   pval_threshold = NA,
#   fdr_threshold = 0.1,
#   Pathway.Subject = NA,#'Metabolic',
#   do_plot = TRUE,
#   mapper_file = features_info,
#   method = "ks",
#   min_member = 2,
#   pathway_col = "SUPER PATHWAY",
#   feature_col = "Metabolite")
# 
# omePath_result <- omePath::omePath(
#   input_data = score_data_filtered,
#   output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HUMAN_CLP_CDT",
#   score_col = 'coef',
#   pval_threshold = NA,
#   fdr_threshold = 0.1,
#   Pathway.Subject = NA,#'Metabolic',
#   do_plot = TRUE,
#   mapper_file = features_info,
#   method = "ks",
#   min_member = 2,
#   pathway_col = "SUB PATHWAY",
#   feature_col = "Metabolite")
# #### HD4
# biomarker_results <- read.delim(
#   "analysis/Maaslin2_HUMAN_HD4_G17_G16/all_results.tsv",
#   sep = '\t',
#   header = T,
#   fill = F,
#   comment.char = "" ,
#   check.names = F,
#   row.names = 1
# )
# score_data_filtered <- biomarker_results[biomarker_results$metadata=="Group" & biomarker_results$value=="17" ,]
# 
# 
# #HD4
# 
# loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX',
#                          type = 'known', sheet = 1, ID = 'Metabolite')
# data <- loaded_data$data
# data <- numeric_dataframe(data)
# features_info <- loaded_data$feature_metadata
# 
# #SUB and SUPER PATHWAYS
# omePath_result <- omePath::omePath(
#   input_data = score_data_filtered,
#   output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_HD4",
#   score_col = 'coef',
#   pval_threshold = NA,
#   fdr_threshold = 0.1,
#   Pathway.Subject = NA,#'Metabolic',
#   do_plot = TRUE,
#   mapper_file = features_info,
#   method = "ks",
#   min_member = 2,
#   pathway_col = "SUPER PATHWAY",
#   feature_col = "Metabolite")
# 
# omePath_result <- omePath::omePath(
#   input_data = score_data_filtered,
#   output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HD4",
#   score_col = 'coef',
#   pval_threshold = NA,
#   fdr_threshold = 0.1,
#   Pathway.Subject = NA,#'Metabolic',
#   do_plot = TRUE,
#   mapper_file = features_info,
#   method = "ks",
#   min_member = 2,
#   pathway_col = "SUB PATHWAY",
#   feature_col = "Metabolite")

## Tweedieverse ##
#CLP#
biomarker_results <- read.delim(
  "analysis/Tweedieverse_HUMAN_CLP_CDT_G17_G16/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="2 months of treatment",]


loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN CLP CDT 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_HUMAN_CLP_CDT_G17_G16",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HUMAN_CLP_CDT_G17_G16",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = 0.1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")
#### HD4
biomarker_results <- read.delim(
  "analysis/Tweedieverse_HUMAN_HD4_G17_G16/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="2 months of treatment" ,]


#HD4
## USE THIS

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP HUMAN HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/tweedie_omePath_SUPER_PATHWAYS_HD4_G17_G16",
  score_col = 'coef',
  pval_threshold =.15,
  fdr_threshold = NA,
  Pathway.Subject = NA,
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_HD4_G17_G16",
  score_col = 'coef',
  pval_threshold = .15,
  fdr_threshold = NA,
  Pathway.Subject = NA,
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

################
## Mouse Data ##
################ 
#Mouse 2 hr
biomarker_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_HD4_G13_G12/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="2 hourse after MEKi" ,]

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/tweedie_omePath_SUPER_PATHWAYS_MOUSE_2hr_HD4_G13_G12",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = .1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/tweedie_omePath_SUB_PATHWAYS_MOUSE_2hr_HD4_G13_G12",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

#Mouse 12 hr
biomarker_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G15_G14/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Group" & biomarker_results$value=="15" ,]

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_MOUSE_12hr_HD4_G15_G14",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_MOUSE_12hr_HD4_G15_G14",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

