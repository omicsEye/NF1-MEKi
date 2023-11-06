library(omePath)
library(omicsArt)
#library(omicsPattern)
fdr_threshold <- 1.00
setwd("~/Library/CloudStorage/Box-Box/NF1_MEKi")
#source('~/Documents/omicsEye/omicsPattern/R/utils.R')
#source('~/Documents/omicsEye/omicsPattern/R/pcl_utils.R')

# NF1 Mice
biomarker_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G11_G9_G10_G8_G7/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
)

#G9 vs G7
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="4_3 weeks of MEKi" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G9_G7",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = 0.1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G9_G7",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = .1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")



#G11 vs G7
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="5_4 - 6 weeks of MEKi" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G11_G7",
  score_col = 'coef',
  pval_threshold = .25,
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
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G11_G7",
  score_col = 'coef',
  pval_threshold = .25,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

#G8 vs G7
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="2_3 weeks of Vehicle" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G8_G7",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = 0.1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUPER PATHWAY",
  feature_col = "Metabolite")

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G8_G7",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = .1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

#G10 vs G7
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="3_4 - 6 weeks of Vehicle" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G10_G7",
  score_col = 'coef',
  pval_threshold = .25,
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
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G10_G7",
  score_col = 'coef',
  pval_threshold = .25,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

##############
# Control Mice
biomarker_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_CDT_G6_G4_G5_G3_G1/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
)

#G4 vs G1
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="4_3 weeks of MEKi",]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G4_G1",
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
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G4_G1",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = .1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

#G6 vs G1
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="5_4 - 6 weeks of MEKi" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G6_G1",
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
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G6_G1",
  score_col = 'coef',
  pval_threshold = NA,
  fdr_threshold = .1,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")

#G11 vs G9
biomarker_results <- read.delim(
  "analysis/Tweedieverse_MOUSE_HD4_G11_G9/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
)
score_data_filtered <- biomarker_results[biomarker_results$metadata=="Treatment" & biomarker_results$value=="4 - 6 weeks of MEKi" ,]
rownames(score_data_filtered) <- score_data_filtered$feature

loaded_data <- load_data(input = 'data/Cleaned_by_Ali_CNMC-01-18VWCLP MOUSE HD4 12-5-19.XLSX',
                         type = 'known', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata

#SUB and SUPER PATHWAYS
omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/enrichment_analysis/omePath_SUPER_PATHWAYS_Mouse_G11_G9",
  score_col = 'coef',
  pval_threshold = .25,
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
  output = "analysis/enrichment_analysis/omePath_SUB_PATHWAYS_Mouse_G11_G9",
  score_col = 'coef',
  pval_threshold = .25,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = features_info,
  method = "ks",
  min_member = 2,
  pathway_col = "SUB PATHWAY",
  feature_col = "Metabolite")
