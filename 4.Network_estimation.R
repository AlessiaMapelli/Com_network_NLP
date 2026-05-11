rm(list=ls(all=TRUE))

#######################
# Step 0 - Configuration step: define the inout and output directories
#######################
# Output directories
OUTPUT_RES_DIR   = "results"


OUTPUT_RES_NET_PATH = paste0(OUTPUT_RES_DIR,"/network_estimation/")
if (!dir.exists(OUTPUT_RES_DIR)) {
  dir.create(OUTPUT_RES_DIR, recursive = TRUE)
}

if (!dir.exists(OUTPUT_RES_NET_PATH)) {
  dir.create(OUTPUT_RES_NET_PATH, recursive = TRUE)
}

#######################
# Step 1 - Network estimation
#######################
# 1.1 Load libraries, GMM functions and inputs
library(readr)
source("functions/ggReg_main_functions.R")
data_path <- file.path(OUTPUT_RES_DIR, "cardiomet_risk_adj_embeddings.csv")
if (!file.exists(data_path)) {
  stop("Data file not found. Please check the path: ", data_path)
}
emb_data <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)
sequence_main_final <- read_csv(file.path(OUTPUT_RES_DIR,"sequence_main_final.csv"))
disease_names <- sort(unique(sequence_main_final$event_name))
conf_matrix <- read_csv(file.path(OUTPUT_RES_DIR,"GGM_prior_clinical_knowledge_input.csv"))

# 1.2 Run GMM model and store the results
results <- GGReg_full_estimation(x = emb_data, # data, nrow: subjects; ncol: features.
                                 scr = TRUE, # TRUE, include prior clinical knowledge
                                 scr_matrix = conf_matrix, #Matrix describing the prior clinical knowledge
                                 lambda_mean_type = "min",
                                 lambda_prec_type = "min",
                                 weight = 1,
                                 asparse = 1,
                                 use_slurm = FALSE,
                                 output_path = OUTPUT_RES_NET_PATH,
                                 verbose = TRUE)

save(results,disease_names, file = paste0(OUTPUT_RES_DIR, "/GGM_adj_main_par_full.RData"))

outcome_df <- results$results$Dic_adj_matrics$Baseline
fix_col <- function(name) {
  name <- gsub("\\.\\.", ", ", name, perl = TRUE)
  name <- gsub("\\.", " ", name, perl = TRUE)
  name <- gsub("\\s+", " ", name, perl = TRUE)
  trimws(name)
}
rownames(outcome_df) <- fix_col(rownames(outcome_df))
colnames(outcome_df) <- fix_col(colnames(outcome_df))

write.csv(outcome_df,
          file = paste0(OUTPUT_RES_DIR, "GGM_adj_main_par_full.csv"),
          row.names = TRUE,
          quote = TRUE)
