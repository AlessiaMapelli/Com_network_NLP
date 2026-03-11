setwd("~/Notebooks/Comorbidities_networks/SIS 2026 - Comorbidity")
rm(list=ls(all=TRUE))


source("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Code_git/4.Method_implementation/ggReg_refactored.R")

# COMORBIDITY NETWORK
data_path <- "phrase_embeddings_main_par_df.csv"

# Check if file exists
if (!file.exists(data_path)) {
  stop("Data file not found. Please check the path: ", data_path)
}

# Load the data
cat("Loading UK Biobank disease embedding data...\n")
ukb_data <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)

library(readr)
sequence_main_final <- read_csv("sequence_main_final.csv")
disease_names <- sort(unique(sequence_main_final$disease_name))

results <- GGReg_full_estimation(x = ukb_data, # expression data, nrow: subjects; ncol: features.
                                 known_ppi = NULL, # previously known PPI 
                                 covariates = NULL,       #covariates to regress on
                                 scr = FALSE,            
                                 lambda_mean_type = "min",
                                 lambda_prec_type = "min",
                                 weight = 1,
                                 asparse = 1,
                                 use_slurm = FALSE,
                                 verbose = TRUE)

save(results,disease_names, file = "results_GGM/GGM_adj_main_par_full.RData")

outcome_df <- results$Dic_Delta_hat$Baseline
fix_col <- function(name) {
  name <- gsub("\\.\\.", ", ", name, perl = TRUE)
  name <- gsub("\\.", " ", name, perl = TRUE)
  name <- gsub("\\s+", " ", name, perl = TRUE)
  trimws(name)
}
rownames(outcome_df) <- fix_col(rownames(outcome_df))
colnames(outcome_df) <- fix_col(colnames(outcome_df))

write.csv(outcome_df,
          file = "results_GGM/GGM_adj_main_par_full.csv",
          row.names = TRUE,
          quote = TRUE)
