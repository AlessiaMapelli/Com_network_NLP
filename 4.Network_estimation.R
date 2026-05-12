rm(list=ls(all=TRUE))

#######################
# Step 0 - Configuration step: define the inout and output directories
#######################
# Output directories
OUTPUT_RES_DIR   = "results"
OUTPUT_PLOT_DIR = "plots"
color_net_plots = "#FFA420"


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
library(dplyr)
library(tidyr)
library(ggplot2)

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

#######################
# Step 2 - Adjacency matrix visualization
#######################
adj_GGM_df <- outcome_df[disease_names, disease_names]
sign_mat <- ifelse(adj_GGM_df == 0, 0, 1)
sign_df <- as.data.frame(as.table(as.matrix(sign_mat)))
colnames(sign_df) <- c("Row", "Col", "Value")
sign_df$Row <- factor(sign_df$Row, levels = rownames(sign_mat))
sign_df$Col <- factor(sign_df$Col, levels = colnames(sign_mat))
sign_df$row_id <- as.numeric(sign_df$Row)
sign_df$col_id <- as.numeric(sign_df$Co)

sign_df <- sign_df %>%
  filter(col_id > row_id) %>%
  select("Row", "Col", "Value")

p_binary <- ggplot(sign_df, aes(x = Col, y = Row, fill = factor(Value))) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c("0" = "lightgray", "1" = color_net_plots),
    labels = c("0" = "Edge not present", "1" = "Edge present"),
    name = ""
  ) +
  coord_fixed() +
  labs(
    title = "Token Similarity Heatmap (GGM) — Lower Triangle",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_binary

ggsave(
  filename = file.path(OUTPUT_PLOT_DIR, "disease_adj_matrix_GGM_heatmap.svg"),
  plot = p_binary,
  width = 12,
  height = 10,
  units = "in"
)

ggsave(
  filename = file.path(OUTPUT_PLOT_DIR, "disease_adj_matrix_GGM_heatmap.png"),
  plot = p_binary,
  width = 12,
  height = 10,
  units = "in"
)
#######################
# Step 2 - Disease transition analysis
#######################
transition_counts <- matrix(
  0,
  nrow = length(disease_names),
  ncol = length(disease_names),
  dimnames = list(disease_names, disease_names)
)

patient_transitions <- matrix(
  0,
  nrow = length(disease_names),
  ncol = length(disease_names),
  dimnames = list(disease_names, disease_names)
)

df_main <- sequence_main_final %>%
  filter(event_name %in% disease_names) %>%
  arrange(ID, age_in_days_since_born)
patient_ids <- unique(df_main$ID)
total_patients <- length(patient_ids)
cat("Processing sequences for", total_patients, "patients...\n")


for (pid in patient_ids) {
  
  seq_patient <- df_main %>%
    filter(ID == pid) %>%
    arrange(age_in_days_since_born) %>%
    pull(event_name)
  
  if (length(seq_patient) < 2) next
  
  patient_has_transition <- matrix(
    FALSE,
    nrow = length(disease_names),
    ncol = length(disease_names),
    dimnames = list(disease_names, disease_names)
  )
  
  for (i in seq_len(length(seq_patient) - 1)) {
      # Consecutive events
      disease_first <- seq_patient[i]
      disease_second <- seq_patient[i+1]
      
      if (!is.na(adj_GGM_df[disease_first, disease_second]) &&
          adj_GGM_df[disease_first, disease_second] != 0) {
        
        transition_counts[disease_first, disease_second] <-
          transition_counts[disease_first, disease_second] + 1
        
        patient_has_transition[disease_first, disease_second] <- TRUE
      }
  }
  
  patient_transitions <- patient_transitions + patient_has_transition
}

transition_percentage <- round(patient_transitions / total_patients * 100, 2)

for (i in disease_names) {
  for (j in disease_names) {
    if (i == j || adj_GGM_df[i, j] == 0) {
      transition_counts[i, j] <- NA
      transition_percentage[i, j] <- NA
    }
  }
}

# 2.1 Transition percentage heatmap

transition_df <- as.data.frame(as.table(transition_percentage))
colnames(transition_df) <- c("Source", "Target", "Percentage")

transition_df$Source <- factor(transition_df$Source, levels = disease_names)
transition_df$Target <- factor(transition_df$Target, levels = disease_names)

p_transition <- ggplot(transition_df, aes(x = Target, y = Source, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    na.value = "white",
    name = "Transition\npercentage"
  ) +
  coord_fixed() +
  labs(
    title = "Transition Percentage Heatmap (GGM)",
    x = "Targets",
    y = "Sources"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_transition

ggsave(
  filename = file.path(OUTPUT_PLOT_DIR, "disease_adj_matrix_GGM_heatmap_transition_percentage.svg"),
  plot = p_transition,
  width = 12,
  height = 10,
  units = "in"
)

ggsave(
  filename = file.path(OUTPUT_PLOT_DIR, "disease_adj_matrix_GGM_heatmap_transition_percentage.png"),
  plot = p_transition,
  width = 12,
  height = 10,
  units = "in"
)

save(
  adj_GGM_df,
  transition_counts,
  transition_percentage,
  file = file.path(OUTPUT_RES_DIR, "GGM_transition_analysis.RData")
)
