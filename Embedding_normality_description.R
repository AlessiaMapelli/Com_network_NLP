# Load required libraries
library(dplyr)
library(ggplot2)
library(corrplot)
library(reshape2)
library(data.table)

# Set working directory and load data
setwd("~/Notebooks/Comorbidities_networks/SIS 2026 - Comorbidity")

############################
# A. Check the embedding for the comorbidity network
############################
data_path <- "phrase_embeddings_main_par_df.csv"

# Check if file exists
if (!file.exists(data_path)) {
  stop("Data file not found. Please check the path: ", data_path)
}

# Load the data
cat("Loading UK Biobank disease embedding data...\n")
ukb_data <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)

cat("Data dimensions:", dim(ukb_data), "\n")
cat("Column names preview:\n")
colnames(ukb_data)


# ===============================================================================
# STEP 1: EXAMINE FIRST EMBEDDINGS - DISTRIBUTION, SUMMARY, MISSING VALUES
# ===============================================================================

first_embedding <- "Female"
cat("\n=== EXAMINING FIRST BIOMARKER:", first_embedding, "===\n")

# Summary statistics
cat("Summary statistics for", first_embedding, ":\n")
print(summary(ukb_data[,first_embedding]))

# Missing values check
missing_count <- sum(is.na(ukb_data[,first_embedding]))
missing_percent <- (missing_count / nrow(ukb_data)) * 100
cat("\nMissing values:", missing_count, "(", round(missing_percent, 2), "%)\n")

# Distribution plot
cat("Creating distribution plot for", first_embedding, "...\n")

# Create a comprehensive plot for the first biomarker
par(mfrow = c(2, 2))

# Histogram
hist(ukb_data[,first_embedding], 
     main = paste("Distribution of", first_embedding),
     xlab = "Embeddings",
     col = "lightblue",
     breaks = 30)

# Box plot
boxplot(ukb_data[,first_embedding],
        main = paste("Boxplot of", first_embedding),
        ylab = "Embeddings",
        col = "lightgreen")

# Q-Q plot for normality check
qqnorm(ukb_data[,first_embedding], 
       main = paste("Q-Q Plot of", first_embedding))
qqline(ukb_data[,first_embedding], col = "red")

# Density plot
if (sum(!is.na(ukb_data[,first_embedding])) > 0) {
  plot(density(ukb_data[,first_embedding], na.rm = TRUE),
       main = paste("Density Plot of", first_embedding),
       xlab = "Embeddings")
}

par(mfrow = c(1, 1))


# ===============================================================================
# STEP 2: OUTLIER DETECTION AND NORMALITY CHECK
# ===============================================================================

cat("\n=== NORMALITY CHECK ===\n")

shapiro_test_results <- rep(NA, ncol(ukb_data))

for(i in 1:ncol(ukb_data)){
  shapiro_test <- shapiro.test(ukb_data[, i])
  cat("Shapiro-Wilk test p-value:", shapiro_test$p.value, "\n")
  shapiro_test_results[i] <- shapiro_test$p.value
}

cat("\n Disease for which the normality assumprions id rejected: \n")
non_normal_scores <- colnames(ukb_data)[shapiro_test_results<0.05]
print(non_normal_scores)

# Create normality plots
# par(mfrow = c(2, 2))
par(mfrow = c(1, 2))
# Histogram after preprocessing
hist(ukb_data[, non_normal_scores[1]], 
     main = paste("Distribution of", non_normal_scores[1]),
     xlab = "Embeddings",
     col = "lightcoral",
     breaks = 30)

# Q-Q plot after preprocessing
qqnorm(ukb_data[, non_normal_scores[1]], 
       main = paste("Q-Q Plot of",  non_normal_scores[1]))
qqline(ukb_data[, non_normal_scores[1]], col = "red")


par(mfrow = c(1, 1))


cat("\n=== OUTLIER DETECTION AND ANALYSIS ===\n")

# Function to detect outliers using multiple methods
detect_outliers <- function(x, method = "IQR", k = 1.5) {
  if (all(is.na(x))) return(rep(FALSE, length(x)))
  
  if (method == "IQR") {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(x < lower_bound | x > upper_bound)
  } else if (method == "z_score") {
    z_scores <- abs(scale(x))
    return(z_scores > k & !is.na(z_scores))
  } else if (method == "modified_z") {
    median_x <- median(x, na.rm = TRUE)
    mad_x <- mad(x, na.rm = TRUE)
    if (mad_x == 0) return(rep(FALSE, length(x)))
    modified_z <- 0.6745 * (x - median_x) / mad_x
    return(abs(modified_z) > k)
  }
}

# Initialize outlier tracking
outlier_summary <- data.frame(
  Disease = character(),
  Total_Values = integer(),
  IQR_Outliers = integer(),
  IQR_Percent = numeric(),
  ZScore_Outliers = integer(),
  ZScore_Percent = numeric(),
  ModZ_Outliers = integer(),
  ModZ_Percent = numeric(),
  Combined_Outliers = integer(),
  Combined_Percent= numeric(),
  stringsAsFactors = FALSE
)

# Analyze outliers for each metabolite
cat("Analyzing outliers across all disease...\n")

outlier_matrices <- list(
  IQR = matrix(FALSE, nrow = nrow(ukb_data), ncol = ncol(ukb_data)),
  ZScore = matrix(FALSE, nrow = nrow(ukb_data), ncol = ncol(ukb_data)),
  ModifiedZ = matrix(FALSE, nrow = nrow(ukb_data), ncol = ncol(ukb_data)),
  Combined = matrix(FALSE, nrow = nrow(ukb_data), ncol = ncol(ukb_data))
)

# Set column names for outlier matrices
for (method in names(outlier_matrices)) {
  colnames(outlier_matrices[[method]]) <- colnames(ukb_data)
}

for (i in 1:ncol(ukb_data)) {
  metabolite_name <- colnames(ukb_data)[i]
  values <- ukb_data[, i]
  
  # Skip if all values are NA
  if (all(is.na(values))) {
    outlier_summary <- rbind(outlier_summary, data.frame(
      Disease = metabolite_name,
      Total_Values = 0,
      IQR_Outliers = 0, IQR_Percent = 0,
      ZScore_Outliers = 0, ZScore_Percent = 0,
      ModZ_Outliers = 0, ModZ_Percent = 0,
      Combined_Outliers = 0, Combined_Percent=0,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # Detect outliers using different methods
  iqr_outliers <- detect_outliers(values, "IQR", 1.5)
  z_outliers <- detect_outliers(values, "z_score", 3)
  modz_outliers <- detect_outliers(values, "modified_z", 3.5)
  combined_outliers <- modz_outliers
  
  # Store outlier information
  outlier_matrices$IQR[, i] <- iqr_outliers
  outlier_matrices$ZScore[, i] <- z_outliers
  outlier_matrices$ModifiedZ[, i] <- modz_outliers
  for(j in 1:nrow(ukb_data)){
    combined_outliers[j] <- all(
      c(iqr_outliers[j], z_outliers[j], modz_outliers[j]))
    outlier_matrices$Combined[j, i] <- combined_outliers[j]
  }
  
  
  # Calculate percentages
  valid_count <- sum(!is.na(values))
  iqr_count <- sum(iqr_outliers, na.rm = TRUE)
  z_count <- sum(z_outliers, na.rm = TRUE)
  modz_count <- sum(modz_outliers, na.rm = TRUE)
  combined_count <- sum(combined_outliers, na.rm = TRUE)
  
  # Add to summary
  outlier_summary <- rbind(outlier_summary, data.frame(
    Disease = metabolite_name,
    Total_Values = valid_count,
    IQR_Outliers = iqr_count,
    IQR_Percent = ifelse(valid_count > 0, round(iqr_count/valid_count*100, 2), 0),
    ZScore_Outliers = z_count,
    ZScore_Percent = ifelse(valid_count > 0, round(z_count/valid_count*100, 2), 0),
    ModZ_Outliers = modz_count,
    ModZ_Percent = ifelse(valid_count > 0, round(modz_count/valid_count*100, 2), 0),
    Combined_Outliers = combined_count,
    Combined_Percent = ifelse(valid_count > 0, round(combined_count/valid_count*100, 2), 0),
    stringsAsFactors = FALSE
  ))
}


# Report outlier statistics
cat("=== OUTLIER DETECTION SUMMARY ===\n")
cat("Total disease analyzed:", nrow(outlier_summary), "\n")

# Overall outlier statistics
total_observations <- sum(outlier_summary$Total_Values)
total_iqr_outliers <- sum(outlier_summary$IQR_Outliers)
total_z_outliers <- sum(outlier_summary$ZScore_Outliers)
total_modz_outliers <- sum(outlier_summary$ModZ_Outliers)
total_comb_outliers <- sum(outlier_summary$Combined_Outliers)

cat("Overall outlier percentages:\n")
cat("- IQR method (k=1.5):", round(total_iqr_outliers/total_observations*100, 2), "%\n")
cat("- Z-Score method (|z|>3):", round(total_z_outliers/total_observations*100, 2), "%\n")
cat("- Modified Z-Score method (|mz|>3.5):", round(total_modz_outliers/total_observations*100, 2), "%\n")
cat("- Combined method (|mz|>3.5 & |z|>3 & k=1.5):", round(total_comb_outliers/total_observations*100, 2), "%\n")

# Find metabolites with high outlier rates
high_outlier_metabolites <- outlier_summary[outlier_summary$IQR_Outliers > 0, ]
if (nrow(high_outlier_metabolites) > 0) {
  cat("\nDisease with outliers (IQR methods):\n")
  print(high_outlier_metabolites[order(high_outlier_metabolites$Combined_Percent, decreasing = TRUE), 
                                 c("Disease", "IQR_Outliers")][1:min(10, nrow(high_outlier_metabolites)), ])
}

sample_outlier_counts <- rowSums(outlier_matrices$IQR, na.rm = TRUE)
high_outlier_samples <- (sample_outlier_counts> 0.1*ncol(ukb_data))
if (sum(high_outlier_samples) > 0) {
  cat(paste0("\n", sum(high_outlier_samples), " samples with  outliers in >10% of embeddings (Combining the methods)\n"))
}

oulier_non_normal_scores <- high_outlier_metabolites[high_outlier_metabolites %in% non_normal_scores] 
if (nrow(oulier_non_normal_scores) > 0) {
  cat("\nDisease with outliers (IQR methods) and non normal:\n")
  print(oulier_non_normal_scores)
}

# ===============================================================================
# STEP 9: CORRELATION ANALYSIS
# ===============================================================================

cat("\n=== CORRELATION ANALYSIS ===\n")

# Calculate correlation matrix
cat("Calculating correlation matrix...\n")
correlation_matrix <- cor(ukb_data, use = "complete.obs")

# Basic correlation statistics
cat("Correlation matrix dimensions:", dim(correlation_matrix), "\n")
cat("Range of correlations (excluding diagonal):\n")
off_diagonal_cors <- correlation_matrix[upper.tri(correlation_matrix)]
cat("Min correlation:", round(min(off_diagonal_cors, na.rm = TRUE), 3), "\n")
cat("Max correlation:", round(max(off_diagonal_cors, na.rm = TRUE), 3), "\n")
cat("Mean absolute correlation:", round(mean(abs(off_diagonal_cors), na.rm = TRUE), 3), "\n")

# Count high correlations
high_cors_08 <- sum(abs(off_diagonal_cors) > 0.8, na.rm = TRUE)
high_cors_09 <- sum(abs(off_diagonal_cors) > 0.9, na.rm = TRUE)
high_cors_095 <- sum(abs(off_diagonal_cors) > 0.95, na.rm = TRUE)

total_pairs <- length(off_diagonal_cors)
cat("Pairs with |correlation| > 0.8:", high_cors_08, "(", round(high_cors_08/total_pairs*100, 2), "%)\n")
cat("Pairs with |correlation| > 0.9:", high_cors_09, "(", round(high_cors_09/total_pairs*100, 2), "%)\n")
cat("Pairs with |correlation| > 0.95:", high_cors_095, "(", round(high_cors_095/total_pairs*100, 2), "%)\n")

# ===============================================================================
# CORRELATION VISUALIZATIONS
# ===============================================================================

cat("\nCreating correlation visualizations...\n")

# 1. Correlation heatmap
cat("Creating full correlation heatmap...\n")

# Create heatmap
corrplot(correlation_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust",
         tl.cex = 0.6, 
         tl.col = "black",
         title = "Metabolite Correlation Matrix",
         mar = c(0,0,2,0))


# 2. Distribution of correlation coefficients
par(mfrow = c(2, 2))

hist(off_diagonal_cors, 
     breaks = 50, 
     main = "Distribution of Correlation Coefficients", 
     xlab = "Correlation Coefficient", 
     col = "lightblue",
     freq = FALSE)
lines(density(off_diagonal_cors, na.rm = TRUE), col = "red", lwd = 2)

# 3. Absolute correlations
hist(abs(off_diagonal_cors), 
     breaks = 50, 
     main = "Distribution of Absolute Correlations", 
     xlab = "Absolute Correlation", 
     col = "lightgreen",
     freq = FALSE)
lines(density(abs(off_diagonal_cors), na.rm = TRUE), col = "red", lwd = 2)

# 4. Correlation vs Distance plot (if hierarchical clustering)

# Hierarchical clustering
dist_matrix <- as.dist(1 - abs(correlation_matrix))
hc <- hclust(dist_matrix, method = "average")

plot(hc, 
     main = "Hierarchical Clustering of Metabolites", 
     xlab = "", 
     sub = "",
     cex = 0.6)


# 5. Network visualization for highly correlated variables
if (!require(igraph, quietly = TRUE)) {
  install.packages("igraph")
  library(igraph)
}

# Create network of high correlations
high_cor_adj <- abs(correlation_matrix) > 0.8
diag(high_cor_adj) <- FALSE

if (sum(high_cor_adj) > 0) {
  # Create igraph object
  cor_network <- graph_from_adjacency_matrix(high_cor_adj, mode = "undirected")
  
  plot(cor_network, 
       vertex.size = 5, 
       vertex.label.cex = 0.6,
       edge.width = 1,
       main = "High Correlation Network (|r| > 0.8)",
       layout = layout_with_fr)
} else {
  plot.new()
  text(0.5, 0.5, "No high correlations (|r| > 0.8) found", cex = 1.2, adj = 0.5)
}

par(mfrow = c(1, 1))
high_cor_adj <- abs(correlation_matrix) > 0.8
diag(high_cor_adj) <- FALSE

if (sum(high_cor_adj) > 0) {
  # Create igraph object
  cor_network <- graph_from_adjacency_matrix(high_cor_adj, mode = "undirected")
  
  plot(cor_network, 
       vertex.size = 5, 
       vertex.label.cex = 0.6,
       edge.width = 1,
       main = "High Correlation Network (|r| > 0.8)",
       layout = layout_with_fr)
} else {
  plot.new()
  text(0.5, 0.5, "No high correlations (|r| > 0.8) found", cex = 1.2, adj = 0.5)
}
write_csv(as.data.frame(correlation_matrix), "corr_matrix_par.csv")

