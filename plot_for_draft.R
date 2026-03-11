### 1. DISEASE INCIDENCE DISTRIBUTION ####
data <- read.csv("/group/diangelantonio/users/alessia_mapelli/Disease_network/Comorb_network/sequence_main_final.csv")
head(data)


library(data.table)
library(ggplot2)

# Convert to data.table if not already
setDT(data)

# 1. Convert age to years and get first incidence per patient per disease
age_incidence <- data[, .(Years = round(age_in_days_since_born / 365, 0)), 
                      by = .(f.eid, disease_name, event_id)][
                        order(f.eid, Years), .SD[1], by = .(f.eid, disease_name)]

# 2. Count events by disease and age
age_dist <- age_incidence[, .N, by = .(disease_name, Years)]

# 3. Calculate total events per disease
total_incidence <- age_dist[, .(total_events = sum(N)), by = disease_name]

# 4. Merge and calculate percentages
age_analysis <- age_dist[total_incidence, on = "disease_name"]
age_analysis[, percentage := round(N/total_events*100, 2)]

# 5. Create the plot (6 rows x 4 columns = 24 facets)
p <- ggplot(age_analysis, aes(x = factor(Years), y = percentage)) +
  geom_point(color = "blue3", position = position_jitter(width = 0.1), size = 0.8) +
  labs(title = "Disease incidence by age - All 24 diseases",
       x = "Age (Years)",
       y = "% of total disease events") +
  facet_wrap(~ disease_name, scales = "free_y", ncol = 4, nrow = 6) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 10)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.5, "lines"))

print(p)

# Optional: print summary statistics
cat("Number of unique diseases:", length(unique(age_analysis$disease_name)), "\n")
cat("Disease names:", paste(unique(age_analysis$disease_name), collapse = ", "), "\n")

# Optional: Create disease name mapping for cleaner labels
disease_mapping <- c(
  "Kidney Diseases, Chronic, Stage 3" = "Chronic Kidney Disease Stage 3",
  "Kidney Diseases, Chronic, Stage 4" = "Chronic Kidney Disease Stage 4",
  "Kidney Diseases, Chronic, Stage 5" = "Chronic Kidney Disease Stage 5",
  "Cardiomyopathy, Dilated" = "Dilated Cardiomyopathy",
  "Pacemaker, Artificial" = "Artificial Pacemaker",
  "Heart Failure, Endstage" = "Heart Failure Endstage",
  "Stroke, Hemorrhagic" = "Hemorrhagic Stroke",
  "Kidney Failure, Chronic, Endstage" = "Chronic Kidney Failure Endstage",
  "Tachycardia, Ventricular" = "Ventricular Tachycardia"
)

# Apply mapping (if name not in mapping, keep original)
age_analysis[, disease_name_clean := ifelse(disease_name %in% names(disease_mapping),
                                            disease_mapping[disease_name],
                                            disease_name)]

# p2 with clean facet labels (no grey boxes)
p2 <- ggplot(age_analysis, aes(x = factor(Years), y = N)) +
  geom_point(color = "blue3", position = position_jitter(width = 0.1), size = 0.8) +
  labs(title = "Disease incidence by age",
       x = "Age (Years)",
       y = "Number of events") +
  facet_wrap(~ disease_name_clean, scales = "free_y", ncol = 4, nrow = 6) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 10)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),  # Rimuove i quadrati grigi
        panel.spacing = unit(0.5, "lines"))

print(p2)

ggsave("disease_incidence_by_age.png", plot = p2, width = 9, height = 11, dpi = 400, bg = "white")

ggsave("disease_incidence_by_age.svg", plot = p2, width = 9, height = 11, device = "svg")



### 2. TOTAL EVENTS PRIMADRY DISEASE ####
data <- read.csv("/group/diangelantonio/users/alessia_mapelli/Disease_network/Comorb_network/sequence_main_final.csv")
setDT(data)
# Convert all events to years (not just first incidence)
all_events_age <- data[, .(Years = round(age_in_days_since_born / 365, 0))]

# Count total events by age
total_events_by_age <- all_events_age[, .N, by = Years]

# Create histogram
p3 <- ggplot(total_events_by_age, aes(x = Years, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  labs(title = "Distribution of disease events by age (primary cardiometabolic diseases)",
       x = "Age (Years)",
       y = "Number of events") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# Optional: save the plot
ggsave("total_events_by_age.png", plot = p3, 
       width = 9, height = 6, dpi = 300, bg = "white")
ggsave("total_events_by_age.svg", plot = p3, 
       width = 9, height = 6, device = "svg")
