#Library----
library(terra)
library(pROC)
library(caret)
library(sf)
library(purrr)
library(tools)
library(tidyverse)
library(rsample)
library(furrr)
library(pbapply)
library(ranger)
library(tidymodels)
library(e1071)

#Functions----
setwd("/media/r_projects/phd_ludovico/2021")
source("R/clean_name.R") #shortcut for :
# clean_name <- function(file) {
#   var_name <- tools::file_path_sans_ext(basename(file))
#   var_name <- gsub("[^[:alnum:]_]", "", var_name)
#   if (grepl("^[0-9]", var_name)) {
#     var_name <- paste0("x", var_name)
#   }
#   var_name
# }

#Load data----
pathfolder<-"ortho server 2021 05.03"
#ortho-----
ortho <- list.files(path=pathfolder,
                    pattern = "Ortho", 
                    full.names = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))

#buffers-----
buffers <- list.files(pattern = "buffer.*shp",
                      full.names = TRUE )%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))

ortho <- ortho[order(names(ortho))]
buffers <- buffers[order(names(buffers))]

#crop and mask----
# Overlay orthomosaic with buffers
#ortho_intersect <- map2(ortho, map(buffers, ~ ext(.x)), ~ intersect(.x, .y))
#WARNING.
ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))
ortho
#Flowers and Grass-----

# Flowers and Grass combined
#Sys.setenv("SHAPE_RESTORE_SHX" = "YES")
fiori_erba <- list.files(pattern = "fiori e erba.*shp", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.)) %>% 
  map(~.[!st_is_empty(.), ])

fiori <- fiori_erba %>% 
  map(~ . %>% filter(id == 1))  # Use filter to select elements where id is 1

erba <- fiori_erba %>% 
  map(~ . %>% filter(id == 2))


#Sort names to match
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

#Values extraction----

# Manually convert a single element for testing
fiori <- map(fiori_erba, ~{
  # Check if the 'id' column exists
  if ("id" %in% names(.x)) {
    # Filter only if there are rows matching the criterion
    filtered <- .x %>% filter(id == 1)
    # Convert to 'sf' only if the result is not empty
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Return NULL for empty data frames
    }
  } else {
    NULL  # Return NULL if 'id' column doesn't exist
  }
})

erba <- map(fiori_erba, ~{
  # Check if the 'id' column exists
  if ("id" %in% names(.x)) {
    # Filter only if there are rows matching the criterion
    filtered <- .x %>% filter(id == 2)
    # Convert to 'sf' only if the result is not empty
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Return NULL for empty data frames
    }
  } else {
    NULL  # Return NULL if 'id' column doesn't exist
  }
})

#Check for NULL areas
# Create a logical vector indicating which elements in 'fiori' are NULL
are_null_fiori <- sapply(fiori, is.null)

# Print indices of NULL elements in 'fiori'
null_indices_fiori <- which(are_null_fiori)
print(null_indices_fiori)

# Repeat the same process for 'erba'
are_null_erba <- sapply(erba, is.null)
null_indices_erba <- which(are_null_erba)
print(null_indices_erba)

#data extraction----
df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))

#Column names for each dataframe
df_fiori <- map(df_fiori, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})

df_erba <- map(df_erba, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})

#labels
df_fiori <- map(df_fiori, ~{
  .$label <- "flower"
  .
})

df_erba <- map(df_erba, ~{
  .$label <- "grass"
  .
})

#Bind dataframes----

# Merge all dataframes into a list
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba,.id = "origin")

set.seed(123) # For reproducibility
size_of_fiori <- nrow(df_fiori_unified)
# Number of rows to sample from each group

# Create a balanced sample
df_erba_unified <- df_erba_unified %>%
  group_by(origin) %>%
  sample_n(size = size_of_fiori, replace = TRUE) %>%
  ungroup()

df_tot <- rbind(df_fiori_unified, df_erba_unified) %>% 
  rowid_to_column(.) %>%
  rename(ID_pixel = rowid)

# Create a 'long' dataframe #ONLY IF WE SEE NA
# df_long <- df_tot %>%
#   tidyr::pivot_longer(cols =c("banda_1", "banda_2", "banda_3"),
#                       names_to = "banda",
#                       values_to = "value")
# df_long %>%
#   filter(is.na(value)) %>%
#   view()

#Splitting into train and test----
#Select the bands

set.seed(123)

# Create a new dataframe with only the columns you want to use
df_selected <- dplyr::select(df_tot, banda_1, banda_2, banda_3, label)

# Split the selected dataframe into training and test sets
data_split <- initial_split(df_selected, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convert to factor-----
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

# GBM Model ----
train_control_gbm <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

tuneGrid_gbm <- expand.grid(
  interaction.depth = c(1, 3, 5), 
  n.trees = c(50, 100, 150),
  shrinkage = c(0.01, 0.1),
  n.minobsinnode = c(10, 20)
)

# Remove any rows with NA values
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

set.seed(123)
GBM <- train(
  label ~ .,
  data = training_data,
  method = "gbm",
  trControl = train_control_gbm,
  tuneGrid = tuneGrid_gbm,
  verbose = FALSE,
  metric = "ROC"
)

print(GBM)

# GBM Model Test ----
predictions_gbm <- predict(GBM, newdata = test_data)
test_accuracy_gbm <- mean(predictions_gbm == test_data$label)
print(paste("GBM Test Accuracy:", test_accuracy_gbm))

confusion_matrix_gbm <- confusionMatrix(as.factor(predictions_gbm), test_data$label, positive = "flower")
print(confusion_matrix_gbm)

# Other evaluation metrics can be calculated as before

# GBM Model Test ----
predictions_gbm <- predict(GBM, newdata = test_data)
test_accuracy_gbm <- mean(predictions_gbm == test_data$label)
print(paste("GBM Test Accuracy:", test_accuracy_gbm))

confusion_matrix_gbm <- confusionMatrix(as.factor(predictions_gbm), test_data$label, positive = "flower")
print(confusion_matrix_gbm)

accuracy_gbm <- confusion_matrix_gbm$overall['Accuracy']
recall_gbm <- confusion_matrix_gbm$byClass["Sensitivity"]
f1_gbm <- confusion_matrix_gbm$byClass["F1"]
auc_roc_gbm <- roc(test_data$label, as.numeric(predictions_gbm))
auc_gbm <- auc(auc_roc_gbm)
precision_gbm <- confusion_matrix_gbm$byClass["Pos Pred Value"]

print(paste("GBM Accuracy: ", accuracy_gbm))
print(paste("GBM Recall: ", recall_gbm))
print(paste("GBM F1-Score: ", f1_gbm))
print(paste("GBM AUC-ROC: ", auc_gbm))
print(paste("GBM Precision: ", precision_gbm))

#WARNING FROM HERE------------------------

#PLOTS----
matrix <- as.data.frame(confusion_matrix_gbm$table)

ggplot(data = matrix, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 1.5) +
  scale_fill_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for the SVM Model", x = "True Label", y = "Prediction")

metrics <- data.frame(
  Metric = c("Accuracy", "Recall", "F1-Score", "AUC-ROC", "Precision"),
  Value = c(accuracy_gbm, recall_gbm, f1_gbm, auc_gbm, precision_gbm)
)

ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Evaluation Metrics for the SVM Model", x = "", y = "Value")

#Attempt flower count----

# Assuming that `ortho` is your list of rasters
# And that the GBM model is saved in the variable `GBM`

# Function to convert raster to dataframe and make predictions
process_and_predict_GBM <- function(raster) {
  # Convert the raster to a dataframe, including only the bands used for training
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rename columns to ensure they match those used in the SVM model
  # Make sure this mapping matches the band names in your raster and how they were named in the training set
  names(pixel_data) <- c('x', 'y', 'banda_1', 'banda_2', 'banda_3')
  
  # Select only the columns of interest (remove 'x', 'y' if not used for prediction)
  pixel_data <- dplyr::select(pixel_data, banda_1, banda_2, banda_3)
  
  # Apply the GBM model to make predictions on the pixels
  predictions <- predict(GBM, newdata = pixel_data)
  
  # Return the predictions
  return(predictions)
}

# Apply the function to each raster in `ortho`
prediction_results_GBM <- lapply(ortho, process_and_predict_GBM)

# Aggregate prediction results for each raster
prediction_counts_GBM <- lapply(prediction_results_GBM, table)

# Print counts for each raster
lapply(prediction_counts_GBM, print)

# Set the area covered by each pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Assuming `prediction_counts` is your list with pixel counts
prediction_counts_cm2 <- lapply(prediction_counts_GBM, function(x) x * area_per_pixel_cm2)

# Print results converted to cm2
lapply(prediction_counts_cm2, print)

#GBM RESULTS----
# $OrthoEL1
# 
# grass    flower 
# 652534.2   1118.0 
# 
# $OrthoEY1
# 
# grass   flower 
# 1213066   40154 
# 
# $OrthoEY3
# 
# grass      flower 
# 1506984.25   10627.75 
# 
# ... (continue listing all results) ...

# Creating a dataframe from results----
drone_results_GBM <- data.frame(
  area = c("OrthoEL1", "OrthoEY1", "OrthoEY3", "OrthoFR1", "OrthoFR2", "OrthoFR4",
           "OrthoGE3", "OrthoGU1", "OrthoGU3", "OrthoGU4", "OrthoMA2", "OrthoMA3",
           "OrthoPA2", "OrthoRA4", "OrthoSG3", "OrthoSI2", "OrthoVA1", "OrthoVA3",
           "OrthoVA4", "OrthoWA1", "OrthoWY2", "OrthoWY3"),
  grass_cm2 = c(652484.00, 1209483.25, 1505023.25, 1445707.00, 2431047.00, 1880813.00,
                1680257.00, 823906.00, 741951.75, 1591141.50, 4654935.25, 1474259.75,
                1271039.50, 931434.50, 1741892.80, 1044875.50, 2240874.75, 2075758.00,
                1437392.20, 836837.50, 1022498.50, 1055802.25),
  flower_cm2 = c(1168.25, 43736.75, 12588.75, 24844.75, 5461.50, 2048.75,
                 5440.00, 964.25, 24265.75, 21682.75, 2785.75, 16940.25,
                 2555.75, 2430.00, 4323.50, 5344.50, 31321.75, 9461.25,
                 63745.50, 4348.75, 21485.25, 5750.25)
)

# Print dataframe to verify data
print(drone_results_GBM)

# Save results to CSV
#write.csv(drone_results_GBM, "Drone_Area_Results_2021.csv", row.names = FALSE)

#Correlations----

# Load field data from CSV file
#field_data <- read.csv("/media/data/r_projects/phd_ludovico/2021/field_data_2021.csv", sep = ";", header = TRUE)

# Sample data as text with each line separated by a newline and each value separated by a comma
field_data_text <- "area,flower_shannon_area,flower_shannon_number,... (and other columns)"
# Convert text data to a dataframe using read.csv and textConnection
field_data <- read.csv(textConnection(field_data_text), header = TRUE)
# Remove rows with specified areas
field_data <- field_data[!field_data$area %in% c("EY2", "FR3", "WY1", "SI1","EY4", "GU2", "VA2", "WY4"), ]

#TRY TO INCREASE R2----
field_data <- field_data[!field_data$area %in% c( "GU2", "VA4"), ]
drone_results_GBM <- drone_results_GBM[!drone_results_GBM$area %in% c( "GU2", "VA4"), ]

#Select the correct areas
field_data <- field_data[, c("area", "flower_shannon_area", "flower_shannon_number", "flower_simpson_area", 
                             "flower_simpson_number", "flower_surface_cm2", "flower_species_richness", "bee_shannon", 
                             "bee_simpson", "bee_species_richness", "bee_abundance")]

# Display the first rows of the created dataframe
head(field_data)

# Remove prefix "Ortho" from the area column in drone results if necessary
drone_results_GBM$area <- gsub("Ortho", "", drone_results_GBM$area)

# Join the columns using dplyr
merged_data <- left_join(field_data, drone_results_GBM[, c("area", "flower_cm2")], by = "area")

# Verify results
head(merged_data)

#FLOWER AB./FLOWER AB. GBM-----

# Calculate the correlation between the two flower area estimates in cm²
correlation <- cor(merged_data$flower_surface_cm2, merged_data$flower_cm2)

# Print the correlation coefficient
print(correlation)

# Build the linear model
linear_model <- lm(flower_cm2 ~ flower_surface_cm2, data = merged_data)

# Summary of the model to obtain R^2 and p-values
model_summary <- summary(linear_model)

# Extract the R^2 value
r_squared <- model_summary$r.squared

# Extract the p-value for the slope (flower_surface_cm2)
p_value <- model_summary$coefficients[2,4]

# Print the R^2 and p-value
print(paste("R^2:", r_squared))
print(paste("p-value:", p_value))

# Format the p-value to display three decimal places
formatted_p_value <- formatC(p_value, format = "f", digits = 6)

annotation_text <- paste("R^2 = ", round(r_squared, digits = 3),
                         "\np-value = ", formatted_p_value)

# Create the improved plot with ggplot2 and add the annotation
ggplot(merged_data, aes(x = flower_surface_cm2, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE,nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Field Data (cm²)", 
       y = "Drone Data GBM (cm²)", 
       title = "Flower Abundance (Field) vs Flower Abundance UAV (GBM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")

#BEE ABUNDANCE/FLOWER AB GBM----

# Calculate correlation
correlation_bee_abundance <- cor(merged_data$bee_abundance, merged_data$flower_cm2)
print(paste("Correlation:", correlation_bee_abundance))

# Build the linear model
linear_model_bee_abundance <- lm(flower_cm2 ~ bee_abundance, data = merged_data)

# Summary of the model to obtain R^2 and p-value
model_summary_bee_abundance <- summary(linear_model_bee_abundance)

# Extract R^2 value
r_squared_bee_abundance <- model_summary_bee_abundance$r.squared
print(paste("R^2:", r_squared_bee_abundance))

# Extract p-value for the slope (bee_abundance)
p_value_bee_abundance <- model_summary_bee_abundance$coefficients[2, 4]
print(paste("p-value:", p_value_bee_abundance))

# Format p-value to display three decimal places
formatted_p_value_bee_abundance <- formatC(p_value_bee_abundance, format = "f", digits = 3)

# Create annotation text
annotation_text_bee_abundance <- paste("R^2 = ", round(r_squared_bee_abundance, digits = 3),
                                       "\np-value = ", formatted_p_value_bee_abundance)

# Visualize the relationship with ggplot2
ggplot(merged_data, aes(x = bee_abundance, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Bee Abundance", 
       y = "Flower Area (cm²)", 
       title = "Correlation between Bee Abundance and Flower Area (GBM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_bee_abundance, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")

#BEE SPECIES RICHNESS / FLOWER AB. GBM----

# Calculate correlation between bee_species_richness and flower_cm2
correlation_richness <- cor(merged_data$bee_species_richness, merged_data$flower_cm2)

# Build linear model for bee_species_richness vs flower_cm2
linear_model_richness <- lm(flower_cm2 ~ bee_species_richness, data = merged_data)

# Summary of the model to obtain R^2 and p-value
model_summary_richness <- summary(linear_model_richness)

# Extract R^2 value
r_squared_richness <- model_summary_richness$r.squared

# Extract p-value for the slope (bee_species_richness)
p_value_richness <- coef(summary(linear_model_richness))["bee_species_richness", "Pr(>|t|)"]

# Format p-value
formatted_p_value_richness <- formatC(p_value_richness, format = "f", digits = 3)

# Create annotation text
annotation_text_richness <- paste("R^2 = ", round(r_squared_richness, digits = 3),
                                  "\np-value = ", formatted_p_value_richness)

# Create improved plot for bee species richness
ggplot(merged_data, aes(x = bee_species_richness, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Bee Species Richness", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Bee Species Richness and Flower Area (GBM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_richness, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")

#SHANNON / FLOW AB GBM----

# Build the linear model for bee_shannon
linear_model_shannon <- lm(flower_cm2 ~ bee_shannon, data = merged_data)

# Summary of the model to obtain R^2 and p-value
model_summary_shannon <- summary(linear_model_shannon)
r_squared_shannon <- model_summary_shannon$r.squared
p_value_shannon <- coef(summary(linear_model_shannon))["bee_shannon", "Pr(>|t|)"]

# Format p-value
formatted_p_value_shannon <- formatC(p_value_shannon, format = "f", digits = 3)

# Create annotation text
annotation_text_shannon <- paste("R^2 = ", round(r_squared_shannon, digits = 3),
                                 "\np-value = ", formatted_p_value_shannon)

# Create improved plot for Shannon index of bees
ggplot(merged_data, aes(x = bee_shannon, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Shannon Index", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Shannon Index and Flower Area (GBM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_shannon, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")

#saveRDS(GBM, file = "/media/r_projects/phd_ludovico/2021/GBM_model_22.rds")





