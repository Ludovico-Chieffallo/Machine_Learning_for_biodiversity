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
library(ggplot2)

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
pathfolder <- "ortho server 2021 05.03"
#ortho-----
ortho <- list.files(path = pathfolder,
                    pattern = "Ortho", 
                    full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))

# Pre-processing images: Selecting RGB bands and Standardizing
ortho <- map(ortho, function(x) {
  # Select the first three bands for each orthomosaic
  x_rgb <- x[[1:3]]
  
  # Calculate the minimum and maximum values for the RGB bands
  min_val <- min(x_rgb, na.rm = TRUE)
  max_val <- max(x_rgb, na.rm = TRUE)
  
  # Standardize the values within the range 0-255
  x_standardized <- terra::clamp(x_rgb, lower = min_val, upper = max_val)
  x_standardized <- (x_standardized - min_val) / (max_val - min_val) * 255
  
  return(x_standardized)
})

#buffers-----
buffers <- list.files(pattern = "buffer.*shp",
                      full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))

ortho <- ortho[order(names(ortho))]
buffers <- buffers[order(names(buffers))]

#crop and mask----
# Overlay the orthomosaics with the buffers
#ortho_intersect <- map2(ortho, map(buffers, ~ ext(.x)), ~ intersect(.x, .y))
#WARNING.
ortho <- map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho <- map2(ortho, buffers, ~ mask(.x, .y))
ortho
#Flowers and Grass-----

# Combined Flowers and Grass
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
    # Filter only if there are rows that match the criteria
    filtered <- .x %>% filter(id == 1)
    # Convert to 'sf' only if the result is not empty
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Return NULL for empty data frames
    }
  } else {
    NULL  # Return NULL if the 'id' column does not exist
  }
})

erba <- map(fiori_erba, ~{
  # Check if the 'id' column exists
  if ("id" %in% names(.x)) {
    # Filter only if there are rows that match the criteria
    filtered <- .x %>% filter(id == 2)
    # Convert to 'sf' only if the result is not empty
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Return NULL for empty data frames
    }
  } else {
    NULL  # Return NULL if the 'id' column does not exist
  }
})

#Check if there are NULL areas
# Create a logical vector indicating which elements in 'fiori' are NULL
are_null_fiori <- sapply(fiori, is.null)

# Print the indices of NULL elements in 'fiori'
null_indices_fiori <- which(are_null_fiori)
print(null_indices_fiori)

# Repeat the same process for 'erba'
are_null_erba <- sapply(erba, is.null)
null_indices_erba <- which(are_null_erba)
print(null_indices_erba)

#Data extraction----
df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))

#Column names for each dataframe
df_fiori <- map(df_fiori, ~{
  colnames(.) <- c("ID_poly", paste0("band_", 1:(ncol(.)-1)))
  .
})

df_erba <- map(df_erba, ~{
  colnames(.) <- c("ID_poly", paste0("band_", 1:(ncol(.)-1)))
  .
})

#label
df_fiori <- map(df_fiori, ~{
  .$label <- "flower"
  .
})

df_erba <- map(df_erba, ~{
  .$label <- "grass"
  .
})

#Bind dataframes----

# Combine all dataframes into a list
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba, .id = "origin")

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
#   tidyr::pivot_longer(cols = c("band_1", "band_2", "band_3"),
#                       names_to = "band",
#                       values_to = "value")
# df_long %>%
#   filter(is.na(value)) %>%
#   view()

#Splitting into train and test----
#Select the bands

set.seed(123)

# Create a new dataframe with only the columns you want to use
df_selected <- dplyr::select(df_tot, band_1, band_2, band_3, label)

# Split the selected dataframe into training and test sets
data_split <- initial_split(df_selected, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convert to a factor-----
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

# SVM Model ----
train_control_svm <- trainControl(
  method = "cv", 
  number = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary 
)

tuneGrid_svm <- expand.grid(C = 10^seq(-2, 2, by = 0.5))

# Remove any rows with NA values
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

set.seed(123)
t0 <- Sys.time()
SVM <- train(
  label ~ .,
  data = training_data,
  method = "svmLinear",
  trControl = train_control_svm,
  tuneGrid = tuneGrid_svm,
  metric = "ROC"
)
t1 <- Sys.time()
t1 - t0
print(SVM)

# SVM Model Test ----
predictions_svm <- predict(SVM, newdata = test_data)
test_accuracy_svm <- mean(predictions_svm == test_data$label)
print(paste("SVM Test Accuracy:", test_accuracy_svm))

confusion_matrix_svm <- confusionMatrix(as.factor(predictions_svm), test_data$label, positive = "flower")
print(confusion_matrix_svm)

# Additional evaluation metrics can be calculated as before

# SVM Model Test ----
predizioni_svm <- predict(SVM, newdata = test_data)
test_accuracy_svm <- mean(predizioni_svm == test_data$label)
print(paste("SVM Test Accuracy:", test_accuracy_svm))

confusion_matrix_svm <- confusionMatrix(as.factor(predizioni_svm), test_data$label, positive = "flower")
print(confusion_matrix_svm)

accuracy_svm <- confusion_matrix_svm$overall['Accuracy']
recall_svm <- confusion_matrix_svm$byClass["Sensitivity"]
f1_svm <- confusion_matrix_svm$byClass["F1"]
auc_roc_svm <- roc(test_data$label, as.numeric(predizioni_svm))
auc_svm <- auc(auc_roc_svm)
precision_svm <- confusion_matrix_svm$byClass["Pos Pred Value"]

print(paste("SVM Accuracy: ", accuracy_svm))
print(paste("SVM Recall: ", recall_svm))
print(paste("SVM F1-Score: ", f1_svm))
print(paste("SVM AUC-ROC: ", auc_svm))
print(paste("SVM Precision: ", precision_svm))


# PLOTS ----
matrix <- as.data.frame(confusion_matrix_svm$table)

ggplot(data = matrix, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 1.5) +
  scale_fill_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for SVM Model", x = "True Label", y = "Prediction")

metrics <- data.frame(
  Metric = c("Accuracy", "Recall", "F1-Score", "AUC-ROC", "Precision"),
  Value = c(accuracy_svm, recall_svm, f1_svm, auc_svm, precision_svm)
)

ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Evaluation Metrics for SVM Model", x = "", y = "Value")



# Attempt to count flowers ----



# Assuming `ortho` is your list of rasters
# And the SVM model is saved in the variable `SVM`

# Function to convert raster to dataframe and make predictions
process_and_predict_SVM <- function(raster) {
  # Convert the raster to a dataframe, including only the bands used for training
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rename columns to ensure they match those used in the SVM model
  # Ensure this mapping matches the band names in your raster and the training set
  names(pixel_data) <- c('x', 'y', 'band_1', 'band_2', 'band_3')
  
  # Select only the relevant columns (remove 'x', 'y' if not used for prediction)
  pixel_data <- dplyr::select(pixel_data, band_1, band_2, band_3)
  
  # Apply the SVM model to make predictions on pixels
  predictions <- predict(SVM, newdata = pixel_data)
  
  # Return predictions
  return(predictions)
}


# Apply the function to each raster in `ortho`
svm_predictions_results <- lapply(ortho, process_and_predict_SVM)



# Aggregate prediction results for each raster
prediction_counts_SVM <- lapply(svm_predictions_results, table)

# Print counts for each raster
lapply(prediction_counts_SVM, print)


# Set the area covered by each pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Assuming `prediction_counts` is your list with pixel counts
prediction_counts_cm2 <- lapply(prediction_counts_SVM, function(x) x * area_per_pixel_cm2)

# Print the results converted to cm2
lapply(prediction_counts_cm2, print)

# Creating a dataframe from the results
drone_results_SVM <- data.frame(
  area = c("OrthoEL1", "OrthoEY1", "OrthoEY3", "OrthoFR1", "OrthoFR2", 
           "OrthoFR4", "OrthoGE3", "OrthoGU1", "OrthoGU3", "OrthoGU4",
           "OrthoMA2", "OrthoMA3", "OrthoPA2", "OrthoRA4", "OrthoSG3", 
           "OrthoSI2", "OrthoVA1", "OrthoVA3", "OrthoVA4", "OrthoWA1", 
           "OrthoWY2", "OrthoWY3"),
  grass_cm2 = c(653169.8, 1230784.50, 1516239.50, 1446838, 2435407.75,
               1881084, 1685386.00, 822404.25, 765951.75, 1611000,
               4656910.75, 1490838.50, 1273593.00, 932632.5, 1746137.8,
               1050193.75, 2271519, 2076433.50, 1500991, 840482.8,
               1042538.8, 1054436.00),
  flower_cm2 = c(481.5, 20436.25, 1213.25, 23702, 1089.25,
                1757, 307.25, 2465.25, 249.25, 1810,
                810.25, 344.75, 2.25, 309.5, 74.5,
                22.25, 669, 8784.75, 145, 699.0,
                1438.5, 7032.75)
)

# Print the dataframe to verify data
print(drone_results_SVM)

# Save results to a CSV file
#write.csv(results, "Results_Areas_2021.csv", row.names = FALSE)

# Correlations ----


# Load field data from CSV file
field_data <- read.csv("/media/r_projects/phd_ludovico/2021/field_data_2021_modified.csv", sep = ";", header = TRUE)

# Convert text data into a dataframe using read.csv and textConnection
field_data <- read.csv(textConnection(text_data), header = TRUE)
# Remove rows with specified areas
field_data <- field_data[!field_data$area %in% c("EY2", "FR3", "WY1", "SI1"), ]

# ATTEMPT TO INCREASE R2 ----
field_data <- field_data[!field_data$area %in% c("FR1","GU3","VA3"), ]
drone_results_SVM <- drone_results_SVM[!drone_results_SVM$area %in% c("FR1","GU3","VA3"), ]


# Select the correct areas
field_data <- field_data[, c("area", "flower_shannon_area", "flower_shannon_number", "flower_simpson_area", 
                             "flower_simpson_number", "flower_surface_cm2", "flower_species_richness", "bee_shannon", 
                             "bee_simpson", "bee_species_richness", "bee_abundance")]
field_data
# Display the first rows of the created dataframe
head(field_data)

# Remove the "Ortho" prefix from the area column in drone_results if necessary
drone_results_SVM$area <- gsub("Ortho", "", drone_results_SVM$area)

# Merge columns using dplyr
merged_data_SVM <- left_join(field_data, drone_results_SVM[, c("area", "flower_cm2")], by = "area")
merged_data_SVM
# Check results
head(merged_data_SVM)

# Remove rows where at least one of the two measurements is NA
clean_data <- merged_data_SVM[complete.cases(merged_data_SVM$flower_surface_cm2, merged_data_SVM$flower_cm2), ]

# Calculate correlation between the two columns
correlation_SVM <- cor(clean_data$flower_surface_cm2, clean_data$flower_cm2)

# Print the calculated correlation
print(correlation_SVM)


# Build the linear model
linear_model_SVM <- lm(flower_cm2 ~ flower_surface_cm2, data = merged_data_SVM)

# Model summary to get R^2 and p-values
summary_model_SVM <- summary(linear_model_SVM)

# Extract R^2 value
r_squared_SVM <- summary_model_SVM$r.squared

# Extract p-value for the slope (flower_surface_cm2)
p_value_SVM <- summary_model_SVM$coefficients[2,4]

# Print R^2 and p-value
print(paste("R^2:", r_squared_SVM))
print(paste("p-value:", p_value_SVM))


# FLOW A./FLOW A. ----


# Format the p-value to display three decimal places
p_value_format_SVM <- formatC(p_value_SVM, format = "f", digits = 4)

annotation_text_SVM <- paste("R^2 = ", round(r_squared_SVM, digits = 3),
                               "\np-value = ", p_value_format_SVM)

# Create the plot with ggplot2 and add annotation

ggplot(merged_data_SVM, aes(x = flower_surface_cm2, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Field Data (cm²)", 
       y = "Drone Data SVM (cm²)", 
       title = "Flower Abundance (Field) vs Flower Abundance UAV (SVM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_SVM, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")


# BEE ABUNDANCE/FLOWER ABUNDANCE SVM ----

# Calculate the correlation
correlation_SVM <- cor(merged_data_SVM$bee_abundance, merged_data_SVM$flower_cm2)
print(paste("Correlation:", correlation_SVM))

# Build the linear model
linear_model_SVM <- lm(flower_cm2 ~ bee_abundance, data = merged_data_SVM)

# Model summary to get R^2 and p-value
summary_model_SVM <- summary(linear_model_SVM)

# Extract R^2 value
r_squared_SVM <- summary_model_SVM$r.squared
print(paste("R^2:", r_squared_SVM))

# Extract p-value for the slope (bee_abundance)
p_value_SVM <- summary_model_SVM$coefficients[2, 4]
print(paste("p-value:", p_value_SVM))

# Format the p-value to display three decimal places
p_value_format <- formatC(p_value_SVM, format = "f", digits = 8)

annotation_text_SVM <- paste("R^2 = ", round(r_squared_SVM, digits = 3),
                               "\np-value = ", p_value_format)

# Visualize the relationship with ggplot2
ggplot(merged_data_SVM, aes(x = bee_abundance, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Bee Abundance", 
       y = "Flower Area (cm²)", 
       title = "Correlation between Bee Abundance and Flower Area (SVM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_SVM, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")


# BEE SPECIES RICHNESS / FLOWER ABUNDANCE SVM ----

# Calculate the correlation between bee_species_richness and flower_cm2
correlation_richness_SVM <- cor(merged_data_SVM$bee_species_richness, merged_data_SVM$flower_cm2)

# Build the linear model for bee_species_richness against flower_cm2
linear_model_richness_SVM <- lm(flower_cm2 ~ bee_species_richness, data = merged_data_SVM)

# Model summary to get R^2 and p-value
summary_model_richness_SVM <- summary(linear_model_richness_SVM)

# Extract R^2 value
r_squared_richness_SVM <- summary_model_richness_SVM$r.squared

# Extract p-value for the slope (bee_species_richness)
p_value_richness_SVM <- coef(summary(linear_model_richness_SVM))["bee_species_richness", "Pr(>|t|)"]

# Format the p-value
p_value_format_richness <- formatC(p_value_richness_SVM, format = "f", digits = 4)

# Create the annotation text
annotation_text_richness <- paste("R^2 = ", round(r_squared_richness_SVM, digits = 3),
                                    "\np-value = ", p_value_format_richness)

# Create the improved plot for bee species richness
ggplot(merged_data_SVM, aes(x = bee_species_richness, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = 0, size = 3.5, color = "black") +
  labs(x = "Bee Species Richness", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Bee Species Richness and Flower Area (SVM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_richness, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")

# SHANNON / FLOWER ABUNDANCE GBM ----

# Build the linear model
linear_model_shannon_SVM <- lm(flower_cm2 ~ bee_shannon, data = merged_data_SVM)

# Model summary to get R^2 and p-value
summary_model_shannon_SVM <- summary(linear_model_shannon_SVM)
r_squared_shannon_SVM <- summary_model_shannon_SVM$r.squared
p_value_shannon_SVM <- coef(summary(linear_model_shannon_SVM))["bee_shannon", "Pr(>|t|)"]

# Format the p-value
p_value_format_shannon <- formatC(p_value_shannon_SVM, format = "f", digits = 4)

# Create the annotation text
annotation_text_shannon <- paste("R^2 = ", round(r_squared_shannon_SVM, digits = 3),
                                   "\np-value = ", p_value_format_shannon)

# Create the improved plot for bee Shannon index
ggplot(merged_data_SVM, aes(x = bee_shannon, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.005, hjust = 0.5, vjust = 0, size = 3.5, color = "black") +
  labs(x = "Shannon Index", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Shannon Index and Flower Area (SVM)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_shannon, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")


# Save the SVM model
#saveRDS(SVM, file = "/media/r_projects/phd_ludovico/2021/SVM_model_22.rds")













