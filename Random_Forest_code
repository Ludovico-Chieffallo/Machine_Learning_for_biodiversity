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
library(doParallel)

#Functions----
setwd("/media/r_projects/phd_ludovico/2021")
source("R/clean_name.R") #shortcut for:
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

# Image preprocessing: Selection of RGB bands and Standardization
cl<- makePSOCKcluster(4)
registerDoParallel(cl)
t0<-Sys.time()
ortho <- map(ortho, function(x) {
  # Select the first three bands for each orthomosaic
  x_rgb <- x[[1:3]]
  
  # Calculate the minimum and maximum values for the RGB bands
  min_val <- min(x_rgb, na.rm = TRUE)
  max_val <- max(x_rgb, na.rm = TRUE)
  
  # Standardize values within the 0-255 range
  x_standardized <- terra::clamp(x_rgb, lower = min_val, upper = max_val)
  x_standardized <- (x_standardized - min_val) / (max_val - min_val) * 255
  
  return(x_standardized)
})
t1<-Sys.time()
t1-t0
#Time difference of 1.263211 hours

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
t0<-Sys.time()
ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))
ortho
t1<-Sys.time()
t1-t0
#Time difference of 25.15504 mins

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


# Sort names to match
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

#Values extraction----

# Manually convert a single element for the test
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

# Check for NULL areas
# Create a logical vector indicating which elements in 'fiori' are NULL
are_null_fiori <- sapply(fiori, is.null)

# Print the indices of NULL elements in 'fiori'
null_indices_fiori <- which(are_null_fiori)
print(null_indices_fiori)

# Repeat the same process for 'erba'
are_null_erba <- sapply(erba, is.null)
null_indices_erba <- which(are_null_erba)
print(null_indices_erba)

# Data extraction----
df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))


# Column names of each dataframe
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
#   tidyr::pivot_longer(cols =c("band_1", "band_2", "band_3"),
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

#RF Model----
train_control_RF <- trainControl(
  method = "cv",             # Cross-validation
  number = 10,               # Number of folds for cross-validation
  savePredictions = "final", # Save predictions for each model
  classProbs = TRUE,         # Calculate class probabilities
  summaryFunction = twoClassSummary # Summary function for binary classification
)

max_mtry <- max(2, floor(sqrt(ncol(training_data) - 1)))

tuneGrid_RF <- expand.grid(
  mtry = seq(2, max_mtry, by=1),
  splitrule = c("gini", "extratrees"),
  min.node.size = c(1, 5, 10, 15)
)

print(tuneGrid_RF)


## Train the model using ranger and the tuning grid
# Calculate class weights
#classes <- unique(training_data$label)
#class_weights <- c("grass" = 0.1, "flower" = 0.9)

# Remove any rows with NA values
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

# Train the model using class weights
set.seed(123)
t0<-Sys.time()
RF <- train(
  label ~ ., 
  data = training_data,
  method = "ranger", 
  trControl = train_control_RF, 
  tuneGrid = tuneGrid_RF,
  metric = "ROC"
)
t1<-Sys.time()
t1-t0
#10.78502 mins
print(RF)
print(RF$bestTune)
stopCluster(cl)

# If you want to train an SVM model, you could use method = "svmLinear"

#RF model test----
predictions_RF <- predict(RF, newdata= test_data)
test_accuracy_RF <- mean(predictions_RF == test_data$label)
print(paste("Test Accuracy:", test_accuracy_RF))

confusion_matrix_RF <- confusionMatrix(as.factor(predictions_RF), test_data$label, positive = "flower")
print(confusion_matrix_RF)

accuracy_RF <- confusion_matrix_RF$overall['Accuracy']
recall_RF <- confusion_matrix_RF$byClass["Sensitivity"]
f1_RF <- confusion_matrix_RF$byClass["F1"]
auc_roc_RF <- roc(test_data$label, as.numeric(predictions_RF))
auc_RF <- auc(auc_roc_RF)
precision_RF <- confusion_matrix_RF$byClass["Pos Pred Value"]

print(paste("Accuracy: ", accuracy_RF))
print(paste("Recall: ", recall_RF))
print(paste("F1-Score: ", f1_RF))
print(paste("AUC-ROC: ", auc_RF))
print(paste("Precision: ", precision_RF))

#saveRDS(RF, file = "/media/r_projects/phd_ludovico/2021/RF_model_22.rds")

#WARNING FROM HERE------------------------

#GRAPHS----
matrix <- as.data.frame(confusion_matrix_RF$table)

ggplot(data = matrix, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 1.5) +
  scale_fill_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for the RF Model", x = "Actual Label", y = "Prediction")

metrics <- data.frame(
  Metric = c("Accuracy", "Recall", "F1-Score", "AUC-ROC", "Precision"),
  Value = c(accuracy_RF, recall_RF, f1_RF, auc_RF, precision_RF)
)

ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Evaluation Metrics for the RF Model", x = "", y = "Value")

#Attempt flower count----

# Assuming that `ortho` is your list of rasters
# And that the RF model is saved in the variable `RF`

# Function to convert raster to dataframe and make predictions
process_and_predict_RF <- function(raster) {
  # Convert the raster into a dataframe, including only the bands used for training
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rename the columns to ensure they match those used in the RF model
  # Make sure this mapping matches the band names in your raster and the names used in the training set
  names(pixel_data) <- c('x', 'y', 'band_1', 'band_2', 'band_3')
  
  # Select only the columns of interest (remove 'x', 'y' if not used for prediction)
  pixel_data <- dplyr::select(pixel_data, band_1, band_2, band_3)
  
  # Apply the RF model to make predictions on the pixels
  predictions <- predict(RF, newdata = pixel_data)
  
  # Return the predictions
  return(predictions)
}

# Apply the function to each raster in `ortho`
prediction_results_RF <- lapply(ortho, process_and_predict_RF)

# Aggregate prediction results for each raster
prediction_counts_RF <- lapply(prediction_results_RF, table)

# Print the counts for each raster
lapply(prediction_counts_RF, print)

# Set the area covered by each pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Assuming `prediction_counts` is your list with pixel counts
prediction_counts_cm2 <- lapply(prediction_counts_RF, function(x) x * area_per_pixel_cm2)

# Print the results converted to cm²
lapply(prediction_counts_cm2, print)

#RF RESULTS-----
#----------------------------------------------------------------------------------------------------

drone_results_RF <- data.frame(
  area = c("OrthoEL1", "OrthoEY1", "OrthoEY3", "OrthoFR1", "OrthoFR2", "OrthoFR4",
           "OrthoGE3", "OrthoGU1", "OrthoGU3", "OrthoGU4", "OrthoMA2", "OrthoMA3",
           "OrthoPA2", "OrthoRA4", "OrthoSG3", "OrthoSI2", "OrthoVA1", "OrthoVA3",
           "OrthoVA4", "OrthoWA1", "OrthoWY2", "OrthoWY3"),
  grass_cm2 = c(2606838, 4828012, 6002191, 5768441, 9737990, 7447449, 6727721, 
                3283969, 2988731, 6406325, 18623283, 5905001, 5010614, 3721440,
                6978128, 4177196, 8961403, 8302703, 5914980, 3359116, 4143781, 
                4208069),
  flower_cm2 = c(7767, 176871, 67620, 113717, 7998, 83917, 15052, 15509, 76073,
                 44917, 7601, 59732, 83767, 10328, 6721, 23668, 127348, 38170, 
                 89565, 5611, 32128, 37806)
)

# Print the dataframe to verify the data
print(drone_results_RF)

# Save the results to a CSV file
#write.csv(drone_results_RF, "Results_Areas_2021.csv", row.names = FALSE)

#Correlations----

# Load field data from the CSV file
field_data <- read.csv("/media/r_projects/phd_ludovico/2021/field_data_2021.csv", sep = ";", header = TRUE)

# Join the columns using dplyr
combined_data <- left_join(field_data, drone_results_RF[, c("area", "flower_cm2")], by = "area")

#FLOWER AB./FLOWER AB. GBM-----

# Calculate the correlation between the two estimates of flower cm²
correlation <- cor(combined_data$flower_surface_cm2, combined_data$flower_cm2)

# Print the correlation coefficient
print(correlation)

# Build the linear model
linear_model <- lm(flower_cm2 ~ flower_surface_cm2, data = combined_data)

# Model summary to obtain R^2 and p-value
model_summary <- summary(linear_model)

# Extract R^2 value
r_squared <- model_summary$r.squared

# Extract p-value for the slope (flower_surface_cm2)
p_value <- model_summary$coefficients[2,4]

# Print R^2 and p-value
print(paste("R^2:", r_squared))
print(paste("p-value:", p_value))

# Format the p-value to display three decimal places
formatted_p_value <- formatC(p_value, format = "f", digits = 6)

annotation_text <- paste("R^2 = ", round(r_squared, digits = 3),
                         "\np-value = ", formatted_p_value)

# Create the enhanced plot with ggplot2 and add the annotation
ggplot(combined_data, aes(x = flower_surface_cm2, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
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

