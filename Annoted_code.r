# Load necessary libraries
library(terra)  # For spatial data manipulation
library(pROC)  # For ROC curve analysis
library(caret)  # For machine learning models
library(sf)  # For working with shapefiles
library(purrr)  # Functional programming utilities
library(tools)  # Tools for handling file paths
library(tidyverse)  # Data wrangling and visualization
library(rsample)  # Data splitting utilities
library(furrr)  # Parallel processing utilities
library(pbapply)  # Progress bars for apply functions
library(ranger)  # Random forest implementation
library(tidymodels)  # Machine learning framework
library(e1071)  # Support Vector Machines and additional ML tools

# Set working directory and load helper functions
setwd("/media/r_projects/phd_ludovico/2021")
source("R/clean_name.R")  # Function to sanitize file names

# Define the path for orthomosaic images
dataset_path <- "ortho server 2021 05.03"

# Load orthomosaic images and apply a cleaning function to their names
ortho <- list.files(path = dataset_path, pattern = "Ortho", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))  # Convert file paths to raster objects

# Load buffer shapefiles and clean their names
buffers <- list.files(pattern = "buffer.*shp", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))  # Read shapefile as an sf object

# Sort orthomosaic and buffer names alphabetically to ensure alignment
ortho <- ortho[order(names(ortho))]
buffers <- buffers[order(names(buffers))]

# Crop and mask orthomosaic images using buffer extents
ortho <- map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho <- map2(ortho, buffers, ~ mask(.x, .y))

# Load and clean shapefiles for flowers and grass
fiori_erba <- list.files(pattern = "fiori e erba.*shp", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.)) %>%  # Read shapefiles
  map(~.[!st_is_empty(.), ])  # Remove empty geometries

# Separate flowers and grass data
fiori <- fiori_erba %>% map(~ . %>% filter(id == 1))
erba <- fiori_erba %>% map(~ . %>% filter(id == 2))

# Ensure orthomosaic, flowers, and grass datasets have matching names
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

# Extract values from orthomosaic images using flower and grass locations
df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))

# Assign meaningful column names
df_fiori <- map(df_fiori, ~{ colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1))); . })
df_erba <- map(df_erba, ~{ colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1))); . })

# Add class labels
df_fiori <- map(df_fiori, ~{ .$label <- "flower"; . })
df_erba <- map(df_erba, ~{ .$label <- "grass"; . })

# Merge all dataframes into a single dataset
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba, .id = "origin")

# Balance the dataset by oversampling grass data to match flower count
set.seed(123)
size_of_fiori <- nrow(df_fiori_unified)
df_erba_unified <- df_erba_unified %>%
  group_by(origin) %>%
  sample_n(size = size_of_fiori, replace = TRUE) %>%
  ungroup()

df_tot <- rbind(df_fiori_unified, df_erba_unified) %>% 
  rowid_to_column(.) %>%
  rename(ID_pixel = rowid)

# Split data into training (70%) and testing (30%)
set.seed(123)
data_split <- initial_split(df_tot, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convert labels to factors for classification tasks
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

# Define parameters for GBM model training
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

# Train the GBM model
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

# Evaluate the GBM model
predictions_gbm <- predict(GBM, newdata = test_data)
test_accuracy_gbm <- mean(predictions_gbm == test_data$label)
print(paste("GBM Test Accuracy:", test_accuracy_gbm))

confusion_matrix_gbm <- confusionMatrix(as.factor(predictions_gbm), test_data$label, positive = "flower")
print(confusion_matrix_gbm)

# Extract and display key performance metrics
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

# Save the trained GBM model
saveRDS(GBM, file = "/media/r_projects/phd_ludovico/2021/GBM_model_22.rds")
