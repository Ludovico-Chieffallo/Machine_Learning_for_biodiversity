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
setwd("/media/data/r_projects/phd_ludovico/2021")
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
pathfolder <- "ortho server 2021 05.03"
#ortho-----
ortho <- list.files(path = pathfolder,
                    pattern = "Ortho", 
                    full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))


# Pre-processing of images: Selection of RGB bands and Standardization
t0 <- Sys.time()
ortho <- map(ortho, function(x) {
  # Select the first three bands for each orthomosaic
  x_rgb <- x[[1:3]]
  
  # Calculate the minimum and maximum values for the RGB bands
  min_val <- min(x_rgb, na.rm = TRUE)
  max_val <- max(x_rgb, na.rm = TRUE)
  
  # Standardize values within the range 0-255
  x_standardized <- terra::clamp(x_rgb, lower = min_val, upper = max_val)
  x_standardized <- (x_standardized - min_val) / (max_val - min_val) * 255
  
  return(x_standardized)
})
t1 <- Sys.time()
t1 - t0

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
t0 <- Sys.time()
ortho <- map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho <- map2(ortho, buffers, ~ mask(.x, .y))
ortho
t1 <- Sys.time()
t1 - t0

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

# Manually convert a single element for testing
fiori <- map(fiori_erba, ~{
  # Check if the 'id' column exists
  if ("id" %in% names(.x)) {
    # Filter only if there are rows matching the criteria
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
    # Filter only if there are rows matching the criteria
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

# Merge all dataframes into a list
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

training_data <- training_data %>% drop_na(band_1, band_2, band_3)
test_data <- test_data %>% drop_na(band_1, band_2, band_3)



# Neural network ----


# Define cross-validation strategy

train_control_NNET <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final",
  classProbs = TRUE,
  summaryFunction = twoClassSummary)


# Nnet hyperparameters grid

nnet_grid <- expand.grid(
  decay = c(0.5, 1e-2, 1e-4),
  size = c(2, 4, 8))


# Nnet training

NNET <- train(label ~ .,
              data = training_data,
              method = "nnet",
              trControl = train_control_NNET,
              tuneGrid = nnet_grid,
              allowParallel = TRUE,
              metric = "ROC",
              importance = 'impurity',
              weights = ifelse(training_data$label == "flower", 10, 1))


# NNET Model Test ----

predictions_NNET <- predict(NNET, newdata = test_data)

pred_char <- as.character(predictions_NNET)
predictions_NNET <- factor(pred_char, levels = c("grass", "flower"))

test_accuracy_NNET <- mean(predictions_NNET == test_data$label)
test_accuracy_NNET

print(paste("NNET Test Accuracy:", test_accuracy_NNET))

confusion_matrix_NNET <- confusionMatrix(as.factor(predictions_NNET), test_data$label, positive = "flower")
print(confusion_matrix_NNET)

# Additional evaluation metrics can be calculated as before


accuracy_NNET <- confusion_matrix_NNET$overall['Accuracy']
recall_NNET <- confusion_matrix_NNET$byClass["Sensitivity"]
f1_NNET <- confusion_matrix_NNET$byClass["F1"]
precision_NNET <- confusion_matrix_NNET$byClass["Pos Pred Value"]



test_data$label <- factor(test_data$label, levels = c("grass", "flower"))
predictions_NNET_prob <- predict(NNET, newdata = test_data, type = "prob")

# Assumendo che "flower" sia la classe positiva, estrai la probabilità corrispondente
prob_flower <- predictions_NNET_prob$flower

# Calcola la ROC specificando i livelli (assicurati che "grass" sia il livello negativo e "flower" quello positivo)
auc_roc_NNET <- pROC::roc(response = test_data$label, predictor = prob_flower, levels = c("grass", "flower"))
auc_NNET <- auc(auc_roc_NNET)

print(auc_NNET)



ci_auc_boot <- pROC::ci.auc(auc_roc_NNET, conf.level = 0.95, boot.n = 2000)
print(ci_auc_boot)









print(paste("_NNET Accuracy: ", accuracy_NNET))
print(paste("_NNET Recall: ", recall_NNET))
print(paste("_NNET F1-Score: ", f1_NNET))
print(paste("_NNET AUC-ROC: ", auc_NNET))
print(paste("_NNET Precision: ", precision_NNET))


# Attempt flower count ----

# Assuming `ortho` is your list of rasters
# And the GBM model is saved in the `GBM` variable

# Function to convert raster to dataframe and make predictions
process_and_predict_NNET <- function(raster) {
  # Convert the raster to a dataframe, including only the bands used for training
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rename columns to match those used in the NNET model
  # Ensure that this mapping matches the band names in your raster and how they were named in the training set
  names(pixel_data) <- c('x', 'y', 'band_1', 'band_2', 'band_3')
  
  # Select only the relevant columns (remove 'x', 'y' if not used for prediction)
  pixel_data <- dplyr::select(pixel_data, band_1, band_2, band_3)
  
  # Apply the NNET model to make predictions on the pixels
  predictions <- predict(NNET, newdata = pixel_data)
  
  # Return the predictions
  return(predictions)
}

# Apply the function to each raster in `ortho`
prediction_results_NNET <- lapply(ortho, process_and_predict_NNET)

# Aggregate prediction results for each raster
prediction_counts_NNET <- lapply(prediction_results_NNET, table)

# Print counts for each raster
lapply(prediction_counts_NNET, print)

# Set the area covered by each pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Assuming `prediction_counts` is your list with pixel counts
prediction_counts_cm2_NNET <- lapply(prediction_counts_NNET, function(x) x * area_per_pixel_cm2)

# Print the results converted to cm²
lapply(prediction_counts_cm2_NNET, print)


# Creating a dataframe from the results ----
drone_results_NNET <- data.frame(
  area = c("OrthoEL1", "OrthoEY1", "OrthoEY3", "OrthoFR1", "OrthoFR2", 
           "OrthoFR4", "OrthoGE3", "OrthoGU1", "OrthoGU3", "OrthoGU4",
           "OrthoMA2", "OrthoMA3", "OrthoPA2", "OrthoRA4", "OrthoSG3", 
           "OrthoSI2", "OrthoVA1", "OrthoVA3", "OrthoVA4", "OrthoWA1", 
           "OrthoWY2", "OrthoWY3"),
  grass_cm2 = c(640018.8, 709945.2, 1256118.0, 1378716, 2412519.25,
                1685528.0, 1632932, 793031.25, 579951, 1539840,
                4582537.25, 1108933, 916572.2, 895130.5, 1732467.8,
                959881, 1875582.2, 2068833.8, 1393572.5, 836679.00,
                984420.50, 1018754),
  flower_cm2 = c(13632.5, 541275.5, 261334.8, 91823, 23977.75,
                 197313.5, 52761, 31838.25, 186250, 72970,
                 75183.75, 382250, 357023.0, 37811.5, 13744.5,
                 90335, 396605.5, 16384.5, 107563.8, 4502.75,
                 59556.75, 42715)
)

# Print the dataframe to verify data
print(drone_results_NNET)

# Saving results to a CSV file
#write.csv(results, "Results_Areas_2021.csv", row.names = FALSE)

# Correlations ----

# Load field data from CSV file
#field_data <- read.csv("/media/data/r_projects/phd_ludovico/2021/field_data_2021.csv", sep = ";", header = TRUE)

# Data text, with each row separated by a newline and each value separated by a comma
data_text <- ",area,flower_shannon_area,flower_shannon_number,flower_simpson_area,flower_simpson_number,flower_surface_cm2,flower_species_richness,bee_shannon,bee_simpson,bee_species_richness,bee_abundance,flowers_rf_original,flowers_rf_1cm,flowers_rf_2cm,flowers_rf_5cm,flowers_svm_original,flowers_svm_1cm,flowers_svm_2cm,flowers_svm_5cm,flowers_nnet_original,flowers_nnet_1cm,flowers_nnet_2cm,flowers_nnet_5cm
1,EL1,0.910721783,0.783767869,0.564777891,0.461355529,1221.1,3,0,0,1,1,0.495117792,1.15071909,1.46230725,2.477797359,3.344074039,3.247341441,3.778030047,5.348648297,0.338154148,0.31,0.601650519,3.4803638
...
"

# Convert data text into a dataframe using read.csv and textConnection
field_data <- read.csv(textConnection(data_text), header = TRUE)
# Remove rows with specific areas
field_data <- field_data[!field_data$area %in% c("EY2", "FR3", "WY1", "SI1"), ]

# Try to increase R2 ----
#field_data <- field_data[!field_data$area %in% c( "MA3", "VA1", "PA2","EY1"), ]
#drone_results_NNET <- drone_results_NNET[!drone_results_NNET$area %in% c( "MA3", "VA1", "PA2","EY1"), ]

# Select the correct areas
field_data_NNET <- field_data[, c("area", "flower_shannon_area", "flower_shannon_number", "flower_simpson_area", 
                             "flower_simpson_number", "flower_surface_cm2", "flower_species_richness", "bee_shannon", 
                             "bee_simpson", "bee_species_richness", "bee_abundance")]

# Show the first rows of the created dataframe
head(field_data_NNET)

# Remove the "Ortho" prefix from the area column in drone results if needed
drone_results_NNET$area <- gsub("Ortho", "", drone_results_NNET$area)

# Join columns using dplyr
joined_data_NNET <- left_join(field_data, drone_results_NNET[, c("area", "flower_cm2")], by = "area")

# Verify the results
head(joined_data_NNET)

# FLOWER AB./FLOWER AB. GBM ----

# Remove rows where at least one of the two measurements is NA
joined_data_NNET <- joined_data_NNET[complete.cases(joined_data_NNET$flower_surface_cm2, joined_data_NNET$flower_cm2), ]

# Calculate the correlation between the two columns
correlation_NNET <- cor(joined_data_NNET$flower_surface_cm2, joined_data_NNET$flower_cm2)

# Print the calculated correlation
print(correlation_NNET)

# Build the linear model
linear_model_NNET <- lm(flower_cm2 ~ flower_surface_cm2, data = joined_data_NNET)

# Model summary to obtain R^2 and p-value
summary_model_NNET <- summary(linear_model_NNET)

# Extract R^2 value
r_squared_NNET <- summary_model_NNET$r.squared

# Extract p-value for the slope (flower_surface_cm2)
p_value_NNET <- summary_model_NNET$coefficients[2, 4]

# Print R^2 and p-value
print(paste("R^2:", r_squared_NNET))
print(paste("p-value:", p_value_NNET))

# Format the p-value to display three decimal places
p_value_format <- formatC(p_value_NNET, format = "f", digits = 4)

# Create annotation text
annotation_text_NNET <- paste("R^2 = ", round(r_squared_NNET, digits = 3),
                              "\np-value = ", p_value_format)

# Create the plot with ggplot2 and add the annotation
ggplot(joined_data_NNET, aes(x = flower_surface_cm2, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Field Data (cm²)", 
       y = "Drone Data NNET (cm²)", 
       title = "Flower Abundance (Field) vs Flower Abundance UAV (NNET)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_NNET, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")


# BEE ABUNDANCE / FLOWER ABUNDANCE NNET ----

# Calculate correlation between bee_abundance and flower_cm2
correlation_NNET <- cor(joined_data_NNET$bee_abundance, joined_data_NNET$flower_cm2)
print(paste("Correlation:", correlation_NNET))

# Build the linear model
linear_model_NNET <- lm(flower_cm2 ~ bee_abundance, data = joined_data_NNET)

# Model summary to obtain R^2 and p-value
summary_model_NNET <- summary(linear_model_NNET)

# Extract R^2 value
r_squared_NNET <- summary_model_NNET$r.squared
print(paste("R^2:", r_squared_NNET))

# Extract p-value for the slope (bee_abundance)
p_value_NNET <- summary_model_NNET$coefficients[2, 4]
print(paste("p-value:", p_value_NNET))

# Format p-value to display three decimal places
p_value_format <- formatC(p_value_NNET, format = "f", digits = 8)

# Create annotation text
annotation_text_NNET <- paste("R^2 = ", round(r_squared_NNET, digits = 3),
                              "\np-value = ", p_value_format)

# Visualize the relationship with ggplot2
ggplot(joined_data_NNET, aes(x = bee_abundance, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Bee Abundance", 
       y = "Flower Area (cm²)", 
       title = "Correlation between Bee Abundance and Flower Area (NNET)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_NNET, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")


# BEE SPECIES RICHNESS / FLOWER ABUNDANCE NNET ----

# Calculate the correlation between bee_species_richness and flower_cm2
correlation_richness_NNET <- cor(joined_data_NNET$bee_species_richness, joined_data_NNET$flower_cm2)

# Build the linear model for bee_species_richness vs flower_cm2
linear_model_richness_NNET <- lm(flower_cm2 ~ bee_species_richness, data = joined_data_NNET)

# Model summary to obtain R^2 and p-value
summary_model_richness_NNET <- summary(linear_model_richness_NNET)

# Extract R^2 value
r_squared_richness_NNET <- summary_model_richness_NNET$r.squared

# Extract p-value for the slope (bee_species_richness)
p_value_richness_NNET <- coef(summary(linear_model_richness_NNET))["bee_species_richness", "Pr(>|t|)"]

# Format p-value
p_value_format_richness <- formatC(p_value_richness_NNET, format = "f", digits = 4)

# Create annotation text
annotation_text_richness <- paste("R^2 = ", round(r_squared_richness_NNET, digits = 3),
                                  "\np-value = ", p_value_format_richness)

# Create improved plot for bee species richness
ggplot(joined_data_NNET, aes(x = bee_species_richness, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Bee Species Richness", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Bee Species Richness and Flower Area (NNET)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotation_text_richness, hjust = 1.8, vjust = 1, size = 5, colour = "steelblue")


# SHANNON / FLOWER ABUNDANCE NNET ----

# Build the linear model
linear_model_shannon_NNET <- lm(flower_cm2 ~ bee_shannon, data = joined_data_NNET)

# Model summary to obtain R^2 and p-value
summary_model_shannon_NNET <- summary(linear_model_shannon_NNET)
r_squared_shannon_NNET <- summary_model_shannon_NNET$r.squared
p_value_shannon_NNET <- coef(summary(linear_model_shannon_NNET))["bee_shannon", "Pr(>|t|)"]

# Format p-value
p_value_format_shannon <- formatC(p_value_shannon_NNET, format = "f", digits = 4)

# Create annotation text
annotation_text_shannon <- paste("R^2 = ", round(r_squared_shannon_NNET, digits = 3),
                                 "\np-value = ", p_value_format_shannon)

# Create improved plot for Shannon index
ggplot(joined_data_NNET, aes(x = bee_shannon, y = flower_cm2, label = area)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", col = "steelblue", linetype = "dashed", size = 1) +
  geom_text(check_overlap = TRUE, nudge_y = 5000, hjust = 0.5, vjust = 0, size = 3.5, color = "black") +
  labs(x = "Shannon Index", 
       y = "Flower Area (cm²)", 
       title = "Relationship between Shannon Index and Flower Area (NNET)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90")) +
  annotate("text", x = Inf, y = Inf, label = annotazione_testo_shannon, hjust = 1.1, vjust = 1, size = 5, colour = "steelblue")

# Save the NNET model ----
#saveRDS(NNET, file = "/media/r_projects/phd_ludovico/2021/NNET_model_22.rds")














