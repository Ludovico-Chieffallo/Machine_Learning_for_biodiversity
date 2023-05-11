#Library----
library(raster)
library(randomForest)
library(sf)
library(tidyr)
library(terra)
library(dplyr)
library(caret)
library(pROC)

#Setwd----
setwd("D:/Desktop/Prova phd ey1")
getwd()

#Load data orthomosaic----
EY1_ortho<- rast("Ortho EY1.tif")
EY1_ortho
#Crop
area_buffer<-st_read("Buffer EY1.shp")

EY1_ortho<-crop(EY1_ortho, extent(area_buffer))
EY1_ortho<-mask(EY1_ortho, area_buffer)

#Load Flowers and Grass----
poligoni_fiori <- st_read("fiori EY1.shp")
poligoni_erba <- st_read("erba EY1.shp")

nrow(poligoni_erba)
nrow(poligoni_fiori)

#Remove empty geometries----
poligoni_fiori <- poligoni_fiori[!st_is_empty(poligoni_fiori), ]
nrow(poligoni_fiori)

poligoni_erba<- poligoni_erba[!st_is_empty(poligoni_erba), ]
nrow(poligoni_erba)

#Extraction for pixel values---- 
valori_fiori <- extract(EY1_ortho, poligoni_fiori)
valori_erba <- extract(EY1_ortho, poligoni_erba)
class(valori_erba)

#Trasform in dataframe----
df_fiori <- do.call(rbind, lapply(valori_fiori, as.data.frame))
df_erba <- do.call(rbind, lapply(valori_erba, as.data.frame))

colnames(df_fiori) <- paste0("banda_", 1:ncol(df_fiori))
colnames(df_erba) <- paste0("banda_", 1:ncol(df_erba))

#Add label (fiori e erba)----
df_fiori$label <- "fiore"
df_erba$label <- "erba"

# Imposta il numero di righe da estrarre
num_righe_campione <- 1000
num_righe_campione2<- 300000
# Estrai un campione casuale di righe dal dataframe
df_fiori <- df_fiori %>% sample_n(num_righe_campione)
df_erba <- df_erba %>% sample_n(num_righe_campione2)

#Creation of single dataset----
dataset <- rbind(df_fiori, df_erba)
dataset$label <- as.factor(dataset$label)


#Training data and testing data----
# Scegli un sottoinsieme casuale dei dati
set.seed(123) # Imposta un seed per la riproducibilità
indici_subset <- sample(nrow(dataset), size = floor(1 * nrow(dataset)))

# Crea il subset
subset_data <- dataset[indici_subset, ]

# Creazione del set di dati di training e test per il subset
set.seed(123) # Imposta un seed per la riproducibilità
indici_training_subset <- sample(nrow(subset_data), size = floor(0.8 * nrow(subset_data)))
training_data_subset <- subset_data[indici_training_subset, ]
training_data_subset <- na.omit(training_data_subset)
training_data_subset
test_data_subset <- subset_data[-indici_training_subset, ]

# Addestramento del modello sul subset
RF<- randomForest(label ~ ., data = training_data_subset)


#Evaluation of model----
predizioni <- predict(RF, test_data_subset)
matrice_confusione <- table(test_data_subset$label, predizioni)
accuracy <- sum(diag(matrice_confusione)) / sum(matrice_confusione)
print(accuracy)

#Rename columns for matches with dataset----
names(EY1_ortho) <- colnames(dataset)[1:3]

#Application of model----
classificazione_RF <- predict(EY1_ortho, RF)
plot(classificazione_RF)


# Salva la classificazione come un file GeoTIFF
output_fil.ename <- "classificazione_RF.tif"
writeRaster(classificazione_RF, output_filename, overwrite = TRUE)

########################################################



# Definizione della griglia di parametri per la cross-validation
param_grid <- expand.grid(mtry = 1)

# Configurazione del controllo della cross-validation
control <- trainControl(method = "cv", number = 5)

# Addestramento del modello con cross-validation
model_cv <- train(label ~ ., data = training_data_subset, method = "rf", trControl = control, tuneGrid = param_grid)

# Stampa dei risultati della cross-validation
print(model_cv)

# Utilizzo del modello ottimizzato per la valutazione
predizioni_cv <- predict(model_cv, test_data_subset)

print(length(test_data_subset$label))
print(length(predizioni_cv))
test_data_subset <- na.omit(test_data_subset)


matrice_confusione_cv <- table(test_data_subset$label, predizioni_cv)
accuracy_cv <- sum(diag(matrice_confusione_cv)) / sum(matrice_confusione_cv)
print(accuracy_cv)

###################################

# Calcola le probabilità di previsione
prob_pred <- predict(model_cv, newdata = test_data_subset, type = "prob")

# Calcola la previsione delle classi
class_pred <- predict(model_cv, newdata = test_data_subset)

# Calcola la confusione Matrix
cm <- confusionMatrix(class_pred, test_data_subset$label)

# Ottieni precisione, recall e F1-score
precision <- posPredValue(class_pred, test_data_subset$label)
recall <- sensitivity(class_pred, test_data_subset$label)
F1 <- (2 * precision * recall) / (precision + recall)

# Calcola l'AUC-ROC
roc_obj <- roc(test_data_subset$label, prob_pred[,2]) # assumendo che la seconda colonna corrisponda alla classe "positiva"
auc <- auc(roc_obj)

print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1-score:", F1))
print(paste("AUC-ROC:", auc))

##########################

install.packages("h2o")
install.packages("bit64")
library(bit64)
library(h2o)
# Inizializza un'istanza di H2O
h2o.init(nthreads = -1) # Usa tutti i core disponibili

# Converti il tuo training_data_subset e test_data_subset in H2OFrame
training_hf <- as.h2o(training_data_subset)
testing_hf <- as.h2o(test_data_subset)

# Identifica i predittori e la risposta
predictors <- setdiff(names(training_data_subset), "label")
response <- "label"

# Esegui AutoML per 20 modelli (o il massimo possibile in 30 minuti)
aml <- h2o.automl(x = predictors,
                  y = response,
                  training_frame = training_hf,
                  leaderboard_frame = testing_hf,
                  max_models = 20,
                  max_runtime_secs = 1800,
                  seed = 1)

# Visualizza il leaderboard
lb <- aml@leaderboard
print(lb)

# Previsioni sul set di test
pred <- h2o.predict(aml@leader, newdata = testing_hf)

# Calcola l'accuratezza
accuracy_automl <- sum(pred$predict == testing_hf$label) / nrow(testing_hf)
print(paste("Accuracy AutoML:", accuracy_automl))

# Fermare l'istanza H2O
h2o.shutdown(prompt = FALSE)
