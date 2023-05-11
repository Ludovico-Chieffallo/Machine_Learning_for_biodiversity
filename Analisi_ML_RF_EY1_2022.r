#Library----
library(raster)
library(randomForest)
library(sf)
library(tidyr)
library(terra)
library(dplyr)

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

