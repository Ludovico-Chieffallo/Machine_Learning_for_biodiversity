install.packages(c("sf", "raster", "caret", "randomForest"))

library(sf)
library(raster)
library(caret)
library(randomForest)



# Se hai un file CSV
poligoni_df <- read.csv("esempio machine learning.csv")

# Se hai un shapefile
poligoni_sf <- st_read("poligoni.shp")


prato_tif <- raster("prato.tif")


poligoni_sf_raster <- raster::extract(prato_tif, poligoni_sf)
colnames(poligoni_sf)[colnames(poligoni_sf) == "id"] <- "Etichetta"
poligoni_sf$Etichetta <- as.factor(poligoni_sf$Etichetta)
str(poligoni_sf)




valori_raster_lista <- raster::extract(prato_tif, poligoni_sf)
valori_raster <- lapply(valori_raster_lista, function(x) mean(x, na.rm = TRUE))
poligoni_sf$RasterValori <- unlist(valori_raster)


set.seed(123)
indici_addestramento <- createDataPartition(poligoni_sf$Etichetta, p = 0.75, list = FALSE)
addestramento <- poligoni_sf[indici_addestramento, ]
validazione <- poligoni_sf[-indici_addestramento, ]


controllo <- trainControl(method = "cv", number = 5)
modello_rf <- train(Etichetta ~ RasterValori, data = addestramento, method = "rf",
                    trControl = controllo, tuneGrid = data.frame(mtry = 2), ntree = 500)


predizioni <- predict(modello_rf, validazione)

livelli_etichetta <- unique(addestramento$Etichetta)
validazione$Etichetta <- factor(validazione$Etichetta, levels = livelli_etichetta)
predizioni <- factor(predizioni, levels = livelli_etichetta)


matrice_confusione <- confusionMatrix(predizioni, validazione$Etichetta)
print(matrice_confusione)


valori_raster <- as.data.frame(prato_tif, xy = TRUE)
names(valori_raster) <- c("x", "y", "RasterValori")
valori_raster


valori_raster$predizioni <- predict(modello_rf, valori_raster)


prato_classificato <- raster::rasterFromXYZ(valori_raster[, c("x", "y", "predizioni")], crs = prato_tif@crs)
