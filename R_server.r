#SAREBBE UTILE SCRIVERE DEL CODICE PER VERIFICARE GLI NA

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




source("R/clean_name.R") #scorciatoia per :
# clean_name <- function(file) {
#   var_name <- tools::file_path_sans_ext(basename(file))
#   var_name <- gsub("[^[:alnum:]_]", "", var_name)
#   if (grepl("^[0-9]", var_name)) {
#     var_name <- paste0("x", var_name)
#   }
#   var_name
# }

ortho <- list.files(pattern = "Ortho", full.names = TRUE, recursive = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))

buffers <- list.files(pattern = "Buffer.*shp",
                      full.names = TRUE, recursive = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))
# remove EL2, EY3, RA2 (2, 5, 10)
buffers[[10]] <- NULL
buffers[[5]] <- NULL
buffers[[2]] <- NULL

ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))

#Vogliamo sapere se lo sfondo del nostro transetto sia composto da NA per 
#escluderli successivamente

outside_na <- map2(ortho, buffers, function(ortho, buffer) {
  # Converti il buffer in SpatVector
  buffer <- terra::vect(buffer)
  
  # Crea un raster di stesse dimensioni dell'immagine originale e riempilo con 1
  base_raster <- ortho * 0 + 1
  
  # Applica la maschera all'immagine creata
  mask_inverse <- terra::mask(base_raster, buffer, inverse = TRUE)
  
  # Estrae i valori al di fuori del buffer
  outside_values <- terra::values(mask_inverse)
  
  # Verifica se questi valori sono tutti NA
  all(is.na(outside_values))
})

# Stampa i risultati
print(outside_na)


#Fiori ed Erba


fiori<- list.files(pattern = "fiori.*shp", full.names = TRUE, recursive = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))%>% 
  map(~.[!st_is_empty(.), ])

erba<-list.files(pattern = "erba.*shp", full.names = TRUE, recursive = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))%>% 
  map(~.[!st_is_empty(.), ])

#Ordiniamo i nomi in modo che coincidano
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))


#Nomi delle colonne di ogni dataframe
df_fiori <- map(df_fiori, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})

df_erba <- map(df_erba, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})


#label
df_fiori <- map(df_fiori, ~{
  .$label <- "fiore"
  .
})

df_erba <- map(df_erba, ~{
  .$label <- "erba"
  .
})

# Unisci tutti i dataframe in una lista
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba,.id = "origin")


# Numero di righe da campionare da ciascun gruppo
sample_size <- 90000

# Crea un campione bilanciato
df_erba_unified <- df_erba_unified %>%
  group_by(origin) %>%
  sample_n(size = sample_size, replace = TRUE) %>%
  ungroup()


df_tot <- rbind(df_fiori_unified, df_erba_unified) %>% 
  rowid_to_column(.) %>%
  rename(ID_pixel = rowid)

# Crea un dataframe 'lungo' #SOLO SE VEDIAMO NA
# df_long <- df_tot %>%
#   tidyr::pivot_longer(cols =c("banda_1", "banda_2", "banda_3"),
#                       names_to = "banda",
#                       values_to = "valore")
# df_long %>%
#   filter(is.na(valore)) %>%
#   view()


#Splitting into train and test----
#Seleziono le bande

set.seed(123)

# Crea un nuovo dataframe con solo le colonne che vuoi usare
df_selected <- df_tot %>%
  select(banda_1,banda_2 , banda_3, label)

# Dividi il dataframe selezionato in set di addestramento e di test
data_split <- initial_split(df_selected, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convertiamo in un fattore
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

#Model----

tuneGrid <- expand.grid(.mtry = c(1, 2, 3),
                        .splitrule = "hellinger",
                        .min.node.size = c(1, 3, 5))

# tuneGrid <- expand.grid(.mtry = c(1, 2, 3),
#                         .splitrule ="gini",
#                         .min.node.size = c(1, 3, 5))

# Define cross-validation strategy
train_control <- trainControl(method = "cv", number = 10,
                              savePredictions = "final",
                              verboseIter = TRUE)

## Train the model using ranger and the tuning grid
# Calcola i pesi delle classi
classi <- unique(training_data$label)
class_weights <- c("erba" = 0.1, "fiore" = 0.9)

# Addestra il modello utilizzando i pesi delle classi
RF <- train(label ~ ., data = training_data, 
            method = "ranger", 
            num.trees = 500,
            trControl = train_control,
            tuneGrid = tuneGrid,
            class.weights = class_weights)

# t0<-Sys.time()
# RF <- train(label ~ ., data = training_data, 
#             method = "ranger", 
#             num.trees = 500,
#             trControl = train_control,
#             tuneGrid = tuneGrid)
# t1<-Sys.time()
# t1-t0

print(RF)

#se volessi addestrare un modello SVM, potrei usare method = "svmLinear"

predizioni <- predict(RF, test_data)
matrice_confusione <- confusionMatrix(predizioni, test_data$label, positive = "fiore")
print(matrice_confusione)
accuracy <- matrice_confusione$overall['Accuracy']
recall <- matrice_confusione$byClass["Sensitivity"]
f1 <- matrice_confusione$byClass["F1"]
auc_roc<- roc(test_data$label, as.numeric(predizioni))
auc <- auc(auc_roc)
precision <- matrice_confusione$byClass["Pos Pred Value"]


print(matrice_confusione)
print(paste("Accuracy: ", accuracy))
print(paste("Recall: ", recall))
print(paste("F1-Score: ", f1))
print(paste("AUC-ROC: ", auc))
print(paste("Precision: ", precision))


#Funzione a tutti gli elementi in ortho
ortho <- map(ortho, function(x) {
  names(x) <- colnames(df_selected)[1:3]
  return(x)
})


#Application of model----
t0<-Sys.time()
classificazione_RF_list <- pbapply::pblapply(names(ortho), function(name) {
  terra::predict(ortho[[name]], RF, na.rm=TRUE) #gli NA sono dello sfondo, quindi possiamo escluderli
})
t1<-Sys.time()
t1-t0



#classificazione_RF_list <- future_map(ortho[[1]], ~predict(., RF),.progress = TRUE)


# Salva la classificazione come un file GeoTIFF
#output_filename <- "classificazione_RF.tif"
#writeRaster(classificazione_RF, output_filename, overwrite = TRUE)


names(classificazione_RF_list) <- names(ortho)  # Assicuriamoci che i nomi corrispondano


walk(names(classificazione_RF_list), function(name) {
  terra::writeRaster(classificazione_RF_list[[name]], 
                     filename = paste0(name, ".tif"), 
                     overwrite=TRUE)
})




