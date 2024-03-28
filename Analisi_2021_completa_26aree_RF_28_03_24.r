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
setwd("/media/data/r_projects/phd_ludovico/2021")
source("R/clean_name.R") #scorciatoia per :
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

#crop e maschera----
# Sovrapponi gli orthomosaic con i buffer
#ortho_intersect <- map2(ortho, map(buffers, ~ ext(.x)), ~ intersect(.x, .y))
#ATTENZIONE.
ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))
ortho
#Fiori ed Erba-----

# Fiori ed Erba combinati
#Sys.setenv("SHAPE_RESTORE_SHX" = "YES")
fiori_erba <- list.files(pattern = "fiori e erba.*shp", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.)) %>% 
  map(~.[!st_is_empty(.), ])

fiori <- fiori_erba %>% 
  map(~ . %>% filter(id == 1))  # Usa filter per selezionare gli elementi dove id è 1

erba <- fiori_erba %>% 
  map(~ . %>% filter(id == 2))


#Ordiniamo i nomi in modo che coincidano
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

#Values extraction----

# Converti manualmente un singolo elemento per il test
fiori <- map(fiori_erba, ~{
  # Verifica che la colonna 'id' esista
  if ("id" %in% names(.x)) {
    # Filtra solo se ci sono righe che corrispondono al criterio
    filtered <- .x %>% filter(id == 1)
    # Converte in 'sf' solo se il risultato non è vuoto
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Restituisce NULL per frame dati vuoti
    }
  } else {
    NULL  # Restituisce NULL se non esiste la colonna 'id'
  }
})

erba <- map(fiori_erba, ~{
  # Verifica che la colonna 'id' esista
  if ("id" %in% names(.x)) {
    # Filtra solo se ci sono righe che corrispondono al criterio
    filtered <- .x %>% filter(id == 2)
    # Converte in 'sf' solo se il risultato non è vuoto
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Restituisce NULL per frame dati vuoti
    }
  } else {
    NULL  # Restituisce NULL se non esiste la colonna 'id'
  }
})

#controllo per capire se ci sono aree NULL
# Crea un vettore logico che indica quali elementi in 'fiori' sono NULL
are_null_fiori <- sapply(fiori, is.null)

# Stampa gli indici degli elementi NULL in 'fiori'
null_indices_fiori <- which(are_null_fiori)
print(null_indices_fiori)

# Ripeti lo stesso processo per 'erba'
are_null_erba <- sapply(erba, is.null)
null_indices_erba <- which(are_null_erba)
print(null_indices_erba)

#estrazione dati----
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

#Bind dataframes----

# Unisci tutti i dataframe in una lista
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba,.id = "origin")


set.seed(123) # Per riproducibilitÃ 
size_of_fiori <- nrow(df_fiori_unified)
# Numero di righe da campionare da ciascun gruppo

# Crea un campione bilanciato
df_erba_unified <- df_erba_unified %>%
  group_by(origin) %>%
  sample_n(size = size_of_fiori, replace = TRUE) %>%
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
df_selected <- dplyr::select(df_tot, banda_1, banda_2, banda_3, label)


# Dividi il dataframe selezionato in set di addestramento e di test
data_split <- initial_split(df_selected, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convertiamo in un fattore-----
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

#RF Model----
train_control <- trainControl(
  method = "cv",             # Cross-validation
  number = 10,               # Numero di fold per la cross-validation
  savePredictions = "final", # Salva le predizioni per ogni modello
  classProbs = TRUE,         # Calcola le probabilitÃ  di appartenenza alle classi
  summaryFunction = twoClassSummary # Funzione di riepilogo per classificazione binaria
)

tuneGrid <- expand.grid(
  mtry = c(1,2,3),       # Numero di variabili considerate a ogni split
  splitrule = c("gini", "extratrees"),
  min.node.size = c(1, 3, 5)
)

## Train the model using ranger and the tuning grid
# Calcola i pesi delle classi
#classi <- unique(training_data$label)
#class_weights <- c("erba" = 0.1, "fiore" = 0.9)

# Rimuovi eventuali righe con valori NA
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

# Addestra il modello utilizzando i pesi delle classi
set.seed(123)
RF <- train(
  label ~ .,
  data = training_data,
  method = "ranger",
  trControl = train_control,
  tuneGrid = tuneGrid,
  metric = "ROC",
  importance = 'impurity', # Qui puoi impostare 'importance' per 'ranger'
  weights = ifelse(training_data$label == "fiore", 10, 1) # Pesi per affrontare lo squilibrio dei dati
)

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

#RF model test----
predizioni <- predict(RF, newdata= test_data)
test_accuracy <- mean(predizioni== test_data$label)
print(paste("Test Accuracy:", test_accuracy))

matrice_confusione <- confusionMatrix(as.factor(predizioni), test_data$label, positive = "fiore")
print(matrice_confusione)

accuracy <- matrice_confusione$overall['Accuracy']
recall <- matrice_confusione$byClass["Sensitivity"]
f1 <- matrice_confusione$byClass["F1"]
auc_roc<- roc(test_data$label, as.numeric(predizioni))
auc <- auc(auc_roc)
precision <- matrice_confusione$byClass["Pos Pred Value"]

print(paste("Accuracy: ", accuracy))
print(paste("Recall: ", recall))
print(paste("F1-Score: ", f1))
print(paste("AUC-ROC: ", auc))
print(paste("Precision: ", precision))



#Valutiamo l'overfitting----

#Confronto tra Prestazioni su Addestramento e Test:
# Calcola le prestazioni sul set di addestramento
train_predictions <- predict(RF, newdata = training_data)
train_conf_matrix <- confusionMatrix(as.factor(train_predictions), training_data$label, positive = "fiore")
print(train_conf_matrix)



#Cross-Validation Esterna
# Esegui k-fold cross-validation

df_selected <- na.omit(df_selected)

set.seed(123)
cv_results <- train(
  label ~ .,
  data = df_selected, # Usa l'intero dataset prima della divisione
  method = "ranger",
  trControl = trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  ),
  tuneGrid = tuneGrid,
  metric = "ROC"
)

print(cv_results)

#Verifichiamo la robustezza----
# Ecco un esempio di come ottenere l'importanza delle variabili per un modello 'ranger'

importance <- varImp(RF)
print(importance)


#Confronto delle Prestazioni su Set Diversi:
# Hai giÃ  le metriche di addestramento e di test; ora confrontale con la cross-validation
print(cv_results$results)

#Analisi della Curva ROC:
# Calcola e traccia la curva ROC per il set di test
test_probabilities <- predict(RF, test_data, type = "prob")
roc_test <- roc(response = test_data$label, test_probabilities[, "fiore"])

# Calcola l'AUC
auc(roc_test)

# Aggiungi l'AUC al grafico
plot(roc_test, main = paste("ROC Curve (AUC =", auc(roc_test), ")"))
